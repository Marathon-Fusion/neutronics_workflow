#!/usr/bin/env python3
"""
Nuclear Reaction Chain Analyzer

A tool for finding nuclear reaction and decay chains that produce specific nuclides.
Supports filtering by half-life, reaction types, and chain length.
"""

import math
import re
from typing import Dict, List, Set, Tuple, Optional, Union
from dataclasses import dataclass, field
from enum import Enum
from collections import defaultdict, deque
import heapq
import numpy as np

# Import centralized element registry
from element_registry import get_symbol, get_display_name


class ReactionType(Enum):
    """Enumeration of nuclear reaction types"""
    NEUTRON_CAPTURE = "n,gamma"
    NEUTRON_PROTON = "n,p"
    NEUTRON_ALPHA = "n,alpha"
    NEUTRON_DEUTERON = "n,d"
    NEUTRON_2N = "n,2n"
    NEUTRON_3N = "n,3n"
    NEUTRON_FISSION = "n,f"
    BETA_MINUS = "beta-"
    BETA_PLUS = "beta+"
    ELECTRON_CAPTURE = "EC"
    ALPHA_DECAY = "alpha"
    SPONTANEOUS_FISSION = "SF"
    ISOMERIC_TRANSITION = "IT"
    NEUTRON_EMISSION = "n"
    PROTON_EMISSION = "p"
    NONELASTIC_SCATTER = 'n,non'
    ELASTIC_SCATTER = 'n,el'
    NEUTRON_TRITON = 'n,t'


@dataclass(frozen=True)
class Nuclide:
    """
    Represents a nuclear isotope with identity based only on (z, a, m)
    
    The nuclide identity (equality, hashing) depends ONLY on (z, a, m).
    Symbol is derived automatically from the centralized element registry
    and serves as display metadata only.
    
    Attributes:
        z: Atomic number (number of protons)
        a: Mass number (protons + neutrons) 
        m: Metastable state (0 for ground state, 1+ for excited states)
        is_stable: True if nuclide is stable (no spontaneous decay), False if radioactive
        _symbol: Cached symbol for display (automatically generated)
    """
    z: int
    a: int
    m: int = 0
    is_stable: Optional[bool] = None
    # Symbol is now derived automatically and cached
    _symbol: str = field(default="", init=False, compare=False, repr=False)
    
    def __post_init__(self):
        """Initialize derived fields after creation"""
        # Set symbol from centralized registry (not part of identity)
        object.__setattr__(self, '_symbol', get_symbol(self.z))
    
    @property
    def symbol(self) -> str:
        """Get element symbol (display only, not part of identity)"""
        if not self._symbol:
            # Lazy initialization if not set
            object.__setattr__(self, '_symbol', get_symbol(self.z))
        return self._symbol
    
    def __eq__(self, other) -> bool:
        """Equality based ONLY on (z, a, m) - symbol ignored"""
        if not isinstance(other, Nuclide):
            return False
        return (self.z, self.a, self.m) == (other.z, other.a, other.m)
    
    def __hash__(self) -> int:
        """Hash based ONLY on (z, a, m) - symbol ignored"""
        return hash((self.z, self.a, self.m))
    
    def __str__(self) -> str:
        """String representation using centralized display formatting"""
        base_name = get_display_name(self.z, self.a, self.m)
        stability_suffix = ""
        if self.is_stable is not None:
            stability_suffix = " (stable)" if self.is_stable else " (radioactive)"
        return f"{base_name}{stability_suffix}"
    
    def __repr__(self) -> str:
        """Representation showing core identity (z, a, m) and stability"""
        return f"Nuclide(z={self.z}, a={self.a}, m={self.m}, is_stable={self.is_stable})"
    
    def get_identity_key(self) -> Tuple[int, int, int]:
        """Get the identity key tuple (z, a, m) for database operations"""
        return (self.z, self.a, self.m)
    
    @classmethod
    def from_identity(cls, z: int, a: int, m: int = 0, is_stable: Optional[bool] = None) -> 'Nuclide':
        """Create nuclide from identity components (recommended constructor)"""
        return cls(z=z, a=a, m=m, is_stable=is_stable)
    
    @classmethod 
    def from_legacy(cls, z: int, a: int, m: int = 0, symbol: str = "", is_stable: Optional[bool] = None) -> 'Nuclide':
        """Create nuclide with legacy symbol parameter (for backward compatibility)"""
        # Symbol parameter is ignored - always derived from registry
        return cls(z=z, a=a, m=m, is_stable=is_stable)



@dataclass
class Reaction:
    """
    Represents a nuclear reaction or decay process
    
    Attributes:
        reactant: Initial nuclide
        product: Final nuclide
        reaction_type: Type of nuclear process
        cross_section: Cross section in barns (deprecated - use tabulated_xs)
        half_life: Half-life in seconds (for decay reactions)
        branching_ratio: Fraction of reactions of this type (0-1)
        threshold_energy: Minimum neutron energy in eV (for neutron reactions)
        q_value: Energy released/absorbed in eV (JENDL units - will be converted to MeV for display)
        tabulated_xs: Raw tabulated cross section data by temperature
        has_tabulated_data: Whether tabulated cross section data is available
    """
    reactant: Nuclide
    product: Nuclide
    reaction_type: ReactionType
    cross_section: float = 0.0  # barns (deprecated - use tabulated_xs)
    half_life: Optional[float] = None  # seconds
    branching_ratio: float = 1.0  # fraction
    threshold_energy: float = 0.0  # eV
    q_value: float = 0.0  # eV (JENDL native units)
    tabulated_xs: Optional[Dict[str, Tuple[List[float], List[float]]]] = None  # {temp: ([energies], [xs_values])}
    has_tabulated_data: bool = False
    
    @property
    def thermal_cross_section(self) -> float:
        """Cross section at thermal energy (0.0253 eV)"""
        return self.cross_section
    
    @property
    def fast_cross_section(self) -> float:
        """Fast neutron cross section - not available in simplified format"""
        return 0.0
    
    @property
    def q_value_mev(self) -> float:
        """Q-value in MeV (converted from eV)"""
        return self.q_value * 1e-6  # Convert eV to MeV
    
    def is_decay(self) -> bool:
        """Check if this is a decay reaction"""
        decay_types = {ReactionType.BETA_MINUS, ReactionType.BETA_PLUS, 
                      ReactionType.ELECTRON_CAPTURE, ReactionType.ALPHA_DECAY,
                      ReactionType.SPONTANEOUS_FISSION, ReactionType.ISOMERIC_TRANSITION,
                      ReactionType.NEUTRON_EMISSION, ReactionType.PROTON_EMISSION}
        return self.reaction_type in decay_types
    
    def is_neutron_induced(self) -> bool:
        """Check if this is a neutron-induced reaction"""
        neutron_types = {ReactionType.NEUTRON_CAPTURE, ReactionType.NEUTRON_PROTON,
                        ReactionType.NEUTRON_ALPHA, ReactionType.NEUTRON_DEUTERON,
                        ReactionType.NEUTRON_2N, ReactionType.NEUTRON_3N, 
                        ReactionType.NEUTRON_FISSION, ReactionType.INELASTIC_SCATTER, ReactionType.ELASTIC_SCATTER}
        return self.reaction_type in neutron_types
    
    def __str__(self) -> str:
        half_life_str = f", t½={self.half_life:.2e}s" if self.half_life else ""
        q_value_str = f", Q={self.q_value_mev:.2f}MeV" if self.q_value != 0.0 else ""
        xs_str = ""
        if self.cross_section > 0:
            xs_str = f", σ={self.cross_section:.2e}b"
        return f"{self.reactant} --{self.reaction_type.value}--> {self.product}{half_life_str}{q_value_str}{xs_str}"
    
    @property
    def energy_released_mev(self) -> float:
        """Energy released in MeV (positive for exothermic)"""
        return self.q_value_mev
    
    @property  
    def energy_released_joules(self) -> float:
        """Energy released in Joules"""
        return self.q_value_mev * 1.602176634e-13  # MeV to Joules
    
    @property
    def is_exothermic(self) -> bool:
        """Check if reaction releases energy"""
        return self.q_value_mev > 0.0
    
    @property
    def is_endothermic(self) -> bool:
        """Check if reaction requires energy input"""
        return self.q_value_mev < 0.0


@dataclass
class ReactionChain:
    """
    Represents a complete nuclear reaction chain
    
    Attributes:
        reactions: List of reactions in sequence
        total_probability: Combined probability of entire chain
        total_time: Estimated time for chain completion
    """
    reactions: List[Reaction]
    total_probability: float = 1.0
    total_time: float = 0.0
    
    def __post_init__(self):
        """Calculate chain statistics"""
        self._calculate_probability()
        self._calculate_time()
    
    def _calculate_probability(self):
        """Calculate total chain probability"""
        prob = 1.0
        for reaction in self.reactions:
            prob *= reaction.branching_ratio
        self.total_probability = prob
    
    def _calculate_time(self):
        """Estimate total chain time (simplified)"""
        time = 0.0
        for reaction in self.reactions:
            if reaction.is_decay() and reaction.half_life:
                # Use half-life as rough time estimate
                time += reaction.half_life
            elif reaction.is_neutron_induced():
                # Assume instant for neutron reactions (flux dependent)
                time += 0.0
        self.total_time = time
    
    @property
    def initial_nuclide(self) -> Nuclide:
        """Get the starting nuclide"""
        return self.reactions[0].reactant if self.reactions else None
    
    @property
    def final_nuclide(self) -> Nuclide:
        """Get the ending nuclide"""
        return self.reactions[-1].product if self.reactions else None
    
    @property
    def length(self) -> int:
        """Get chain length"""
        return len(self.reactions)
    
    @property
    def total_energy_released_mev(self) -> float:
        """Total energy released by the entire chain in MeV (sum of all reaction Q-values)"""
        return sum(reaction.q_value for reaction in self.reactions)
    
    @property  
    def total_energy_released_joules(self) -> float:
        """Total energy released by the entire chain in Joules"""
        return self.total_energy_released_mev * 1.602176634e-13  # MeV to Joules conversion
    
    @property
    def is_net_exothermic(self) -> bool:
        """Check if the entire chain releases net energy (positive total Q-value)"""
        return self.total_energy_released_mev > 0.0
    
    @property
    def is_net_endothermic(self) -> bool:
        """Check if the entire chain requires net energy input (negative total Q-value)"""
        return self.total_energy_released_mev < 0.0
    
    @property
    def energy_classification(self) -> str:
        """Return 'Exothermic' if net energy releasing, 'Endothermic' if net energy consuming, or 'Neutral' if zero"""
        total_energy = self.total_energy_released_mev
        if total_energy > 0.0:
            return "Exothermic"
        elif total_energy < 0.0:
            return "Endothermic"
        else:
            return "Neutral"
    
    @property
    def energy_per_step(self) -> float:
        """Average energy released per reaction step in MeV"""
        return self.total_energy_released_mev / len(self.reactions) if self.reactions else 0.0

    def __str__(self) -> str:
        chain_str = " -> ".join([str(r.reactant) for r in self.reactions])
        if self.reactions:
            chain_str += f" -> {self.reactions[-1].product}"
        return f"Chain: {chain_str} (P={self.total_probability:.2e}, t={self.total_time:.2e}s, ΔE={self.total_energy_released_mev:.3f}MeV)"


class HalfLifeFilter:
    """Filter reactions based on half-life constraints"""
    
    def __init__(self, min_half_life: Optional[float] = None, 
                 max_half_life: Optional[float] = None):
        """
        Initialize half-life filter
        
        Args:
            min_half_life: Minimum half-life in seconds
            max_half_life: Maximum half-life in seconds
        """
        self.min_half_life = min_half_life
        self.max_half_life = max_half_life
    
    def passes_filter(self, reaction: Reaction) -> bool:
        """Check if reaction passes half-life filter"""
        if not reaction.is_decay() or reaction.half_life is None:
            return True  # Non-decay reactions always pass
        
        if self.min_half_life is not None and reaction.half_life < self.min_half_life:
            return False
        
        if self.max_half_life is not None and reaction.half_life > self.max_half_life:
            return False
        
        return True


# Test the core data structures
def test_nuclide():
    """Test Nuclide class functionality with new identity-based system"""
    print("Testing Nuclide class...")
    
    # Test basic creation using new system
    u235 = Nuclide.from_identity(92, 235, 0)
    assert u235.z == 92
    assert u235.a == 235
    assert u235.m == 0
    assert u235.symbol == 'U'  # Symbol derived automatically
    print(f"✓ Basic creation: {u235}")
    
    # Test metastable state
    tc99m = Nuclide.from_identity(43, 99, 1)
    assert tc99m.m == 1
    assert tc99m.symbol == 'Tc'  # Symbol derived automatically
    print(f"✓ Metastable state: {tc99m}")
    
    # Test string representation
    assert str(u235) == "U-235"
    assert str(tc99m) == "Tc-99m1"
    print("✓ String representation works")
    
    # Test NEW BEHAVIOR: equality ignores symbol differences
    u235_copy = Nuclide.from_identity(92, 235, 0)
    u235_legacy = Nuclide.from_legacy(92, 235, 0, 'URANIUM')  # Different symbol ignored
    assert u235 == u235_copy
    assert u235 == u235_legacy  # Same (z,a,m) = equal despite different symbol param
    assert hash(u235) == hash(u235_copy)
    assert hash(u235) == hash(u235_legacy)
    print("✓ New identity-based equality works (symbol ignored)")
    
    # Test identity key
    assert u235.get_identity_key() == (92, 235, 0)
    print("✓ Identity key extraction works")
    
    # Test that symbol is always consistent for same Z
    u238 = Nuclide.from_identity(92, 238, 0)
    assert u235.symbol == u238.symbol == 'U'  # Both uranium isotopes have same symbol
    print("✓ Symbol consistency across isotopes")
    
    # Test set behavior - should deduplicate based on (z,a,m) only
    nuclide_set = {u235, u235_copy, u235_legacy}
    assert len(nuclide_set) == 1  # All three are the same nuclide
    print("✓ Set deduplication works with new identity system")
    
    print("Nuclide tests passed! ✅\n")


def test_reaction():
    """Test Reaction class functionality"""
    print("Testing Reaction class...")
    
    # Test neutron capture
    u238 = Nuclide.from_identity(92, 238, 0)
    u239 = Nuclide.from_identity(92, 239, 0)
    capture = Reaction(u238, u239, ReactionType.NEUTRON_CAPTURE, 
                      cross_section=2.7, branching_ratio=1.0)
    
    assert capture.is_neutron_induced()
    assert not capture.is_decay()
    print(f"✓ Neutron capture: {capture}")
    
    # Test decay reaction
    np239 = Nuclide.from_identity(93, 239, 0)
    pu239 = Nuclide.from_identity(94, 239, 0)
    decay = Reaction(np239, pu239, ReactionType.BETA_MINUS,
                    half_life=2.356e5, branching_ratio=1.0)  # ~2.35 days
    
    assert decay.is_decay()
    assert not decay.is_neutron_induced()
    print(f"✓ Beta decay: {decay}")
    
    print("Reaction tests passed!\n")


def test_reaction_chain():
    """Test ReactionChain class functionality"""
    print("Testing ReactionChain class...")
    
    # Create a simple 2-step chain: U-238 + n -> U-239 -> Np-239 -> Pu-239
    u238 = Nuclide.from_identity(92, 238, 0)
    u239 = Nuclide.from_identity(92, 239, 0)
    np239 = Nuclide.from_identity(93, 239, 0)
    pu239 = Nuclide.from_identity(94, 239, 0)
    
    reactions = [
        Reaction(u238, u239, ReactionType.NEUTRON_CAPTURE, 
                cross_section=2.7, branching_ratio=1.0),
        Reaction(u239, np239, ReactionType.BETA_MINUS,
                half_life=1.408e3, branching_ratio=1.0),  # ~23 min
        Reaction(np239, pu239, ReactionType.BETA_MINUS,
                half_life=2.356e5, branching_ratio=1.0)   # ~2.35 days
    ]
    
    chain = ReactionChain(reactions)
    
    assert chain.length == 3
    assert chain.initial_nuclide == u238
    assert chain.final_nuclide == pu239
    assert chain.total_probability == 1.0  # All branching ratios = 1.0
    print(f"✓ Chain creation: {chain}")
    
    print("ReactionChain tests passed!\n")


def test_half_life_filter():
    """Test HalfLifeFilter functionality"""
    print("Testing HalfLifeFilter class...")
    
    # Create test reactions with different half-lives
    h1 = Nuclide.from_identity(1, 1)
    h2 = Nuclide.from_identity(1, 2)
    short_decay = Reaction(h1, h2, ReactionType.BETA_MINUS,
                          half_life=1e2)  # 100 seconds
    medium_decay = Reaction(h1, h2, ReactionType.BETA_MINUS,
                           half_life=1e6)  # ~11.6 days
    long_decay = Reaction(h1, h2, ReactionType.BETA_MINUS,
                         half_life=1e10)  # ~317 years
    neutron_reaction = Reaction(h1, h2, ReactionType.NEUTRON_CAPTURE)
    
    # Test minimum half-life filter
    min_filter = HalfLifeFilter(min_half_life=1e5)  # > ~1.16 days
    assert not min_filter.passes_filter(short_decay)
    assert min_filter.passes_filter(medium_decay)
    assert min_filter.passes_filter(long_decay)
    assert min_filter.passes_filter(neutron_reaction)  # Non-decay always passes
    print("✓ Minimum half-life filter works")
    
    # Test maximum half-life filter
    max_filter = HalfLifeFilter(max_half_life=1e7)  # < ~116 days
    assert max_filter.passes_filter(short_decay)
    assert max_filter.passes_filter(medium_decay)
    assert not max_filter.passes_filter(long_decay)
    assert max_filter.passes_filter(neutron_reaction)
    print("✓ Maximum half-life filter works")
    
    # Test range filter
    range_filter = HalfLifeFilter(min_half_life=1e5, max_half_life=1e7)
    assert not range_filter.passes_filter(short_decay)
    assert range_filter.passes_filter(medium_decay)
    assert not range_filter.passes_filter(long_decay)
    print("✓ Range half-life filter works")
    
    print("HalfLifeFilter tests passed!\n")


if __name__ == "__main__":
    # Run all tests
    test_nuclide()
    test_reaction()
    test_reaction_chain()
    test_half_life_filter()
    print("All core data structure tests passed! ✅")