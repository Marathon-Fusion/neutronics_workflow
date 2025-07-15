#!/usr/bin/env python3
"""
JENDL Data Loader for Nuclear Chain Analyzer

Loads neutron cross-sections and decay data from JENDL files using OpenMC's data API.
Provides comprehensive MT reaction extraction with energy information.
"""

import os
import re
import math
from typing import Dict, List, Set, Tuple, Optional, Union
from dataclasses import dataclass, field
from pathlib import Path
try:
    import openmc.data
    OPENMC_AVAILABLE = True
except ImportError:
    raise ImportError(
        "OpenMC is required for JENDL data processing. Please install OpenMC:\n"
        "  pip install openmc\n"
        "or use conda:\n"
        "  conda install -c conda-forge openmc\n\n"
        "For analysis of existing .pkl files, use simple_pickle_loader.py instead."
    )

from nuclear_chain_analyzer import Nuclide, Reaction, ReactionType
import numpy as np


# Comprehensive MT number to reaction type mapping
MT_TO_REACTION_TYPE = {
    # Basic neutron reactions
    102: ReactionType.NEUTRON_CAPTURE,     # (n,gamma)
    103: ReactionType.NEUTRON_PROTON,      # (n,p)
    107: ReactionType.NEUTRON_ALPHA,       # (n,alpha)
    104: ReactionType.NEUTRON_DEUTERON,    # (n,d)
    105: ReactionType.NEUTRON_TRITON,      # (n,t)
    106: "NEUTRON_HELIUM3",                # (n,3He)
    
    # Multiple neutron emission
    16: ReactionType.NEUTRON_2N,           # (n,2n)
    17: ReactionType.NEUTRON_3N,           # (n,3n)
    37: "NEUTRON_4N",                      # (n,4n)
    152: "NEUTRON_5N",                     # (n,5n)
    153: "NEUTRON_6N",                     # (n,6n)
    160: "NEUTRON_7N",                     # (n,7n)
    161: "NEUTRON_8N",                     # (n,8n)
    
    # Multiple particle production
    28: "NEUTRON_NP",                      # (n,np)
    22: "NEUTRON_NA",                      # (n,na)
    32: "NEUTRON_ND",                      # (n,nd)
    33: "NEUTRON_NT",                      # (n,nt)
    34: "NEUTRON_N3HE",                    # (n,n3He)
    
    # Multiple charged particles
    111: "NEUTRON_2P",                     # (n,2p)
    108: "NEUTRON_2A",                     # (n,2a)
    112: "NEUTRON_PA",                     # (n,pa)
    115: "NEUTRON_PD",                     # (n,pd)
    116: "NEUTRON_PT",                     # (n,pt)
    117: "NEUTRON_DA",                     # (n,da)
    
    # Complex reactions
    41: "NEUTRON_2NP",                     # (n,2np)
    42: "NEUTRON_3NP",                     # (n,3np)
    44: "NEUTRON_N2P",                     # (n,n2p)
    45: "NEUTRON_NPA",                     # (n,npa)
}

# Extended ReactionType enum values (these would be added to the enum)
EXTENDED_REACTION_TYPES = {
    "NEUTRON_PROTON": "n,p",
    "NEUTRON_ALPHA": "n,alpha", 
    "NEUTRON_DEUTERON": "n,d",
    "NEUTRON_TRITON": "n,t",
    "NEUTRON_HELIUM3": "n,3He",
    "NEUTRON_4N": "n,4n",
    "NEUTRON_5N": "n,5n",
    "NEUTRON_6N": "n,6n",
    "NEUTRON_7N": "n,7n",
    "NEUTRON_8N": "n,8n",
    "NEUTRON_NP": "n,np",
    "NEUTRON_NA": "n,na",
    "NEUTRON_ND": "n,nd",
    "NEUTRON_NT": "n,nt",
    "NEUTRON_N3HE": "n,n3He",
    "NEUTRON_2P": "n,2p",
    "NEUTRON_2A": "n,2a",
    "NEUTRON_PA": "n,pa",
    "NEUTRON_PD": "n,pd",
    "NEUTRON_PT": "n,pt",
    "NEUTRON_DA": "n,da",
    "NEUTRON_2NP": "n,2np",
    "NEUTRON_3NP": "n,3np",
    "NEUTRON_N2P": "n,n2p",
    "NEUTRON_NPA": "n,npa",
}

# ENDF decay mode to ReactionType mapping
DECAY_MODE_MAPPING = {
    'beta-': ReactionType.BETA_MINUS,
    'beta+': ReactionType.BETA_PLUS,
    'ec': ReactionType.ELECTRON_CAPTURE,
    'alpha': ReactionType.ALPHA_DECAY,
    'sf': ReactionType.SPONTANEOUS_FISSION,
    'it': ReactionType.ISOMERIC_TRANSITION,
    'IT': ReactionType.ISOMERIC_TRANSITION,     # Uppercase variant
    'n': ReactionType.NEUTRON_EMISSION,
    'p': ReactionType.PROTON_EMISSION,
    # Also add uppercase variants for common decay modes
    'BETA-': ReactionType.BETA_MINUS,
    'BETA+': ReactionType.BETA_PLUS,
    'EC': ReactionType.ELECTRON_CAPTURE,
    'ALPHA': ReactionType.ALPHA_DECAY,
    'SF': ReactionType.SPONTANEOUS_FISSION,
    'N': ReactionType.NEUTRON_EMISSION,
    'P': ReactionType.PROTON_EMISSION,
    # Combined decay modes
    'ec/beta+': ReactionType.ELECTRON_CAPTURE,  # Treat as EC (dominant mode)
    'beta+/ec': ReactionType.BETA_PLUS,         # Treat as beta+ (dominant mode)
    'beta-/alpha': ReactionType.BETA_MINUS,     # Treat as beta- (dominant mode)
    'alpha/beta-': ReactionType.ALPHA_DECAY,    # Treat as alpha (dominant mode)
    'beta-/n': ReactionType.BETA_MINUS,         # Treat as beta- (dominant mode)
    'n/beta-': ReactionType.NEUTRON_EMISSION,   # Treat as neutron emission (dominant mode)
    'beta+/p': ReactionType.BETA_PLUS,          # Treat as beta+ (dominant mode)
    'p/beta+': ReactionType.PROTON_EMISSION,    # Treat as proton emission (dominant mode)
}


@dataclass
class JENDLReaction(Reaction):
    """Enhanced reaction class with JENDL data"""
    mt_number: int = 0                     # ENDF MT reaction number
    daughter_products: List[Nuclide] = field(default_factory=list)
    
    def __post_init__(self):
        """Calculate energy information after initialization"""
        # Note: is_exothermic is a property, no need to set it explicitly
        pass
    
    @property
    def energy_released_mev(self) -> float:
        """Energy released in MeV (positive for exothermic)"""
        return self.q_value_mev  # Use the inherited property that converts eV to MeV
    
    @property
    def energy_released_joules(self) -> float:
        """Energy released in Joules"""
        return self.q_value_mev * 1.602176634e-13  # MeV to Joules
    
    def __str__(self) -> str:
        base_str = super().__str__()
        mt_str = f", MT={self.mt_number}" if self.mt_number != 0 else ""
        return base_str + mt_str


def determine_products_from_mt(mt_number: int, parent_nuclide: Nuclide, available_nuclide_objects: Dict[Tuple[int, int, int], Nuclide] = None) -> List[Nuclide]:
    """
    Determine daughter products from MT number and parent nuclide
    Only references products that exist in JENDL source files - NO NUCLIDE CREATION
    
    Args:
        mt_number: ENDF MT reaction number
        parent_nuclide: Parent nuclide undergoing reaction
        available_nuclide_objects: Dict mapping (z, a, m) -> existing Nuclide objects from JENDL files
        
    Returns:
        List of daughter nuclides (only existing ones from JENDL files - no creation)
    """
    z, a = parent_nuclide.z, parent_nuclide.a
    
    def find_existing_nuclide(new_z, new_a):
        """Helper to find existing nuclide from JENDL files - NO CREATION"""
        product_identity = (new_z, new_a, 0)
        
        # Only return nuclide if it exists in JENDL files
        if available_nuclide_objects is not None and product_identity in available_nuclide_objects:
            return available_nuclide_objects[product_identity]
            
        return None
    
    # Product determination based on MT number - ONLY find existing nuclides
    products = []
    
    if mt_number == 102:  # (n,gamma)
        product = find_existing_nuclide(z, a + 1)
        if product: products = [product]
    elif mt_number == 103:  # (n,p)
        product = find_existing_nuclide(z - 1, a)
        if product: products = [product]
    elif mt_number == 107:  # (n,alpha)
        product = find_existing_nuclide(z - 2, a - 3)
        if product: products = [product]
    elif mt_number == 104:  # (n,d)
        product = find_existing_nuclide(z - 1, a - 1)
        if product: products = [product]
    elif mt_number == 105:  # (n,t)
        product = find_existing_nuclide(z - 1, a - 2)
        if product: products = [product]
    elif mt_number == 106:  # (n,3He)
        product = find_existing_nuclide(z - 2, a - 2)
        if product: products = [product]
    elif mt_number == 16:  # (n,2n)
        product = find_existing_nuclide(z, a - 1)
        if product: products = [product]
    elif mt_number == 17:  # (n,3n)
        product = find_existing_nuclide(z, a - 2)
        if product: products = [product]
    elif mt_number == 37:  # (n,4n)
        product = find_existing_nuclide(z, a - 3)
        if product: products = [product]
    elif mt_number == 152:  # (n,5n)
        product = find_existing_nuclide(z, a - 4)
        if product: products = [product]
    elif mt_number == 153:  # (n,6n)
        product = find_existing_nuclide(z, a - 5)
        if product: products = [product]
    elif mt_number == 28:  # (n,np)
        product = find_existing_nuclide(z - 1, a - 1)
        if product: products = [product]
    elif mt_number == 22:  # (n,na)
        product = find_existing_nuclide(z - 2, a - 4)
        if product: products = [product]
    elif mt_number == 32:  # (n,nd)
        product = find_existing_nuclide(z - 1, a - 2)
        if product: products = [product]
    elif mt_number == 33:  # (n,nt)
        product = find_existing_nuclide(z - 1, a - 3)
        if product: products = [product]
    elif mt_number == 34:  # (n,n3He)
        product = find_existing_nuclide(z - 2, a - 3)
        if product: products = [product]
    elif mt_number == 111:  # (n,2p)
        product = find_existing_nuclide(z - 2, a - 1)
        if product: products = [product]
    elif mt_number == 108:  # (n,2a)
        product = find_existing_nuclide(z - 4, a - 7)
        if product: products = [product]
    elif mt_number == 112:  # (n,pa)
        product = find_existing_nuclide(z - 3, a - 4)
        if product: products = [product]
    elif mt_number == 115:  # (n,pd)
        product = find_existing_nuclide(z - 2, a - 2)
        if product: products = [product]
    elif mt_number == 116:  # (n,pt)
        product = find_existing_nuclide(z - 2, a - 3)
        if product: products = [product]
    elif mt_number == 117:  # (n,da)
        product = find_existing_nuclide(z - 3, a - 5)
        if product: products = [product]
    else:
        # Unknown or complex reaction - return empty list
        return []
    
    # Return only existing nuclides (no None values since we only add existing products)
    return products


class JENDLDataLoader:
    """
    Loads nuclear data from JENDL files using OpenMC's data API
    """
    
    def __init__(self, jendl_path: str = "/opt/marathon/workspace/jendl/"):
        """
        Initialize JENDL data loader
        
        Args:
            jendl_path: Path to JENDL data directory
        """
        self.jendl_path = Path(jendl_path)
        self.neutron_path = self.jendl_path / "jendl5-n" / "jendl5-n"
        self.decay_path = self.jendl_path / "jendl5-dec_upd5" / "jendl5-dec_upd5"
        
        # Cache for loaded data
        self._neutron_data_cache = {}
        self._decay_data_cache = {}
        self._available_nuclides = None
        
        # Tracking for cross section data extraction
        self.non_tabulated_reactions = []  # Track reactions without tabulated data
        self.xs_extraction_stats = {
            'total_reactions': 0,
            'tabulated_reactions': 0,
            'non_tabulated_reactions': 0,
            'extraction_errors': 0
        }
        
        # Element symbols now handled by centralized element_registry module
        
        print(f"JENDL Data Loader initialized:")
        print(f"  Neutron data path: {self.neutron_path}")
        print(f"  Decay data path: {self.decay_path}")
    
    def parse_filename(self, filename: str) -> Optional[Tuple[int, int, int]]:
        """
        Parse JENDL filename to extract Z, A, and metastable state
        
        Args:
            filename: JENDL data filename (e.g., 'n_092-U-238.dat' or 'dec-092-U-238m1.dat')
            
        Returns:
            Tuple of (Z, A, metastable_state) or None if parsing fails
        """
        # Pattern for neutron files: n_ZZZ-Symbol-AAA.dat or n_ZZZ-Symbol-AAAm?.dat
        neutron_pattern = r'n_(\d{3})-([A-Z][a-z]?)-(\d{3})(m\d+)?\.dat'
        
        # Pattern for decay files: dec-ZZZ-Symbol-AAA.dat or dec-ZZZ-Symbol-AAAm?.dat  
        decay_pattern = r'dec-(\d{3})-([A-Z][a-z]?)-(\d{3})(m\d+)?\.dat'
        
        for pattern in [neutron_pattern, decay_pattern]:
            match = re.match(pattern, filename)
            if match:
                z = int(match.group(1))
                symbol = match.group(2)
                a = int(match.group(3))
                metastable = match.group(4)
                
                # Parse metastable state
                m = 0
                if metastable:
                    m = int(metastable[1:])  # Extract number after 'm'
                
                return z, a, m
        
        return None
    
    def get_available_nuclides(self) -> Set[Nuclide]:
        """
        Discover all available nuclides from JENDL files
        
        Returns:
            Set of available nuclides
        """
        if self._available_nuclides is not None:
            return self._available_nuclides
        
        nuclides = set()
        
        # Scan neutron data files
        if self.neutron_path.exists():
            for filename in os.listdir(self.neutron_path):
                if filename.endswith('.dat'):
                    parsed = self.parse_filename(filename)
                    if parsed:
                        z, a, m = parsed
                        # Skip neutron (Z=0) - not a real nuclide
                        if z == 0:
                            continue
                        # Create nuclide using centralized system (symbol derived automatically)
                        nuclides.add(Nuclide.from_identity(z, a, m, is_stable=None))
        
        # Scan decay data files  
        if self.decay_path.exists():
            for filename in os.listdir(self.decay_path):
                if filename.endswith('.dat'):
                    parsed = self.parse_filename(filename)
                    if parsed:
                        z, a, m = parsed
                        # Skip neutron (Z=0) - not a real nuclide
                        if z == 0:
                            continue
                        # Create nuclide using centralized system (symbol derived automatically)
                        nuclides.add(Nuclide.from_identity(z, a, m, is_stable=None))
        
        self._available_nuclides = nuclides
        print(f"Found {len(nuclides)} available nuclides in JENDL data")
        return nuclides
    
    def load_neutron_data(self, nuclide: Nuclide) -> Optional[openmc.data.IncidentNeutron]:
        """
        Load neutron cross-section data for a nuclide
        
        Args:
            nuclide: Target nuclide
            
        Returns:
            OpenMC IncidentNeutron data or None if not available
        """
        # Check cache first
        cache_key = str(nuclide)
        if cache_key in self._neutron_data_cache:
            return self._neutron_data_cache[cache_key]
        
        # Construct filename
        metastable_suffix = f"m{nuclide.m}" if nuclide.m > 0 else ""
        filename = f"n_{nuclide.z:03d}-{nuclide.symbol}-{nuclide.a:03d}{metastable_suffix}.dat"
        filepath = self.neutron_path / filename
        
        if not filepath.exists():
            print(f"Warning: Neutron data file not found: {filename}")
            self._neutron_data_cache[cache_key] = None
            return None
        
        try:
            # Load using OpenMC ENDF reader
            print(f"Loading neutron data for {nuclide} from {filename}")
            neutron_data = openmc.data.IncidentNeutron.from_endf(str(filepath))
            self._neutron_data_cache[cache_key] = neutron_data
            return neutron_data
        except Exception as e:
            print(f"Error loading neutron data for {nuclide}: {e}")
            # Log more details for debugging
            import traceback
            print(f"  Full error: {traceback.format_exc()}")
            self._neutron_data_cache[cache_key] = None
            return None
    
    def load_decay_data(self, nuclide: Nuclide) -> Optional[openmc.data.Decay]:
        """
        Load decay data for a nuclide
        
        Args:
            nuclide: Target nuclide
            
        Returns:
            OpenMC Decay data or None if not available
        """
        # Check cache first
        cache_key = str(nuclide)
        if cache_key in self._decay_data_cache:
            return self._decay_data_cache[cache_key]
        
        # Construct filename
        metastable_suffix = f"m{nuclide.m}" if nuclide.m > 0 else ""
        filename = f"dec-{nuclide.z:03d}-{nuclide.symbol}-{nuclide.a:03d}{metastable_suffix}.dat"
        filepath = self.decay_path / filename
        
        if not filepath.exists():
            # Not all nuclides have decay data (stable nuclides)
            self._decay_data_cache[cache_key] = None
            return None
        
        try:
            # Load using OpenMC ENDF reader
            print(f"Loading decay data for {nuclide} from {filename}")
            decay_data = openmc.data.Decay.from_endf(str(filepath))
            self._decay_data_cache[cache_key] = decay_data
            return decay_data
        except Exception as e:
            print(f"Error loading decay data for {nuclide}: {e}")
            # Log more details for debugging
            import traceback
            print(f"  Full error: {traceback.format_exc()}")
            self._decay_data_cache[cache_key] = None
            return None
    
    
    def get_neutron_reactions(self, nuclide: Nuclide) -> List[JENDLReaction]:
        """
        Extract all neutron reactions for a nuclide
        
        Args:
            nuclide: Target nuclide
            
        Returns:
            List of JENDL reactions
        """
        reactions = []
        neutron_data = self.load_neutron_data(nuclide)
        
        if neutron_data is None:
            return reactions
        
        # Extract reactions from OpenMC data
        print(f"Extracting neutron reactions for {nuclide}")
        
        # Extract reactions using OpenMC data (reactions is a dict keyed by MT number)
        for mt, reaction in neutron_data.reactions.items():
            if mt in MT_TO_REACTION_TYPE:
                try:
                    # Extract raw tabulated cross section data
                    tabulated_xs = None
                    has_tabulated_data = False
                    
                    self.xs_extraction_stats['total_reactions'] += 1
                    
                    if hasattr(reaction, 'xs') and reaction.xs:
                        try:
                            tabulated_xs = {}
                            for temp_key, xs_func in reaction.xs.items():
                                # Check if this is tabulated data (has .x and .y attributes)
                                if hasattr(xs_func, 'x') and hasattr(xs_func, 'y'):
                                    # Extract raw tabulated data points
                                    energy_points = xs_func.x.tolist() if hasattr(xs_func.x, 'tolist') else list(xs_func.x)
                                    xs_values = xs_func.y.tolist() if hasattr(xs_func.y, 'tolist') else list(xs_func.y)
                                    tabulated_xs[temp_key] = (energy_points, xs_values)
                                    has_tabulated_data = True
                                else:
                                    # Non-tabulated data - track for reporting
                                    self.non_tabulated_reactions.append({
                                        'nuclide': str(nuclide),
                                        'mt_number': mt,
                                        'reaction_type': MT_TO_REACTION_TYPE[mt].value if mt in MT_TO_REACTION_TYPE else f"MT-{mt}",
                                        'temp_key': temp_key,
                                        'xs_func_type': type(xs_func).__name__,
                                        'reason': f'Non-tabulated cross section function: {type(xs_func).__name__}'
                                    })
                            
                            if has_tabulated_data:
                                self.xs_extraction_stats['tabulated_reactions'] += 1
                            else:
                                self.xs_extraction_stats['non_tabulated_reactions'] += 1
                                
                        except Exception as e:
                            self.xs_extraction_stats['extraction_errors'] += 1
                            self.non_tabulated_reactions.append({
                                'nuclide': str(nuclide),
                                'mt_number': mt,
                                'reaction_type': MT_TO_REACTION_TYPE[mt].value if mt in MT_TO_REACTION_TYPE else f"MT-{mt}",
                                'temp_key': 'unknown',
                                'xs_func_type': 'unknown',
                                'reason': f'Extraction error: {str(e)}'
                            })
                    else:
                        # No cross section data available
                        self.xs_extraction_stats['non_tabulated_reactions'] += 1
                        self.non_tabulated_reactions.append({
                            'nuclide': str(nuclide),
                            'mt_number': mt,
                            'reaction_type': MT_TO_REACTION_TYPE[mt].value if mt in MT_TO_REACTION_TYPE else f"MT-{mt}",
                            'temp_key': 'N/A',
                            'xs_func_type': 'N/A',
                            'reason': 'No cross section data available'
                        })
                    
                    # Extract Q-value (in eV - JENDL native units) and threshold
                    q_value = getattr(reaction, 'q_value', 0.0)  # eV
                    threshold = 0.0
                    if hasattr(reaction, 'threshold'):
                        threshold = reaction.threshold
                    
                    # Determine products (only reference existing nuclides from JENDL files)
                    available_nuclides = self.get_available_nuclides()
                    available_nuclide_objects = {(n.z, n.a, n.m): n for n in available_nuclides}
                    products = determine_products_from_mt(mt, nuclide, available_nuclide_objects)
                    if products:
                        reaction_type = MT_TO_REACTION_TYPE[mt]
                        
                        # Skip reactions that need extended enum for now
                        if isinstance(reaction_type, str):
                            print(f"  Skipping {reaction_type} (MT {mt}) - needs extended enum")
                            continue
                        
                        jendl_reaction = JENDLReaction(
                            reactant=nuclide,
                            product=products[0],
                            reaction_type=reaction_type,
                            mt_number=mt,
                            cross_section=0.0,  # Deprecated - use tabulated_xs
                            q_value=q_value,  # eV (JENDL native units)
                            threshold_energy=threshold,
                            daughter_products=products,
                            tabulated_xs=tabulated_xs,
                            has_tabulated_data=has_tabulated_data
                        )
                        reactions.append(jendl_reaction)
                        
                        # Basic logging
                        q_mev = q_value * 1e-6  # Convert to MeV for display
                        xs_info = f"tabulated_xs: {len(tabulated_xs)} temps" if has_tabulated_data else "no tabulated data"
                        print(f"  Added {reaction_type.value} (MT {mt}): Q={q_mev:.2f}MeV, {xs_info}")
                        
                except Exception as e:
                    print(f"  Error processing MT {mt}: {e}")
                    continue
        
        return reactions
    
    def get_decay_reactions(self, nuclide: Nuclide) -> List[JENDLReaction]:
        """
        Extract decay reactions for a nuclide
        
        Args:
            nuclide: Target nuclide
            
        Returns:
            List of decay reactions
        """
        reactions = []
        decay_data = self.load_decay_data(nuclide)
        
        if decay_data is None:
            return reactions
        
        # Extract decay modes from OpenMC data
        print(f"Extracting decay reactions for {nuclide}")
        
        # Create available nuclide lookup ONCE (not inside the loop)
        available_nuclides = self.get_available_nuclides()
        available_nuclide_objects = {(n.z, n.a, n.m): n for n in available_nuclides}
        
        # Extract half-life
        half_life = None
        if hasattr(decay_data, 'half_life'):
            if hasattr(decay_data.half_life, 'nominal_value'):
                half_life = decay_data.half_life.nominal_value
            else:
                half_life = decay_data.half_life
        
        # Extract decay modes
        if hasattr(decay_data, 'modes'):
            for mode in decay_data.modes:
                try:
                    # Get decay types from mode.modes list
                    if not hasattr(mode, 'modes') or not mode.modes:
                        continue
                    
                    # Get branching ratio
                    branching_ratio = 1.0
                    if hasattr(mode, 'branching_ratio'):
                        if hasattr(mode.branching_ratio, 'nominal_value'):
                            branching_ratio = mode.branching_ratio.nominal_value
                        else:
                            branching_ratio = mode.branching_ratio
                    
                    # Get Q-value from energy
                    q_value = 0.0
                    if hasattr(mode, 'energy'):
                        if hasattr(mode.energy, 'nominal_value'):
                            q_value = mode.energy.nominal_value
                        else:
                            q_value = mode.energy
                    
                    # Process each decay type in the mode
                    for decay_type in mode.modes:
                        try:
                            # Map decay type to reaction type
                            if decay_type in DECAY_MODE_MAPPING:
                                reaction_type = DECAY_MODE_MAPPING[decay_type]
                                
                                # Determine decay product
                                if decay_type in ['beta-', 'beta-/alpha', 'beta-/n']:
                                    product_z = nuclide.z + 1
                                    product_a = nuclide.a
                                elif decay_type in ['beta+', 'ec', 'ec/beta+', 'beta+/ec', 'beta+/p']:
                                    product_z = nuclide.z - 1
                                    product_a = nuclide.a
                                elif decay_type in ['alpha', 'alpha/beta-']:
                                    product_z = nuclide.z - 2
                                    product_a = nuclide.a - 4
                                elif decay_type in ['it', 'IT']:
                                    # For isomeric transition, product is ground state of same nuclide
                                    product_z = nuclide.z
                                    product_a = nuclide.a
                                    # Note: Should transition from metastable to ground state
                                elif decay_type in ['n', 'n/beta-']:
                                    # Neutron emission
                                    product_z = nuclide.z
                                    product_a = nuclide.a - 1
                                elif decay_type in ['p', 'p/beta+']:
                                    # Proton emission
                                    product_z = nuclide.z - 1
                                    product_a = nuclide.a - 1
                                elif decay_type == 'sf':
                                    # Spontaneous fission - complex, skip for now
                                    print(f"  Skipping spontaneous fission - complex decay products")
                                    continue
                                else:
                                    print(f"  Unsupported decay type: {decay_type}")
                                    continue
                                
                                # Create product nuclide using centralized system
                                # For isomeric transition, handle metastable state correctly
                                if decay_type in ['it', 'IT'] and nuclide.m > 0:
                                    # Transition from excited state to ground state or lower excited state
                                    product_m = max(0, nuclide.m - 1)
                                else:
                                    product_m = 0
                                # Only reference existing nuclides from JENDL files - NO CREATION
                                product_identity = (product_z, product_a, product_m)
                                
                                if product_identity not in available_nuclide_objects:
                                    # Skip decay reaction if product doesn't exist in JENDL files
                                    print(f"    Skipping {decay_type}: product {product_identity} not in available nuclides")
                                    continue
                                    
                                product = available_nuclide_objects[product_identity]
                                
                                decay_reaction = JENDLReaction(
                                    reactant=nuclide,
                                    product=product,
                                    reaction_type=reaction_type,
                                    half_life=half_life,
                                    branching_ratio=branching_ratio,
                                    q_value=q_value
                                )
                                reactions.append(decay_reaction)
                                print(f"  Added {decay_type} decay: t½={half_life:.2e}s, BR={branching_ratio:.3f}, Q={q_value:.2f}eV")
                            else:
                                print(f"  Skipping unknown decay type: {decay_type}")
                        
                        except Exception as e:
                            print(f"  Error processing decay type {decay_type}: {e}")
                            continue
                        
                except Exception as e:
                    print(f"  Error processing decay mode: {e}")
                    continue
        
        return reactions
    
    def get_all_reactions(self, nuclide: Nuclide) -> List[JENDLReaction]:
        """
        Get all reactions (neutron + decay) for a nuclide
        
        Args:
            nuclide: Target nuclide
            
        Returns:
            List of all available reactions
        """
        reactions = []
        reactions.extend(self.get_neutron_reactions(nuclide))
        reactions.extend(self.get_decay_reactions(nuclide))
        return reactions
    
    def determine_nuclide_stability(self, nuclide: Nuclide) -> bool:
        """
        Determine if a nuclide is stable based on available decay data
        
        Args:
            nuclide: Target nuclide
            
        Returns:
            True if stable (no decay data), False if radioactive
        """
        decay_data = self.load_decay_data(nuclide)
        return decay_data is None
    
    def get_nuclide_with_stability(self, z: int, a: int, m: int = 0, symbol: str = "") -> Nuclide:
        """
        Create a nuclide with determined stability
        
        Args:
            z: Atomic number
            a: Mass number
            m: Metastable state
            symbol: Element symbol (ignored - symbol derived from registry)
            
        Returns:
            Nuclide with stability determined
        """
        # Create temporary nuclide to check stability (symbol parameter ignored)
        temp_nuclide = Nuclide.from_identity(z, a, m, is_stable=None)
        is_stable = self.determine_nuclide_stability(temp_nuclide)
        
        return Nuclide.from_identity(z, a, m, is_stable=is_stable)
    
    def get_xs_extraction_stats(self) -> Dict:
        """
        Get statistics about cross section data extraction
        
        Returns:
            Dictionary with extraction statistics
        """
        return self.xs_extraction_stats.copy()
    
    def get_non_tabulated_reactions(self) -> List[Dict]:
        """
        Get list of reactions without tabulated cross section data
        
        Returns:
            List of dictionaries with reaction details
        """
        return self.non_tabulated_reactions.copy()
    
    def print_xs_extraction_summary(self):
        """Print summary of cross section extraction results"""
        stats = self.xs_extraction_stats
        total = stats['total_reactions']
        tabulated = stats['tabulated_reactions']
        non_tabulated = stats['non_tabulated_reactions']
        errors = stats['extraction_errors']
        
        print(f"\n=== Cross Section Data Extraction Summary ===")
        print(f"Total reactions processed: {total}")
        print(f"Reactions with tabulated data: {tabulated} ({100*tabulated/total:.1f}%)")
        print(f"Reactions without tabulated data: {non_tabulated} ({100*non_tabulated/total:.1f}%)")
        print(f"Extraction errors: {errors} ({100*errors/total:.1f}%)")
        print(f"=============================================")


# Test functions
def test_jendl_data_loader():
    """Test JENDL data loader with real files"""
    print("Testing JENDL Data Loader...")
    
    loader = JENDLDataLoader()
    
    # Test filename parsing
    test_files = [
        "n_092-U-238.dat",
        "n_094-Pu-239.dat", 
        "dec-092-U-238.dat",
        "dec-043-Tc-099m1.dat"
    ]
    
    for filename in test_files:
        parsed = loader.parse_filename(filename)
        print(f"  {filename} -> {parsed}")
    
    # Test available nuclides discovery
    nuclides = loader.get_available_nuclides()
    print(f"Found {len(nuclides)} available nuclides")
    
    # Test loading specific nuclides
    test_nuclides = [
        Nuclide.from_identity(92, 238, 0),
    ]
    
    for nuclide in test_nuclides:
        print(f"\nTesting {nuclide}:")
        reactions = loader.get_all_reactions(nuclide)
        print(f"  Found {len(reactions)} reactions")
        for reaction in reactions:
            print(f"    {reaction}")
    
    print("JENDL Data Loader test completed!")


def test_mt_product_determination():
    """Test MT number to product determination"""
    print("\nTesting MT product determination...")
    
    u238 = Nuclide.from_identity(92, 238, 0)
    
    test_mts = [102, 103, 107, 16, 17, 28, 22, 32]
    
    for mt in test_mts:
        products = determine_products_from_mt(mt, u238)
        mt_name = MT_TO_REACTION_TYPE.get(mt, f"MT{mt}")
        print(f"  {u238} + n -> {mt_name} -> {products}")
    
    print("MT product determination test completed!")


if __name__ == "__main__":
    test_jendl_data_loader()
    test_mt_product_determination()