#!/usr/bin/env python3
"""
Simple Pickle Loader - No OpenMC Dependencies

Loads JENDL nuclear database from pickle files without requiring OpenMC.
Perfect for local analysis after generating the pkl file on HPC.
"""

import pickle
import time
from typing import Dict, List, Set, Optional
from collections import defaultdict
from pathlib import Path
import os

# No OpenMC dependencies needed for pickle loading

from nuclear_chain_analyzer import Nuclide, Reaction, ReactionType
from element_registry import get_symbol, get_atomic_number, normalize_symbol


class SimplePickleNuclearDatabase:
    """
    Nuclear database that loads from pickle files without OpenMC dependencies
    Completely independent implementation - no inheritance from OpenMC-dependent classes
    """
    
    def __init__(self, cache_dir: str = r"C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp\My crap", pkl_filename: str = "jendl_nuclear_data_appended.pkl"):
        """
        Initialize database from pickle file
        
        Args:
            cache_dir: Directory containing pickle file
            pkl_filename: Name of pickle file
        """
        print("Loading JENDL Nuclear Database from pickle file...")
        start_time = time.time()
        
        # Initialize parent class data structures
        self.reactions = []
        self.nuclides = set()
        self.reactions_from = defaultdict(list)
        self.reactions_to = defaultdict(list)
        
        # Load from pickle file
        cache_path = Path(cache_dir).resolve()  # Get absolute path
        pickle_file = cache_path / pkl_filename
        
        if not pickle_file.exists():
            available_pkl_files = list(cache_path.glob('*.pkl'))
            error_msg = f"Pickle file not found: {pickle_file}\n"
            error_msg += f"Current directory: {Path.cwd()}\n"
            if available_pkl_files:
                error_msg += f"Available pkl files: {[f.name for f in available_pkl_files]}"
            else:
                error_msg += "No pkl files found in directory"
            raise FileNotFoundError(error_msg)
        
        print(f"Loading from: {pickle_file}")
        
        try:
            with open(pickle_file, 'rb') as f:
                cache_data = pickle.load(f)
            
            self.nuclides = set(cache_data['nuclides'])
            self.reactions = cache_data['reactions']
            
            print(f"Loaded {len(self.nuclides)} nuclides and {len(self.reactions)} reactions")
            
            # OPTION 1: Validate stability information (should all be determined by loader)
            unknown_stability = [n for n in self.nuclides if n.is_stable is None]
            if unknown_stability:
                print(f"WARNING: {len(unknown_stability)} nuclides have unknown stability!")
                print("This suggests an older database format - applying stability determination...")
                self._determine_stability()
            else:
                stable_count = sum(1 for n in self.nuclides if n.is_stable is True)
                radioactive_count = sum(1 for n in self.nuclides if n.is_stable is False)
                print(f"✓ Stability information complete: {stable_count} stable, {radioactive_count} radioactive")
                
                # OPTION 1: Check for reaction reference consistency
                sample_reactions = self.reactions[:10] if len(self.reactions) > 10 else self.reactions
                reactions_need_update = any(
                    r.reactant.is_stable is None or r.product.is_stable is None 
                    for r in sample_reactions
                )
                if reactions_need_update:
                    print("Reactions have stability inconsistencies - updating references...")
                    self._update_reactions_with_stability()
                else:
                    print("✓ Reaction references are consistent with stability information")
            
        except Exception as e:
            raise RuntimeError(f"Failed to load pickle file: {e}")
        
        # Build internal data structures
        self._build_database_structures()
        
        load_time = time.time() - start_time
        print(f"Database loaded in {load_time:.2f} seconds")
    
    def add_reaction(self, reaction: Reaction):
        """Add a reaction to the database"""
        self.reactions.append(reaction)
        self.reactions_from[reaction.reactant].append(reaction)
        self.reactions_to[reaction.product].append(reaction)
        self.nuclides.add(reaction.reactant)
        self.nuclides.add(reaction.product)
    
    def _determine_stability(self):
        """
        Determine stability for each nuclide based on available decay reactions
        A nuclide is stable if it has no decay reactions
        """
        print("Determining nuclide stability...")
        
        # Track which nuclides have decay reactions
        decaying_nuclides = set()
        for reaction in self.reactions:
            if reaction.is_decay():
                decaying_nuclides.add(reaction.reactant)
        
        # Update stability for all nuclides
        updated_nuclides = set()
        stable_count = 0
        radioactive_count = 0
        
        for nuclide in self.nuclides:
            # If stability is already known, keep it
            if nuclide.is_stable is not None:
                updated_nuclides.add(nuclide)
                if nuclide.is_stable:
                    stable_count += 1
                else:
                    radioactive_count += 1
                continue
            
            # Determine stability based on decay reactions
            is_stable = nuclide not in decaying_nuclides
            
            # Create new nuclide object with stability information
            updated_nuclide = Nuclide(
                z=nuclide.z,
                a=nuclide.a, 
                m=nuclide.m,
                symbol=nuclide.symbol,
                is_stable=is_stable
            )
            updated_nuclides.add(updated_nuclide)
            
            if is_stable:
                stable_count += 1
            else:
                radioactive_count += 1
        
        # Replace nuclides set with updated one
        self.nuclides = updated_nuclides
        
        print(f"Stability determined: {stable_count} stable, {radioactive_count} radioactive nuclides")
        
        # Update reactions to use the new nuclides with stability information
        self._update_reactions_with_stability()
    
    def _update_reactions_with_stability(self):
        """
        Update all reactions to use nuclides with stability information
        This ensures lookup tables work correctly after stability determination
        """
        print("Updating reactions with stability information...")
        
        # Create mapping from (z,a,m) identity to updated nuclide (symbol no longer part of key)
        nuclide_map = {}
        for nuclide in self.nuclides:
            key = nuclide.get_identity_key()  # Returns (z, a, m)
            nuclide_map[key] = nuclide
        
        # Update each reaction's reactant and product
        updated_reactions = []
        for reaction in self.reactions:
            # Find updated reactant using identity key
            reactant_key = (reaction.reactant.z, reaction.reactant.a, reaction.reactant.m)
            updated_reactant = nuclide_map.get(reactant_key, reaction.reactant)
            
            # Find updated product using identity key
            product_key = (reaction.product.z, reaction.product.a, reaction.product.m)
            updated_product = nuclide_map.get(product_key, reaction.product)
            
            # Create new reaction with updated nuclides
            # Preserve all other reaction attributes
            updated_reaction = type(reaction)(
                reactant=updated_reactant,
                product=updated_product,
                reaction_type=reaction.reaction_type,
                cross_section=reaction.cross_section,
                half_life=reaction.half_life,
                branching_ratio=reaction.branching_ratio,
                threshold_energy=reaction.threshold_energy,
                q_value=getattr(reaction, 'q_value', 0.0),
                **{k: v for k, v in reaction.__dict__.items() 
                   if k not in ['reactant', 'product', 'reaction_type', 'cross_section', 
                               'half_life', 'branching_ratio', 'threshold_energy', 'q_value']}
            )
            updated_reactions.append(updated_reaction)
        
        # Replace reactions list
        self.reactions = updated_reactions
        print(f"Updated {len(updated_reactions)} reactions with stability information")
    
    def _build_database_structures(self):
        """Build internal database structures from loaded data"""
        print("Building reaction lookup tables...")
        
        # Build reaction lookup tables
        for reaction in self.reactions:
            self.reactions_from[reaction.reactant].append(reaction)
            self.reactions_to[reaction.product].append(reaction)
    
    def find_nuclide(self, symbol: str, mass_number: int, metastable: int = 0) -> Optional[Nuclide]:
        """
        Find a nuclide by symbol and mass number (now uses atomic number lookup)
        
        Args:
            symbol: Chemical symbol (e.g., 'U', 'Pu')
            mass_number: Mass number
            metastable: Metastable state (default 0)
            
        Returns:
            Nuclide object or None if not found
        """
        # Convert symbol to atomic number for reliable lookup
        from element_registry import get_atomic_number, normalize_symbol
        atomic_number = get_atomic_number(normalize_symbol(symbol))
        if atomic_number is None:
            return None
        
        # Search by (z, a, m) identity - more reliable than symbol matching
        for nuclide in self.nuclides:
            if (nuclide.z == atomic_number and 
                nuclide.a == mass_number and 
                nuclide.m == metastable):
                return nuclide
        
        return None
    
    def get_reactions_from(self, nuclide: Nuclide) -> List[Reaction]:
        """Get all reactions starting from a nuclide"""
        # First try direct lookup (fastest)
        direct_result = self.reactions_from.get(nuclide, [])
        if direct_result:
            return direct_result
        
        # If direct lookup fails, search by Z, A, M (handles stability mismatch)
        for db_nuclide, reactions in self.reactions_from.items():
            if (db_nuclide.z == nuclide.z and 
                db_nuclide.a == nuclide.a and 
                db_nuclide.m == nuclide.m):
                return reactions
        
        return []
    
    def get_reactions_to(self, nuclide: Nuclide) -> List[Reaction]:
        """Get all reactions producing a nuclide"""
        # First try direct lookup (fastest)
        direct_result = self.reactions_to.get(nuclide, [])
        if direct_result:
            return direct_result
        
        # If direct lookup fails, search by Z, A, M (handles stability mismatch)
        for db_nuclide, reactions in self.reactions_to.items():
            if (db_nuclide.z == nuclide.z and 
                db_nuclide.a == nuclide.a and 
                db_nuclide.m == nuclide.m):
                return reactions
        
        return []
    
    def get_all_nuclides(self) -> Set[Nuclide]:
        """Get all available nuclides"""
        return self.nuclides.copy()
    
    def get_database_stats(self) -> Dict:
        """Get database statistics"""
        stats = {
            'total_nuclides': len(self.nuclides),
            'total_reactions': len(self.reactions),
            'neutron_reactions': 0,
            'decay_reactions': 0,
            'reaction_types': defaultdict(int),
            'elements': set(),
            'mass_range': [float('inf'), 0]
        }
        
        # Analyze reactions
        for reaction in self.reactions:
            if reaction.is_neutron_induced():
                stats['neutron_reactions'] += 1
            elif reaction.is_decay():
                stats['decay_reactions'] += 1
            
            stats['reaction_types'][reaction.reaction_type] += 1
        
        # Analyze nuclides
        for nuclide in self.nuclides:
            stats['elements'].add(nuclide.symbol)
            stats['mass_range'][0] = min(stats['mass_range'][0], nuclide.a)
            stats['mass_range'][1] = max(stats['mass_range'][1], nuclide.a)
        
        stats['total_elements'] = len(stats['elements'])
        return stats
    
    def get_all_reactions(self) -> List[Reaction]:
        """Get all reactions in the database"""
        return self.reactions.copy()


if __name__ == "__main__":
    os.chdir(r"C:\Users\lowin\Documents\Marathon\nuclear-parser-temp\nuclear-parser-temp")
    # Test the simple loader
    try:
        db = SimplePickleNuclearDatabase()
        stats = db.get_database_stats()
        print(f"\nDatabase Statistics:")
        print(f"  Total nuclides: {stats['total_nuclides']}")
        print(f"  Total reactions: {stats['total_reactions']}")
        print(f"  Neutron reactions: {stats['neutron_reactions']}")
        print(f"  Decay reactions: {stats['decay_reactions']}")
        print(f"  Elements: {stats['total_elements']}")
        
        # Test finding Au-197
        au197 = db.find_nuclide('Au', 197)
        if au197:
            print(f"\nFound {au197}")
            reactions_to = len(db.get_reactions_to(au197))
            print(f"  {reactions_to} reactions can produce Au-197")
        else:
            print("\nAu-197 not found in database")
            
    except Exception as e:
        print(f"Error: {e}")