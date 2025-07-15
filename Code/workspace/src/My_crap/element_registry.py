#!/usr/bin/env python3
"""
Element Registry - Centralized Periodic Table Management

Provides a single authoritative source for element symbols and atomic numbers.
Eliminates inconsistent symbol mappings across the nuclear database system.
"""

from typing import Dict, Optional, Set
import re


class ElementRegistry:
    """
    Centralized registry for chemical element information
    
    Provides consistent element symbol mapping and validation across
    the entire nuclear database system.
    """
    
    # Complete periodic table mapping (Z -> symbol)
    _ELEMENT_SYMBOLS: Dict[int, str] = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
        11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
        21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
        31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
        41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
        61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
        71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
        81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
        91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
        101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
        111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og',
        # Future superheavy elements (predicted)
        119: 'Uue', 120: 'Ubn', 121: 'Ubu', 122: 'Ubb', 123: 'Ubt', 124: 'Ubq', 125: 'Ubp', 126: 'Ubh'
    }
    
    # Reverse mapping (symbol -> Z) for lookups
    _SYMBOL_TO_Z: Dict[str, int] = {symbol: z for z, symbol in _ELEMENT_SYMBOLS.items()}
    
    @classmethod
    def get_symbol(cls, atomic_number: int) -> str:
        """
        Get element symbol for given atomic number
        
        Args:
            atomic_number: Atomic number (Z)
            
        Returns:
            Element symbol (e.g., 'U', 'Pu') or systematic name for unknown elements
        """
        if atomic_number in cls._ELEMENT_SYMBOLS:
            return cls._ELEMENT_SYMBOLS[atomic_number]
        
        # Generate systematic name for unknown superheavy elements
        # Following IUPAC systematic element naming convention
        return cls._generate_systematic_name(atomic_number)
    
    @classmethod
    def get_atomic_number(cls, symbol: str) -> Optional[int]:
        """
        Get atomic number for given element symbol
        
        Args:
            symbol: Element symbol (e.g., 'U', 'Pu')
            
        Returns:
            Atomic number or None if symbol not found
        """
        # Normalize symbol (proper case)
        normalized_symbol = cls.normalize_symbol(symbol)
        return cls._SYMBOL_TO_Z.get(normalized_symbol)
    
    @classmethod
    def normalize_symbol(cls, symbol: str) -> str:
        """
        Normalize element symbol to standard format
        
        Args:
            symbol: Element symbol in any case
            
        Returns:
            Normalized symbol (first letter uppercase, rest lowercase)
        """
        if not symbol:
            return ""
        
        # Handle systematic names (e.g., "Uue", "Z119")
        if symbol.startswith('Z') and symbol[1:].isdigit():
            # Convert Z notation to systematic name
            z = int(symbol[1:])
            return cls.get_symbol(z)
        
        # Standard normalization: First letter uppercase, rest lowercase
        return symbol[0].upper() + symbol[1:].lower()
    
    @classmethod
    def _generate_systematic_name(cls, atomic_number: int) -> str:
        """
        Generate IUPAC systematic name for superheavy elements
        
        Args:
            atomic_number: Atomic number
            
        Returns:
            Systematic element name (e.g., Uue for 119)
        """
        if atomic_number <= 118:
            return f"Z{atomic_number}"  # Fallback for missing known elements
        
        # IUPAC systematic naming: digits -> latin roots -> symbol
        digit_names = {
            '0': 'nil', '1': 'un', '2': 'bi', '3': 'tri', '4': 'quad',
            '5': 'pent', '6': 'hex', '7': 'sept', '8': 'oct', '9': 'enn'
        }
        
        # Convert atomic number to systematic name
        digits = str(atomic_number)
        roots = [digit_names[digit] for digit in digits]
        systematic_name = ''.join(roots) + 'ium'
        
        # Create symbol from first letters, capitalize first
        symbol_parts = []
        for root in roots:
            symbol_parts.append(root[0])
        
        # Add 'ium' ending indicator if needed for uniqueness
        if len(symbol_parts) == 2:
            symbol_parts.append('m')  # From 'ium'
        elif len(symbol_parts) == 1:
            symbol_parts.extend(['u', 'm'])  # From 'ium'
        
        symbol = ''.join(symbol_parts)
        return symbol[0].upper() + symbol[1:].lower()
    
    @classmethod
    def is_valid_symbol(cls, symbol: str) -> bool:
        """
        Check if symbol is a valid element symbol
        
        Args:
            symbol: Element symbol to validate
            
        Returns:
            True if valid element symbol
        """
        normalized = cls.normalize_symbol(symbol)
        return normalized in cls._SYMBOL_TO_Z or normalized.startswith('Z')
    
    @classmethod
    def get_all_known_elements(cls) -> Set[int]:
        """
        Get set of all known atomic numbers
        
        Returns:
            Set of atomic numbers for all known elements
        """
        return set(cls._ELEMENT_SYMBOLS.keys())
    
    @classmethod
    def validate_element_data(cls, atomic_number: int, symbol: str) -> bool:
        """
        Validate that atomic number and symbol are consistent
        
        Args:
            atomic_number: Atomic number
            symbol: Element symbol
            
        Returns:
            True if atomic number and symbol match
        """
        expected_symbol = cls.get_symbol(atomic_number)
        normalized_symbol = cls.normalize_symbol(symbol)
        
        return expected_symbol == normalized_symbol
    
    @classmethod
    def get_display_name(cls, atomic_number: int, mass_number: int, metastable: int = 0) -> str:
        """
        Get formatted display name for nuclide
        
        Args:
            atomic_number: Atomic number (Z)
            mass_number: Mass number (A)
            metastable: Metastable state (default 0)
            
        Returns:
            Formatted nuclide name (e.g., "U-238", "Tc-99m1")
        """
        symbol = cls.get_symbol(atomic_number)
        metastable_suffix = f"m{metastable}" if metastable > 0 else ""
        return f"{symbol}-{mass_number}{metastable_suffix}"


# Create singleton instance for module-level access
_registry = ElementRegistry()

# Module-level convenience functions
def get_symbol(atomic_number: int) -> str:
    """Get element symbol for atomic number"""
    return _registry.get_symbol(atomic_number)

def get_atomic_number(symbol: str) -> Optional[int]:
    """Get atomic number for element symbol"""
    return _registry.get_atomic_number(symbol)

def normalize_symbol(symbol: str) -> str:
    """Normalize element symbol to standard format"""
    return _registry.normalize_symbol(symbol)

def is_valid_symbol(symbol: str) -> bool:
    """Check if symbol is valid"""
    return _registry.is_valid_symbol(symbol)

def validate_element_data(atomic_number: int, symbol: str) -> bool:
    """Validate atomic number and symbol consistency"""
    return _registry.validate_element_data(atomic_number, symbol)

def get_display_name(atomic_number: int, mass_number: int, metastable: int = 0) -> str:
    """Get formatted display name for nuclide"""
    return _registry.get_display_name(atomic_number, mass_number, metastable)


def test_element_registry():
    """Test element registry functionality"""
    print("Testing Element Registry...")
    
    # Test known elements
    assert get_symbol(1) == 'H'
    assert get_symbol(92) == 'U'
    assert get_symbol(118) == 'Og'
    print("✓ Known element symbols correct")
    
    # Test reverse lookup
    assert get_atomic_number('H') == 1
    assert get_atomic_number('U') == 92
    assert get_atomic_number('u') == 92  # Case insensitive
    print("✓ Reverse lookups work")
    
    # Test normalization
    assert normalize_symbol('uranium') == 'Uranium'  # Not a real symbol, but normalized
    assert normalize_symbol('u') == 'U'
    assert normalize_symbol('Au') == 'Au'
    print("✓ Symbol normalization works")
    
    # Test systematic names for superheavy elements
    symbol_119 = get_symbol(119)
    print(f"✓ Element 119 symbol: {symbol_119}")
    
    # Test validation
    assert validate_element_data(92, 'U') == True
    assert validate_element_data(92, 'Pu') == False
    print("✓ Element data validation works")
    
    # Test display names
    assert get_display_name(92, 238) == 'U-238'
    assert get_display_name(43, 99, 1) == 'Tc-99m1'
    print("✓ Display name formatting works")
    
    print("Element Registry tests passed! ✅")


if __name__ == "__main__":
    test_element_registry()