# msp_parser.py
from __future__ import annotations
from typing import List, Tuple, Optional, Callable, Any
import re
from collections import Counter
from dataclasses import dataclass

MONOISOTOPIC = {
    'H': 1.007825032, 'C': 12.0, 'N': 14.003074, 'O': 15.99491462,
    'P': 30.973762, 'S': 31.972071, 'Cl': 34.968853, 'Br': 78.918338,
    'F': 18.9984032, 'I': 126.90447
}

_INCHI_FORMULA_PATTERN = re.compile(r'InChI=1S?/([^/]+)')


@dataclass
class MSPRecord:
    """Data class for MSP record with type hints."""
    peaks: List[Tuple[float, float]]
    charge: Optional[int] = None
    ion_mode: Optional[str] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    pubmed: Optional[str] = None
    ccs: Optional[float] = None
    col_energy1: Optional[float] = None
    col_energy2: Optional[float] = None
    ms_level: Optional[int] = None
    instrument_type: Optional[str] = None
    compound_name: Optional[str] = None
    precursor_mz: Optional[float] = None
    inchikey: Optional[str] = None
    num_peaks: Optional[int] = None


class MSPParser:
    """Parser for custom MSP-like text files."""

    # Field mappings: field_name -> (parser_function, target_field_name)
    FIELD_PARSERS = {
        'CHARGE': (int, 'charge'),
        'IONMODE': (str.upper, 'ion_mode'),
        'SMILES': (str, 'smiles'),
        'INCHI': (str, 'inchi'),
        'PUBMED': (str, 'pubmed'),
        'CCS': (float, 'ccs'),
        'COL_ENERGY1': (float, 'col_energy1'),
        'COL_ENERGY2': (float, 'col_energy2'),
        'MS_LEVEL': (int, 'ms_level'),
        'INSTRUMENT_TYPE': (str, 'instrument_type'),
        'COMPOUND_NAME': (str, 'compound_name'),
        'PRECURSOR_MZ': (float, 'precursor_mz'),
        'INCHIKEY': (str, 'inchikey'),
        'NUM_PEAKS': (int, 'num_peaks'),
    }

    @classmethod
    def parse_lines(cls, lines: List[str]) -> List[MSPRecord]:
        """Parse MSP lines into structured records."""
        records = []
        current_data = {}
        current_peaks = []

        def flush_record():
            """Flush current record to records list."""
            nonlocal current_data, current_peaks
            if current_data and current_peaks:
                current_data['peaks'] = current_peaks
                try:
                    records.append(MSPRecord(**current_data))
                except TypeError as e:
                    print(f"Error creating MSPRecord with data: {current_data}")
                    print(f"Error: {e}")
            current_data = {}
            current_peaks = []

        for line in lines:
            stripped_line = line.strip()
            if not stripped_line:
                continue

            # Handle field lines
            if ':' in stripped_line:
                field_name, value = map(str.strip, stripped_line.split(':', 1))
                field_name_clean = field_name.replace(' ', '_').upper()

                if field_name_clean == 'CHARGE':
                    flush_record()

                if field_name_clean in cls.FIELD_PARSERS:
                    parser, target_field = cls.FIELD_PARSERS[field_name_clean]
                    try:
                        current_data[target_field] = parser(value)
                    except (ValueError, TypeError) as e:
                        print(f"Warning: Could not parse {field_name}: {value} - {e}")
            else:
                # Handle peak data
                try:
                    mz, intensity = map(float, stripped_line.split())
                    current_peaks.append((mz, intensity))
                except ValueError:
                    # Skip lines that can't be parsed as peaks
                    continue

        # Don't forget the last record
        flush_record()

        return records

    @classmethod
    def parse_file(cls, filepath: str) -> List[MSPRecord]:
        """Parse MSP file directly from file path."""
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        return cls.parse_lines(lines)


def extract_formula_and_mass_from_inchi(inchi: str) -> tuple[Optional[str], Optional[str]]:
    """Extract formula from InChI and compute exact mass.

    Args:
        inchi: InChI string

    Returns:
        Tuple of (formula, exact_mass) or (None, None) if parsing fails
    """
    if not inchi:
        return None, None

    match = _INCHI_FORMULA_PATTERN.match(inchi)
    if not match:
        return None, None

    formula = match.group(1)

    # Parse formula components
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    element_counts = Counter({
        element: int(count) if count else 1
        for element, count in tokens
    })

    # Calculate exact mass
    try:
        exact_mass = sum(
            MONOISOTOPIC[element] * count
            for element, count in element_counts.items()
        )
        return formula, f"{exact_mass:.4f}"
    except KeyError as e:
        # Handle unknown elements
        print(f"Warning: Unknown element in formula: {e}")
        return formula, None