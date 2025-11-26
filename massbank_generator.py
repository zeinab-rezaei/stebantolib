from pathlib import Path
from typing import List, Tuple, Dict, Any

from classifiers.classyfire import get_compound_classes
from parsers.msp_parser import MSPRecord, extract_formula_and_mass_from_inchi
from utils.splash_utils import splash_for_peaks

from pathlib import Path
from typing import List, Tuple, Dict, Any

from classifiers.classyfire import get_compound_classes
from parsers.msp_parser import MSPRecord, extract_formula_and_mass_from_inchi
from utils.splash_utils import splash_for_peaks


class MassBankRecordBuilder:
    """Builds MassBank records from MSP records."""

    def __init__(self, accession_prefix: str, debug_classify: bool = False):
        self.accession_prefix = accession_prefix
        self.debug_classify = debug_classify

    def build_record(self, record: MSPRecord, index: int) -> Tuple[str, str]:
        """Convert one MSPRecord into a MassBank3-compatible text block."""
        accession = self._generate_accession(index)
        title = self._build_title(record)
        formula, exact_mass = extract_formula_and_mass_from_inchi(record.inchi or "")
        classes = self._get_compound_classes(record)

        mb_lines = [
            *self._build_header(accession, title, record),  # Pass record here
            *self._build_compound_section(record, formula, exact_mass, classes),
            *self._build_instrument_section(record),
            *self._build_spectral_section(record),
            "//"
        ]

        return accession, "\n".join(mb_lines)

    def _generate_accession(self, index: int) -> str:
        return f"MSBNK-{self.accession_prefix}{index:06d}"

    def _build_title(self, record: MSPRecord) -> str:
        title_parts = [
            record.compound_name,
            record.instrument_type,
            f"MS{record.ms_level}",
            f"CE:{record.col_energy1}",
            "[M+H]+"
        ]
        return "; ".join(str(part) for part in title_parts if part is not None)

    def _get_compound_classes(self, record: MSPRecord) -> Dict[str, Any]:
        """Get compound classification with optional debugging."""
        classes = get_compound_classes(
            (record.inchikey or "").strip() or None,
            (record.smiles or "").strip() or None
        )

        # if self.debug_classify:
        #     self._log_classification_debug(record, classes)

        return classes or {}

    # def _log_classification_debug(self, record: MSPRecord, classes: Dict[str, Any]) -> None:
    #     """Log classification debug information."""
    #     source = classes.get("source") if classes else "NONE"
    #     print(f"[classify] used={source} | IK={record.inchikey} | NAME={record.compound_name}")
    #     if classes:
    #         print(f"          classes={classes}")

    def _build_header(self, accession: str, title: str, record: MSPRecord) -> List[str]:
        # Build comment with CCS as dictionary-like format
        if record.ccs:
            comment = f"Ccs: {record.ccs}"
        else:
            comment = "Ccs: not available"

        return [
            f"ACCESSION: {accession}",
            f"RECORD_TITLE: {title}",
            "DATE: 2025.08.05",
            "AUTHORS: Helmholtz Centre for Infection Research",
            "LICENSE: CC BY-NC",
            "COPYRIGHT: HZI",
            f"COMMENT: {comment}",
        ]

    def _build_compound_section(self, record: MSPRecord, formula: str, exact_mass: str, classes: Dict[str, Any]) -> \
    List[str]:
        section = [
            f"CH$NAME: {record.compound_name or 'Unknown'}",
            self._build_compound_class_line(classes),
            f"CH$FORMULA: {formula or 'CH4'}",
            f"CH$EXACT_MASS: {exact_mass or '16'}",
            f"CH$SMILES: {record.smiles or 'C'}",
            f"CH$IUPAC: {record.inchi or 'InChI=1S/CH4/h1H4'}",
            f"CH$LINK: INCHIKEY {record.inchikey or 'VNWKTOKETHGBQD-UHFFFAOYSA-N'}",
        ]

        # Add ChemOnt link if classification is available
        chemont_link = self._build_chemont_link(classes)
        if chemont_link:
            section.append(chemont_link)

        return section

    # def _build_compound_class_line(self, classes: Dict[str, Any]) -> str:
    #     if not classes:
    #         return "CH$COMPOUND_CLASS: unavailable"
    #
    #     return (
    #         f"CH$COMPOUND_CLASS: "
    #         f"ontology={classes.get('ontology', '')}; "
    #         f"kingdom={classes.get('kingdom', '')}; "
    #         f"superclass={classes.get('superclass', '')}; "
    #         f"class={classes.get('class', '')}; "
    #         f"subclass={classes.get('subclass', '')}; "
    #         f"direct_parent={classes.get('direct_parent', '')}; "
    #         f"chemont_id={classes.get('ontology_id', '')}; "
    #         f"source={classes.get('source', '')}"
    #     )
    def _build_compound_class_line(self, classes: Dict[str, Any]) -> str:
        """Build compound class line using direct_parent from classification"""
        if not classes:
            return "CH$COMPOUND_CLASS: unavailable"

        # Use direct_parent from classification if available, otherwise use default
        direct_parent = classes.get('direct_parent', 'Natural Product')
        compound_classes = [direct_parent, "Natural Product"]

        return f"CH$COMPOUND_CLASS: {'; '.join(compound_classes)}"

    def _build_chemont_link(self, classes: Dict[str, Any]) -> str:
        """Build ChemOnt link in the format expected by your React frontend."""
        if not classes or not classes.get('ontology_id'):
            return ""

        ontology_id = classes.get('ontology_id', '')
        kingdom = classes.get('kingdom', '')
        superclass = classes.get('superclass', '')
        compound_class = classes.get('class', '')
        subclass = classes.get('subclass', '')
        direct_parent = classes.get('direct_parent', '')

        # Build the identifier string in the format your React code expects
        identifier_parts = [ontology_id]
        if kingdom:
            identifier_parts.append(kingdom)
        if superclass:
            identifier_parts.append(superclass)
        if compound_class:
            identifier_parts.append(compound_class)
        if subclass:
            identifier_parts.append(subclass)
        if direct_parent:
            identifier_parts.append(direct_parent)

        identifier = "; ".join(identifier_parts)

        return f"CH$LINK: CHEMONT {identifier}"

    def _build_instrument_section(self, record: MSPRecord) -> List[str]:
        return [
            f"AC$INSTRUMENT: { 'Bruker timsTOF'}",
            f"AC$INSTRUMENT_TYPE: { 'timsTOF'}",
            f"AC$MASS_SPECTROMETRY: MS_TYPE MS{record.ms_level or '2'}",
            f"AC$MASS_SPECTROMETRY: ION_MODE {record.ion_mode or 'Positive'}",
            f"AC$MASS_SPECTROMETRY: COLLISION_ENERGY {record.col_energy1 or '0'}",
            f"MS$FOCUSED_ION: PRECURSOR_M/Z {record.precursor_mz or '0'}",
            "MS$DATA_PROCESSING: Converted with convert_msp_to_massbank.py",
        ]

    def _build_spectral_section(self, record: MSPRecord) -> List[str]:
        section = [
            f"PK$SPLASH: {splash_for_peaks(record.peaks)}",
            f"PK$NUM_PEAK: {record.num_peaks or len(record.peaks)}",
            "PK$PEAK: m/z int. rel.int.",
        ]

        # Add peaks with relative intensities
        if record.peaks:
            max_intensity = max(i for _, i in record.peaks)
            for mz, intensity in record.peaks:
                rel_int = int(round(intensity / max_intensity * 999))
                section.append(f"  {mz:.5f} {intensity:.6g} {rel_int}")

        return section


def process_records(records: List[MSPRecord], outdir: Path, accession_prefix: str,
                    debug_classify: bool = False) -> None:
    """Process all records and write MassBank files with progress logging."""
    builder = MassBankRecordBuilder(accession_prefix, debug_classify)

    total = len(records)
    print(f"Found {total} MSP records. Starting Process...\n")

    for idx, record in enumerate(records, start=1):
        name = record.compound_name or "Unknown"
        print(f"[{idx}/{total}] Writing MassBank record for: {name!r}")
        accession, mb_text = builder.build_record(record, idx)
        outfile = outdir / f"{accession}.txt"
        outfile.write_text(mb_text, encoding="utf-8")

    print(f"\nDone. Wrote {total} MassBank records to {outdir.resolve()}")

