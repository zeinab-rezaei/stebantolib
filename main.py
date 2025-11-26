import argparse
from pathlib import Path

from parsers.msp_parser import MSPParser
from massbank_generator import process_records


def main():
    """Main entry point for the MSP to MassBank converter."""
    args = parse_arguments()

    # Setup output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse and process records → writes individual files into outdir
    records = MSPParser.parse_file(args.input)
    process_records(records, outdir, args.accession_prefix, debug_classify=args.debug_classify)

    # --- NEW PART: create output.txt with all records ---
    combined_path = Path("output.txt")   # or choose another path/name
    with combined_path.open("w", encoding="utf-8") as fout:
        for txt_file in sorted(outdir.glob("*.txt")):
            text = txt_file.read_text(encoding="utf-8")
            fout.write(text)
            # ensure a newline between records
            if not text.endswith("\n"):
                fout.write("\n")

    print(f"Wrote {len(records)} MassBank files to {outdir.resolve()}")
    print(f"Wrote combined MassBank file to {combined_path.resolve()}")

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Convert MSP to individual MassBank files")

    parser.add_argument(
        "--input",
        default="input.msp",
        help="Input MSP file"
    )
    parser.add_argument(
        "--outdir",
        default="massbank_records",
        help="Directory to write *.txt files"
    )
    parser.add_argument(
        "--accession-prefix",
        default="HZI-CBIO-AA-",
        help="ACCESSION prefix"
    )
    parser.add_argument(
        "--debug-classify",
        action="store_true",
        help="Log ClassyFire/NPClassifier lookups and save raw JSON"
    )

    return parser.parse_args()



if __name__ == "__main__":
    main()