#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO

def embl_to_fasta(input_path):
    """
    Convert an EMBL file to a FASTA file with the same base name in the same directory.
    """
    input_file = Path(input_path)
    if not input_file.exists() or not input_file.is_file():
        print(f"Error: File not found - {input_path}")
        return

    output_file = input_file.with_suffix(".fasta")

    try:
        count = SeqIO.convert(str(input_file), "embl", str(output_file), "fasta")
        print(f"Converted {count} record(s) to {output_file.name}")
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert an EMBL file to a FASTA file (same name, .fasta extension)"
    )
    parser.add_argument(
        "input", help="Path to the input EMBL file"
    )
    args = parser.parse_args()
    embl_to_fasta(args.input)

if __name__ == "__main__":
    main()
