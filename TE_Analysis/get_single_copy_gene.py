import os
import glob
from typing import Tuple, Optional
from Bio.Seq import Seq
import pysam
import argparse


def parse_first_single_copy_busco(full_table_path: str) -> Optional[str]:
    with open(full_table_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            busco_id, status = fields[0], fields[1]
            if status.lower() == "complete":
                return busco_id
    return None


def find_gff_file(busco_run_path: str) -> str:
    gff_files = glob.glob(os.path.join(busco_run_path, "**", "*.gff"), recursive=True)
    if not gff_files:
        raise FileNotFoundError("No GFF file found in BUSCO run directory.")
    return gff_files[0]


def get_coordinates_from_gff(gff_file: str, busco_id: str) -> Tuple[str, int, int, str]:
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or busco_id not in line:
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            return chrom, start, end, strand
    raise ValueError(f"BUSCO ID {busco_id} not found in GFF file.")


def extract_sequence(reference_fasta: str, chrom: str, start: int, end: int, strand: str) -> str:
    fasta = pysam.FastaFile(reference_fasta)
    try:
        sequence = fasta.fetch(chrom, start - 1, end)
    except Exception as e:
        raise RuntimeError(f"Could not extract {chrom}:{start}-{end} from FASTA: {e}")
    return str(Seq(sequence).reverse_complement()) if strand == "-" else sequence


def write_fasta(out_file: str, header: str, sequence: str):
    with open(out_file, "w") as f:
        f.write(f">{header}\n{sequence}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract the first complete BUSCO gene (nucleotide sequence) from a reference genome."
    )
    parser.add_argument("busco_dir", help="Top-level BUSCO output directory (e.g. BUSCO_Analysis/)")
    parser.add_argument("reference_fasta", help="Path to reference genome FASTA file")
    parser.add_argument(
        "--busco_set",
        default="run_drosophilidae_odb12",
        help="Name of the BUSCO run folder inside the BUSCO directory (default: run_drosophilidae_odb12)"
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Prefix for the output FASTA file (default: empty, uses folder name)"
    )

    args = parser.parse_args()

    busco_run_path = os.path.join(args.busco_dir, args.busco_set)
    full_table_path = os.path.join(busco_run_path, "full_table.tsv")

    if not os.path.isfile(full_table_path):
        parser.error(f"full_table.tsv not found in {busco_run_path}")

    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference FASTA not found: {args.reference_fasta}")

    # Step 1: BUSCO ID
    busco_id = parse_first_single_copy_busco(full_table_path)
    if not busco_id:
        parser.error("No complete BUSCO gene found in full_table.tsv")
    print(f"Found BUSCO ID: {busco_id}")

    # Step 2: Locate GFF
    gff_file = find_gff_file(busco_run_path)
    print(f"Using GFF file: {gff_file}")

    # Step 3: Coordinates
    chrom, start, end, strand = get_coordinates_from_gff(gff_file, busco_id)
    print(f"Coordinates: {chrom}:{start}-{end} ({strand})")

    # Step 4: Extract sequence
    seq = extract_sequence(args.reference_fasta, chrom, start, end, strand)

    # get folder name
    folder_name = os.path.basename(os.path.dirname(args.busco_dir))

    # If prefix is provided, use it; otherwise, use the folder name
    if args.prefix:
        folder_name = args.prefix   

    # Step 5: Output
    out_file = f"{folder_name}_{busco_id}_nt.fa"
    write_fasta(out_file, f"{busco_id}_{chrom}_{start}_{end}_{strand}", seq)
    print(f"Sequence written to: {out_file}")


if __name__ == "__main__":
    main()
