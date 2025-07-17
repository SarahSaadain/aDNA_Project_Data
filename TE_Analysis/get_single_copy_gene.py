import os
import glob
from typing import Tuple, Optional, List
from Bio.Seq import Seq
import pysam
import argparse


def parse_single_copy_buscos(full_table_path: str, num_genes: int) -> List[str]:
    busco_ids = []
    with open(full_table_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            busco_id, status = fields[0], fields[1]
            if status.lower() == "complete":
                busco_ids.append(busco_id)
            if len(busco_ids) >= num_genes:
                break
    return busco_ids


def get_coordinates_from_gff(gff_file: str, busco_id: str) -> Tuple[str, int, int, str]:
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            return chrom, start, end, strand
    raise ValueError(f"No valid feature found in {gff_file} for BUSCO ID {busco_id}")


def extract_sequence(reference_fasta: str, chrom: str, start: int, end: int, strand: str) -> str:
    fasta = pysam.FastaFile(reference_fasta)
    try:
        sequence = fasta.fetch(chrom, start - 1, end)
    except Exception as e:
        raise RuntimeError(f"Could not extract {chrom}:{start}-{end} from FASTA: {e}")
    return str(Seq(sequence).reverse_complement()) if strand == "-" else sequence


def write_multi_fasta(out_file: str, entries: List[Tuple[str, str]]):
    with open(out_file, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract up to N complete single-copy BUSCO gene sequences from a reference genome."
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
    parser.add_argument(
        "--num_genes",
        type=int,
        default=3,
        help="Number of single-copy complete BUSCO genes to extract (default: 3)"
    )

    args = parser.parse_args()

    busco_run_path = os.path.join(args.busco_dir, args.busco_set)
    full_table_path = os.path.join(busco_run_path, "full_table.tsv")
    gff_dir = os.path.join(busco_run_path, "busco_sequences","single_copy_busco_sequences")

    if not os.path.isfile(full_table_path):
        parser.error(f"full_table.tsv not found in {busco_run_path}")

    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference FASTA not found: {args.reference_fasta}")

    if not os.path.isdir(gff_dir):
        parser.error(f"GFF directory not found: {gff_dir}")

    # Step 1: BUSCO IDs
    busco_ids = parse_single_copy_buscos(full_table_path, args.num_genes)
    if not busco_ids:
        parser.error("No complete BUSCO genes found in full_table.tsv")
    print(f"Found BUSCO IDs: {', '.join(busco_ids)}")

    # Step 2–4: Extract sequences from per-BUSCO GFFs
    fasta_entries = []
    for busco_id in busco_ids:
        gff_file = os.path.join(gff_dir, f"{busco_id}.gff")
        if not os.path.isfile(gff_file):
            print(f"Warning: GFF file not found for {busco_id}, skipping.")
            continue
        try:
            chrom, start, end, strand = get_coordinates_from_gff(gff_file, busco_id)
            print(f"{busco_id}: {chrom}:{start}-{end} ({strand})")
            seq = extract_sequence(args.reference_fasta, chrom, start, end, strand)
            header = f"{busco_id}_{chrom}_{start}_{end}_{strand}"
            fasta_entries.append((header, seq))
        except Exception as e:
            print(f"Error processing {busco_id}: {e}")

    if not fasta_entries:
        parser.error("No valid BUSCO sequences were extracted.")

    # Step 5: Output
    folder_name = os.path.basename(os.path.dirname(args.busco_dir))
    out_prefix = args.prefix if args.prefix else folder_name
    out_file = f"{out_prefix}_scg.fasta"

    write_multi_fasta(out_file, fasta_entries)
    print(f"Sequences written to: {out_file}")


if __name__ == "__main__":
    main()
