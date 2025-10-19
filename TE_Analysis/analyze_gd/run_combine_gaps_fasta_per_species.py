import os
import glob
import subprocess
from collections import defaultdict
import sys
from Bio import SeqIO


def process_species(species, file_list, output_base_dir):
    """Combine FASTA files and run deduplication + deviaTE for a species."""
    species_dir = os.path.join(output_base_dir, species)
    os.makedirs(species_dir, exist_ok=True)

    combined_path = os.path.join(species_dir, f"{species}_combined_gaps.fasta")

    print(f"\n🔬 Processing species: {species} ({len(file_list)} files)")

    for folder_name, fasta_path in file_list:
        if not os.path.isfile(fasta_path):
            sys.exit(f"❌ Missing input file: {fasta_path}")

    # Combine input files if needed
    if os.path.exists(combined_path):
        print(f"  ⏩ Skipping combine: {combined_path} already exists.")
    else:
        with open(combined_path, "w") as outfile:
            for folder_name, fasta_path in file_list:
                with open(fasta_path, "r") as infile:
                    for line in infile:
                        if line.startswith(">"):
                            header = line.strip()[1:]
                            new_header = f">{folder_name}|{header}"
                            outfile.write(new_header + "\n")
                        else:
                            outfile.write(line)
        print(f"  ✔ Combined FASTA saved to: {combined_path}")

def main():
    input_base_dir = "GenomeDeltaResult"
    output_base_dir = "GenomeDeltaResult"
    os.makedirs(output_base_dir, exist_ok=True)

    if not os.path.isdir(input_base_dir):
        sys.exit(f"❌ Input directory not found: {input_base_dir}")

    fasta_files = glob.glob(os.path.join(input_base_dir, "*", "file-GD.fasta"))
    if not fasta_files:
        sys.exit("❌ No candidate FASTA files found.")

    print(f"✅ Found {len(fasta_files)} FASTA files.")

    species_groups = defaultdict(list)
    for filepath in fasta_files:
        folder_name = os.path.basename(os.path.dirname(filepath))
        species_prefix = folder_name[:4]
        species_groups[species_prefix].append((folder_name, filepath))

    for species, file_list in species_groups.items():
        process_species(species, file_list, output_base_dir)

    print("✅ All processing completed.")


if __name__ == "__main__":
    main()
