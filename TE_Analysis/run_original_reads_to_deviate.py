import os
import glob
import subprocess
from collections import defaultdict
import sys
import argparse
from Bio import SeqIO

def run_deviate(fastq_path, work_dir, library_file_species, scg_names):
    """Call deviaTE on the deduplicated FASTQ inside the species-specific folder, with SCG genes."""
    print(f"  ▶ Running deviaTE on {os.path.basename(fastq_path)}")

    deviate_cmd = [
        "deviaTE",
        "--input", fastq_path
    ]

    if library_file_species:
        deviate_cmd.append("--library")
        deviate_cmd.append(os.path.basename(library_file_species))

    if scg_names:
        deviate_cmd.append("--single_copy_genes")
        deviate_cmd.extend(scg_names)
    else:
        deviate_cmd.append("--rpm")

    print(f"Deviate command: {' '.join(deviate_cmd)}")

    try:
        subprocess.run(deviate_cmd, check=True, cwd=work_dir)
        print(f"  ✔ deviaTE completed in: {work_dir}")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ deviaTE failed in {work_dir}: {e}")


def process_fastq_file(fastq_path, output_base_dir, scg_dir, library_file=None):
    """Process a single .fastq.gz file with deviaTE."""
    filename = os.path.basename(fastq_path)
    #take everything before the first dot as the individual
    individual = filename.split(".")[0]

    #tfirst 4 letter are the species
    species = individual[:4]

    individual_dir = os.path.join(output_base_dir, individual)
    os.makedirs(individual_dir, exist_ok=True)

    print(f"\n🔬 Processing FASTQ.GZ: {filename}")

    # 🆕 Prepare SCG file
    single_copy_gene_file = os.path.join(scg_dir, f"{species}_scg.fasta")
    if os.path.isfile(single_copy_gene_file):
        print(f"🔍 Using SCG file: {single_copy_gene_file}")
    else:
        print(f"⚠️ No SCG file found for {species}, using default SCGs.")
        single_copy_gene_file = None

    # 🆕 Prepare SCG names and library
    scg_names = []
    library_file_species = os.path.join(individual_dir, f"{individual}_te_library.fasta")

    if single_copy_gene_file:
        with open(single_copy_gene_file, "r") as scg_file:
            i = 0
            for line in scg_file:
                if line.startswith(">"):
                    i += 1
                    scg_names.append(f"{individual}_SCG_{i}")

        if not os.path.exists(library_file_species):
            print(f"  🧬 Creating TE+SCG library at {library_file_species}...")
            with open(library_file_species, "w") as outfile:
                with open(library_file, "r") as lib_file:
                    for line in lib_file:
                        line = line.rstrip()
                        if not line:
                            continue
                        outfile.write(line.upper() + "\n")

                with open(single_copy_gene_file, "r") as scg_file:
                    i = 0
                    for line in scg_file:
                        line = line.rstrip()
                        if not line:
                            continue
                        if line.startswith(">"):
                            i += 1
                            outfile.write(f">{individual}_SCG_{i}\n")
                        else:
                            outfile.write(line.upper() + "\n")
    else:
        scg_names = None
        library_file_species = library_file

    if library_file_species and not os.path.exists(library_file_species):
        library_file_species = None

    if scg_names:
        print(f"  🧬 Using single-copy genes: {', '.join(scg_names)}")
    else:
        print("  🧬 No single-copy genes found.")

    run_deviate(fastq_path, individual_dir, library_file_species, scg_names)


def main():
    parser = argparse.ArgumentParser(description="Run deviaTE on FASTQ files in a given folder.")
    parser.add_argument("--read_dir", type=str, help="Directory containing .fastq.gz files")  # 🆕
    parser.add_argument("--output_dir", type=str, default="DeviaTE_Analysis_original", help="Output base directory.")
    parser.add_argument("--scg_dir", type=str, help="Directory containing SCG fasta files")
    parser.add_argument("--library_base", type=str, help="TE library fasta file to use")

    args = parser.parse_args()

    if not args.read_dir or not os.path.isdir(args.read_dir):
        sys.exit("❌ You must provide a valid --read_dir path.")

    os.makedirs(args.output_dir, exist_ok=True)
    scg_dir = args.scg_dir if args.scg_dir else args.read_dir

    fastq_files = glob.glob(os.path.join(args.read_dir, "*.fastq.gz"))
    if not fastq_files:
        sys.exit("❌ No .fastq.gz files found in the specified --read_dir.")

    print(f"✅ Found {len(fastq_files)} FASTQ files.")

    for fastq_path in fastq_files:
        process_fastq_file(fastq_path, args.output_dir, scg_dir, args.library_base)

    print("\n✅ All FASTQ files processed.")


if __name__ == "__main__":
    main()
