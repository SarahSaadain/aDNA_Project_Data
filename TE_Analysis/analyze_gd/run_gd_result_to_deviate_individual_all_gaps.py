import os
import glob
import subprocess
from collections import defaultdict
import sys
import argparse
from Bio import SeqIO

def run_fastp(input_fasta, output_fastq, species, work_dir):
    """Run fastp deduplication inside the species-specific folder."""
    print(f"  ▶ Running fastp for {species}...")

    fastp_cmd = [
        "fastp",
        "--in1", os.path.basename(input_fasta),
        "--out1", os.path.basename(output_fastq),
        "--dedup",
        "--disable_adapter_trimming",
        "--thread", "4",
        "--report_title", f"{species} Deduplication Report",
        "--html", f"{species}_fastp_report.html",
        "--json", f"{species}_fastp_report.json"
    ]

    try:
        subprocess.run(fastp_cmd, check=True, cwd=work_dir)
        print(f"  ✔ Deduplicated FASTQ saved to: {output_fastq}")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ fastp failed for {species}: {e}")
        return False

    return True


def run_deviate(fastq_path, work_dir, library_file_species, scg_names):
    """Call deviaTE on the deduplicated FASTQ inside the species-specific folder, with SCG genes."""
    print(f"  ▶ Running deviaTE on {os.path.basename(fastq_path)}")

    deviate_cmd = [
        "deviaTE",
        "--input", os.path.basename(fastq_path)
    ]

    if library_file_species:
        deviate_cmd.append("--library")
        deviate_cmd.append(os.path.basename(library_file_species))

    if scg_names:
        deviate_cmd.append("--single_copy_genes")
        deviate_cmd.extend(scg_names)

    try:
        subprocess.run(deviate_cmd, check=True, cwd=work_dir)
        print(f"  ✔ deviaTE completed in: {work_dir}")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ deviaTE failed in {work_dir}: {e}")


def make_dummy_quality(seq):
    """Return dummy quality scores for each base."""
    return [40] * len(seq)  # 40 corresponds to high quality "I"


def process_species(individual, file_list, output_base_dir, single_copy_gene_file=None, library_file=None):
    """Combine FASTA files and run deduplication + deviaTE for a species."""
    individual_dir = os.path.join(output_base_dir, individual)
    os.makedirs(individual_dir, exist_ok=True)

    combined_path = os.path.join(individual_dir, f"{individual}_combined_candidates.fasta")
    combined_fastq = os.path.join(individual_dir, f"{individual}_combined_candidates.fastq")
    dedup_fastq = os.path.join(individual_dir, f"{individual}_deduplicated_candidates.fastq")

    print(f"\n🔬 Processing species: {individual} ({len(file_list)} files)")

    for folder_name, fasta_path in file_list:
        if not os.path.isfile(fasta_path):
            sys.exit(f"❌ Missing input file: {fasta_path}")

    library_file_species = os.path.join(individual_dir, f"{individual}_te_library.fasta")

    scg_names = []
    if single_copy_gene_file:
        with open(single_copy_gene_file, "r") as scg_file:

            i = 0
            for line in scg_file:
                if line.startswith(">"):
                    i += 1
                    scg_names.append(f"{individual}_SCG_{i}")

        if not os.path.exists(library_file_species):
            print (f"  🧬 Creating TE library for {individual} at {library_file_species}...")
            with open(library_file_species, "w") as outfile:
                # Write the first file (TE library)
                with open(library_file, "r") as lib_file:
                    for line in lib_file:
                        line = line.rstrip()  # Remove trailing whitespace including \n
                        
                        if not line:
                            continue  # skip empty lines

                        outfile.write(line.upper() + "\n")

                # Write the second file (SCG file)
                with open(single_copy_gene_file, "r") as scg_file:
                    i = 0
                    for line in scg_file:
                        line = line.rstrip()  # Remove trailing whitespace including \n

                        if not line:
                            continue  # skip empty lines

                        if line.startswith(">"):
                            i += 1
                            outfile.write(f">{individual}_SCG_{i}\n")
                        else:
                            outfile.write(line.upper() + "\n")

    else:
        scg_names = ["Dmel_rpl32", "Dmel_piwi"]
        library_file_species = None

    if library_file_species != None and not os.path.exists(library_file_species):
        library_file_species = None

    print(f"  🧬 Using single-copy genes: {', '.join(scg_names)}")

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
                            outfile.write(line.upper())
            if single_copy_gene_file:
                with open(single_copy_gene_file, "r") as scg_file:
                    i = 0
                    for line in scg_file:
                        if line.startswith(">"):
                            i += 1
                            outfile.write(f">{individual}_SCG_{i}\n")
                        else:
                            outfile.write(line.upper())

        print(f"  ✔ Combined FASTA saved to: {combined_path}")

    if not os.path.exists(combined_path):
        sys.exit(f"❌ Combined FASTA file not found: {combined_path}")

    with open(combined_fastq, "w") as out_handle:
        for record in SeqIO.parse(combined_path, "fasta"):
            record.letter_annotations["phred_quality"] = make_dummy_quality(record.seq)
            SeqIO.write(record, out_handle, "fastq")

    if os.path.exists(dedup_fastq):
        print(f"  ⏩ Skipping deduplication: {dedup_fastq} already exists.")
    else:
        success = run_fastp(combined_fastq, dedup_fastq, individual, individual_dir)
        if not success:
            return

    run_deviate(dedup_fastq, individual_dir, library_file_species, scg_names)


def main():
    parser = argparse.ArgumentParser(description="Process candidate FASTA files with fastp and deviaTE.")
    parser.add_argument(
        "--input_dir", type=str, default="GenomeDeltaResult",
        help="Input base directory containing species folders."
    )
    parser.add_argument(
        "--output_dir", type=str, default="DeviaTE_Analysis",
        help="Output base directory."
    )
    parser.add_argument(
        "--scg_dir", type=str, default=None,
        help="Directory containing SCG fasta files named like <species>_scg_<id>_nt.fa."
    )
    parser.add_argument(
        "--library_base", type=str, default=None,
        help="Base name for the TE library file."
    )

    parser.add_argument(
        "--folder_prefix", type=str, default=None,
        help="Optional prefix to filter input folders by name (e.g., Dmel, Dana)."
    )

    args = parser.parse_args()

    input_base_dir = args.input_dir
    output_base_dir = args.output_dir
    scg_dir = args.scg_dir if args.scg_dir else input_base_dir

    os.makedirs(output_base_dir, exist_ok=True)

    if not os.path.isdir(input_base_dir):
        sys.exit(f"❌ Input directory not found: {input_base_dir}")

    all_candidate_paths = glob.glob(os.path.join(input_base_dir, "*", "file-GD.fasta"))

    # Filter folders by prefix if provided
    if args.folder_prefix:
        fasta_files = [
            path for path in all_candidate_paths
            if os.path.basename(os.path.dirname(path)).startswith(args.folder_prefix)
        ]
    else:
        fasta_files = all_candidate_paths
    if not fasta_files:
        sys.exit("❌ No candidate FASTA files found.")

    print(f"✅ Found {len(fasta_files)} FASTA files.")

    species_groups = defaultdict(list)
    for filepath in fasta_files:
        folder_name = os.path.basename(os.path.dirname(filepath))
        species_groups[folder_name].append((folder_name, filepath))

    for individual, file_list in species_groups.items():
        species = individual[:4]  # Use first 4 characters as species name
        single_copy_gene_file = os.path.join(scg_dir, f"{species}_scg.fasta")

        if os.path.isfile(single_copy_gene_file):
            print(f"🔍 Using SCG file for individual {individual}: {single_copy_gene_file}")
        else:
            single_copy_gene_file = None
            print(f"⚠️ No SCG file found for individual {individual}, using default SCG genes.")

        process_species(individual, file_list, output_base_dir, single_copy_gene_file, args.library_base)

    print("\n✅ All processing completed.")


if __name__ == "__main__":
    main()
