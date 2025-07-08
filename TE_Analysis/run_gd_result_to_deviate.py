import os
import glob
import subprocess
from collections import defaultdict
import sys

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


def run_deviate(fastq_path, work_dir):
    """Call deviaTE on the deduplicated FASTQ inside the species-specific folder."""
    print(f"  ▶ Running deviaTE on {os.path.basename(fastq_path)}")
    try:
        subprocess.run(["deviaTE", "--input", os.path.basename(fastq_path)], check=True, cwd=work_dir)
        print(f"  ✔ deviaTE completed in: {work_dir}")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ deviaTE failed in {work_dir}: {e}")

def process_species(species, file_list, output_base_dir):
    """Combine FASTA files and run deduplication + deviaTE for a species."""
    species_dir = os.path.join(output_base_dir, species)
    os.makedirs(species_dir, exist_ok=True)

    combined_path = os.path.join(species_dir, f"{species}_combined_candidates.fasta")
    dedup_fastq = os.path.join(species_dir, f"{species}_deduplicated_candidates.fastq")

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

    # Deduplicate if needed
    if os.path.exists(dedup_fastq):
        print(f"  ⏩ Skipping deduplication: {dedup_fastq} already exists.")
    else:
        success = run_fastp(combined_path, dedup_fastq, species, species_dir)
        if not success:
            return  # Skip deviaTE if fastp failed

    # Run deviaTE
    run_deviate(dedup_fastq, species_dir)


def main():
    input_base_dir = "GenomeDeltaResult"
    output_base_dir = "DeviaTE_Analysis"
    os.makedirs(output_base_dir, exist_ok=True)

    if not os.path.isdir(input_base_dir):
        sys.exit(f"❌ Input directory not found: {input_base_dir}")

    fasta_files = glob.glob(os.path.join(input_base_dir, "*", "file-GD-candidates.fasta"))
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

    print("\n✅ All processing completed.")


if __name__ == "__main__":
    main()
