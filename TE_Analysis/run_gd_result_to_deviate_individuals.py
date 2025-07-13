import os
import glob
import subprocess
from collections import defaultdict
import sys
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


def run_deviate(fastq_path, work_dir):
    """Call deviaTE on the deduplicated FASTQ inside the species-specific folder, with SCG file."""
    print(f"  ▶ Running deviaTE on {os.path.basename(fastq_path)} with SCG")

    # If you are analyzing TEs in Drosophila specifying a --library or --annotation of 
    # reference sequences is optional. By default deviaTE automatically downloads and 
    # uses the TE library from https://github.com/bergmanlab/drosophila-transposons 
    # if no library and annotation are given.
    # For single-copy gene normalization in Drosophila five genes are automatically 
    # added to the library (Dmel_rpl32, Dmel_piwi, Dmel_Act5C, Dmel_RpII140 and Dmel_p53), 
    # which can be used for normalisation

    deviate_cmd = [
        "deviaTE",
        "--input", os.path.basename(fastq_path),
        "--single_copy_genes", "Dmel_rpl32", "Dmel_piwi"
    ]

    try:
        subprocess.run(deviate_cmd, check=True, cwd=work_dir)
        print(f"  ✔ deviaTE completed in: {work_dir}")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ deviaTE failed in {work_dir}: {e}")


# Dummy quality score (same length as sequence, using "I" = high quality)
def make_dummy_quality(seq):
    return [40] * len(seq)  # 40 corresponds to ASCII "I" in Phred+33

def process_species(species, file_list, output_base_dir, single_copy_gene_file=None):
    """Combine FASTA files and run deduplication + deviaTE for a species."""
    species_dir = os.path.join(output_base_dir, species)
    os.makedirs(species_dir, exist_ok=True)

    combined_path = os.path.join(species_dir, f"{species}_combined_candidates.fasta")
    combined_fastq = os.path.join(species_dir, f"{species}_combined_candidates.fastq")
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
            if single_copy_gene_file:
                with open(single_copy_gene_file, "r") as scg_file:
                    for line in scg_file:
                        outfile.write(line)
                        

        print(f"  ✔ Combined FASTA saved to: {combined_path}")


    if not os.path.exists(combined_path):
        sys.exit(f"❌ Combined FASTA file not found: {combined_path}")
    else:
        # Read FASTA and write to FASTQ with fake quality
        with open(combined_fastq, "w") as out_handle:
            for record in SeqIO.parse(combined_path, "fasta"):
                record.letter_annotations["phred_quality"] = make_dummy_quality(record.seq)
                SeqIO.write(record, out_handle, "fastq")

    # Deduplicate if needed
    if os.path.exists(dedup_fastq):
        print(f"  ⏩ Skipping deduplication: {dedup_fastq} already exists.")
    else:
        success = run_fastp(combined_fastq, dedup_fastq, species, species_dir)
        if not success:
            return  # Skip deviaTE if fastp failed

    # Run deviaTE
    run_deviate(dedup_fastq, species_dir)


def main():
    single_copy_gene_file = "dmel_single_copy_gene.fasta"  # Path to single-copy gene file if needed
    input_base_dir = "GenomeDeltaResult"
    output_base_dir = "DeviaTE_Analysis_candidates"
    os.makedirs(output_base_dir, exist_ok=True)

    if not os.path.isfile(single_copy_gene_file):
        print(f"⚠️ Warning: Single-copy gene file not found: {single_copy_gene_file}. Proceeding without it.")
        single_copy_gene_file = None

    if not os.path.isdir(input_base_dir):
        sys.exit(f"❌ Input directory not found: {input_base_dir}")

    fasta_files = glob.glob(os.path.join(input_base_dir, "*", "file-GD.fasta"))
    if not fasta_files:
        sys.exit("❌ No candidate FASTA files found.")

    print(f"✅ Found {len(fasta_files)} FASTA files.")

    for species, file_list in species_groups.items():
        process_species(species, file_list, output_base_dir, single_copy_gene_file)

    print("\n✅ All processing completed.")


if __name__ == "__main__":
    main()


