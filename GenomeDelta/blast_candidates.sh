#!/bin/bash

# Make sure blast_results/ exists
mkdir -p blast_results

# Loop over all matching FASTA files
for fasta in */file-GD-candidates.fasta; do
    # Extract the folder name (e.g., Dbus02)
    folder=$(dirname "$fasta")
    
    # Set the output directory
    output_dir_blast="blast_results/$folder"
    # Make the output directory if it doesn't exist
    mkdir -p "$output_dir_blast"

    echo "Running BLAST for $fasta..."
    python blast_fasta_against_ncbi.py "$fasta" "$output_dir_blast"

    # Set the output directory
    output_dir_dfam="dfam_results/$folder"
    # Make the output directory if it doesn't exist
    mkdir -p "$output_dir_dfam"

    echo "Running DFAM for $fasta..."
    python blast_fasta_against_dfam.py "$fasta" "$output_dir_dfam"
done
