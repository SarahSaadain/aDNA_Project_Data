#!/bin/bash

# Make sure blast_results/ exists
mkdir -p blast_results

# Loop over all matching FASTA files
for fasta in */file-GD-candidates.fasta; do
    # Extract the folder name (e.g., Dbus02)
    folder=$(dirname "$fasta")
    
    # Set the output directory
    output_dir="blast_results/$folder"

    # Make the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Run the BLAST script
    echo "Running BLAST for $fasta..."
    python blast_fasta_against_ncbi.py "$fasta" "$output_dir"
done
