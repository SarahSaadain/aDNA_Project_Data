#!/bin/bash

# Make sure folders/ exists
mkdir -p NCBI_Analysis
mkdir -p DFAM_Analysis

# Loop over all matching FASTA files
for fasta in GenomeDeltaResult/*/file-GD-candidates.fasta; do
    # Extract the folder name without GenomeDeltaResult/

    folder=$(dirname "$fasta") 
    folder=${folder#GenomeDeltaResult/}
    
    # Set the output directory
    output_dir_blast="NCBI_Analysis/$folder"
    # Make the output directory if it doesn't exist
    mkdir -p "$output_dir_blast"

    echo "Running BLAST for $fasta..."
    python -u check_fasta_against_ncbi.py "$fasta" "$output_dir_blast"

    # Set the output directory
    output_dir_dfam="DFAM_Analysis/$folder"
    # Make the output directory if it doesn't exist
    mkdir -p "$output_dir_dfam"

    echo "Running DFAM for $fasta..."
    python -u check_fasta_against_dfam.py "$fasta" "$output_dir_dfam"
done
