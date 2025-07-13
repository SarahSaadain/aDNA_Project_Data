#!/bin/bash

# Create base output directory
mkdir -p /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis

# Loop over all input FASTA files
for fasta in /mnt/data5/sarah/aDNA/D*/raw/ref_genome/*.fna; do
    # Get species identifier: extract species from the full path
    # Example: Dhis/raw/ref_genome/genome.fna → Dhis/raw/ref_genome
    species=$(dirname "${fasta#/mnt/data5/sarah/aDNA/}")

    # Remove the "raw/ref_genome" part to get just the species identifier
    species=${species%/raw/ref_genome}

    # Define output path inside BUSCO_Analysis
    species_output_dir=/mnt/data5/sarah/TE_Analysis/BUSCO_Analysis/$species
    output_dir_busco=$species_output_dir/${species}_busco
    
    # Create the output directories
    mkdir -p "$species_output_dir"
    mkdir -p "$output_dir_busco"

    dataset_path=/mnt/data5/sarah/TE_Analysis/BUSCO_Analysis/busco_downloads/lineages/drosophilidae_odb12

    # log path
    log_path=$species_output_dir/busco_${species}.log

    echo "Running BUSCO on $fasta → $output_dir_busco"

    # Note: busco will create temporary files in the dataset directory. if the processes are running in parallel, they compete for the same files.
    # this is a problem and will cause BUSCO to fail / not create the output files.
    # to avoid this, we run BUSCO sequentially for each species.
    python -u /home/sarah/miniconda3/bin/busco  --in $fasta --out $output_dir_busco -l $dataset_path -m genome -c 16 -f 2>&1 &

done
