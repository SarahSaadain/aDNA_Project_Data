#!/bin/bash

# Default threshold
min_threshold=0

# Parse arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <species_folder1> [<species_folder2> ...] [--min <threshold>]"
    exit 1
fi

species_dirs=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --min)
            min_threshold="$2"
            shift 2
            ;;
        *)
            species_dirs+=("$1")
            shift
            ;;
    esac
done

# Output header
echo -e "Species\tCluster\tUnique_Individuals"

# Process each species
for species_dir in "${species_dirs[@]}"; do
    cluster_dir=$(find "$species_dir" -maxdepth 1 -type d -name "*-GD-clusters" | head -n 1)

    if [[ -z "$cluster_dir" ]]; then
        echo "Warning: No cluster directory found in $species_dir" >&2
        continue
    fi

    species_name=$(basename "$species_dir")

    # Process each fasta file (cluster)
    for fasta in "$cluster_dir"/*.fasta; do
        cluster_name=$(basename "$fasta")

        # Extract unique individual IDs
        count=$(grep '^>' "$fasta" | cut -d'|' -f1 | sed 's/^>//' | sort -u | wc -l)

        if [[ "$count" -ge "$min_threshold" ]]; then
            echo -e "${species_name}\t${cluster_name}\t${count}"
        fi
    done
done
