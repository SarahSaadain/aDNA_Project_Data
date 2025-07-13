#!/bin/bash

# Default parameters
min_bitscore=1000
refine_d=2500
remove_temp=0
prefix="file"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gap_fasta) gap_fasta="$2"; shift ;;
        --of) mapped_folder="$2"; shift ;;
        --prefix) prefix="$2"; shift ;;
        --min_bitscore) min_bitscore="$2"; shift ;;
        --refine) refine_set=1 ;;
        --refine_d) refine_d="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$gap_fasta" || -z "$mapped_folder" ]]; then
    echo "Usage: $0 --gap_fasta <gap_fasta_file> --of <mapped_folder> [--prefix <prefix>] [--min_bitscore <value>] [--refine] [--refine_d <value>]"
    exit 1
fi

# Prepare variables
filename="$prefix"
current_dir=$(dirname "$(readlink -f "$0")")

mkdir -p "$mapped_folder"

echo "Starting GenomeDelta from GAPS FASTA file..."
echo "Input Gaps FASTA: $gap_fasta"
echo "Mapped Folder: $mapped_folder"
echo "Prefix: $prefix"
echo "Minimum Bitscore: $min_bitscore"
[[ $refine_set -eq 1 ]] && echo "Refine Mode: ENABLED with distance $refine_d"

# BLAST analysis
blast_out="${mapped_folder}/${filename}-GD.blast"
if [[ ! -f "$blast_out" ]]; then
    echo "Running BLAST to find repetitive gaps..."
    blastn -query "$gap_fasta" -subject "$gap_fasta" -out "${blast_out}.tmp" -outfmt 6
    awk -v min="$min_bitscore" '($12)>=min' "${blast_out}.tmp" > "$blast_out"

    [[ "$remove_temp" -eq 1 ]] && rm "${blast_out}.tmp"
else
    echo "BLAST file $blast_out already exists, skipping."
fi

# Check BLAST output
if [[ ! -s "$blast_out" ]]; then
    echo "ERROR: No repetitive gaps found. Check input or reduce stringency."
    exit 1
fi

# Extract non-repetitive sequences
non_rep_fa="${mapped_folder}/${filename}-GD-non_rep.fasta"
clusters_dir="${mapped_folder}/${filename}-GD-clusters/"

if [[ ! -f "$non_rep_fa" ]]; then
    mkdir -p "$clusters_dir"
    python "$current_dir/scripts/blast2clusters.py" "$blast_out" "$gap_fasta" "$clusters_dir"
    mv "$clusters_dir/non_rep.fasta" "$non_rep_fa"
else
    echo "Non-repetitive fasta file already exists: $non_rep_fa"
fi

# Generate consensus sequences
echo "Extracting consensus sequences..."

count=$(ls "$clusters_dir"/*.fasta 2>/dev/null | wc -l)
if [[ "$count" -eq 0 ]]; then
    echo "No repetitive cluster FASTA files found. Exiting."
    [[ "$remove_temp" -eq 1 ]] && rm -r "$clusters_dir"
    exit 1
fi

i=0
for fasta in "$clusters_dir"/*.fasta; do
    i=$((i+1))
    echo "Processing ($i/$count): $fasta"

    msa="${fasta%.fasta}.MSA"
    consensus="${fasta%.fasta}.consensus"

    if [[ -f "$consensus" ]]; then
        echo "Consensus already exists for $fasta"
        continue
    fi

    echo "Aligning with MAFFT..."
    mafft --thread 4 --auto "$fasta" > "$msa"

    if [[ ! -s "$msa" ]]; then
        echo "MAFFT failed or produced empty alignment. Skipping."
        rm -f "$msa"
        continue
    fi

    echo "Generating consensus..."
    python "$current_dir/scripts/MSA2consensus.py" "$msa" "$consensus"

    if [[ "$remove_temp" -eq 1 ]]; then
        rm -f "$msa"
    fi
done

# Optional refinement
if [[ "$refine_set" -eq 1 ]]; then
    refined_dir="${mapped_folder}/${filename}-GD-clusters-refined"
    mkdir -p "$refined_dir"
    bash "$current_dir/scripts/find-coupled-clusters.sh" "$clusters_dir" "$refined_dir" "$refine_d"
fi

# Concatenate all consensus sequences
candidates_fa="${mapped_folder}/${filename}-GD-candidates.fasta"
if [[ ! -f "$candidates_fa" ]]; then
    echo "Concatenating consensus sequences..."
    cat "$clusters_dir"/*.consensus > "$candidates_fa"
    
    if [[ "$refine_set" -eq 1 ]]; then
        find "$refined_dir" -name "*.consensus" -exec cat {} + >> "$candidates_fa"
    fi
fi

# Indexing
echo "Indexing final FASTA files..."
samtools faidx "$candidates_fa"
samtools faidx "$non_rep_fa"

# Visualization
echo "Running visualization..."
Rscript "$current_dir/scripts/visualization.R" \
    --rep "${candidates_fa}.fai" \
    --nonrep "${non_rep_fa}.fai" \
    --output1 "${mapped_folder}/${filename}-GD-candidates.png" \
    --output2 "${mapped_folder}/${filename}-GD-n_
