#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import pysam

def bam_to_bed_with_coverage(bam_file, fai_file):
    # Prepare output filename
    base_name = os.path.splitext(os.path.basename(bam_file))[0]
    out_file = f"{base_name}_mapped_reads_coverage.tsv"

    # Load TE lengths from .fai
    te_lengths = {}
    with open(fai_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            read_name = parts[0]
            read_length = int(parts[1])
            te_lengths[read_name] = read_length

    # Open BAM file with pysam
    bam = pysam.AlignmentFile(bam_file, "rb")

    records = []
    for read in bam.fetch(until_eof=True):
        chrom = bam.get_reference_name(read.reference_id)
        start = read.reference_start
        end = read.reference_end
        strand = "-" if read.is_reverse else "+"
        mapq = read.mapping_quality
        read_name = read.query_name
        aln_len = end - start
        read_length = te_lengths.get(read_name, None)
        percent_cov = (aln_len / read_length * 100) if read_length else None

        # Alignment type
        if read.is_secondary:
            aln_type = "secondary"
        elif read.is_supplementary:
            aln_type = "supplementary"
        else:
            aln_type = "primary"

        # Add read coordinates
        read_start = read.query_alignment_start   # 0-based start on the read
        read_end = read.query_alignment_end       # 0-based end (exclusive) on the read

        records.append({
            "Chromosome": chrom,
            "Start": start,
            "End": end,
            "MAPQ": mapq,
            "Strand": strand, 
            "Read_name": read_name,
            "Alignment_type": aln_type,
            "read_length": read_length if read_length else "NA",
            "Align_len": aln_len,
            "Read_start": read_start,
            "Read_end": read_end,
            "Percent_coverage": percent_cov if percent_cov else "NA"
        })


    bam.close()

    # Convert to DataFrame
    df = pd.DataFrame(records)
    # Sort by Percent_coverage descending
    df['Percent_coverage'] = pd.to_numeric(df['Percent_coverage'], errors='coerce')
    #df.sort_values(by='Percent_coverage', ascending=False, inplace=True)
    # Save as TSV
    df.to_csv(out_file, sep="\t", index=False)
    print(f"Output written to: {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Calculate TE coverage from BAM alignments and TE FASTA index")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (sorted and indexed)")
    parser.add_argument("-f", "--fai", required=True, help="TE FASTA index file (.fai)")

    args = parser.parse_args()
    bam_to_bed_with_coverage(args.bam, args.fai)


if __name__ == "__main__":
    main()
