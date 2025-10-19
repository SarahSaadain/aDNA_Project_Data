#!/usr/bin/env python
# Roberts_normalization_script_modified.py
# This script processes SAM files to compute summary statistics and normalized coverage for transposable elements (TEs) and single-copy genes (SCGs).
# Author: Robert Kofler (modified)

import re
import argparse
import collections
import os

# Function to calculate the weight of a read based on its CIGAR string and read length
# Only 'M' (match) operations contribute to the weight; other operations are ignored
# Returns: match length / read length

def get_weight(cigar,readlen):
    pe=[]

    # Parse the CIGAR string into (length, operation) tuples
    for fi in re.finditer(r"(\d+)([HSIDMN])", cigar):
        num=int(fi.group(1))
        id=fi.group(2)
        pe.append((num,id))

    matchsum=0

    # Sum the lengths of 'M' operations
    for num,id in pe:
        if id=="M":
            matchsum+=num
        elif id=="D" or id=="N" or id=="I" or id=="S" or id=="H":
            pass
        else:
            raise Exception("unknown cigar"+id)
    return float(matchsum)/float(readlen)

# Function to read a .fai index file and extract TE and SCG sequence names and lengths
# Returns: dict of {name: length}, list of SCG names, list of TE names

def readfai(fai):
     te_seqence_names_list=[]
     scg_sequence_names_list=[]

     te_and_scg_length_directory={}
     for line in open(fai):
          # Example line: LTR65_te\t669
          a = line.rstrip("\n").split("\t")
          
          sequence_length = int(a[1])

          if a[0].endswith("_te"):
               te_sequence_name=a[0][:-3]
               te_seqence_names_list.append(te_sequence_name)
               te_and_scg_length_directory[te_sequence_name] = sequence_length
          elif a[0].endswith("_scg"):
               scg_sequence_name=a[0][:-4]
               scg_sequence_names_list.append(scg_sequence_name)
               te_and_scg_length_directory[scg_sequence_name]=sequence_length
          else:
               # Only _te and _scg suffixes are supported
               raise Exception("Unknown reference (supported: _te, _scg); found: " + a[0])
     
     return(te_and_scg_length_directory, scg_sequence_names_list, te_seqence_names_list)

# Argument parser setup for command-line usage
parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None, dest="sam", required=True, help="A sam file")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="min mapping quality")
parser.add_argument("--fai", type=str, required=True, dest="fai", default=None, help="fai fasta-index (samtools faidx) of the TE database")
parser.add_argument('--summary', type=str, required=False, dest="summary", help="Output file for summary statistics")
parser.add_argument('--coverage', type=str, required=False, dest="coverage", help="Output file for normalized coverage results")
args = parser.parse_args()
min_mapping_quality_threshhold = args.minmq

sam_basename = os.path.splitext(os.path.basename(args.sam.name))[0]
summary_filename = args.summary if args.summary else f"{sam_basename}_summary.tsv"
coverage_filename = args.coverage if args.coverage else f"{sam_basename}_normalized.tsv"

print(f"Starting normalization")
print(f"SAM file: {args.sam.name}")
print(f"FAI file: {args.fai}")
print(f"Summary output file: {summary_filename}")
print(f"Coverage output file: {coverage_filename}")
print(f"Minimum mapping quality threshold: {min_mapping_quality_threshhold}")

#check if input files exist
if not os.path.isfile(args.fai):
    raise Exception(f"FAI file does not exist: {args.fai}")
if not os.path.isfile(args.sam.name):
    raise Exception(f"SAM file does not exist: {args.sam.name}")

# check if output files already exist to prevent overwriting
if os.path.isfile(summary_filename):
    raise Exception(f"Summary output file already exists: {summary_filename}")
if os.path.isfile(coverage_filename):
    raise Exception(f"Coverage output file already exists: {coverage_filename}")

# Read sequence lengths and lists from the .fai file
te_and_scg_length_directory, scg_names_list, te_names_list = readfai(args.fai)

# Initialize counters for TE and SCG weights
te_weights_col  = collections.defaultdict(lambda:0.0)
scg_weights_col = collections.defaultdict(lambda:0.0)

# Initialize summary statistics counters
count_reads=0                      # Total reads
count_mapped_reads=0               # Mapped reads
count_reads_passed_filter=0        # Reads passing mapping quality filter
count_reads_passed_te = 0          # Reads passing mapping to TEs
count_reads_passed_scg = 0         # Reads passing mapping to SCGs
weighted_sum_passed_reads=0.0      # Weighted sum of reads passing filters

sum_te_weights=0.0          # Weighted sum mapping to TEs
sum_scg_weights=0.0         # Weighted sum mapping to SCGs

# Main loop: process each line in the SAM file
print("Processing SAM file and calculating statistics...")
for line in args.sam:
     """
     Example SAM line:
     r1\t16\tM14653_te\t172\t70\t23M\t*\t0\t0\tATGTCGAGTTTCGTGCCGAATAA\tFFFFFFFFFFFFFFFFFFBBBBB\tPG:Z:novoalign\tAS:i:0\tUQ:i:0\tNM:i:0\tMD:Z:23
     """
     if line.startswith("@"):
          continue  # Skip header lines

     line=line.rstrip("\n")
     line_components = line.split("\t")

     flag = int(line_components[1])

     count_reads += 1  # Count all reads

     if flag & 0x004 > 0:   # Remove unmapped reads
          continue

     count_mapped_reads+=1 # Count mapped reads

     read_mapping_quality = int(line_components[4])

     if read_mapping_quality < min_mapping_quality_threshhold:          # Remove reads below min mapping quality
          continue

     count_reads_passed_filter+=1 # Count reads passing mapping quality

     read_weight = get_weight(line_components[5],len(line_components[9])) # Compute weight of the read
     
     weighted_sum_passed_reads+=read_weight
     
     ref_sequence_name = line_components[2]

     # Assign read to TE or SCG, accumulate weighted counts
     if ref_sequence_name.endswith("_te"):
          te_sequence_name = ref_sequence_name[:-3]
          te_weights_col[te_sequence_name] += read_weight
          sum_te_weights += read_weight
          count_reads_passed_te += 1
     elif ref_sequence_name.endswith("_scg"):
          scg_sequence_name = ref_sequence_name[:-4]
          scg_weights_col[scg_sequence_name] += read_weight
          sum_scg_weights += read_weight
          count_reads_passed_scg += 1
     else:
          raise Exception("unknown reference (supported: _te, _scg); found: " + ref_sequence_name)

# Calculate mean SCG coverage for normalization
mean_scg_coverage=0.0
coverage_scg=[]

for scg, weight in scg_weights_col.items():
     length = te_and_scg_length_directory[scg]
     coverage = float(weight)/float(length)
     coverage_scg.append(coverage)
     mean_scg_coverage+=coverage

mean_scg_coverage=mean_scg_coverage/float(len(coverage_scg))

print("Writing summary statistics to file...")
# Print summary statistics to the specified summary file
with open(summary_filename, 'w') as summary_file:
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","all_reads",count_reads))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","mapped_reads",count_mapped_reads))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","reads_with_mapq",count_reads_passed_filter))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","weighted_reads_with_mapq",round(weighted_sum_passed_reads,2)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","mapping_to_tes",count_reads_passed_te))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","mapping_to_scgs",count_reads_passed_scg))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","mapping_to_te_weighted",round(sum_te_weights,2)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary","mapping_to_scg_weighted",round(sum_scg_weights,2)))

print(f"All reads: {count_reads}")
print(f"Mapped reads: {count_mapped_reads}")
print(f"Reads with sufficient mapping quality: {count_reads_passed_filter}")
print(f"Weighted reads with sufficient mapping quality: {weighted_sum_passed_reads}")
print(f"Mapped reads to TEs: {count_reads_passed_te}")
print(f"Mapped reads to SCGs: {count_reads_passed_scg}")
print(f"Sum weighted reads to TEs: {sum_te_weights}")
print(f"Sum weighted reads to SCGs: {sum_scg_weights}")

print("Summary statistics written.")

print("Writing normalized coverage results to file...")
with open(coverage_filename, 'w') as coverage_file:
    coverage_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("type","name","length","weighted_mapped_reads","normalized_coverage"))
    # Normalized coverage for each SCG
    for scg in scg_names_list:
        weight = scg_weights_col[scg]
        length = te_and_scg_length_directory[scg]
        coverage = float(weight)/float(length)
        normalized_coverage = coverage / mean_scg_coverage
        coverage_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("scg",scg,length,round(weight,2),round(normalized_coverage,2)))
    # Normalized coverage for each TE
    for te in te_names_list:
        length = te_and_scg_length_directory[te]
        weight = te_weights_col[te]
        coverage = float(weight)/float(length)
        normalized_coverage = coverage/mean_scg_coverage
        coverage_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("te",te,length,round(weight,2),round(normalized_coverage,2)))
print("Normalized coverage results written.")
print("Script completed successfully.")