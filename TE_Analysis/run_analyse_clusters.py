#!/usr/bin/env python3

import argparse
import os
import glob
import json
import re
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Count unique individuals in cluster FASTA files and check for DFAM hits."
    )
    parser.add_argument(
        'species_dirs',
        nargs='+',
        help='One or more species directories (e.g., GenomeDelta/Dbus01)'
    )
    parser.add_argument(
        '--min_individuals',
        dest='min_threshold',
        type=int,
        default=0,
        help='Minimum number of unique individuals to include'
    )
    parser.add_argument(
        '--min_ratio',
        dest='min_ratio',
        type=float,
        default=0.00,
        help='Minimum ratio of sequences to unique individuals to include'
    )
    parser.add_argument(
        '--dfam_dir',
        default='DFAM_Analysis',
        help='Base directory containing DFAM analysis results'
    )
    parser.add_argument(
        '--only-dfam-hits',
        action='store_true',
        help='Only include clusters with DFAM hits'
    )
    return parser.parse_args()

def find_cluster_dir(species_dir):
    for root, dirs, _ in os.walk(species_dir):
        for d in dirs:
            if d.endswith('-GD-clusters'):
                return os.path.join(root, d)
        break  # only check top-level
    return None

import re
from collections import defaultdict

def analyze_fasta_all_regions(fasta_path):
    """
    Parses a FASTA file to:
    - Count unique individuals and total sequences
    - Extract and merge coordinate regions for all reference regions

    Returns:
    - unique_individual_count: int
    - total_sequence_count: int
    - unique_ids: set of str
    - merged_regions_str: str (e.g., 'NC_046604.1:1000-2000,NC_046604.1:2500-3000')
    """
    unique_ids = set()
    num_sequences = 0
    has_pipe = False

    # Matches any region like NC_XXXXXX.X:START-END
    pattern = re.compile(r"(?P<ref>[^\s:]+):(?P<start>\d+)-(?P<end>\d+)")

    regions_by_ref = defaultdict(list)

    with open(fasta_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                num_sequences += 1
                header = line.strip().lstrip('>')

                # Track unique individuals
                if '|' in header:
                    has_pipe = True
                    id_part = header.split('|')[0]
                    unique_ids.add(id_part)

                    header = header.split('|')[1]  # Use the part after the pipe

                # Extract region
                match = pattern.search(header)
                if match:
                    ref = match.group('ref')
                    start = int(match.group('start'))
                    end = int(match.group('end'))
                    regions_by_ref[ref].append((start, end))

    unique_count = len(unique_ids) if has_pipe else 1

    # Merge overlapping regions per reference
    merged_region_list = []
    for ref, regions in regions_by_ref.items():
        merged = []
        regions.sort()
        current_start, current_end = regions[0]
        for start, end in regions[1:]:
            if start <= current_end:  # Overlapping or adjacent
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))

        # Format each merged region as ref:start-end
        for start, end in merged:
            merged_region_list.append(f"{ref}:{start}-{end}")

    return unique_count, num_sequences, unique_ids, merged_region_list

import os
import glob
import json

def check_dfam_hit(species_name, cluster_name, dfam_base_dir):
    cluster_id = os.path.splitext(cluster_name)[0]  # remove .fasta
    dfam_species_dir = os.path.join(dfam_base_dir, species_name)
    pattern = os.path.join(dfam_species_dir, f"{cluster_id}.consensus*.json")
    matches = glob.glob(pattern)

    if not matches:
        return "error", "error", []

    for json_file in matches:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

                if isinstance(data, dict) and 'results' in data:
                    any_hits = False
                    any_tandem = False
                    hit_info_list = []

                    for result in data['results']:
                        hits = result.get('hits', [])

                        if hits:
                            any_hits = True
                            for hit in hits:
                                description = f"\"{hit.get('description', '')}\""
                                # Collect coordinates as a tuple (start, end)
                                

                                hit_info_list.append({
                                    "description": description,
                                    "bit_score": hit.get('bit_score', 0),
                                    "start": hit.get("ali_start"),
                                    "end": hit.get("ali_end"),
                                    "strand": hit.get("strand", "?")
                                })

                        if result.get('tandem_repeats'):
                            any_tandem = True

                    return (
                        "yes" if any_hits else "no",
                        "yes" if any_tandem else "no",
                        hit_info_list
                    )
                else:
                    return "error", "error", []
        except Exception:
            return "error", "error", []

    return "error", "error", []


def main():
    args = parse_arguments()
    print("Species\tCluster\tSequences\tUnique_Individuals\tRatio_Seq_Indiv\tDFAM_Hit\tDFAM_Tandem\tIndividuals\tRegions\tTE_Descriptions\tTE_Coordinates\tTE_Strand\tTE_Bit_Score")

    for species_dir in args.species_dirs:
        cluster_dir = find_cluster_dir(species_dir)
        if not cluster_dir:
            print(f"Warning: No cluster directory found in {species_dir}", file=os.sys.stderr)
            continue

        species_name = os.path.basename(os.path.normpath(species_dir))
        fasta_files = glob.glob(os.path.join(cluster_dir, '*.fasta'))

        for fasta_file in fasta_files:
            cluster_name = os.path.basename(fasta_file)
            count_individuals, count_sequences, unique_ids, regions = analyze_fasta_all_regions(fasta_file)

            if count_individuals < args.min_threshold:
                continue

            ratio_seq_indiv = count_sequences / count_individuals if count_individuals > 0 else 0

            if ratio_seq_indiv < args.min_ratio:
                continue

            # sort unique IDs for consistent output
            unique_ids = sorted(unique_ids)

            dfam_hit, dfam_tandem, hit_info_list = check_dfam_hit(species_name, cluster_name, args.dfam_dir)
            if not args.only_dfam_hits or dfam_hit == "yes" or dfam_hit == "error":

                if len(hit_info_list) == 0:
                    print(f"{species_name}\t{cluster_name}\t{count_sequences}\t{count_individuals}\t{ratio_seq_indiv}\t{dfam_hit}\t{dfam_tandem}\t{','.join(unique_ids)}\t{','.join(regions)}\tNo DFAM hits")
                else:
                    for hit_info in hit_info_list:
                        description = hit_info['description']
                        bit_score = hit_info['bit_score']
                        start = hit_info['start']
                        end = hit_info['end']
                        strand = hit_info['strand']

                        if strand == "-":
                            start = hit_info['end']
                            end = hit_info['start']

                        print(f"{species_name}\t{cluster_name}\t{count_sequences}\t{count_individuals}\t{ratio_seq_indiv}\t{dfam_hit}\t{dfam_tandem}\t{','.join(unique_ids)}\t{','.join(regions)}\t{description}\t{start}-{end}\t{strand}\t{bit_score}")

if __name__ == '__main__':
    main()
