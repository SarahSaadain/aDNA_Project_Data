#!/usr/bin/env python3

import argparse
import os
import glob
import json

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

def analyze_fasta_all_regions(fasta_path):
    """
    Parses a FASTA file to:
    - Count unique individuals and total sequences
    - Extract and merge coordinate regions for all reference regions

    Returns:
    - unique_individual_count: int
    - total_sequence_count: int
    - unique_ids: set of str
    - merged_regions: dict {region_prefix: list of merged (start, end) tuples}
    """
    unique_ids = set()
    num_sequences = 0
    has_pipe = False

    # Matches any region like NC_XXXXXX.X:START-END
    pattern = re.compile(r"(?P<ref>NC_\d+\.\d+):(?P<start>\d+)-(?P<end>\d+)")

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

                # Extract regions
                match = pattern.search(header)
                if match:
                    ref = match.group('ref')
                    start = int(match.group('start'))
                    end = int(match.group('end'))
                    regions_by_ref[ref].append((start, end))

    unique_count = len(unique_ids) if has_pipe else 1

    # Merge overlapping regions per reference
    merged_regions = {}
    for ref, regions in regions_by_ref.items():
        merged = []
        regions.sort()
        current_start, current_end = regions[0]
        for start, end in regions[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))
        merged_regions[ref] = merged

    return unique_count, num_sequences, unique_ids, merged_regions

def count_individuals_and_sequences(fasta_path):
    unique_ids = set()
    num_sequences = 0
    has_pipe = False

    with open(fasta_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                num_sequences += 1
                header = line.strip().lstrip('>')
                if '|' in header:
                    has_pipe = True
                    id_part = header.split('|')[0]
                    unique_ids.add(id_part)

    unique_count = len(unique_ids) if has_pipe else 1
    return unique_count, num_sequences, unique_ids

def check_dfam_hit(species_name, cluster_name, dfam_base_dir):
    cluster_id = os.path.splitext(cluster_name)[0]  # remove .fasta
    dfam_species_dir = os.path.join(dfam_base_dir, species_name)
    pattern = os.path.join(dfam_species_dir, f"{cluster_id}.consensus*.json")
    matches = glob.glob(pattern)

    if not matches:
        return "error", "error"

    for json_file in matches:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

                if isinstance(data, dict) and 'results' in data:
                    for result in data['results']:
                        has_hits = bool(result.get('hits'))
                        has_tandem = bool(result.get('tandem_repeats'))
                        return ("yes" if has_hits else "no", 
                                "yes" if has_tandem else "no")
                else:
                    return "error", "error"
        except Exception:
            return "error", "error"

    return "error", "error"

def main():
    args = parse_arguments()
    print("Species\tCluster\tSequences\tUnique_Individuals\tRatio_Seq_Indiv\tDFAM_Hit\tDFAM_Tandem\tIndividuals")

    for species_dir in args.species_dirs:
        cluster_dir = find_cluster_dir(species_dir)
        if not cluster_dir:
            print(f"Warning: No cluster directory found in {species_dir}", file=os.sys.stderr)
            continue

        species_name = os.path.basename(os.path.normpath(species_dir))
        fasta_files = glob.glob(os.path.join(cluster_dir, '*.fasta'))

        for fasta_file in fasta_files:
            cluster_name = os.path.basename(fasta_file)
            count_individuals, count_sequences, unique_ids = count_individuals_and_sequences(fasta_file)

            if count_individuals < args.min_threshold:
                continue

            ratio_seq_indiv = count_sequences / count_individuals if count_individuals > 0 else 0

            if ratio_seq_indiv < args.min_ratio:
                continue

            # sort unique IDs for consistent output
            unique_ids = sorted(unique_ids)

            dfam_hit, dfam_tandem = check_dfam_hit(species_name, cluster_name, args.dfam_dir)
            if not args.only_dfam_hits or dfam_hit == "yes" or dfam_hit == "error":
                print(f"{species_name}\t{cluster_name}\t{count_sequences}\t{count_individuals}\t{ratio_seq_indiv}\t{dfam_hit}\t{dfam_tandem}\t{','.join(unique_ids)}")

if __name__ == '__main__':
    main()
