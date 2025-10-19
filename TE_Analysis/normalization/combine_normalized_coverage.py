#!/usr/bin/env python3
import os
import glob
import csv

def get_individual_name(filename):
    # Extract the individual name from the filename (before _normalized.tsv)
    base = os.path.basename(filename)
    return base.split('_normalized.tsv')[0]

# Find all _normalized.tsv files (recursively)
normalized_files = glob.glob('**/*_normalized.tsv', recursive=True)

te_coverage = {}  # {te_name: {individual: normalized_coverage}}
individuals = []

for file in normalized_files:
    individual = get_individual_name(file)
    individuals.append(individual)
    with open(file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['type'] == 'te':
                te = row['name']
                coverage = row['normalized_coverage']
                if te not in te_coverage:
                    te_coverage[te] = {}
                te_coverage[te][individual] = coverage

# Sort individuals and TEs
individuals = sorted(set(individuals))
te_names = sorted(te_coverage.keys())

# Write combined output
output_file = 'combined_normalized_coverage.tsv'
with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['te'] + individuals)
    for te in te_names:
        row = [te]
        for ind in individuals:
            row.append(te_coverage[te].get(ind, 'NA'))
        writer.writerow(row)

print(f"Combined normalized coverage written to {output_file}")
