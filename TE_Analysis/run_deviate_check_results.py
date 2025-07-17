import os
import glob

# Directory containing all subfolders
base_dir = "/mnt/data5/sarah/TE_Analysis/DeviaTE_Analysis_candidates"

# Columns we're interested in
columns_of_interest = ["A", "C", "G", "T", "cov", "hq_cov"]

def has_nonzero_values(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) < 10:
                continue  # skip malformed lines

            # Extract relevant columns by index
            try:
                values = list(map(float, parts[4:10]))  # A, C, G, T, cov, hq_cov
            except ValueError:
                continue  # skip lines with non-numeric values

            if any(val != 0.0 for val in values):
                return True
    return False

def find_candidates(base_dir):
    pattern = os.path.join(base_dir, "**", "*.deviate")
    all_files = glob.glob(pattern, recursive=True)
    for file_path in all_files:
        if has_nonzero_values(file_path):
            print(file_path)

if __name__ == "__main__":
    find_candidates(base_dir)

