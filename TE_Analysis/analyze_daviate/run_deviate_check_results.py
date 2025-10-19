import os
import glob
import argparse

def calculate_coverage(file_path):
    total_lines = 0
    nonzero_lines = 0

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) < 10:
                continue  # skip malformed lines

            total_lines += 1

            try:
                values = list(map(float, parts[4:10]))  # A, C, G, T, cov, hq_cov
            except ValueError:
                continue

            if any(val != 0.0 for val in values):
                nonzero_lines += 1

    coverage = (nonzero_lines / total_lines * 100.0) if total_lines > 0 else 0.0
    return coverage, nonzero_lines, total_lines

def find_candidates(base_dir, folder_prefix="", include_scg=False,
                    check_coverage=True, coverage_minimum=0.0):
    pattern = os.path.join(base_dir, "**", "*.deviate")
    all_files = glob.glob(pattern, recursive=True)

    all_files.sort()

    # Print header
    print("file_path\tnonzero_lines\ttotal_lines\tcoverage_percent")

    for file_path in all_files:
        relative_path = os.path.relpath(file_path, base_dir)
        parts = relative_path.split(os.sep)
        if not parts or not parts[0].startswith(folder_prefix):
            continue

        if not include_scg and "_SCG_" in os.path.basename(file_path):
            continue

        if check_coverage:
            coverage, nonzero, total = calculate_coverage(file_path)
            if coverage >= coverage_minimum:
                print(f"{relative_path}\t{nonzero}\t{total}\t{coverage:.2f}")
        else:
            if has_nonzero_values(file_path):
                print(f"{relative_path}\tNA\tNA\tNA\tyes")

def has_nonzero_values(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) < 10:
                continue
            try:
                values = list(map(float, parts[4:10]))
            except ValueError:
                continue
            if any(val != 0.0 for val in values):
                return True
    return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scan .deviate files for nonzero entries.")
    parser.add_argument("base_dir", help="Base directory containing subfolders with .deviate files.")
    parser.add_argument("--folder_prefix", default="", help="Only analyze folders that start with this prefix.")
    parser.add_argument("--include_scg", action="store_true", help="Include files with '_SCG_' in the name.")
    parser.add_argument("--check_coverage", action=argparse.BooleanOptionalAction, default=True,
                        help="Calculate percent of non-zero coverage (default: True).")
    parser.add_argument("--coverage_minimum", type=float, default=0.0,
                        help="Minimum percent coverage required to include a file in output.")

    args = parser.parse_args()

    find_candidates(
        args.base_dir,
        folder_prefix=args.folder_prefix,
        include_scg=args.include_scg,
        check_coverage=args.check_coverage,
        coverage_minimum=args.coverage_minimum
    )
