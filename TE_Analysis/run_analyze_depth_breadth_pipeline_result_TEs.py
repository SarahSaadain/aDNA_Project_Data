import pandas as pd
import argparse
import sys

def read_input_file(file_path: str) -> pd.DataFrame:
    """Reads the CSV file into a pandas DataFrame."""
    try:
        df = pd.read_csv(file_path)
        return df
    except Exception as e:
        print(f"Error reading file '{file_path}': {e}")
        sys.exit(1)

def pivot_percent_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Pivots the DataFrame to have scaffold as rows and Filename as columns with percent_covered as values."""
    if 'scaffold' not in df or 'Filename' not in df or 'percent_covered' not in df:
        print("Input file must contain 'scaffold', 'Filename', and 'percent_covered' columns.")
        sys.exit(1)
        
    pivot_df = df.pivot_table(
        index='scaffold',
        columns='Filename',
        values='percent_covered'
    )
    return pivot_df.round(4)

def save_output(df: pd.DataFrame, output_file: str):
    """Saves the pivoted DataFrame to a CSV file."""
    df.to_csv(output_file)
    print(f"Pivoted data saved to: {output_file}")

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Pivot scaffold percent coverage from a CSV file.")
    parser.add_argument("input_file", help="Path to the input CSV file.")
    parser.add_argument("-o", "--output", default="pivoted_output.csv", help="Path to save the output CSV file.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    df = read_input_file(args.input_file)
    pivoted_df = pivot_percent_coverage(df)
    save_output(pivoted_df, args.output)

if __name__ == "__main__":
    main()
