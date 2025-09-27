#!/usr/bin/env python3

# Import required libraries
import argparse    # For handling command-line arguments
import csv         # For reading tab- or comma-delimited files
import sys         # For exiting the program if errors occur

def read_samples(path):
    """
    Read the first input (one-column file with sample IDs).
    Returns a list of sample IDs.
    """
    samples = []                 # Initialize empty list
    with open(path, "r") as f:   # Open file in read mode
        for line in f:           # Iterate over lines
            line = line.strip()  # Remove whitespace and newline
            if line:             # Skip empty lines
                samples.append(line)  # Add sample to list
    return samples               # Return list

def read_values(path, col_index):
    """
    Read the second input file (multi-column).
    Ignore column 1, use column 2 for sample IDs,
    and a user-specified column for values.
    Returns a dictionary mapping sample ID -> value.
    """
    values = {}
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter=",")  # force tab-delimited
        for parts in reader:
            if len(parts) <= col_index:
                continue
            sample_id = parts[1]
            value = parts[col_index]
            values[sample_id] = value
    return values

def main():
    parser = argparse.ArgumentParser(description="Assign values to samples")
    parser.add_argument("--file1", required=True, help="One-column file of sample IDs")
    parser.add_argument("--file2", required=True, help="Multi-column file with values")
    parser.add_argument("--col", required=True, type=int, help="Column index (>2) for values, 1-based")
    args = parser.parse_args()

    col_index = args.col - 1
    if col_index < 2:  # Ensure column >=3
        sys.exit("Error: column index must be >2 (3 onwards).")

    samples = read_samples(args.file1)      # Read sample IDs
    values = read_values(args.file2, col_index)  # Read sample->value mapping

    for sample in samples:
        if sample in values:                 # If sample exists
            val = values[sample]            # Assign value
        else:
            val = "-9"                       # Else assign -9
        print(f"{sample}\t{sample}\t{val}")  # Output tab-separated

if __name__ == "__main__":
    main()
