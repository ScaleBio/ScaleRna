#!/usr/bin/env python
"""
Pre-process samples.csv to simplify parsing downstream.
Rename deprecriated column names and fill in defaults
Basic validation

Writes the normalized samples.csv to stdout
"""
import argparse
import csv
import sys

def validateName(name):
    """
    Check sample and library names for invalid characters
    Print error and exit for invalid names

    Args:
        name: sample/library name to check
    """
    for n in name:
        if not (n.isalnum() or n in "-."):
            print(f"Name should only contain [a-z],[A-Z],[0-9], dash (-) or dot (.): {name}",
                  file=sys.stderr)
            sys.exit(1)

def main(samplesCsv):
    """
    Writes normalized samples.csv to stdout
    Args:
        samplesCsv: Path to samples.csv file
    """
    rows = {} # sample -> samples.csv row
    with open(samplesCsv) as csvfile:
        samples = csv.reader(csvfile)
        cols = next(samples)
        # Trim empty columns
        while not cols[-1]: 
            cols = cols[:-1]

        if cols[0] != "sample":
            print("First column should be 'sample' (the sample name)", file=sys.stderr)
            sys.exit(1)

        for row in samples:
            if not row or len(row[0].strip()) == 0:
                continue
            # Trim empty columns
            while not row[-1] and len(row) > len(cols):
                row = row[:-1]

            if len(row) != len(cols):
                print("Unexpected number of columns:", row, sep="\n", file=sys.stderr)
                sys.exit(1)
            sample = row[0].strip()
            validateName(sample)
            if sample in rows:
                print(f"Duplicate sample name: {sample}", file=sys.stderr)
                sys.exit(1)
            row[0] = sample
            rows[sample] = row

    if "libName" not in cols and "fastqName" in cols:
        cols[cols.index('fastqName')] = "libName"
    if "libIndex" not in cols and "index" in cols:
        cols[cols.index('index')] = "libIndex"
    if "libIndex" not in cols and "fastqIndex" in cols:
        cols[cols.index('fastqIndex')] = "libIndex"
    # Normalize index sequences
    if "libIndex" in cols:
        libIndexInd = cols.index("libIndex")
        for s,r in rows.items():
            index = "".join(r[libIndexInd].split()) # Remove whitespace
            index = index.upper()
            index = ';'.join(sorted(index.split(';'))) # sort sequences
            r[libIndexInd] = index
        # Default libName based on index
        if "libName" not in cols:
            cols.insert(libIndexInd, "libName")
            for r in rows.values():
                name = r[libIndexInd].split(';')[0]
                r.insert(libIndexInd, name)
    elif "libName" not in cols: # Only one (unnamed) library; Use default name
        cols.insert(1, "libName")
        for r in rows.values():
            r.insert(1, "ScaleRna")
    
    w = csv.writer(sys.stdout)
    w.writerow(cols)
    w.writerows(rows.values())   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Standardize samples.csv for the workflow (defaults, column names, etc.')
    parser.add_argument('samples', metavar='SAMPLES.csv',
                    help='CSV with sample names and information for the workflow run')
    args = parser.parse_args()

    main(args.samples)
