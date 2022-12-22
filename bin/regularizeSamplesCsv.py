#!/usr/bin/env python

# Pre-process samples.csv to simplify parsing downstream.
# Rename deprecriated column names and fill in defaults

import argparse
import csv
import sys

def main(samplesCsv):
    rows = []
    with open(samplesCsv) as csvfile:
        samples = csv.reader(csvfile)
        cols = next(samples)
        for row in samples:
            if not row or len(row[0].strip()) == 0:
                continue
            if len(row) != len(cols):
                print("Unexpected number of columns:",row, sep="\n",file=sys.stderr)
                sys.exit(1)
            rows.append(row)

    if "libName" not in cols and "fastqName" in cols:
        cols[cols.index('fastqName')] = "libName"
    if "libIndex" not in cols and "index" in cols:
        cols[cols.index('index')] = "libIndex"
    if "libIndex" not in cols and "fastqIndex" in cols:
        cols[cols.index('fastqIndex')] = "libIndex"
    
    if "libIndex" in cols: # Normalize index sequences
        libIndexInd = cols.index("libIndex")
        for r in rows:
            index = "".join(r[libIndexInd].split()) # Remove whitespace
            index = index.upper()
            index = ';'.join(sorted(index.split(';'))) # sort sequences
            r[libIndexInd] = index
        if "libName" not in cols: # Default libName based on index
            cols.insert(libIndexInd, "libName")
            for r in rows:
                name = r[libIndexInd].split(';')[0]
                r.insert(libIndexInd, name)
    elif "libName" not in cols: # Only one (unnamed) library; Use default name
        cols.insert(1, "libName")
        for r in rows:
            r.insert(1, "ScaleBio")
    
    w = csv.writer(sys.stdout)
    w.writerow(cols)
    w.writerows(rows)   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Standardize samples.csv for the workflow (defaults, column names, etc.')
    parser.add_argument('samples', metavar='SAMPLES.csv',
                    help='CSV with sample names and information for the workflow run')
    args = parser.parse_args()

    main(args.samples)
