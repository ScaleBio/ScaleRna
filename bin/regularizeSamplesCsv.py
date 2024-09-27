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
import json
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List

def validateName(name):
    """
    Check sample name or ID for invalid characters
    Print error and exit for invalid names

    Args:
        name: sample name or ID
    """
    for n in name:
        if not (n.isalnum() or n in "-."):
            print(f"Name should only contain [a-z],[A-Z],[0-9], dash (-) or dot (.): {name}",
                  file=sys.stderr)
            sys.exit(1)
    if not name[0].isalpha():
        print(f"Name should start with a letter: {name}",
                file=sys.stderr)
        sys.exit(1)

def makeUniqueSampleIds(cols, rows):
    """
    Add an 'id' column with a unique value for each row

    The same sampleName can occur once per library. (e.g. multiple PCR plates).
    Hence we can use "sampleName.libName" as unique id
    (if multiple libraries or libName is not specified Scale.RNA is used as the default libName).
    """
    sampleInd = cols.index("sample")
    libInd = cols.index("libName") if "libName" in cols else None
    for r in rows:
        if libInd is None:
            r.insert(0, f"{r[sampleInd]}.ScaleRNA")
        else:
            r.insert(0, f"{r[sampleInd]}.{r[libInd]}")
    cols.insert(0, "id")

def check_whether_barcode_ranges_overlap(all_samples_barcode_range: Dict, libraryStruct: Path):
    """
    Check whether the user supplied barcode ranges overlap amongst samples
    Throw exception if there is an overlap
    
    Args:
        all_samples_barcode_range: Dictionary of barcode ranges for all samples with libName as key
        libraryStruct: Path to the library structure json file
    """
    lib_struct_dir = libraryStruct.parent
    lib_struct = json.load(open(libraryStruct))
    # Get the txt file that corresponds to the sample_barcode in the lib json
    for barcode_info in lib_struct['barcodes']:
        if lib_struct['sample_barcode'] == barcode_info['name']:
            fname = barcode_info['sequences']
            break
    barcode_whitelist = pd.read_csv(lib_struct_dir / fname, sep="\t", names=['barcode', 'well'])
    for libName in all_samples_barcode_range:
        # Needs to be rest for each libName
        verbose_barcodes_list = []
        for sample_barcodes in all_samples_barcode_range[libName]:
            # String can be a single barcode or a range of barcodes separated by semi-colon
            for semi_colon_separated in sample_barcodes.split(";"):
                if "-" in semi_colon_separated:
                    # Get well coordinate that corresponds to starting of the barcodes range for a sample
                    starting = barcode_whitelist.index[barcode_whitelist['well']==semi_colon_separated.split("-")[0]][0]
                    # Get well coordinate that corresponds to end of the barcodes range for a sample
                    end = barcode_whitelist.index[barcode_whitelist['well']==semi_colon_separated.split("-")[1]][0]
                    # Retrieve all the well coordinates that correspond to the barcodes range for a sample
                    all_barcodes = barcode_whitelist.loc[starting:end, 'well'].tolist()
                    # extend because all_barcodes is a list
                    verbose_barcodes_list.extend(all_barcodes)
                else:
                    verbose_barcodes_list.append(semi_colon_separated)
        # Check whether the barcode ranges overlap amongst individual samples
        if len(verbose_barcodes_list) != len(set(verbose_barcodes_list)):
            print(f"The barcodes range mentioned for each sample overlap amongst individual samples", file=sys.stderr)
            sys.exit(1)

def main(samplesCsv:Path, splitFastq:bool, reporting:bool, libraryStruct:Path, singleLibrary:Optional[str], resultDir:Optional[str]):
    """
    Writes normalized samples.csv to stdout
    Args:
        samplesCsv: Path to samples.csv file
        reporting: Generate a samples.csv for a merging/reporting only workflow
        resultDir: Path or URL to previous workflow outputs
    """
    # Load column headers and all lines from the csv into lists
    rows: List[List[str]] = [] # Each row of the CSV, split by column
    # utf8-sig ignores 'BOM' characters at the beginning of the file, which some users might
    # accidentally add and which will break downstream parsing of samples.csv column headers.
    with open(samplesCsv, encoding='utf-8-sig') as csvfile:
        samples = csv.reader(csvfile)
        cols = next(samples) # CSV column header (in same order as values in rows)
        # Trim empty columns
        while not cols[-1]: 
            cols = cols[:-1]
        if 'sample' not in cols:
            print("'sample' column missing (the sample name)", file=sys.stderr)
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
            rows.append(row)

    # Process csv (add defaults, rename columns etc.)
    # Fix legacy column names
    if "libName" not in cols and "fastqName" in cols:
        cols[cols.index('fastqName')] = "libName"
    if "libIndex" not in cols and "index" in cols:
        cols[cols.index('index')] = "libIndex"
    if "libIndex" not in cols and "fastqIndex" in cols:
        cols[cols.index('fastqIndex')] = "libIndex"

    if "libIndex" in cols:
        libIndexInd = cols.index("libIndex")
        # Default libName based on index
        if "libName" not in cols:
            cols.insert(libIndexInd+1, "libName")
            for r in rows:
                name = r[libIndexInd].strip()
                if ';' in name:
                    raise ValueError("libName is required if libIndex lists multiple index sequences")
                r.insert(libIndexInd+1, name)
        #Normalize index sequences
        for r in rows:
            index = "".join(r[libIndexInd].split()) # Remove whitespace
            index = index.upper()
            index = ';'.join(sorted(index.split(';'))) # sort sequences
            r[libIndexInd] = index

    # Resolve sample names and IDs
    if cols[0] != "id":
        if "id" in cols:
            print("id has to be the first column", file=sys.stderr)
            sys.exit(1)
        else:
            makeUniqueSampleIds(cols, rows)
    ids = set()
    for r in rows:
        id = r[0]
        validateName(id)
        if id in ids:
            print(f"Duplicate sample name: {id}", file=sys.stderr)
            sys.exit(1)
        else:
            ids.add(id)
        validateName(r[cols.index("sample")])

    all_samples_barcode_range = {}

    if not reporting:
        if "libName" not in cols: # Only one (unnamed) library; Use default name
            cols.insert(2, "libName")
            for r in rows:
                r.insert(2, "ScaleRNA")
        libNameIndex = cols.index("libName")
        if singleLibrary: #Overwrite libName with single fixed value
            for r in rows:
                r[libNameIndex] = singleLibrary
        for r in rows:
            all_samples_barcode_range[r[libNameIndex]] = []
        if "barcodes" not in cols:
            cols.append("barcodes")
            for r in rows:
                r.append("1A-12H")
                all_samples_barcode_range[r[libNameIndex]].append("1A-12H")
        else: # check whether each entry for barcodes is empty or not
            barcodesIndex = cols.index("barcodes")
            for r in rows:
                if r[barcodesIndex].strip() == "":
                    r[barcodesIndex]="1A-12H"
                    all_samples_barcode_range[r[libNameIndex]].append("1A-12H")
                else:
                    r[barcodesIndex] = r[barcodesIndex].replace(" ", "")
                    all_samples_barcode_range[r[libNameIndex]].append(r[barcodesIndex])
        check_whether_barcode_ranges_overlap(all_samples_barcode_range, libraryStruct)
        if splitFastq:
            if "split" not in cols:
                cols.append("split")
                for r in rows:
                    r.append("true")
    # "group" is used for merging samples with --merge; default is by "sample" (e.g. different libraries of the same sample)
    if "group" not in cols:
        nameInd = cols.index("sample")
        cols.append("group")
        for r in rows:
            r.append(r[nameInd])
    if resultDir and "resultDir" not in cols: # Set a single resultDir for all samples
        cols.append("resultDir")
        for r in rows:
            r.append(resultDir)
    w = csv.writer(sys.stdout)
    w.writerow(cols)
    w.writerows(rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Standardize samples.csv for the workflow (defaults, column names, etc.')
    parser.add_argument('samples', metavar='SAMPLES.csv', help='CSV with sample names and information for the workflow run')
    parser.add_argument('--splitSample', help='Flag that denotes whether bcParser will split per RT well to ensure parallelization', action="store_true", default=False)
    parser.add_argument('--reporting', help='set for use in merging/reporting only workflow', action="store_true", default=False)
    parser.add_argument('--resultDir', help="Previous pipeline output directory (used as default for 'resultDir column)")
    parser.add_argument('--libraryStruct', help="Library structure json file", type=Path)
    parser.add_argument('--singleLibrary', help="Set libName for all samples to this fixed value")
    args = parser.parse_args()
    main(args.samples, args.splitSample, args.reporting, args.libraryStruct, args.singleLibrary, args.resultDir)
