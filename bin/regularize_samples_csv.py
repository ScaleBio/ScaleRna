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
import pandas as pd
from pathlib import Path
from typing import Optional, Dict
from scale_utils.validation import validateName
from scale_utils.lib_json_parser import LibJsonParser


COL_NAMES_WHITELIST = [
    "id",
    "sample",
    "barcodes",
    "libName",
    "libIndex",
    "libIndex2",
    "scalePlexLibName",
    "scalePlexLibIndex",
    "scalePlexLibIndex2",
    "scalePlexBarcodes",
    "split",
    "group",
    "resultDir",
    "rnaId",
    "expectedCells"
]


def drop_hash_cols(cols, rows):
    """
    If scalePlexLibName and scalePlexLibIndex are present, drop them
    """
    for hash_col in ["scalePlexLibName", "scalePlexLibIndex", "scalePlexLibIndex2"]:
        if hash_col in cols:
            for r in rows:
                r.pop(cols.index(hash_col))
            cols.remove(hash_col)


def makeUniqueSampleIds(cols, rows, quantum):
    """
    Add an 'id' column with a unique value for each row

    The same sampleName can occur more than once per library. (e.g. multiple PCR plates).
    Hence we can use "sampleName.libName" as unique id
    """
    for r in rows:
        validateName(r[cols.index("sample")], display_name="sample")
        validateName(r[cols.index("libName")], display_name="libName", other_chars="-")
        if quantum and "libIndex2" in cols:
            r.insert(0, f"{r[cols.index('sample')]}.{r[cols.index('libIndex2')]}")
        else:
            r.insert(0, f"{r[cols.index('sample')]}.{r[cols.index('libName')]}")
    # add column name at end so indices are correct in line above
    cols.insert(0, "id")


def fail(msg: str, sep=" "):
    """
    Helper function to print an error message and exit

    Args:
        msg: Error message
        sep: Separator to use when printing the error message
    """
    print(msg, file=sys.stderr, sep=sep)
    # Custom exit code to prevent retries
    # Used in nextflow config
    sys.exit(42)


def unwrapSemiColonSeparatedValues(cols: list, rows: list, index: int, scale_plex_col_name: str) -> list:
    """
    Unwrap each entry separated by a semi-colon into its own row

    Args:
        cols: Column headers
        rows: List of rows from the csv
        index: Index of the column to unwrap
        scale_plex_col_name: Name of the scalePlex column

    Returns:
        final_rows: List of rows with unwrapped values
    """
    final_rows = []
    for r in rows:
        # Unwrap each entry separated by a semi-colon into its own row
        if ";" in r[index]:
            split_libs = r[index].split(";")
            split_libs = [libname.strip() for libname in split_libs]
            split_scaleplexLibs = []
            if scale_plex_col_name in cols:
                split_scaleplexLibs = r[cols.index(scale_plex_col_name)].split(";")
                split_scaleplexLibs = [libname.strip() for libname in split_scaleplexLibs]
                if len(split_libs) != len(split_scaleplexLibs):
                    fail(
                        f"Number of scalePlex entries {r[cols.index(scale_plex_col_name)]} "
                        f"must match number of rna entries {r[index]}"
                    )
            for idx, split_libname in enumerate(split_libs):
                unwrapped_row = r[:]
                unwrapped_row[index] = split_libname
                if split_scaleplexLibs:
                    unwrapped_row[cols.index(scale_plex_col_name)] = split_scaleplexLibs[idx]
                final_rows.append(unwrapped_row)
        else:
            final_rows.append(r)
    return final_rows


def readSamplesCsv(samplesCsv: Path, quantum: bool, fastq_input: bool, reporting: bool) -> tuple[list, list[list[str]]]:
    """
    Load samples.csv into a list.

    Perform basic validation and cleanup.
    Args:
        samplesCsv: Path to samples.csv file
        quantum: Flag to indicate analysis is for data generated on a ScaleQuantum kit
        fastq_input: Flag to indicate that fastq files will be used as input to the workflow
        reporting: Flag to indicate this is a reporting run

    Returns:
        cols: Column headers
        rows: List of rows from the csv
    """
    rows = []
    # Create a mapping of lowercase column names to their original names
    col_names_mapping = {col.lower(): col for col in COL_NAMES_WHITELIST}
    # utf8-sig ignores 'BOM' characters at the beginning of the file, which some users might
    # accidentally add and which will break downstream parsing of samples.csv column headers.
    with open(samplesCsv, encoding="utf-8-sig") as csvfile:
        samples = csv.reader(csvfile)
        cols = next(samples)  # CSV column header (in same order as values in rows)
        # Normalize column names
        cols = [
            col_names_mapping[col.lower()] if col.lower() in col_names_mapping else fail(f"Invalid column name {col}")
            for col in cols
        ]
        # Trim empty columns
        while not cols[-1]:
            cols = cols[:-1]
        if len(cols) != len(set(cols)):
            fail("Duplicate column names, check samplesheet")
        if "sample" not in cols:
            fail("'sample' column missing (the sample name)")
        for row in samples:
            if not row or len(row[0].strip()) == 0:
                continue
            # Remove leading and trailing spaces from all values in the row
            row = [row_element.strip() for row_element in row]
            # Trim empty columns
            while not row[-1] and len(row) > len(cols):
                row = row[:-1]

            if len(row) != len(cols):
                fail(f"Unexpected number of columns: {row}", sep="\n")
            rows.append(row)
    if "libName" not in cols and "libIndex" in cols:
        for r in rows:
            if ";" in r[cols.index("libIndex")]:
                fail("libName is required if libIndex lists multiple index sequences")
    if "libName" in cols and "libIndex" in cols:
        for r in rows:
            if ";" in r[cols.index("libIndex")] and ";" in r[cols.index("libName")]:
                fail("Cannot list multiple libIndex if listing multiple libName")
    if quantum and not reporting:
        if "libIndex2" in cols and "libName" in cols:
            fail("Please provide a value for only the libIndex2 column or the libName column, not both")
        if not fastq_input:
            if "libName" in cols:
                fail("When starting from bcl files, do not provide a value for the libName column")
    # Normalize index sequences
    for idx_col in ["libIndex", "libIndex2", "scalePlexLibIndex", "scalePlexLibIndex2"]:
        if idx_col in cols:
            for r in rows:
                index = "".join(r[cols.index(idx_col)].split())  # Remove whitespace
                index = index.upper()
                index = ";".join(sorted(index.split(";")))  # sort sequences
                r[cols.index(idx_col)] = index
    final_rows = []
    if quantum and "libName" in cols:
        final_rows = unwrapSemiColonSeparatedValues(cols, rows, cols.index("libName"), "scalePlexLibName")
    else:
        final_rows = rows
    if quantum and "libIndex2" in cols:
        final_rows = unwrapSemiColonSeparatedValues(cols, final_rows, cols.index("libIndex2"), "scalePlexLibIndex2")
    else:
        final_rows = final_rows

    return cols, final_rows


def check_whether_barcode_ranges_overlap(all_samples_barcode_range: Dict, sample_barcode_fname: Path):
    """
    Check whether the user supplied barcode ranges overlap amongst samples
    Throw exception if there is an overlap

    Args:
        all_samples_barcode_range: Dictionary of barcode ranges for all samples with libName as key
        sample_barcode_fname: Path to sample barcode whitelist file
    """
    barcode_whitelist = pd.read_csv(sample_barcode_fname, sep="\t", names=["well", "barcode"])
    for libName in all_samples_barcode_range:
        # Needs to be rest for each libName
        verbose_barcodes_list = []
        for sample_barcodes in all_samples_barcode_range[libName]:
            # String can be a single barcode or a range of barcodes separated by semi-colon
            for semi_colon_separated in sample_barcodes.split(";"):
                if "-" in semi_colon_separated:
                    # Get well coordinate that corresponds to starting of the barcodes range for a sample
                    starting = barcode_whitelist.index[barcode_whitelist["well"] == semi_colon_separated.split("-")[0]][
                        0
                    ]
                    # Get well coordinate that corresponds to end of the barcodes range for a sample
                    end = barcode_whitelist.index[barcode_whitelist["well"] == semi_colon_separated.split("-")[1]][0]
                    # Retrieve all the well coordinates that correspond to the barcodes range for a sample
                    all_barcodes = barcode_whitelist.loc[starting:end, "well"].tolist()
                    # extend because all_barcodes is a list
                    verbose_barcodes_list.extend(all_barcodes)
                else:
                    verbose_barcodes_list.append(semi_colon_separated)
        # Check whether the barcode ranges overlap amongst individual samples
        if len(verbose_barcodes_list) != len(set(verbose_barcodes_list)):
            fail("The barcodes range mentioned for each sample overlap amongst individual samples")


def get_first_and_last_entry(file: Path, sep: str = "\t") -> tuple[str, str]:
    """
    Get the first and last entry of a file
    Args:
        file: Path to file
        sep: Separator used in the file
    Returns:
        (first_entry, last_entry): Tuple of first and last entry
    """
    with open(file) as f:
        lines = f.readlines()
        first_entry = lines[0].strip().split(sep)[0]
        last_entry = lines[-1].strip().split(sep)[0]

    return first_entry, last_entry


def write_unwrapped_barcodes_range_samples_csv(cols: list, rows: list, sample_barcode_fname: Path):
    """
    Write a csv with a singular well coordinate for the barcodes column on each line with same id,
    sample name and library name
    Relevant for analyzing ultima data where each well coordinate maps to a single file

    Args:
        cols: Column headers
        rows: Rows from the csv
        sample_barcode_fname: Path to sample barcode whitelist file
    """
    unwrapped_rows = []
    barcode_whitelist_wells = pd.read_csv(sample_barcode_fname, sep="\t", names=["well", "barcode"])["well"].to_list()
    for r in rows:
        unwrapped_rows.append(
            unwrap_barcodes_range(
                r[cols.index("barcodes")], r[cols.index("id")], r[cols.index("libName")], barcode_whitelist_wells
            )
        )
    flattened_unwrapped_rows = [
        item for sublist in unwrapped_rows for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    with open("barcode_range_unwrapped_samples.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["id", "libName", "barcode"])
        writer.writerows(flattened_unwrapped_rows)


def get_all_barcodes_in_range(barcode_range: str, barcode_whitelist_wells):
    """
    Get all barcodes in the range defined by @barcode_range

    Args:
        barcode_range: Range of barcodes
        barcode_whitelist_wells: All barcode well coordinates in whitelist file

    Returns:
        List of all barcodes in the range
    """
    starting_barcode = barcode_range.split("-")[0]
    ending_barcode = barcode_range.split("-")[1]
    starting_index = barcode_whitelist_wells.index(starting_barcode)
    ending_index = barcode_whitelist_wells.index(ending_barcode)
    return barcode_whitelist_wells[starting_index : ending_index + 1]


def unwrap_barcodes_range(
    barcodes_range: str, sample_id: str, libName: str, barcode_whitelist_wells: list
) -> list[tuple[str, str]]:
    """
    Unwrap the barcodes range into individual barcodes

    Args:
        barcodes_range: Range of barcodes
        sample_id: Sample id
        libName: Library name
        barcode_whitelist_wells: All barcode well coordinates in whitelist file

    Returns:
        List of tuples containing sample_id, libName, barcode
    """
    if ";" in barcodes_range:
        unwrapped_barcodes = []
        split_barcodes = barcodes_range.split(";")
        for split_barcode in split_barcodes:
            if "-" in split_barcode:
                unwrapped_barcodes.extend(get_all_barcodes_in_range(split_barcode, barcode_whitelist_wells))
            else:
                unwrapped_barcodes.append(split_barcode)
    elif "-" in barcodes_range:
        unwrapped_barcodes = get_all_barcodes_in_range(barcodes_range, barcode_whitelist_wells)
    else:
        return (sample_id, libName, barcodes_range)
    return [(sample_id, libName, barcode) for barcode in unwrapped_barcodes]


def add_all_alias_for_index(
    cols: list[str],
    rows: list[list[str]],
    all_aliases: list[str],
    col_name: str,
    scaleplex_col_name: str,
    scaleplex_to_rna_mapping: dict[str:str],
) -> list[list[str]]:
    """
    For a particular barcode (corresponds to @read), add all aliases from the whitelist file for each sample

    Args:
        cols: Column headers
        rows: Rows from the csv
        all_aliases: All aliases from whitelist file for the given read
        col_name: Column name to add aliases for
        libraryStruct: Path to library structure json
        scaleplex_col_name: Column name to add scalePlex aliases for
        scaleplex_to_rna_mapping: RNA PCR to ScalePlex PCR mapping

    Returns:
        List of rows with all aliases added
    """
    rows_with_all_aliases = []

    for r in rows:
        for alias in all_aliases:
            r_copy = r[:]
            r_copy.insert(cols.index(col_name), alias)
            if scaleplex_col_name:
                scaleplex_alias = rna_to_scaleplex_mapping.get(alias)
                if scaleplex_alias:
                    r_copy.insert(cols.index(scaleplex_col_name), scaleplex_alias)
                else:
                    fail(f"ScalePlex alias not found for {alias} in {scaleplex_to_rna_mapping}")
            rows_with_all_aliases.append(r_copy)

    return rows_with_all_aliases


def main(
    samplesCsv: Path,
    splitFastq: bool,
    scalePlex: bool,
    reporting: bool,
    libraryStruct: Path,
    rna_to_scaleplex_mapping: dict[str:str],
    quantum: bool,
    resultDir: Optional[str],
    ultima: bool,
    samples_csv_fname: str,
    fastq_input: str,
):
    """
    Writes normalized samples.csv to stdout
    Args:
        samplesCsv: Path to samples.csv file
        splitFastq: Split fastq files per RT well
        scalePlex: Add matched hash libraries for each sample to samples.csv
        reporting: Generate a samples.csv for a merging/reporting only workflow
        quantum: Analysis is for data generated on a ScaleQuantum kit
        resultDir: Path or URL to previous workflow outputs
        libraryStruct: Path to library structure json
        rna_to_scaleplex_mapping: RNA PCR to ScalePlex PCR mapping
        ultima: Analysis is for data generated on an ultima sequencer
        samples_csv_fname: Name of the samples.csv file
        fastq_input: Flag to indicate that fastq files will be used as input to the workflow
    """
    lib_json_obj = LibJsonParser(libraryStruct)
    sample_barcode_whitelist_first_entry, sample_barcode_whitelist_second_entry = lib_json_obj.get_barcodes_range_of_whitelist_file(lib_json_obj.sample_barcode_fname)
    barcodes_range = f"{sample_barcode_whitelist_first_entry}-{sample_barcode_whitelist_second_entry}"
    # Load column headers and all lines from the csv into lists
    cols, rows = readSamplesCsv(samplesCsv, quantum, fastq_input, reporting)
    # Process csv (add defaults, rename columns etc.)
    # Default libName based on index
    if "libName" not in cols:
        if quantum:
            if "libIndex2" not in cols:
                cols.insert(len(cols), "libIndex2")
                if scalePlex:
                    if "scalePlexLibIndex2" not in cols:
                        cols.insert(len(cols), "scalePlexLibIndex2")
                    else:
                        fail("Can not specify scalePlexLibIndex2 without libIndex2 in samples.csv")
                rows = add_all_alias_for_index(
                    cols,
                    rows,
                    lib_json_obj.all_whitelist_contents_by_read["index2"],
                    "libIndex2",
                    "scalePlexLibIndex2" if scalePlex else None,
                    rna_to_scaleplex_mapping,
                )
                if scalePlex:
                    # Default scalePlexLibName based on scalePlexLibIndex2 value added in add_all_alias_for_index
                    if "scalePlexLibName" in cols:
                        fail("Can not specify scalePlexLibName without libName in samples.csv")
                    idx = cols.index("scalePlexLibIndex2")
                    for r in rows:
                        r.insert(idx, r[cols.index("scalePlexLibIndex2")])
                    cols.insert(idx, "scalePlexLibName")
        # Add libName column after libIndex or at end of cols
        idx = cols.index("libIndex") if "libIndex" in cols else len(cols)
        for r in rows:
            if "libIndex" in cols:
                name = r[cols.index("libIndex")]
            elif "libIndex2" in cols and quantum:
                name = r[cols.index("libIndex2")]
            else:
                name = "ScaleRNA"
            r.insert(idx, name)
        # insert column name at idx after row updates
        cols.insert(idx, "libName")

    # Resolve sample names and IDs
    if cols[0] != "id":
        if "id" in cols:
            fail("id has to be the first column")
        else:
            makeUniqueSampleIds(cols, rows, quantum)
    ids = set()
    for r in rows:
        id = r[0]
        if id in ids:
            fail(f"Duplicate sample name: {id}")
        else:
            ids.add(id)

    all_samples_barcode_range = {}

    if not reporting:
        libNameIndex = cols.index("libName")
        for r in rows:
            all_samples_barcode_range[r[libNameIndex]] = []
        if "barcodes" not in cols:
            cols.append("barcodes")
            for r in rows:
                r.append(barcodes_range)
                all_samples_barcode_range[r[libNameIndex]].append(barcodes_range)
        else:  # check whether each entry for barcodes is empty or not
            barcodesIndex = cols.index("barcodes")
            for r in rows:
                if r[barcodesIndex].strip() == "":
                    r[barcodesIndex] = barcodes_range
                    all_samples_barcode_range[r[libNameIndex]].append(barcodes_range)
                else:
                    r[barcodesIndex] = r[barcodesIndex].replace(" ", "")
                    all_samples_barcode_range[r[libNameIndex]].append(r[barcodesIndex])
        check_whether_barcode_ranges_overlap(all_samples_barcode_range, lib_json_obj.sample_barcode_fname)
        if splitFastq:
            if "split" not in cols:
                cols.append("split")
                for r in rows:
                    r.append("true")
    # "group" is used for merging samples with --merge; default is by "sample"
    # (e.g. different libraries of the same sample)
    if "group" not in cols:
        nameInd = cols.index("sample")
        cols.append("group")
        for r in rows:
            r.append(r[nameInd])
    if resultDir and "resultDir" not in cols:  # Set a single resultDir for all samples
        cols.append("resultDir")
        for r in rows:
            r.append(resultDir)
    # Add hash libraries for scalePlex
    if scalePlex:
        hash_rows = []
        # Store a reference to the RNA sample-lib combination
        cols.append("rnaId")
        for r in rows:
            # For each original row add a corresponding hash library row to samples.csv
            hash_row = r[:]
            if "scalePlexLibName" in cols:
                hash_lib_name = r[cols.index("scalePlexLibName")]
                if hash_lib_name == r[cols.index("libName")]:
                    fail(f"ScalePlex library name must be distinct from RNA library name: {hash_lib_name}")
                hash_row[cols.index("libName")] = hash_lib_name
            else:
                if quantum and "libIndex2" in cols:
                    hash_row[cols.index("libName")] = rna_to_scaleplex_mapping[r[cols.index("libIndex2")]]
                else:
                    hash_row[cols.index("libName")] = f"{r[cols.index('libName')]}-ScalePlex"
            if "scalePlexLibIndex" in cols or "scalePlexLibIndex2" in cols:
                if quantum:
                    hash_row[cols.index("libIndex2")] = r[cols.index("scalePlexLibIndex2")]
                else:
                    hash_row[cols.index("libIndex")] = r[cols.index("scalePlexLibIndex")]
            else:
                if quantum and "libIndex2" in cols:
                    hash_row[cols.index("libIndex2")] = rna_to_scaleplex_mapping[r[cols.index("libIndex2")]]
                if "libIndex" in cols:
                    libIndex = r[cols.index("libIndex")]
                    if all(char in "ACGTacgt;" for char in libIndex):
                        # Custom libIndex specified for RNA samples
                        # Set to empty string so that hash index1 seqs will be used
                        hash_row[cols.index("libIndex")] = ""
                    else:
                        # set to corresponding hash library named index seq
                        hash_row[cols.index("libIndex")] = libIndex.replace("RNA", "ScalePlex")
            # Add rnaId
            hash_row.append(r[cols.index("id")])
            # drop id, new one will be created
            hash_row = hash_row[1:]
            hash_rows.append(hash_row)
        makeUniqueSampleIds(cols[1:], hash_rows, quantum)
        rows.extend(hash_rows)
        drop_hash_cols(cols, rows)
    if ultima:
        write_unwrapped_barcodes_range_samples_csv(cols, rows, lib_json_obj.sample_barcode_fname)
    with open(samples_csv_fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Standardize samples.csv for the workflow (defaults, column names, etc."
    )
    parser.add_argument(
        "samples", metavar="SAMPLES.csv", help="CSV with sample names and information for the workflow run"
    )
    parser.add_argument(
        "--splitSample",
        help="Flag that denotes whether bcParser will split per RT well to ensure parallelization",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--reporting", help="set for use in merging/reporting only workflow", action="store_true", default=False
    )
    parser.add_argument(
        "--resultDir", help="Previous pipeline output directory (used as default for 'resultDir column)"
    )
    parser.add_argument("--libraryStruct", help="Library structure json file", type=Path)
    parser.add_argument("--scalePlexToRnaMapping", type=Path, help="Path to scaleplex to rna mapping file")
    parser.add_argument(
        "--scalePlex", help="Flag to indicate samples have matched hash libraries", action="store_true", default=False
    )
    parser.add_argument(
        "--quantum",
        help="Flag to indicate analysis is for data generated on a ScaleQuantum kit",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--ultima",
        help="Flag to indicate analysis is for data generated on an ultima sequencer",
        action="store_true",
        default=False,
    )
    parser.add_argument("--samples_csv_fname", help="Name of the samples.csv file", type=str, default="samples.csv")
    parser.add_argument(
        "--fastq",
        help="Flag to indicate that fastq files will be used as input to the workflow",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    rna_to_scaleplex_mapping = {}
    if args.scalePlexToRnaMapping:
        df = pd.read_csv(args.scalePlexToRnaMapping, sep="\t", header=None)
        rna_to_scaleplex_mapping = dict(zip(df[1], df[0]))

    main(
        samplesCsv=args.samples,
        splitFastq=args.splitSample,
        scalePlex=args.scalePlex,
        reporting=args.reporting,
        libraryStruct=args.libraryStruct,
        rna_to_scaleplex_mapping=rna_to_scaleplex_mapping,
        quantum=args.quantum,
        resultDir=args.resultDir,
        ultima=args.ultima,
        samples_csv_fname=args.samples_csv_fname,
        fastq_input=args.fastq,
    )
