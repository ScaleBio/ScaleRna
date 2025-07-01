#!/usr/bin/env python
"""
Detect library pool from index2 fastq files or unaligned ultima cram files

pool_mapping.csv will be created in the current directory with the following columns:
    file: The input file name
    pool: The detected pool name or "Unknown" if no match is found or "Ambiguous" if multiple pools are detected
"""
import argparse
from pathlib import Path
from collections import Counter
import csv
import gzip
import json
import sys
import pandas as pd
import pysam

NUM_READS = 1000
THRESHOLD = 0.1  # % of reads that must match known pool
COMPLEMENT = str.maketrans(dict(zip("ACGT", "TGCA")))


def load_pool_whitelist(lib_json: dict[str, dict]):
    """
    Find the PCR whitelist from the library structure JSON file
    and return a dictionary of the whitelist mapping alias to sequences
    """
    pool_mapping = {}
    for lib in lib_json:
        pcr_pool_whitelist = None
        for bc in lib_json[lib]["barcodes"]:
            if bc.get("name", "").upper() == "PCR" or bc.get("alias", "").upper() == "PCR":
                pcr_pool_whitelist = lib.parent / bc["sequences"]
                substr = slice(bc.get("start", 0), bc.get("start", 0) + bc.get("length", 0))

        if pcr_pool_whitelist is None:
            raise ValueError("PCR whitelist not found in library structure")

        with open(pcr_pool_whitelist) as f:
            for line in f:
                split_line = line.strip().split("\t")
                # key is the alias, value is a tuple whose elements are a list of sequences and the substr of read
                pool_mapping[split_line[0]] = (split_line[1:], substr)
    return pool_mapping


def get_barcode_entry_idx_in_lib_json(lib_json: dict[str, dict]):
    """
    Get the index of the PCR and RT cell barcodes in the library structure JSON file
    """
    pcr_idx = []
    rt_idx = []
    for lib in lib_json:
        sample_barcode = lib_json[lib]["sample_barcode"]
        for idx, bc in enumerate(lib_json[lib]["barcodes"]):
            if bc.get("name", "").upper() == "PCR" or bc.get("alias", "").upper() == "PCR":
                pcr_idx.append(idx)
            if bc.get("name", "").upper() == sample_barcode or bc.get("alias", "").upper() == sample_barcode:
                rt_idx.append(idx)

    if all(x == pcr_idx[0] for x in pcr_idx) and all(x == rt_idx[0] for x in rt_idx):
        return pcr_idx[0], rt_idx[0]
    else:
        raise ValueError(
            "The position of the PCR/RT cell barcode is not the same in the RNA lib json vs the ScalePlex lib json"
        )


def parse_fastq(fastq_file: Path):
    with gzip.open(fastq_file, "rt") as f:
        while True:
            # Read the header
            if not f.readline():  # End of file
                break
            seq = f.readline().strip()
            f.readline()  # The '+' line
            f.readline()  # The quality line

            yield seq


def get_hamming_distance(observed_seq: str, candidate: str) -> int:
    """Compute hamming distance between two sequences

    Args:
        observed_seq: index read
        candidate: sequence from the whitelist

    Returns:
        Hamming distance, representative of the number of mismatches
    """
    return sum(1 for a, b in zip(observed_seq, candidate) if a != b)


def match_seq_to_get_orientation(observed_seq: str, whitelist: dict) -> bool:
    """Check for matches in the whitelist while prioritizing exact matches to get the orientation

    Args:
        observed_seq: index read
        whitelist: alias as key and list of sequences as first element of value, substr as second element of value

    Returns:
        Boolean indicating if the sequence is reverse complemented
    """
    dist_fwd = []
    dist_rev = []
    for _alias, seqs in whitelist.items():
        for candidate in seqs[0]:
            dist_fwd.append(get_hamming_distance(observed_seq[seqs[1]], candidate))
            dist_rev.append(get_hamming_distance(observed_seq[seqs[1]], candidate[::-1].translate(COMPLEMENT)))

    # 0 indicates exact match, 1 indicates mismatch
    if dist_fwd.count(0) > 1 or dist_rev.count(0) > 1:
        raise ValueError(f"{observed_seq} has multiple exact matches in whitelist")
    if dist_fwd.count(0) == 1 and dist_rev.count(0) == 1:
        print(f"{observed_seq} has exact matches in both orientations in whitelist", file=sys.stderr)
        return None
    if dist_fwd.count(1) >= 1 and dist_rev.count(1) >= 1:
        print(f"{observed_seq} has multiple mismatched matches in both orientations in whitelist", file=sys.stderr)
        return None

    if any([i <= 1 for i in dist_fwd]):
        return False
    elif any([i <= 1 for i in dist_rev]):
        return True
    else:
        return None


def match_seq(observed_seq: str, whitelist: dict, is_rc: bool) -> str:
    """Check for matches in the whitelist

    Args:
        observed_seq: index read
        whitelist: alias as key and list of sequences as first element of value, substr as second element of value
        is_rc: orientation of the fastq file

    Returns:
        String containing the alias

    """
    mismatched_alias = []
    for alias, seqs in whitelist.items():
        for candidate in seqs[0]:
            # Match according to orientation
            if is_rc:
                candidate = candidate[::-1].translate(COMPLEMENT)
            hamming_distance = get_hamming_distance(observed_seq[seqs[1]], candidate)
            if hamming_distance == 0:
                return alias
            elif hamming_distance == 1:
                mismatched_alias.append(alias)
    if len(mismatched_alias) == 1:
        return mismatched_alias[0]

    return None


def sniff_fastq(fastq_dir: Path, whitelist: dict, allow_ambiguous: bool):
    """
    Detect library pool from fastq files

    This function reads the first 1000 reads from each fastq file and detects the library pool.
    It is assumed the correct fastq file (index2) is present in the directory.

    Args:
        fastq_dir: Directory containing fastq files
        whitelist: Dictionary of pool whitelist
        allow_ambiguous: If True, don't error if fastq has multiple pools
    """
    fq_to_pool_mapping = []
    for fq_file in fastq_dir.glob("*.fastq.gz"):
        rc_info = []
        is_rc = None
        for i, seq in enumerate(parse_fastq(fq_file)):
            if i >= NUM_READS:
                break
            # Get orientation of the sequence
            rc_info.append(match_seq_to_get_orientation(seq, whitelist))
        # Indicates no matches
        if all(is_rc is None for is_rc in rc_info):
            aliases = [None]
        else:
            aliases = []
            rc_counts = Counter(rc_info)
            # Determine orientation of the fastq file
            is_rc = rc_counts[True] > rc_counts[False]
            # Match reads again with determined orientation
            for i, seq in enumerate(parse_fastq(fq_file)):
                if i >= NUM_READS:
                    break
                aliases.append(match_seq(seq, whitelist, is_rc))
        match = pd.Series(aliases).value_counts(normalize=True)
        fq_to_pool_mapping.append(return_mapping_entry(match, fq_file, is_rc))

    write_to_csv(fq_to_pool_mapping, "pool_mapping.csv")
    raise_ambiguous_error(fq_to_pool_mapping, allow_ambiguous)


def write_to_csv(mapping: list, filename: str, header: list = ["file", "pool", "is_rc"]):
    """Write the mapping information to a CSV file

    Args:
        mapping: List of mappings
        filename: Name of the output CSV file
        header: Header for the CSV file
    """
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(mapping)


def return_mapping_entry(match: pd.Series, fname: Path, is_rc: bool):
    """Return the mapping entry for a given match

    Args:
        match: Normalized match series
        fname: File name
        is_rc: Orientation of the fastq file

    Returns:
        The mapping entry for the given match
    """
    if len(match) == 0:
        return [fname, "Unknown", None]
    if len(match) > 1:
        return [fname, "Ambiguous", None]
    if len(match) == 1:
        if match.iloc[0] >= THRESHOLD:
            return [fname, match.index[0], is_rc]
        else:
            return [fname, "Unknown", None]


def raise_ambiguous_error(mapping: list, allow_ambiguous: bool, ultima: bool = False):
    """Raise an error if ambiguous pools are detected

    Args:
        mapping: List of mappings
        allow_ambiguous: If True, don't error if fastq has multiple pools
    """
    if not allow_ambiguous:
        ambiguous = [
            str(x[0]) for x in mapping if x[1] == "Ambiguous" or (ultima and x[2] == "Ambiguous")
        ]  # check whether rt is ambiguous
        if ambiguous:
            raise ValueError(f"Ambiguous pools detected in the following files: {ambiguous}")


def sniff_cram(ultima_cram_dir: Path, pcr_idx: int, rt_idx: int, allow_ambiguous: bool, whitelist: dict):
    """
    Read 1000 entries and detect pcr barcode alias and rt barcode alias from cram files

    Args:
        ultima_cram_dir: Directory containing cram files
        pcr_idx: The entry number of the PCR cell barcode in the lib json.
            Corresponds to the index of the cell barcode in the read
        rt_idx: The entry number of the RT cell barcode in the lib json.
            Corresponds to the index of the cell barcode in the read
        allow_ambiguous: If True, don't error if cram has multiple pools
        whitelist: Dictionary of pool whitelist
    """
    cram_to_pool_mapping = []
    for cram in ultima_cram_dir.glob("*.cram"):
        pcr_match = get_cram_match_info_for_cb(cram, pcr_idx, whitelist, "PCR")
        rt_match = get_cram_match_info_for_cb(cram, rt_idx, whitelist, "RT")
        pcr_entry = return_mapping_entry(pcr_match, cram, None)
        rt_entry = return_mapping_entry(rt_match, cram, None)
        # fname, pcr_pool, rt, is_rc
        entry = [pcr_entry[0], pcr_entry[1], rt_entry[1]]
        cram_to_pool_mapping.append(entry)

    write_to_csv(cram_to_pool_mapping, "pool_mapping.csv", ["file", "pcr_pool", "rt"])
    raise_ambiguous_error(cram_to_pool_mapping, allow_ambiguous, True)


def match_cram_pcr_alias(cb_alias: str, whitelist: dict) -> str:
    """Match the PCR alias to the whitelist

    Args:
        cb_alias: The cell barcode alias
        whitelist: Dictionary of pool whitelist

    Returns:
        The matched alias or None
    """
    for alias in whitelist:
        if cb_alias == alias:
            return alias
    return None


def get_cram_match_info_for_cb(
    cram: Path, idx: int, whitelist: dict, barcode_alias: str, tag: str = "CB", separator: str = "+"
):
    """Normalize counts for barcode matches in a cram file"""
    aliases = []
    with pysam.AlignmentFile(cram, "rb", check_sq=False) as cram_file:
        for i, read in enumerate(cram_file.fetch(until_eof=True)):
            if i >= NUM_READS:
                break
            try:
                cb = read.get_tag(tag).split(separator)
                if barcode_alias == "PCR":
                    aliases.append(match_cram_pcr_alias(cb[idx], whitelist))  # need to match to whitelist
                else:
                    aliases.append(cb[idx])  # for rt, no matching, so just append the alias
            except KeyError:
                aliases.append(None)  # CB tag not found
    return pd.Series(aliases).value_counts(normalize=True)


def main():
    parser = argparse.ArgumentParser("Read first 1000 reads of fastq files and detect library pool")
    parser.add_argument("--fastqDir", type=Path, help="Directory containing fastq files")
    parser.add_argument("--ultimaCramDir", type=Path, help="Directory containing Ultima unaligned CRAM files")
    parser.add_argument("--libraryStruct", required=True, type=Path)
    parser.add_argument("--scalePlexLibraryStruct", required=False, type=Path)
    parser.add_argument("--allowAmbiguous", action="store_true", help="Don't error if fastq has multiple pools")
    args = parser.parse_args()

    lib_json = {args.libraryStruct: json.load(open(args.libraryStruct))}
    if args.scalePlexLibraryStruct is not None:
        lib_json[args.scalePlexLibraryStruct] = json.load(open(args.scalePlexLibraryStruct))
    whitelist = load_pool_whitelist(lib_json)
    if args.fastqDir:
        sniff_fastq(
            fastq_dir=args.fastqDir,
            whitelist=whitelist,
            allow_ambiguous=args.allowAmbiguous,
        )
    elif args.ultimaCramDir:
        pcr_idx, rt_idx = get_barcode_entry_idx_in_lib_json(lib_json)
        sniff_cram(
            ultima_cram_dir=args.ultimaCramDir,
            pcr_idx=pcr_idx,
            rt_idx=rt_idx,
            allow_ambiguous=args.allowAmbiguous,
            whitelist=whitelist,
        )
    else:
        raise ValueError("Either --fastqDir or --ultimaCramDir must be provided")


if __name__ == "__main__":
    main()
