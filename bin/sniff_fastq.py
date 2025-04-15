#!/usr/bin/env python
"""
Detect library pool from index2 fastq files

fastq_to_pool_mapping.csv will be created in the current directory with the following columns:
    fastq_file: The fastq file name
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

NUM_READS = 1000
THRESHOLD = 0.1  # % of reads that must match known pool
COMPLEMENT = str.maketrans(dict(zip("ACGT", "TGCA")))


def load_pool_whitelist(lib_struct: Path, scale_plex_lib_struct: Path):
    """
    Find the PCR whitelist from the library structure JSON file
    and return a dictionary of the whitelist mapping alias to sequences
    """
    pool_mapping = {}
    lib_json = [json.load(open(lib_struct))]
    if scale_plex_lib_struct is not None:
        lib_json.append(json.load(open(scale_plex_lib_struct)))
    for lib in lib_json:
        pcr_pool_whitelist = None
        for bc in lib["barcodes"]:
            if bc.get("name", "").upper() == "PCR" or bc.get("alias", "").upper() == "PCR":
                pcr_pool_whitelist = Path(lib_struct.parent / bc["sequences"])
                substr = slice(bc["start"], bc["start"] + bc["length"])
                break

        if pcr_pool_whitelist is None:
            raise ValueError("PCR whitelist not found in library structure")

        with open(pcr_pool_whitelist) as f:
            for line in f:
                split_line = line.strip().split("\t")
                # key is the alias, value is a tuple whose elements are a list of sequences and the substr of read
                pool_mapping[split_line[0]] = (split_line[1:], substr)

    return pool_mapping


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
        if len(match) == 0:
            fq_to_pool_mapping.append([fq_file, "Unknown", None])
        if len(match) > 1:
            fq_to_pool_mapping.append([fq_file, "Ambiguous", None])
        if len(match) == 1:
            if match.iloc[0] >= THRESHOLD:
                fq_to_pool_mapping.append([fq_file, match.index[0], is_rc])
            else:
                fq_to_pool_mapping.append([fq_file, "Unknown", None])

    # Writing to a CSV file
    with open("fastq_to_pool_mapping.csv", "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["fastq_file", "pool", "is_rc"])
        writer.writerows(fq_to_pool_mapping)

    if not allow_ambiguous:
        ambiguous = [str(x[0]) for x in fq_to_pool_mapping if x[1] == "Ambiguous"]
        if ambiguous:
            raise ValueError(f"Ambiguous pools detected in the following fastq files: {ambiguous}")


def main():
    parser = argparse.ArgumentParser("Read first 1000 reads of fastq files and detect library pool")
    parser.add_argument("--fastqDir", required=True, type=Path, help="Directory containing fastq files")
    parser.add_argument("--libraryStruct", required=True, type=Path)
    parser.add_argument("--scalePlexLibraryStruct", required=False, type=Path)
    parser.add_argument("--allowAmbiguous", action="store_true", help="Don't error if fastq has multiple pools")
    args = parser.parse_args()

    whitelist = load_pool_whitelist(args.libraryStruct, args.scalePlexLibraryStruct)

    sniff_fastq(fastq_dir=args.fastqDir, whitelist=whitelist, allow_ambiguous=args.allowAmbiguous)


if __name__ == "__main__":
    main()
