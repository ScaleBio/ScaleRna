#!/usr/bin/env python
"""
Python script to create bcl_convert samplesheet.csv from workflow
samples.csv
Contents of samplesheet.csv are written to stdout
"""

import argparse
import csv
import dataclasses
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Any, List, Dict
import xml.etree.ElementTree as ET
from scale_utils.validation import validateName

# Fixed samplesheet.csv settings
SETTINGS = {
    "CreateFastqForIndexReads": "1",
    "MinimumTrimmedReadLength": "16",
    "MaskShortReads": "16",
}


@dataclass
class Library:
    """Info on index sequences for one library (Illumina samplesheet.csv sample)"""

    name: str
    # List of index1 (i7) barcode sequences
    index1: List[str] = dataclasses.field(default_factory=lambda: [])
    # List of index2 (i5) barcode sequences
    # Orientation as given in samples.csv; i.e. never reverse complemented
    index2: List[str] = dataclasses.field(default_factory=lambda: [])
    # Unique Dual Indexing; i.e. i5&i7 occur in unique combinations
    # Means len(index1) == len(index2) and index1[i] is paired only with index2[i]
    isUniqueDual: bool = False


@dataclass
class ReadInfo:
    """Information about one read (mate) from RunInfo.xml"""

    num: int  # Read number (1 = R1, 2=R2 or I1, etc.)
    bp: int  # read length
    isIndex: bool  # Is index read

    # Is this read in reverse complement; relevant for Index2 (i5) on some platforms
    # Even though all of NextSeq550, NovaSeq 6k, NextSeq2k and NovaSeq X have i5 in reverse orientation, this
    # flag in RunInfo.xml is only set on NextSeq2k, NovaSeqX and later. If bcl_convert sees this flag
    # it rev. comps. the sequences in samplesheet.csv. Hence if this flag is present we need to output
    # the opposite orientation then when it is not; even if the underlying sequencing read is the same
    # (NovaSeq 6k vs. NextSeq2k)
    isRevComp: bool


@dataclass
class IndexInfo:
    """Index setting info to generate a samplesheet.csv fitting all of samples.csv, library.json and runInfo.xml"""

    index1Used: bool  # Include 'index' column in samplesheet
    index2Used: bool  # Include 'index2' column in samplesheet
    revCompIndex2: bool  # Rev. comp. index2 sequences to match sequencing orientation


def revComp(seq: str) -> str:
    """Reverse complement a sequence"""
    return seq.translate(str.maketrans("acgtACGT", "tgcaTGCA"))[::-1]


def load_bcs(tsv: Path) -> Dict[str, str]:
    """
    Load barcode sequences from a .tsv file as used by library.json

    Returns:
        Dict mapping name (e.g. well position) to barcode sequence (or sequence to sequence if no names in the file)
    """
    bcs = {}
    for line in open(tsv):
        if line.startswith("#"):
            continue
        line = line.strip().split()
        # No name given, use sequence as name
        if len(line) == 1:
            seq = line[0]
            name = seq
        else:
            name = line[0]
            seq = line[1:]
        bcs[name] = seq
    return bcs


def load_libDef(jlib: Path) -> Dict[str, Any]:
    """
    Load fastq generation settings from library structure definition json

    Args:
        jlib: Library.json file

    Returns:
        Dict: SettingName -> Value
    """
    lib = json.load(open(jlib))
    bclOpts = lib.get("bcl_convert", {})
    res = {}
    res["adapter"] = bclOpts.get("adapter", None)
    res["libIndexUsed"] = None
    res["library_barcodes"] = {}
    if bclOpts.get("library_barcodes"):
        if isinstance(bclOpts["library_barcodes"], dict):
            res["library_barcodes"] = bclOpts["library_barcodes"]["sequences"]
            res["libIndexUsed"] = bclOpts["library_barcodes"]["libIndexUsed"]
        else:
            res["library_barcodes"] = bclOpts["library_barcodes"]
    res["uniqueDualIndex"] = bclOpts.get("uniqueDualIndex", False)
    res["index2RevComp"] = bclOpts.get("index2RevComp", False)
    res["split_on"] = bclOpts.get("split_on")
    res["kit"] = bclOpts.get("kit")
    res["index1BClen"] = bclOpts.get("index1BClen")
    res["index2BClen"] = bclOpts.get("index2BClen")
    if seqsFn := bclOpts.get("index1Seqs"):
        res["index1Seqs"] = jlib.parent / seqsFn
    if seqsFn := bclOpts.get("index2Seqs"):
        res["index2Seqs"] = jlib.parent / seqsFn
    return res


def load_run(runInfo: Path) -> Dict[str, ReadInfo]:
    """
    Function to load read-lengths from RunInfo.xml in Illumina RunFolder

    Args:
        runInfo (str): Path to RunInfo.xml

    Returns:
        Dict mapping "R1" -> read1Info, "I1" -> index1Info, etc.
    """
    reads = {}
    nextRead = 1
    nextIndex = 1
    xml: Any = ET.parse(open(runInfo))
    for read in xml.getroot().findall("Run/Reads/Read"):
        readInfo = ReadInfo(
            num=read.attrib["Number"],
            bp=read.attrib["NumCycles"],
            isIndex=read.get("IsIndexedRead") == "Y",
            isRevComp=read.get("IsReverseComplement") == "Y",
        )
        if readInfo.isIndex:
            reads[f"I{nextIndex}"] = readInfo
            nextIndex += 1
        else:
            reads[f"R{nextRead}"] = readInfo
            nextRead += 1
    return reads


def load_libraries(samplesCsv: Path, namedIndexSeqs: Dict[str, List[str]], libDef: Dict[str, str]) -> List[Library]:
    """
    Load information about library demux (Illumina / fastq samples) from samples.csv

    Each line in samples.csv is one sample, if multiple samples in the same library
    the library definition is repeated
    """
    libs: Dict[str, Library] = {}  # LibName -> Lib
    with open(samplesCsv) as csvfile:
        for row in csv.DictReader(csvfile):
            # For each library (fastq file set) get the library name and
            # associated index sequences
            # If neither is given, use the sample name (assuming one
            # sample per library)
            name = row["libName"]
            validateName(name, "Library name", other_chars="-")
            lib = Library(name=name)

            # Library (fastq) index sequences
            # ;-separated list of sequences or names from 'libraryIndex'
            # file (which can be multiple sequences per set)
            if seqs := row.get("libIndex"):
                for s in seqs.split(";"):
                    s = s.strip()
                    if s in namedIndexSeqs:
                        lib.index1.extend(namedIndexSeqs[seqs])
                    elif any((n not in "ACTGactg" for n in s)):
                        raise ValueError(f"Unknown library index: {s}")
                    else:
                        lib.index1.append(s)
            else:
                if libDef["libIndexUsed"] != "libIndex2":
                    # default in all library index1 sequences
                    for seqs in namedIndexSeqs.values():
                        lib.index1.extend(seqs)

            if seqs := row.get("libIndex2", "").strip():
                if seqs == "*":
                    lib.isUniqueDual = True
                    lib.index2 = lib.index1
                else:
                    for s in seqs.split(";"):
                        s = s.strip()
                        if s in namedIndexSeqs:
                            lib.index2.extend(namedIndexSeqs[seqs])
                        elif any((n not in "ACTGactg" for n in s)):
                            raise ValueError(f"Unknown library index2: {s}")
                        else:
                            lib.index2.append(s)
            else:
                if libDef["libIndexUsed"] == "libIndex2" or libDef["libIndexUsed"] == "scalePlexLibIndex":
                    # default in all library index2 sequences
                    for seqs in namedIndexSeqs.values():
                        lib.index2.extend(seqs)
            # Check for consistency with previous occurance of this library (from other samples)
            if oldLib := libs.get(name):
                if oldLib.index1 != lib.index1:
                    raise ValueError(f"Mismatched index sequences for library {name}")
                if oldLib.index2 != lib.index2:
                    raise ValueError(f"Mismatched index2 sequences for library {name}")
                if oldLib.isUniqueDual != lib.isUniqueDual:
                    raise ValueError(f"Mismatched UDI information for library {name}")
            libs[name] = lib
    return list(libs.values())


def assign_index(libraries: List[Library], libDef: Dict, reads: Dict[str, ReadInfo]) -> IndexInfo:
    """
    Finalize index sequences for each library (fastq sample)

    If no barcodes are given for demux in samples.csv, use the full list of sequences from
    the library structure .json

    Check with index reads are actually used (used in the library and sequenced in the run)
    and whether Index2 should be reverse complemented
    """
    # lib.json can define a the full list of barcode sequences to use for index1 and/or 2
    if bcFn := libDef.get("index1Seqs"):
        allSeqs = load_bcs(bcFn)
        for lib in libraries:
            if not lib.index1:
                lib.index1 = list(allSeqs.values())
    if bcFn := libDef.get("index2Seqs"):
        allSeqs = load_bcs(bcFn)
        for lib in libraries:
            if not lib.index2:
                lib.index2 = list(allSeqs.values())
    # Index read is used in samplesheet.csv if we have index sequences
    # and the read was actually sequenced in this run
    index1Used = ("I1" in reads) and any(lib.index1 for lib in libraries)
    index2Used = ("I2" in reads) and any(lib.index2 for lib in libraries)
    # If we use an index read for demux, check that we have sequences from all libraries
    for lib in libraries:
        if index1Used and not lib.index1:
            raise ValueError(f"Missing index1 for {lib.name}")
        if index2Used and not lib.index2:
            raise ValueError(f"Missing index2 for {lib.name}")

    print(libDef, file=sys.stderr)
    revCompIndex2 = libDef["index2RevComp"]
    # Rev. comp. index2 if needed (depending on sequencing platform and flags in RunInfo.xml used by bcl_convert)
    if index2Used and reads["I2"].isRevComp:
        revCompIndex2 = not revCompIndex2
    print(revCompIndex2, file=sys.stderr)
    return IndexInfo(index1Used, index2Used, revCompIndex2)


def print_settings(settings: Dict):
    """
    Print the settings section of Samplesheet.csv to stdout

    Args:
        settings: samplesheet settings as key,value pairs
    """
    print("[Settings]")
    for s, v in settings.items():
        print(f"{s},{v}")


def main(samplesCsv: Path, libJson: Path, runInfo: Path, splitFastq: bool, settings: dict[str, str], libDef: dict):
    """
    Prepare the bcl_convert samplesheet.csv and print to stdout

    Args:
        samplesCsv: Path to input samples.csv file
        libJson: Path to library structure definition
        runInfo: Path to RunInfo.xml
        splitFastq: Split sample by PCR index barcode sequence
        settings: Fixed settings for samplesheet.csv
        libDef: Library structure definition
    """
    if splitFastq and not libDef["split_on"]:
        raise Exception(
            "If using splitFastq please provide either index1 or index2"
            "as a value of the split_on field in the library structure json"
        )

    if fqBcs := libDef["library_barcodes"]:
        # Search for file relative to lib.json directory
        fqBcs = libJson.parent.joinpath(fqBcs)
        fqBcs = load_bcs(fqBcs)

    reads = load_run(runInfo)
    # Adapter trimming
    if adapt := libDef["adapter"]:
        settings["AdapterRead1"] = adapt
        settings["AdapterRead2"] = adapt
    # When fastqs are split on PCR index, we do not also split on lane to keep
    # the file number reasonable (e.g. for fastqc)
    if splitFastq:
        settings["NoLaneSplitting"] = "true"
    libs = load_libraries(samplesCsv, fqBcs, libDef)
    reads = load_run(runInfo)
    indexInfo = assign_index(libs, libDef, reads)
    # If splitFastq is set, we split the samples by PCR barcode for parallelization
    # in the pattern ("SampleName_PcrIndex")
    # Those will be merged again later in the workflow
    if libDef["split_on"] == "index1":
        if splitFastq and indexInfo.index1Used and (bcFn := libDef.get("index1Seqs")):
            indexSeqs = load_bcs(bcFn)
            splitSamples = []
            for lib in libs:
                for bc_name, seq in indexSeqs.items():
                    if seq in lib.index1:
                        subLib = Library(name=f"{lib.name}_{bc_name}")
                        subLib.index1 = [seq]
                        subLib.index2 = lib.index2
                        splitSamples.append(subLib)
            libs = splitSamples
    elif libDef["split_on"] == "index2":
        if splitFastq and indexInfo.index2Used and (bcFn := libDef.get("index2Seqs")):
            indexSeqs = load_bcs(bcFn)
            splitSamples = []
            for lib in libs:
                for bc_name, seq in indexSeqs.items():
                    if seq in lib.index2:
                        subLib = Library(name=f"{lib.name}_{bc_name}")
                        subLib.index1 = lib.index1
                        subLib.index2 = [seq]
                        splitSamples.append(subLib)
            libs = splitSamples
    elif libDef["split_on"] == "index1,index2":
        if (
            splitFastq
            and indexInfo.index2Used
            and indexInfo.index1Used
            and (bc1Fn := libDef.get("index1Seqs"))
            and (bc2Fn := libDef.get("index2Seqs"))
        ):
            indexSeqs1 = load_bcs(bc1Fn)
            indexSeqs2 = load_bcs(bc2Fn)
            splitSamples = []
            for lib in libs:
                for bc1_name, seq1 in indexSeqs1.items():
                    for bc2_name, seq2 in indexSeqs2.items():
                        if seq1 in lib.index1 and (all(seq in lib.index2 for seq in seq2) or seq2 in lib.index2):
                            subLib = Library(name=f"{lib.name}_{bc1_name}_{bc2_name}")
                            subLib.index1 = seq1
                            subLib.index2 = seq2
                            splitSamples.append(subLib)
            libs = splitSamples
    else:
        raise Exception("Valid value of split_on not provided. Please provide either index1 or index2")
    # If an index read is not used for library / fastq demux,
    # we need to define it as a 'UMI' in OverrideCycles
    overrideCycles = []
    for readName, readInfo in reads.items():
        if not readInfo.isIndex:
            overrideCycles.append(f"Y{readInfo.bp}")
        else:
            if (indexInfo.index1Used and readName == "I1") or (indexInfo.index2Used and readName == "I2"):
                if readName == "I1" and libDef.get("index1BClen"):
                    overrideCycles.append(f"I{libDef['index1BClen']}U{int(readInfo.bp) - int(libDef['index1BClen'])}")
                elif readName == "I2" and libDef.get("index2BClen"):
                    overrideCycles.append(f"U{int(readInfo.bp) - int(libDef['index2BClen'])}I{libDef['index2BClen']}")
                else:
                    overrideCycles.append(f"I{readInfo.bp}")
            else:
                overrideCycles.append(f"U{readInfo.bp}")
    settings["OverrideCycles"] = ";".join(overrideCycles)
    # Add UMI settings only if at least one read marked as UMI
    if "U" in settings["OverrideCycles"]:
        settings["TrimUMI"] = "0"

    # Needed because two Quantum PCR sequences are less than 3 hamming distance apart
    # Remove once assay team fixes the sequences
    if "Quantum" in libJson.name:
        settings["BarcodeMismatchesIndex2"] = "0"
    print_settings(settings)

    print("[Data]")
    # Write header line
    if indexInfo.index1Used and indexInfo.index2Used:
        print("Sample_ID", "index", "index2", sep=",")
    elif indexInfo.index1Used:
        print("Sample_ID", "index", sep=",")
    elif indexInfo.index2Used:
        print("Sample_ID", "index2", sep=",")
    else:
        print("Sample_ID", sep=",")
    # Write one line per index sequence combination
    # repeat sampleName if multiple indicies / combinations
    for lib in libs:
        i1s = lib.index1 if isinstance(lib.index1, list) else [lib.index1]
        i2s = lib.index2 if isinstance(lib.index2, list) else [lib.index2]
        if indexInfo.revCompIndex2:
            i2s = [revComp(i) for i in i2s]
        if indexInfo.index2Used and lib.isUniqueDual:
            # Only matched i1,i2 pairs
            for i1, i2 in zip(i1s, i2s):
                print(lib.name, i1, i2, sep=",")
        elif indexInfo.index1Used:
            for i1 in i1s:
                if indexInfo.index2Used:
                    # All i1 * i2 combinations
                    for i2 in i2s:
                        print(lib.name, i1, i2, sep=",")
                else:
                    print(lib.name, i1, sep=",")
        elif indexInfo.index2Used:
            for i2 in i2s:
                print(lib.name, i2, sep=",")
        else:
            print(lib.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create bcl_convert samplesheet.csv " "from workflow samples.csv")
    parser.add_argument(
        "samples", metavar="SAMPLES.csv", type=Path, help="CSV with samples and index sequences for scATAC workflow run"
    )
    parser.add_argument("libDef", metavar="LIBRARY.json", type=Path, help="Library structure definition")
    parser.add_argument("runinfo", metavar="RUNINFO.xml", type=Path, help="Sequencer runinfo (in runfolder)")
    parser.add_argument("--splitFastq", action="store_true", help="Split sample by PCR index barcode sequence")
    args = parser.parse_args()

    libDef = load_libDef(args.libDef)

    main(args.samples, args.libDef, args.runinfo, args.splitFastq, SETTINGS, libDef)
