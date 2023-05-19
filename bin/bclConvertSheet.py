#!/usr/bin/env python

"""
Python script to create bcl_convert samplesheet.csv from workflow
samples.csv

Args:
    samples (str): Path to csv file with samples and index sequences
        for scRNA workflow run
    lib (str): Path to json file with library structure definition
    runinfo (str): Path to xml file with sequencer runinfo

Output:
    Contents of samplesheet.csv are written to stdout
"""
import json
import csv
import xml.etree.ElementTree as ET
import argparse
import sys
from pathlib import Path

# Fixed samplesheet.csv settings
SETTINGS = {"CreateFastqForIndexReads": "1",
            "MinimumTrimmedReadLength": "16",
            "MaskShortReads": "16",
            }


def load_libBarcodes(tsv):
    """
    Function to load named library index sequences (e.g. S701P, ...)
    from a .tsv file. Each name can refer to one or multiple sequences
    (when using multiple indices for one sample)

    Args:
        tsv (str): Path to tsv file

    Returns:
        dict name -> [seqs]
    """
    if not tsv:
        return {}

    bcs = {}

    for line in open(tsv):
        if line.startswith("#"):
            continue

        line_split = line.strip().split()

        if len(line_split) <= 1:
            continue

        bcs[line_split[0]] = line_split[1:]

    return bcs


def load_lib(jlib: Path):
    """
    Function to load fastq generation settings from library.json

    Args:
        jlib (str): Path to library json file

    Returns:
        Dictionary with settings
    """
    res = {}
    lib = json.load(open(jlib))
    bclOpts = lib.get("bcl_convert", {})
    res['adapter'] = bclOpts.get('adapter', None)
    if (fqBcs := bclOpts.get('library_barcodes', None)):
                fqBcs = jlib.parent.joinpath(fqBcs)
    res['library_barcodes'] = fqBcs

    for bc in lib.get("barcodes",[]):
        if bc.get("name") == "P7":
            indexSeqsFn = bc.get("sequences")
            if indexSeqsFn: 
                res['indexSeqsFn'] = jlib.parent / indexSeqsFn
    return res

def load_index_seqs(fn: Path):
    """
    Load sequenes for one barcode from a text file specified in library.json

    Args:
        full path to the barcode sequence text file
    
    Returns:
        Dict mapping name (e.g. well position) to barcode sequence (or sequence to sequence if no names in the file)
    """
    res = {}
    for l in open(fn):
        line = l.strip().split()
        seq = line[0]
        name = line[1] if len(line) >= 1 else seq
        res[name] = seq
    return res

def load_run(runInfo):
    """
    Function to load read-lengths from RunInfo.xml in Illumina RunFolder

    Args:
        runInfo (str): Path to RunInfo.xml

    Returns:
        List containing reads
    """
    reads = []
    xml = ET.parse(open(runInfo))

    for read in xml.getroot().findall("Run/Reads/Read"):
        reads.append(
            (read.attrib['NumCycles'],
             read.attrib['IsIndexedRead'] == 'Y'))

    return reads


def load_samples(samplesCsv, bcs):
    """
    Function to load libraries (fastq samples) from samples.csv

    Args:
        samplesCsv (str): Path to samples csv
        bcs (dict):
    Returns:
        Dictionary with the format {library_name: [indexSeqs], ...}
    """
    libs = {}

    with open(samplesCsv) as csvfile:
        for row in csv.DictReader(csvfile):
            # For each library (fastq file set) get the library name and
            # associated index sequences
            # If neither is given, use the sample name (assuming one
            # sample per library)
            name = row["libName"]

            for n in name:
                if not (n.isalnum() or n in "-."):
                    raise ValueError(f"sample name should only contain"
                                     f" [a-z],[A-Z],[0-9], dash (-) or "
                                     f"dot (.): {name}")

            # Library (fastq) index sequences
            # ;-separated list of sequences or names from 'libraryIndex'
            # file (which can be multiple sequences per set)
            outSeqs = []
            indexRead = None

            if (seqs := row.get("libIndex")):
                indexRead = "index1"

                for s in seqs.strip().split(';'):
                    if s in bcs:
                        outSeqs.extend(bcs[seqs])

                    else:
                        for n in s:
                            if n not in 'ACTGactg':
                                raise ValueError(f"Unknown library index "
                                                 f"name / sequence: {s}")
                        outSeqs.append(s)

            if name not in libs:
                libs[name] = outSeqs

            elif libs[name] != outSeqs:
                raise ValueError(f"Mismatched index sequences for library"
                                 f" {name}:\n {libs[name]}\n{outSeqs}")

    return (libs, indexRead)


def print_settings(settings):
    """
    Function to print the settings section of Samplesheet.csv

    Args:
        settings (dictionary): Dictionary containing key value pairs
            responsible for the settings section of samplesheet.csv

    Output:
        Writes to stdout the settings
    """
    print("[Settings]")
    for (s, v) in settings.items():
        print(f"{s},{v}")


def main(samplesCsv: Path, libJson: Path, runInfo: Path, splitFastq: bool, settings):

    lib = load_lib(libJson)

    if (fqBcs := lib['library_barcodes']):
        # Search for file relative to lib.json directory
        fqBcs = load_libBarcodes(fqBcs)

    samples, indexReadUsed = load_samples(samplesCsv, fqBcs)

    # If indexSeqsFn is given in the lib.json, and sampls.csv does not already
    # specify multiple libraries with different PCR index sequences (libIndex)
    # we list all allowed PCR BCs for this library 
    # (i.e. all expected barcodes go in one sample and any others go to 'Undetermined')
    # 
    # If splitFastq is set, we split the samples by PCR barcode for parallelization
    # in the pattern ("SampleName_PcrIndex")
    # Those will be merged again later in the workflow
    if (fn := lib.get("indexSeqsFn")) and not indexReadUsed:
        indexSeqs = load_index_seqs(fn)
        indexReadUsed = "index1"
        if not splitFastq:
            for sample in samples:
                samples[sample] = list(indexSeqs.values())
        else:
            splitSamples = {}
            for sample in samples:
                for name,seq in indexSeqs.items():
                    splitSamples[f"{sample}_{name}"] = [seq]
            samples = splitSamples
    
    reads = load_run(runInfo)
    # Adapter trimming
    if (adapt := lib['adapter']):
        settings['AdapterRead1'] = adapt
        settings['AdapterRead2'] = adapt

    # If an index read is not used for library / fastq demux,
    # we need to define it as a 'UMI' in OverrideCycles
    overrideCycles = []
    nextIndexRead = "index1"
    hasUmi = False

    for (length, index) in reads:
        if not index:
            overrideCycles.append(f"Y{length}")

        else:
            if indexReadUsed == nextIndexRead:
                overrideCycles.append(f"I{length}")
            else:
                overrideCycles.append(f"U{length}")
                hasUmi = True

            nextIndexRead = "index2"

    settings['OverrideCycles'] = ';'.join(overrideCycles)

    # Add UMI settings only if at least one read marked as UMI
    if hasUmi:
        settings['TrimUMI'] = '0'

    print_settings(settings)

    print("[Data]")

    if indexReadUsed:
        print("Sample_ID", "index", sep=',')
    else:
        print("Sample_ID")

    for name, seqs in samples.items():
        if not seqs:
            print(name)
        else:
            for i in seqs:
                # One line per index sequence (with repeated sampleName
                # if multiple indicies)
                print(name, i, sep=',')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create bcl_convert samplesheet.csv '
                    'from workflow samples.csv')
    parser.add_argument('samples', metavar='SAMPLES.csv', type=Path,
        help='CSV with samples and index sequences for scATAC workflow run')
    parser.add_argument('lib', metavar='LIBRARY.json', type=Path,
        help='Library structure definition')
    parser.add_argument('runinfo', metavar='RUNINFO.xml', type=Path,
        help='Sequencer runinfo (in runfolder)')
    parser.add_argument('--splitFastq', action="store_true",
        help="Split sample by PCR index barcode sequence")
    args = parser.parse_args()

    main(args.samples, args.lib, args.runinfo, args.splitFastq, SETTINGS)
