#!/usr/bin/env python
"""
Python script to merge multiple demuxJson.json files
"""
from collections import defaultdict
import argparse
import json
import os

def load_json(fname):
    """
    Function to generate a dictionary from a json file

    Args:
        fname (str): Path to file
    
    Returns:
        Dictionary with json contents
    """
    with open(fname) as f:
        return json.load(f)

def merge(bc_jsons, lib_json):
    """
    Function to merge multiple demuxMetrics.json files

    Args:
        bc_jsons (list): List of demuxMetrics.json files
        lib_json (str): Path to library json file
    
    Output:
        Writes combined demuxMetrics.json to cwd
    """
    # Get contents of all jsons in a single list
    bc_jsons_list = [load_json(bc_json) for bc_json in bc_jsons]

    # Initialize dict that will hold combined information
    master_dict = defaultdict(dict,{ k:{} for k in list(bc_jsons_list[0].keys()) })
    
    # Load library json
    lib_json_dict = load_json(lib_json)
    
    # Get barcode levels
    names_in_lib_json = [x["name"] for x in lib_json_dict["barcodes"]]

    # Initialize dict which will hold combined metrics from "barcodes" section
    summation_metrics = defaultdict(defaultdict, {k:{} for k in names_in_lib_json})
    
    for key in summation_metrics:
        summation_metrics[key] = defaultdict(int)
    
    # Loop over all jsons
    for bc_json in bc_jsons_list:
        for barcode_level in bc_json["barcodes"]:
            for read_metrics in bc_json["barcodes"][barcode_level]:
                summation_metrics[barcode_level][read_metrics] += bc_json["barcodes"][barcode_level][read_metrics][0]

    # Loop to assign metrics from summation_metrics to master_dict
    for barcode_level in summation_metrics:
        master_dict["barcodes"][barcode_level] = {}
        total = sum(summation_metrics[barcode_level].values())

        for key in summation_metrics[barcode_level]:
            master_dict["barcodes"][barcode_level][key] = [summation_metrics[barcode_level][key],
                                                           f"{round(100*(summation_metrics[barcode_level][key]/total), 1)}%"]

    # Initialize dict which will hold combined metrics from "reads" section        
    summation_metrics = dict.fromkeys(bc_jsons_list[0]["reads"].keys(), 0)
    for bc_json in bc_jsons_list:
        for key in bc_json["reads"]:
            summation_metrics[key] += bc_json["reads"][key][0]
    
    total_reads = sum(summation_metrics.values())
    
    # Loop to assign metrics from summation_metrics to master_dict
    for key in summation_metrics:
        master_dict["reads"][key] = [summation_metrics[key], f"{round(100*(summation_metrics[key]/total), 1)}%" ]
    
    # Initialize dict which will hold combined metrics from "samples" section
    summation_metrics = defaultdict(dict,{ k:{} for k in list(bc_jsons_list[0]["samples"].keys()) })
    for key in summation_metrics:
        summation_metrics[key]["reads"] = 0
        summation_metrics[key]["barcodes"] = {}
    
    # Initialize per well metrics in the master dict holding combined information
    for bc_json in bc_jsons_list:
        for sample_name in bc_json["samples"]:
            for well in bc_json["samples"][sample_name]["barcodes"]:
                summation_metrics[sample_name]["barcodes"][well] = {"sequence": "", "reads": 0}

    total = 0
    for bc_json in bc_jsons_list:
        for sample_name in bc_json["samples"]:
            for key in bc_json["samples"][sample_name]:
                # If key is "name" or "number" value will be same in combined dict as individual dicts
                if key == "name" or key == "number":
                    summation_metrics[sample_name][key] = bc_json["samples"][sample_name][key]
                elif key == "reads":
                    # 0 for number of reads and 1 for percent of reads
                    total += bc_json["samples"][sample_name][key][0]
                    summation_metrics[sample_name][key] += bc_json["samples"][sample_name][key][0]
                # Key is "barcodes"
                else:
                    for well in bc_json["samples"][sample_name][key]:
                        # Sum up reads in all demux jsons for that specific barcode
                        summation_metrics[sample_name][key][well]["reads"] += bc_json["samples"][sample_name][key][well]["reads"]
                        summation_metrics[sample_name][key][well]["sequence"] = bc_json["samples"][sample_name][key][well]["sequence"]

    # Construct master dict that will hold final merged information
    for sample_name in summation_metrics:
        master_dict["samples"][sample_name] = {}
        for key in summation_metrics[sample_name]:
            if key == "name" or key == "number" or key == "barcodes":
                master_dict["samples"][sample_name][key] = summation_metrics[sample_name][key]
            else:
                master_dict["samples"][sample_name][key] = [summation_metrics[sample_name][key],
                                                            f"{round(100*(int(summation_metrics[sample_name][key])/total_reads), 1)}%"]

    with open('metrics.json', "w") as f:
        json.dump(master_dict, f)


def main():
    parser = argparse.ArgumentParser(
        description="Merge bc parser output files into one")
    parser.add_argument("--bc_jsons", nargs='+', type=str,
                        help="bc parser json files that need to be concatenated",
                        required=True)
    parser.add_argument("--lib_json", type=str, help="Path to library json file",
                        required=True)
    args = parser.parse_args()
    
    merge(args.bc_jsons, args.lib_json)


if __name__ == "__main__":
    main()
