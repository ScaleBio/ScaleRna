#!/usr/bin/env python
"""
Merge cell metrics csv files
"""
import argparse
from pathlib import Path
import pandas as pd

def merge_metrics(id:str, csv_files:list[Path], rna_ids:list[str]):
    """
    Merge cell metrics
    
    Args:
        id: Sample/lib name
        csv_files: List of cell metrics files to merge
        rna_ids: Associated rna ids
    """
    dfs = []
    for idx, file in enumerate(csv_files):
        df = pd.read_csv(file)
        df["Cell_Barcode"] = df["Cell_Barcode"] + "_" + rna_ids[idx]
        dfs.append(df)
    merged = pd.concat(dfs)
    merged.to_csv(f"{id}.cellMetrics.csv", index=False)

def main():
    parser = argparse.ArgumentParser("Merge individual metrics files for sample")
    parser.add_argument("metrics", nargs="+", type=Path, help="Metrics files to merge")
    parser.add_argument("--id", required=True, help="Sample and lib name")
    parser.add_argument("--rnaId", nargs="+", type=str, required=True, help="Sample and lib name")
    args = parser.parse_args()

    merge_metrics(
        id=args.id,
        csv_files=args.metrics,
        rna_ids=args.rnaId
    )
    
if __name__ == "__main__":
    main()
 