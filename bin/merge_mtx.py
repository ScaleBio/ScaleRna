#!/usr/bin/env python
"""
Merge raw umi count matrices
"""
import argparse
from pathlib import Path
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack

def merge_mtx(id:str, mtx_files:list[Path]):
    """
    Merge raw umi matrices
    
    Sparse coordinate matrices are horizontally concatenated, rows are features and columns are cells.
    
    Args:
        sample: Sample name
        mtx_files: List of sparse coordinate mtx files to merge
    """
    raw_mtx = [mmread(file) for file in mtx_files]
    mmwrite(f"{id}.raw.umi.mtx", hstack(raw_mtx))

def main():
    parser = argparse.ArgumentParser("Merge individual umi count matrices for sample")
    parser.add_argument("raw_mtx", nargs="+", type=Path, help="sparse coordinate mtx files to merge")
    parser.add_argument("--id", required=True, help="Sample and lib name")
    args = parser.parse_args()

    merge_mtx(
        id=args.id,
        mtx_files=args.raw_mtx
    )
    
if __name__ == "__main__":
    main()
 