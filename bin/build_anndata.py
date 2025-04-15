#!/usr/bin/env python

"""
This script takes the filtered.matrix directory and the allCells.csv file
for a sample and creates an HDF5 backed AnnData object.

Args:
    mtx_dir: Path to the filtered.matrix directory
    all_cells: Path to the allCells.csv file
"""

import argparse
import h5py
import os
import gzip
import pandas
import numpy as np
import anndata
from pathlib import Path
from natsort import natsorted


def get_mtx_shape(mtx_file: Path) -> dict[int, int, int]:
    """
    Reads the header of a Matrix Market file to extract the shape of the matrix.

    Args:
        mtx_file: Path to the matrix.mtx file.

    Returns:
        dict: Contains the number of genes, cells, and counts (non-zero values) in the matrix.
    """
    with gzip.open(mtx_file, "rt") as f:
        # Skip headers
        f.readline()
        f.readline()
        n_genes, n_cells, n_counts = map(int, f.readline().split())

    print("Matrix Shape")
    print(f"Genes: {n_genes}, Cells: {n_cells}, Counts: {n_counts}")

    return {"genes": n_genes, "cells": n_cells, "counts": n_counts}


def add_features_to_anndata(feature_file: Path, dest_file: Path) -> None:
    """
    Adds gene features (e.g., gene IDs, symbols, annotations) to the HDF5 AnnData file.

    Args:
        feature_file: Path to features.tsv file.
        dest_file: Path to the output HDF5 file.
    """
    with h5py.File(dest_file, "r+") as hdf:

        cols = ["gene_id", "gene_symbol", "annotation"]

        # Read features.tsv
        gene_dat = pandas.read_csv(filepath_or_buffer=feature_file, sep="\t", names=cols)

        if "var" in hdf:
            # Create "var" group in HDF5 and add attributes
            var = hdf["/var"]
        else:
            var = hdf.create_group("var")

        var.attrs["_index"] = cols[0]
        var.attrs["column-order"] = cols
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"

        for col in cols:
            var.create_dataset(
                col, data=gene_dat[col].tolist(), dtype=h5py.string_dtype(), compression="gzip", chunks=True
            )

    print(f"Added {feature_file} to {dest_file}")


def add_all_cells_to_anndata(all_cells: Path, dest_file: Path, mtx_shape: dict, chunk_size: int = 10000) -> None:
    """
    Adds cell metadata from a CSV file to the HDF5 AnnData file.

    Args:
        all_cells: Path to the allCells.csv file.
        dest_file: Path to the output HDF5 file.
        mtx_shape: Shape information of the matrix.
        chunk_size: Number of rows to read at a time.
    """
    with h5py.File(dest_file, "r+") as hdf:
        # Read initial rows to determine column names and data types
        csv_chunk = pandas.read_csv(
            all_cells, nrows=100, dtype={"flags": "str", "classification": "str", "species": "str"}
        )
        dtypes = dict(csv_chunk.dtypes)
        cols = list(dtypes.keys())
        datasets = {}

        # Create "obs" group in HDF5 and add attributes
        obs = hdf.create_group("obs")
        obs.attrs["_index"] = cols[0]
        obs.attrs["column-order"] = cols
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"

        # Create datasets for each column
        for col in cols:
            if dtypes[col] == "object":
                dtype = h5py.string_dtype()
            else:
                dtype = dtypes[col]

            datasets[col] = obs.create_dataset(
                col, shape=(0,), maxshape=(mtx_shape["cells"],), dtype=dtype, compression="gzip", chunks=True
            )

        # Read the CSV file in chunks and append to datasets
        for chunk in pandas.read_csv(
            all_cells,
            iterator=True,
            chunksize=chunk_size,
            keep_default_na=False,
            dtype={"flags": "str", "classification": "str", "species": "str"},
        ):
            if chunk.shape[0] != 0:
                for col in chunk.columns:
                    values = chunk[col].tolist()
                    current_size = datasets[col].shape[0]
                    datasets[col].resize(current_size + len(values), axis=0)
                    datasets[col][-len(values) :] = values

    print(f"Added {all_cells} to {dest_file}")


def add_mtx_to_anndata(mtx_file: Path, dest_file: Path, mtx_shape: dict, chunk_size: int = 10000) -> None:
    """
    Adds sparse matrix data from the Matrix Market file to the HDF5 AnnData file.

    Args:
        mtx_file: Path to the matrix.mtx file.
        dest_file: Path to the output HDF5 file.
        mtx_shape: Shape information of the matrix.
    """
    n_cells = mtx_shape["cells"]
    n_genes = mtx_shape["genes"]
    n_counts = mtx_shape["counts"]

    with h5py.File(dest_file, "r+") as hdf:
        # Create "X" group and datasets
        hdf_mtx = hdf.create_group("X")
        mtx_data = hdf_mtx.create_dataset(
            "data", dtype=np.float32, shape=(n_counts,), chunks=True, compression="gzip"
        )  # Nonzero values from .mtx
        mtx_ind = hdf_mtx.create_dataset(
            "indices", dtype=np.int32, shape=(n_counts,), chunks=True, compression="gzip"
        )  # Column or Gene Indices for non-zero values from .mtx
        mtx_indptr = hdf_mtx.create_dataset(
            "indptr",
            dtype=np.int64,
            shape=(n_cells + 1,),
            chunks=True,
            compression="gzip",  # Starting position in mtx_data and mtx_indices for each row of the matrix
        )
        mtx_indptr[0] = 0  # Initialize with 0

        hdf_mtx.attrs["encoding-type"] = "csr_matrix"
        hdf_mtx.attrs["encoding-version"] = "0.1.0"
        hdf_mtx.attrs["shape"] = (n_cells, n_genes)

        # Process matrix data
        current_pos = 0
        cell_counts = np.zeros(n_cells + 1, dtype=np.int64)  # Array to track number of non-zero values per cell

        mtx_iter = pandas.read_csv(
            mtx_file,
            sep="\\s+",
            skiprows=3,
            header=None,
            names=["gene", "barcode", "count"],
            dtype={"gene": int, "barcode": int, "count": float},
            iterator=True,
            chunksize=chunk_size,
        )

        for chunk in mtx_iter:
            chunk["gene"] -= 1  # Convert from 1-based to 0-based indexing
            chunk["barcode"] -= 1
            chunk_size = len(chunk)

            mtx_data[current_pos : current_pos + chunk_size] = chunk["count"].values
            mtx_ind[current_pos : current_pos + chunk_size] = chunk["gene"].values

            # Update number of non-zero values for each cell
            for barcode, count in zip(chunk["barcode"].values, chunk["count"].values):
                cell_counts[barcode + 1] += 1

            current_pos += chunk_size

        np.cumsum(cell_counts, out=cell_counts)
        mtx_indptr[:] = cell_counts

    print(f"Added {mtx_file} to {dest_file}")


def convert_strings_to_categoricals(hdf5_file: Path, grps=["obs", "var"]) -> None:
    """
    Converts string datasets in the HDF5 file to categorical datasets where applicable.

    Args:
        hdf5_file: Path to the HDF5 AnnData file.
    """
    with h5py.File(hdf5_file, "r+") as hdf:
        for grp in grps:
            for dataset in hdf[grp].items():
                # Check if dataset consists of strings.
                # If so attempt to convert to categorical.
                if dataset[1].dtype == "object":
                    dat = dataset[1].asstr()[:]
                    dat[dat == ""] = np.nan
                    unique_values = np.unique(dat[~pandas.isna(dat)])
                    if len(unique_values) < len(dat):
                        sorted_cats = natsorted(unique_values)
                        dat = pandas.Categorical(dat, categories=sorted_cats)
                        del hdf[f"{grp}/{dataset[0]}"]
                        new_grp = hdf[grp].create_group(dataset[0])
                        new_grp.create_dataset(
                            "categories", data=dat.categories, dtype=h5py.string_dtype(), compression="gzip"
                        )
                        new_grp.create_dataset("codes", data=dat.codes, compression="gzip")
                        new_grp.attrs["encoding-type"] = "categorical"
                        new_grp.attrs["encoding-version"] = "0.2.0"
                        new_grp.attrs["ordered"] = False


def create_anndata(mtx_dir, all_cells_path, sample_id):

    out_file = Path(f"{sample_id}_anndata.h5ad")
    mtx_path = Path(os.path.join(mtx_dir, "matrix.mtx.gz"))
    feature_path = Path(os.path.join(mtx_dir, "features.tsv.gz"))

    with h5py.File(out_file, "w") as hdf:
        hdf.attrs["encoding-type"] = "anndata"
        hdf.attrs["encoding-version"] = "0.1.0"

        # These groups are technically optional, but the on-disk concatenation fails if they are not present.
        hdf.create_group("obsm")
        hdf.create_group("layers")

        mtx_info = get_mtx_shape(mtx_file=mtx_path)
        add_features_to_anndata(feature_file=feature_path, dest_file=out_file)
        add_all_cells_to_anndata(all_cells=all_cells_path, dest_file=out_file, mtx_shape=mtx_info)
        convert_strings_to_categoricals(hdf5_file=out_file)
        add_mtx_to_anndata(mtx_file=mtx_path, dest_file=out_file, mtx_shape=mtx_info)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--matrix_dir", type=Path, required=True, nargs="+", help="Path to directory with .mtx files.")
    parser.add_argument("--all_cells", type=Path, required=True, nargs="+", help="Path to allCells.csv file.")
    parser.add_argument("--sample_ids", type=str, required=True, nargs="+", help="Sample IDs for each sample.")

    args = parser.parse_args()

    for mtx_path, ac_path, samp_id in zip(args.matrix_dir, args.all_cells, args.sample_ids):
        create_anndata(mtx_dir=mtx_path, all_cells_path=ac_path, sample_id=samp_id)

    if len(args.sample_ids) > 1:
        anndata_list = [f"{samp_id}_anndata.h5ad" for samp_id in args.sample_ids]
        # When trying this command I ran into the issue described on github here:
        # https://github.com/scverse/anndata/issues/1505
        # This is only an issue when using the "merge" parameter.
        # As a result of this issue, the concatenated anndata object does not have the var slot.abs
        # So I add the var data to the anndata object after concatenation.
        # anndata.experimental.concat_on_disk(
        #   in_files=anndata_list, out_file="scale_merged_anndata.h5ad", merge = "same")
        anndata.experimental.concat_on_disk(in_files=anndata_list, out_file="merged_anndata.h5ad")
        # Use features.tsv file from first sample for var data.
        feat_file = Path(os.path.join(args.matrix_dir[0], "features.tsv.gz"))
        add_features_to_anndata(feature_file=feat_file, dest_file="merged_anndata.h5ad")
        convert_strings_to_categoricals(hdf5_file="merged_anndata.h5ad", grps=["var"])


if __name__ == "__main__":
    main()
