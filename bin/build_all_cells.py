#!/usr/bin/env python
"""
Calculates metrics for one sample from STARsolo outputs.
Takes in the path to the STARsolo output directory for one sample and returns an allCells.csv
"""


import argparse
import pandas as pd
import polars as pl
import polars.selectors as cs
import duckdb
from pathlib import Path
from scale_utils import io


def read_genes(sample_specific_file_paths):
    """
    Reads in the features.tsv file from the STARsolo output directory,
    labelling each gene as either mm (mouse) or hs (human).

    Args:
        sample_specific_file_paths (Dict[str, Path]):
            Dictionary where each key is a custom file identifier
            and each value is a path to the identified file in the STARsolo output directory

    Returns:
        A dataframe of genes (DataFrame)
    """
    features = pl.read_csv(
        sample_specific_file_paths["features"], separator="\t", has_header=False, columns=[0], new_columns=["ID"]
    )
    features = features.with_row_index(offset=1)  # add 1-based index column to match entries in mtx file
    features = features.with_columns(
        pl.when(pl.col("ID").str.starts_with("ENSG")).then(pl.lit("human")).otherwise(pl.lit("mouse")).alias("species")
    )
    # check that human and mouse genes are not intermixed
    # the species column should only switch one time from the previous row
    if (features["species"] != features["species"].shift(1)).sum() > 1:
        raise ValueError("Barnyard genomes can not have genes intermixed from different species")
    return features


def read_barcodes(sample_specific_file_paths):
    """
    Reads in the barcodes.tsv file from the STARsolo output directory.

    Args:
        sample_specific_file_paths (Dict[str, Path]):
            Dictionary where each key is a custom file identifier
            and each value is a path to the identified file in the STARsolo output directory

    Returns:
        A dataframe of barcodes (DataFrame)
    """
    barcodes = pl.read_csv(sample_specific_file_paths["barcodes"], has_header=False, new_columns=["cell_id"])
    barcodes = barcodes.with_row_index(
        name="barcode", offset=1
    )  # add 1-based index column to match entries in mtx file
    return barcodes


def count_transcripts(sample_specific_file_paths, features, barcodes, isBarnyard):
    """
    Computes barnyard-specific transcript counts from the raw STARsolo count matrix for each barcode.

    Args:
        sample_specific_file_paths (Dict[str, Path]):
            Dictionary where each key is a custom file identifier
            and each value is a path to the identified file in the STARsolo output directory
        features (DataFrame): A dataframe of genes
        barcodes (DataFrame): A dataframe of barcodes
        isBarnyard (bool): Whether or not this is a barnyard sample (i.e. mixed mouse and human)

    Returns:
        A dataframe of unique transcript counts for mouse and human per barcode (DataFrame)
    """
    barnyard_query = ""
    if isBarnyard:
        for species in ["human", "mouse"]:
            idx = features["index"].filter(features["species"] == species)
            barnyard_query += (
                f"SUM(CASE WHEN gene >= {idx.min()} AND gene <= {idx.max()} "
                f"THEN count ELSE 0 END) AS {species}_counts,\n"
                f"SUM(CASE WHEN gene >= {idx.min()} AND gene <= {idx.max()} THEN 1 ELSE 0 END) AS {species}_genes,\n"
            )
    all_cells = duckdb.sql(
        f"""
    SELECT
        barcode,
        COUNT(count) AS genes,
        SUM(count) as counts,
        {barnyard_query}
    FROM read_csv(
        '{sample_specific_file_paths['mtx']}',
        delim=' ',
        skip=3,
        columns = {{
            'gene': 'UINTEGER',
            'barcode': 'UINTEGER',
            'count': 'FLOAT'
        }})
    WHERE count != 0
    GROUP BY barcode;
    """
    ).pl()

    all_cells = (
        barcodes.join(
            all_cells.with_columns(
                # round numeric columns to nearest integer
                cs.ends_with("counts").round().cast(pl.UInt32),
                cs.ends_with("genes").cast(pl.UInt32),
            ),
            on="barcode",
            how="left",
        )
        .drop(["barcode"])
        .fill_null(0)
    )  # for barcodes with no non-zero counts, fill in 0s

    if isBarnyard:
        # set counts and genes to max of human and mouse counts
        all_cells = all_cells.with_columns(
            counts=pl.max_horizontal(cs.ends_with("_counts")), genes=pl.max_horizontal(cs.ends_with("_genes"))
        )
    return all_cells


def build_allCells(sample_specific_file_paths, genes, barcodes, isBarnyard):
    """
    Builds the allCells.csv file for this sample.

    Args:
        sample_specific_file_paths (Dict[str, Path]):
            Dictionary where each key is a custom file identifier
            and each value is a path to the identified file in the STARsolo output directory
        genes (DataFrame): A dataframe of genes
        barcodes (DataFrame): A dataframe of barcodes
        isBarnyard (bool): Whether or not this is a barnyard sample (i.e. mixed mouse and human)

    Returns:
        An allCells.csv containing metrics computed across all barcodes in this sample (DataFrame)
    """
    all_cells = count_transcripts(sample_specific_file_paths, genes, barcodes, isBarnyard)
    stats_cols = [
        "CB",
        "cbMatch",
        "exonic",
        "intronic",
        "exonicAS",
        "countedU",
        "countedM",
        "mito",
        "genomeU",
        "genomeM",
        "featureU",
        "featureM",
        "nUMIunique",
    ]
    read_stats = pl.read_csv(
        sample_specific_file_paths["stats"], separator="\t", skip_rows_after_header=1, columns=stats_cols
    )  # skip first row of stats file that has row for CBnotInPasslist

    # Compute metrics for each barcode
    # See https://github.com/alexdobin/STAR/discussions/1826#discussioncomment-5596619 for a description of the metrics
    read_stats = read_stats.with_columns(
        (pl.col("genomeU") + pl.col("genomeM")).alias("mappedReads"),
        (pl.col("featureU") + pl.col("featureM")).alias("geneReads"),
        (pl.col("countedU") + pl.col("countedM")).alias("countedReads"),
        (pl.lit(1) - pl.col("nUMIunique") / pl.col("countedU")).round(3).fill_nan(None).alias("Saturation"),
        (pl.col("mito") / (pl.col("genomeU") + pl.col("genomeM"))).round(3).fill_nan(None).alias("mitoProp"),
    )
    read_stats = read_stats.rename(
        {
            "CB": "cell_id",
            "cbMatch": "totalReads",
            "exonic": "exonReads",
            "intronic": "intronReads",
            "exonicAS": "antisenseReads",
            "countedM": "countedMultiGeneReads",
            "mito": "mitoReads",
        }
    )
    # polars 0.20.* does not support drop with strict=False
    # polars 1.* will error if a column is not found
    stats_cols = [col for col in stats_cols if col in read_stats.columns]
    read_stats = read_stats.drop(stats_cols)
    all_cells = all_cells.join(read_stats, on="cell_id", how="left")
    return all_cells.select(
        pl.col("cell_id"),
        cs.ends_with("counts"),
        cs.ends_with("genes"),
        pl.col("totalReads"),
        pl.col("countedReads"),
        pl.col("mappedReads"),
        pl.col("geneReads"),
        pl.col("exonReads"),
        pl.col("intronReads"),
        pl.col("antisenseReads"),
        pl.col("mitoReads"),
        pl.col("countedMultiGeneReads"),
        pl.col("Saturation"),
        pl.col("mitoProp"),
    )


def split_barcodes(lib_struct: Path, all_cells: pl.DataFrame, is_merge: bool) -> pd.DataFrame:
    """
    Splits cell barcodes into constituent aliases and add those columns to allCells.

    Args:
        libStruct: Path to the library structure json for this sample
        allCells: An allCells.csv containing metrics computed across all barcodes in this sample

    Returns:
        allCells DataFrame concatenated with the alias columns
    """
    lib_json = io.readJSON(lib_struct, preserveDictOrder=True)
    barcode_info = lib_json["barcodes"]
    aliases = [bc["alias"] for bc in barcode_info if bc.get("type", None) not in ["library_index", "umi", "target"]]
    barcodes = all_cells["cell_id"]
    # If the data is from samples that were merged.
    # The sample ID needs to be removed from the CB before mapping to an alias.
    if is_merge:
        barcodes = barcodes.str.split("_").list.first()
    barcodes = barcodes.str.split("+").list.to_struct(fields=aliases).struct.unnest()
    bead_bc_cols = ["bead1", "bead2", "bead3"]
    if all([bc in aliases for bc in bead_bc_cols]):
        # create concatenated bead barcode column
        barcodes = barcodes.with_columns(pl.concat_str(bead_bc_cols, separator="+").alias("bead_bc"))
        bead_bc_cols_to_drop = [
            bc.get("alias") for bc in barcode_info if bc.get("alias") in bead_bc_cols and not bc.get("plate")
        ]
        barcodes = barcodes.drop(bead_bc_cols_to_drop)
    return all_cells.with_columns(barcodes)


def main():
    parser = argparse.ArgumentParser()

    # Required and optional arguments for specifying the STARsolo outputs for this sample
    parser.add_argument(
        "--STARsolo_out", type=Path, required=True, help="Path to the STARsolo outputs for this sample."
    )
    parser.add_argument(
        "--feature_type", type=str, required=False, default="GeneFull_Ex50pAS", help="STARsolo feature type used."
    )
    parser.add_argument(
        "--matrix_type",
        type=str,
        required=False,
        default="UniqueAndMult-PropUnique.mtx",
        help="STARsolo matrix type used.",
    )
    parser.add_argument(
        "--isBarnyard",
        default=False,
        action="store_true",
        help="If set, this sample will be interpreted as a barnyard sample (i.e. mixed mouse and human).",
    )
    parser.add_argument("--threads", type=int, required=False, default=1, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="8 GB", help="Memory allocated to task")

    # Optional argument to specify the name of the sample for which cells are being called
    parser.add_argument(
        "--sample", type=str, required=False, default="example", help="Unique string to identify this sample."
    )

    # Optional argument to specify the library structure for this sample
    parser.add_argument(
        "--libStruct", type=Path, required=False, help="Path to the library structure json for this sample."
    )

    # Optional argument to specify whether this is a merged workflow
    parser.add_argument(
        "--isMerge",
        default=False,
        action="store_true",
        help="If set, workflow is merged, and sample names will be extracted from the group name.",
    )

    args = parser.parse_args()
    mem_limit, mem_unit = args.memory.split()
    mem_limit = (
        f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"  # allocate duckdb memory based on number of threads
    )
    duckdb.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    """
    )

    sample_specific_file_paths = io.resolve_sample_specific_file_paths(
        args.STARsolo_out, args.feature_type, args.matrix_type
    )
    genes = read_genes(sample_specific_file_paths)
    barcodes = read_barcodes(sample_specific_file_paths)
    allCells = build_allCells(sample_specific_file_paths, genes, barcodes, args.isBarnyard)

    if args.libStruct is not None:
        allCells = split_barcodes(args.libStruct, allCells, args.isMerge)
    # When this is a merged sample, the sample name is actually a "group" name;
    # in order to get the correct sample name we extract it from the barcode
    if args.isMerge:
        allCells = allCells.with_columns(pl.col("cell_id").str.split("_").list.last().alias("sample"))
    else:
        allCells = allCells.with_columns(pl.lit(args.sample).alias("sample"))

    metricsDir = Path(".", f"{args.sample}_metrics")
    metricsDir.mkdir(parents=True)

    # Write allBarcodes.parquet for this sample
    allCells.write_parquet(f"{metricsDir}/{args.sample}_allBarcodes.parquet")


if __name__ == "__main__":
    main()
