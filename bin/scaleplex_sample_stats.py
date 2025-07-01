#!/usr/bin/env python
"""Create scaleplex_stats.csv and metrics.csv files from the output of ScalePlex assignment"""
import argparse
from pathlib import Path
from dataclasses import dataclass
import polars as pl
from scaleplex_assignment import AssignmentCodes


@dataclass
class Metrics:
    """Overall sample statistics"""

    meanReadsPerCell: int = 0  # number of reads per cell. Should be multiplied by correct for ~usable amount
    nUMIPerCell: int = 0  # median number of HASH UMI molecules per cell
    passingPercent: int = 0  # percent of cells that had at least thresh Guide UMI's detected
    maxFailPercent: int = 0  # percent of cells with maxFail assignment
    enrichFailPercent: int = 0  # percent of cells with enrichFail assignment
    indeterminatePercent: int = 0
    unexpectedPercent: int = 0
    saturation: int = 0  # median saturation of cells
    readsInCells: int = 0  # percent of hash reads that ended up in called cells

    def print(self, outFn: Path):
        with open(outFn, "w") as out:
            print("ReadsPerCell", f"{self.meanReadsPerCell:}", sep=",", file=out)
            print("nUMIPerCell", f"{self.nUMIPerCell:}", sep=",", file=out)
            print("passingPercent", f"{self.passingPercent:.1%}", sep=",", file=out)
            print("maxFailPercent", f"{self.maxFailPercent:.1%}", sep=",", file=out)
            print("enrichFailPercent", f"{self.enrichFailPercent:.1%}", sep=",", file=out)
            print("indeterminatePercent", f"{self.indeterminatePercent:.1%}", sep=",", file=out)
            print("unexpectedPercent", f"{self.unexpectedPercent:.1%}", sep=",", file=out)
            print("saturation", f"{self.saturation:.1%}", sep=",", file=out)
            print("readsInCells", f"{self.readsInCells:.1%}", sep=",", file=out)


codes = AssignmentCodes()


def main(
    all_cells: list[Path],
    cell_stats: list[Path],
    id: str,
    out_dir: Path,
):
    metrics = Metrics()
    rna = pl.scan_csv(all_cells).collect().to_pandas().set_index("cell_id")
    cell_stats_passing = (
        pl.scan_parquet(cell_stats)
        .filter(pl.col("pass"))
        .select(["totalReads", "counts", "Saturation", "assigned_scaleplex"])
        .collect()
        .to_pandas()
    )

    total_reads = pl.scan_parquet(cell_stats).select("totalReads").sum().collect()[0, 0]

    # Create output directories
    out_dir.mkdir(parents=True, exist_ok=True)
    assigned_cells = rna[~rna["assigned_scaleplex"].isin(codes.errors)]  # remove cells with assignment errors
    scaleplex_stats = (
        assigned_cells.groupby("assigned_scaleplex")
        .agg(
            passing_reads=("totalReads", "sum"),
            cells_called=("assigned_scaleplex", "size"),
            mean_passing_reads=("totalReads", "mean"),
            median_utc=("counts", "median"),
            median_genes=("genes", "median"),
        )
        .sort_values(["cells_called", "passing_reads"], ascending=False)
        .astype("int")
    )
    scaleplex_stats.to_csv(out_dir / f"{id}.scaleplex_stats.csv")

    metrics.meanReadsPerCell = round(cell_stats_passing["totalReads"].mean(0), 1)
    metrics.nUMIPerCell = round(cell_stats_passing["counts"].median(0), 1)
    metrics.saturation = round(cell_stats_passing["Saturation"].median(0), 1)
    metrics.readsInCells = cell_stats_passing["totalReads"].sum() / total_reads

    assignment_prop = cell_stats_passing["assigned_scaleplex"].value_counts(normalize=True)
    metrics.passingPercent = assignment_prop.drop(codes.errors, errors="ignore").sum()
    # Report out different assignment error stats
    metrics.indeterminatePercent = assignment_prop[codes.indeterminate] if codes.indeterminate in assignment_prop else 0
    metrics.maxFailPercent = assignment_prop[codes.max_fail] if codes.max_fail in assignment_prop else 0
    metrics.enrichFailPercent = assignment_prop[codes.enrich_fail] if codes.enrich_fail in assignment_prop else 0
    metrics.unexpectedPercent = assignment_prop[codes.unexpected] if codes.unexpected in assignment_prop else 0

    metrics.print(out_dir / "metrics.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create ScalePlex sample statistics")
    parser.add_argument(
        "--all_cells",
        nargs="+",
        type=Path,
        help="Passing RNA cells with ScalePlex assignment",
    )
    parser.add_argument(
        "--cell_stats",
        nargs="+",
        type=Path,
        help="Metrics on ScalePlex cell-barcodes",
    )
    parser.add_argument("--id", required=True, help="Sample and lib name")
    parser.add_argument("--outDir", type=Path, help="Directory for outputs")
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(
        args.all_cells,
        args.cell_stats,
        args.id,
        args.outDir,
    )
