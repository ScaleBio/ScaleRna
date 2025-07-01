#!/usr/bin/env python
"""
Across all cell-barcode combinations, calculate the normalized KL divergence for each bead.
Output the bead barcodes below minimum divergence threshold to a CSV file.
"""
import argparse
import duckdb
from pathlib import Path
import pandas as pd


def filter_beads(barcodes_dir: Path, min_utc: int, sample_barcode: str):
    """Use observed distribution of counts across RT barcodes to calculate KL divergence for each bead.

    The KL divergence (relative entropy) is calculated between the distribution of counts across RT
    for the whole library (Q) and the distribution of counts across RT for each bead (P).

    https://en.wikipedia.org/wiki/Kullback-Leibler_divergence
    KL Divergence: D_KL(P || Q) = Σ P(x) * log(P(x) / Q(x))
    Measures how one probability distribution P diverges from a second, expected distribution Q.

    Low KL divergence for a bead indicates that the distribution of counts across RT barcodes matches
    the overall dataset (more uniform) rather than being concentrated among one or a few RT barcodes.
    """
    duckdb.sql(
        f"""
    -- 1. Get the counts for each bead and {sample_barcode}
    WITH cell_metrics AS (
        SELECT counts, totalReads, {sample_barcode}, bead_bc
        FROM '{barcodes_dir / "*allBarcodes.parquet"}'
    ),
    -- 2. Get the total counts for each {sample_barcode}
    lib_{sample_barcode}_totals AS (
        SELECT {sample_barcode}, SUM(counts) AS total_counts
        FROM cell_metrics
        GROUP BY {sample_barcode}
    ),
    -- 3. Get the total number of {sample_barcode}s
    num_{sample_barcode}s AS (
        SELECT COUNT({sample_barcode}) AS count FROM lib_{sample_barcode}_totals
    ),
    -- 4. Get the total counts for all {sample_barcode}s in library
    lib_counts AS(
        SELECT SUM(total_counts) as total FROM lib_{sample_barcode}_totals
    ),
    -- 5. Calculate the probability of each {sample_barcode} in the library (q in KL divergence)
    lib_{sample_barcode}_probs AS (
        SELECT
          lpt.{sample_barcode},
          lpt.total_counts * 1.0 / lc.total as q
        FROM lib_{sample_barcode}_totals lpt, lib_counts lc
    ),
    -- 6. Get the total counts for each bead with more than one RT barcode above minUTC
    bead_totals AS (
        SELECT bead_bc, SUM(counts) AS total_counts, SUM(totalReads) AS total_reads
        FROM cell_metrics
        GROUP BY bead_bc
        HAVING COUNT(DISTINCT CASE WHEN counts >= {min_utc} THEN {sample_barcode} END) > 1
    ),
    -- 7. Calculate the probability of each {sample_barcode} in the bead (p in KL divergence)
    bead_{sample_barcode}_probs AS (
        SELECT
            cm.bead_bc,
            cm.{sample_barcode},
            cm.counts * 1.0 / bt.total_counts as p
        FROM cell_metrics cm
        JOIN bead_totals bt ON cm.bead_bc = bt.bead_bc
    ),
    -- 8. Calculate the KL divergence for each bead
    kl_divergence AS (
      SELECT
        bpp.bead_bc,
        SUM(CASE
          WHEN bpp.p > 0 AND lpp.q > 0 THEN bpp.p * LOG2(bpp.p / lpp.q)
          ELSE 0
        END) AS kl_div
      FROM bead_{sample_barcode}_probs bpp
      JOIN lib_{sample_barcode}_probs lpp ON bpp.{sample_barcode} = lpp.{sample_barcode}
      GROUP BY bpp.bead_bc
    ),
    -- 9. Normalize the KL divergence by the number of {sample_barcode}s and calculate the ambient score
    score AS(
        SELECT
            kld.bead_bc,
            kld.kl_div / LOG2(num_{sample_barcode}s.count) AS kl_norm,
            bt.total_counts AS bead_total,
            bt.total_reads AS bead_reads
        FROM kl_divergence kld
        JOIN bead_totals bt ON kld.bead_bc = bt.bead_bc
        CROSS JOIN num_{sample_barcode}s
    )
    -- 10. Classify beads as ambient with kl_norm below the threshold if there are at least 8 RT barcodes in library
    SELECT bead_bc, kl_norm, bead_total, bead_reads
    FROM score
    CROSS JOIN num_{sample_barcode}s
    WHERE num_{sample_barcode}s.count >= 8;
    """
    ).df().astype(
        {
            # In case dataframe is empty, need to specify bead_bc is VARCHAR
            "bead_bc": pd.StringDtype(storage="pyarrow"),
        }
    ).to_parquet(
        barcodes_dir / "bead_scores.parquet", index=False
    )


def main():
    parser = argparse.ArgumentParser("Output bead KL divergence scores to a CSV file")
    parser.add_argument(
        "--barcodesDir",
        required=True,
        type=Path,
        help="Directory containing allBarcodes.parquet files for a single library",
    )
    parser.add_argument(
        "--minUTC",
        type=int,
        required=False,
        default=100,
        help="Minimum number of counts to consider a barcode as a potential cell",
    )
    parser.add_argument(
        "--sampleBarcode", type=str, required=False, default="RT", help="Sample barcode name used in assay"
    )
    parser.add_argument("--threads", type=int, required=False, default=1, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="8 GB", help="Memory allocated to task")

    args = parser.parse_args()
    mem_limit, mem_unit = args.memory.split()
    # allocate duckdb memory based on number of threads
    mem_limit = f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"
    duckdb.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    """
    )

    filter_beads(
        barcodes_dir=args.barcodesDir,
        min_utc=args.minUTC,
        sample_barcode=args.sampleBarcode,
    )


if __name__ == "__main__":
    main()
