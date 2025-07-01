#!/usr/bin/env python
"""
Concatenate columns (cells) for gene expression matrices in Matrix Market format, into single output file
"""

import duckdb
import argparse
import gzip
from pathlib import Path
import shutil


def concat_mtx(dirs: list[Path], output_dir: Path, con: duckdb.DuckDBPyConnection) -> None:
    """
    Concatenate sparse matrices in Matrix Market format

    Args:
        files: List of filtered mtx directories to concatenate
        output_dir: Output path
        con: duckdb connection
    """
    nnz = 0
    row_size = []
    col_size = []
    for idx, mtx_dir in enumerate(dirs):
        mtx_file = Path(mtx_dir) / "matrix.mtx.gz"
        # Open mtx file to get dimensions and check if it represents an empty mtx
        with gzip.open(mtx_file, "rt") as f:
            # https://math.nist.gov/MatrixMarket/formats.html
            f.readline()  # %%MatrixMarket matrix coordinate real general
            f.readline()  # %
            dims = f.readline().split()  # 36910 185326 0000262179120
            nnz += int(dims[2])
            row_size.append(int(dims[0]))
            col_size.append(int(dims[1]))
        if col_size[idx] == 0:
            # Skip matrices with no cells
            continue
        con.sql(
            f"""
        INSERT INTO mtx (
            SELECT gene, barcode + {sum(col_size[:idx])} as barcode, count
            FROM read_csv(
                '{mtx_file}',
                delim=' ',
                skip={3},
                columns = {{
                    'gene': 'UINTEGER',
                    'barcode': 'UINTEGER',
                    'count': 'FLOAT'
                }})
        );
        """
        )
    if len(set(row_size)) != 1:
        raise AssertionError("Row dimensions must match to concatenate Matrix Market files by column")
    con.sql(
        f"""
    COPY mtx
    TO '{output_dir / "matrix.mtx.gz"}' (
        FORMAT csv,
        DELIMITER ' ',
        HEADER false,
        PREFIX '%%MatrixMarket matrix coordinate real general\n%\n{row_size[0]} {sum(col_size)} {str(nnz).zfill(13)}\n',
        SUFFIX '\n'
    );
    """
    )

    # copy barcodes and features tsv files
    with gzip.open(output_dir / "barcodes.tsv.gz", "wt") as barcodes_file:
        for mtx_dir in dirs:
            with gzip.open(mtx_dir / "barcodes.tsv.gz", "rt") as f:
                barcodes_file.write(f.read())
    with gzip.open(dirs[0] / "features.tsv.gz", "rb") as f_in, gzip.open(output_dir / "features.tsv.gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def main():
    parser = argparse.ArgumentParser("Combine Matrix Market files into single output file, stack columns")
    parser.add_argument("inputDirs", nargs="+", type=Path, help="Filtered matrix directories to merge")
    parser.add_argument("--id", type=str, help="Sample ID")
    parser.add_argument("--threads", type=int, required=False, default=2, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="7 GB", help="Memory allocated to task")
    args = parser.parse_args()

    output_dir = Path(".") / f"{args.id}_filtered_star_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    mem_limit, mem_unit = args.memory.split()
    # allocate duckdb memory based on number of threads
    mem_limit = f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"
    # create persistent duckdb connection
    con = duckdb.connect("tables.db")
    con.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    CREATE OR REPLACE TABLE mtx (
        gene UINTEGER,
        barcode UINTEGER,
        count FLOAT
    );
    """
    )
    concat_mtx(dirs=args.inputDirs, output_dir=output_dir, con=con)
    con.close()


if __name__ == "__main__":
    main()
