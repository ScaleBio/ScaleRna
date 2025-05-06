#!/usr/bin/env python
"""
Concatenate barcode stats from all libraries into a single file
"""
import argparse
import re
import pandas as pd
from pathlib import Path


def concat_library_metrics(type_level_matches: list, overall_matches: list, lib_names: list, output_fname: str):
    """
    Concatenate barcode stats from all libraries into a single file

    Args:
        type_level_matches: Filenames containing barcode type level match metrics (ambiguous, exact, corrected, etc)
        overall_matches: Filenames containing barcode overall match metrics (pass, barcodeerror, tooshorterror, etc)
        lib_names: Names of the libraries
        output_fname: Output filename that will contain concatenated metrics
    
    Output:
        Writes csv file containing concatenated metrics
    """
    dataframes = []

    for type_level_match, overall_match, lib_name in zip(type_level_matches, overall_matches, lib_names):
        df_type_level_match = pd.read_csv(type_level_match)
        df_type_level_match['Status'] = df_type_level_match.apply(lambda row: f"{row['Barcode']}-{row['Match']}", axis=1)
        df_type_level_match = df_type_level_match.pivot_table(index=None, columns='Status', values='Reads', aggfunc='sum')
        
        df_overall_match = pd.read_csv(overall_match).drop(columns=['Reads']).rename(columns={'Percent': 'Reads'})
        df_overall_match = df_overall_match.pivot_table(index=None, columns='Status', values='Reads', aggfunc='sum')
        
        df_type_level_match['Status'] = 'Reads'
        df_overall_match['Status'] = 'Reads'
        
        df = pd.merge(df_type_level_match, df_overall_match, on='Status').drop(columns=['Status'])
        df.insert(0, 'Library', lib_name)
        
        dataframes.append(df)
    
    pd.concat(dataframes, ignore_index=True).to_csv(output_fname, index=False)


def main():
    parser = argparse.ArgumentParser(description="Generate metrics for library report")
    parser.add_argument(
        "--typeLevelMatches", nargs="+", type=str, required=True, help="Per library type level match metrics that need to be concatenated"
    )
    parser.add_argument(
        "--overallMatches", nargs="+", type=str, required=True, help="Per library overall match metrics that need to be concatenated"
    )
    parser.add_argument("--output_fname", type=Path, default="allLibraries.barcodeStats.csv", help="Output filename that will contain concatenated metrics")

    args = parser.parse_args()
    lib_names = []
    
    for fname in args.typeLevelMatches:
        lib_names.append(re.sub(r'^library_|\.typeLevelMatches\.csv$', '', fname))

    concat_library_metrics(args.typeLevelMatches, args.overallMatches, lib_names, args.output_fname)


if __name__ == "__main__":
    main()
