#!/usr/bin/env python3
import argparse
import pandas as pd
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Left join two tables using specified columns (by index). Supports custom delimiters and optional headers."
    )
    parser.add_argument("path1", help="Path to the first table")
    parser.add_argument("path2", help="Path to the second table")
    parser.add_argument(
        "delim1", help="Delimiter for the first table (e.g. ',' or '\\t')"
    )
    parser.add_argument("delim2", help="Delimiter for the second table")
    parser.add_argument(
        "join_col1", type=int, help="Join column index from table 1 (1-based)"
    )
    parser.add_argument(
        "join_col2", type=int, help="Join column index from table 2 (1-based)"
    )
    parser.add_argument(
        "has_header1", choices=["True", "False"], help="Does table 1 have a header?"
    )
    parser.add_argument(
        "has_header2", choices=["True", "False"], help="Does table 2 have a header?"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    join_col1 = args.join_col1 - 1
    join_col2 = args.join_col2 - 1
    has_header1 = args.has_header1 == "True"
    has_header2 = args.has_header2 == "True"
    header1 = 0 if has_header1 else None
    header2 = 0 if has_header2 else None
    # Read data using python engine to avoid regex warnings
    df1 = pd.read_csv(
        args.path1, delimiter=args.delim1, header=header1, engine='python'
    )
    df2 = pd.read_csv(
        args.path2, delimiter=args.delim2, header=header2, engine='python'
    )
    # If no header, assign unique column names
    if not has_header1:
        df1.columns = [f"X{i+1}" for i in range(len(df1.columns))]
    if not has_header2:
        df2.columns = [f"Y{i+1}" for i in range(len(df2.columns))]
    colname1 = df1.columns[join_col1]
    colname2 = df2.columns[join_col2]
    # Do the join by specifying left_on and right_on
    result = pd.merge(
        df1,
        df2,
        how='left',
        left_on=colname1,
        right_on=colname2,
        suffixes=('', '_right'),
    )
    # If the join columns were different, optionally drop the redundant one
    if colname1 != colname2:
        result = result.drop(columns=[colname2])
    # Output result to stdout with no headers/index
    result.to_csv(sys.stdout, sep='\t', index=False, header=False)


if __name__ == "__main__":
    main()
