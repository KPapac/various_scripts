#!/usr/bin/python3
import numpy as np
import pandas as pd
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description='...')
parser.add_argument('dS_file', help="A 2NG dN or dS file from codeML")
parser.add_argument('strain1', nargs='?', help="First strain")
parser.add_argument('--strain2', help="Second strain (optional)", default=None)
parser.add_argument('--multiple_strains', nargs='+', help="Specify multiple strains")
args = parser.parse_args()
if args.multiple_strains and (args.strain1 or args.strain2):
    parser.error("Cannot use --multiple_strains with strain1 or --strain2.")
if args.strain2 and not args.strain1:
    parser.error("strain1 must be provided if --strain2 is used.")


def read_data(data):
    with open(data) as f:
        lines = f.readlines()[1:]
        names = []
        matrix = []
        for line in lines:
            parts = line.split()
            names.append(parts[0])  # First part is the name
            row = [
                float(x) for x in parts[1:]
            ]  # Remaining parts are the pairwise comparisons
            matrix.append(
                row + [np.nan] * (len(lines) - len(row))
            )  # Pad rows with NaN to form a square matrix
        df = pd.DataFrame(matrix, index=names, columns=names)
        return df


def get_comparisons(df, name, name2=None):
    matches = df.index[df.index.str.contains(name, case=True, regex=True)]
    if len(matches) == 0:
        return f"ERROR: No matches found for pattern: {name}"

    elif len(matches) > 1:
        return f"ERROR: Multiple matches found: {list(matches)}"

    name = matches[0]
    if name2 == None:
        return pd.concat(
            [df.loc[name].dropna().to_frame(), df[name].dropna().to_frame()]
        )

    else:
        matches = df.index[df.index.str.contains(name2, case=True, regex=True)]
        if len(matches) == 0:
            return f"ERROR: No matches found for pattern: {name2}"

        elif len(matches) > 1:
            return f"ERROR: Multiple matches found: {list(matches)}"

        name2 = matches[0]
        if pd.isna(df.loc[name, name2]):
            if pd.isna(df.loc[name2, name]):
                print(f"No value found for comparison {name}-{name2}")
                sys.exit()
            return df.loc[name2, name]

        return df.loc[name, name2]


def main():
    df = read_data(args.dS_file)
    if args.multiple_strains != None:
        pairwise_comparisons = {}
        for i in range(0, len(args.multiple_strains) - 1):
            pairwise_comparisons[
                f"{args.multiple_strains[i]}-{args.multiple_strains[i+1]}"
            ] = get_comparisons(
                df, args.multiple_strains[i], args.multiple_strains[i + 1]
            )
        for key in pairwise_comparisons:
            print(f"{args.dS_file} {key} {pairwise_comparisons[key]}")
    else:
        pairwise_comparisons = get_comparisons(df, args.strain1, args.strain2)
        print(pairwise_comparisons)


if __name__ == "__main__":
    main()
