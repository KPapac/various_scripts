#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def create_dotplot(
    tsv_file,
    x_column,
    y_column,
    output_file,
    delimiter="\t",
    log_x=False,
    log_y=False,
    x_range=None,
    y_range=None,
):
    # Load the file into a pandas DataFrame
    try:
        # Enable reading large files in chunks if necessary
        chunk_size = 10 ** 6  # Define a reasonable chunk size
        data_chunks = pd.read_csv(tsv_file, sep=delimiter, chunksize=chunk_size)
        data = pd.concat(data_chunks, ignore_index=True)
    except Exception as e:
        print(f"Error reading the file: {e}")
        return

    # Check if specified columns exist
    if x_column not in data.columns or y_column not in data.columns:
        print(
            f"Error: One or both columns '{x_column}' and '{y_column}' are not in the file."
        )
        print(f"Available columns: {', '.join(data.columns)}")
        print(
            "Hint: Use a tool like `head` or `pandas.read_csv` with a small sample to inspect the file and check available columns."
        )
        return

    # Apply logarithmic transformation if requested
    if log_x:
        data[x_column] = (
            data[x_column].apply( lambda x: pd.NA if x <= 0 else x).dropna().apply(
                lambda x: np.log(x)
            )
        )
    if log_y:
        data[y_column] = (
            data[y_column].apply( lambda y: pd.NA if y <= 0 else y).dropna().apply(
                lambda y: np.log(y)
            )
        )
    # Create the dot plot
    try:
        plt.figure(figsize=(8, 6))
        plt.scatter(data[x_column], data[y_column], alpha=0.7)
        plt.xlabel(f"{'log(' + x_column + ')' if log_x else x_column}")
        plt.ylabel(f"{'log(' + y_column + ')' if log_y else y_column}")
        plt.title("Dot Plot")
        plt.grid(True, linestyle="--", alpha=0.6)
        # Set axis ranges if provided
        if x_range:
            plt.xlim(x_range)
        if y_range:
            plt.ylim(y_range)
        # Save the plot
        plt.savefig(output_file)
        print(f"Dot plot saved as '{output_file}'")
        plt.show()
    except Exception as e:
        print(f"Error creating the plot: {e}")


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Create a dot plot from a file.")
    parser.add_argument("tsv_file", type=str, help="Path to the input file.")
    parser.add_argument("x_column", type=str, help="Column name to use for the x-axis.")
    parser.add_argument("y_column", type=str, help="Column name to use for the y-axis.")
    parser.add_argument(
        "output_file",
        type=str,
        nargs="?",
        default="dotplot.png",
        help="Path to save the output plot (default: 'dotplot.png').",
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default="\t",
        help="Delimiter used in the input file (default: '\t').",
    )
    parser.add_argument(
        "--log_x",
        action="store_true",
        help="Apply a logarithmic scale to the x-axis column.",
    )
    parser.add_argument(
        "--log_y",
        action="store_true",
        help="Apply a logarithmic scale to the y-axis column.",
    )
    parser.add_argument(
        "--x_range",
        type=float,
        nargs=2,
        metavar=("X_MIN", "X_MAX"),
        help="Specify the range of the x-axis as two float values.",
    )
    parser.add_argument(
        "--y_range",
        type=float,
        nargs=2,
        metavar=("Y_MIN", "Y_MAX"),
        help="Specify the range of the y-axis as two float values.",
    )
    args = parser.parse_args()
    # Call the function with parsed arguments
    create_dotplot(
        args.tsv_file,
        args.x_column,
        args.y_column,
        args.output_file,
        args.delimiter,
        args.log_x,
        args.log_y,
        args.x_range,
        args.y_range,
    )
