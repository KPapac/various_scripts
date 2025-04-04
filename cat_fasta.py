#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def concatenate_fasta_sequences(
    input_handle, output_file=None, fasta_header="concatenated_sequence"
):
    # Read all sequences from the input FASTA file
    sequences = [str(record.seq) for record in SeqIO.parse(input_handle, "fasta")]
    # Concatenate all sequences into one
    concatenated_sequence = ''.join(sequences)
    # Format the sequence to 60 characters per line
    formatted_sequence = '\n'.join(
        [
            concatenated_sequence[i: i + 60]
            for i in range(0, len(concatenated_sequence), 60)
        ]
    )
    # Output to file or standard output
    if output_file:
        with open(output_file, "w") as output_handle:
            output_handle.write(f">{fasta_header}\n" + formatted_sequence + "\n")
    else:
        print(f">{fasta_header}\n" + formatted_sequence)


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate sequences from a multi-FASTA file and format the output with 60 characters per line."
    )
    parser.add_argument(
        "input_fasta",
        nargs="?",
        default="-",
        help="Path to the input FASTA file (or '-' for stdin)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the output file (optional, prints to stdout by default)",
    )
    parser.add_argument(
        "-f",
        "--fasta-header",
        default="concatenated_sequence",
        help="Custom FASTA header (default: 'concatenated_sequence')",
    )
    args = parser.parse_args()
    # Use stdin if input_fasta is '-'
    input_handle = sys.stdin if args.input_fasta == "-" else args.input_fasta
    concatenate_fasta_sequences(input_handle, args.output, args.fasta_header)


if __name__ == "__main__":
    main()
