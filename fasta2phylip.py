#!/usr/bin/env python3
import argparse
from Bio import AlignIO


def convert_fasta_to_phylip(input_file, output_file, strict=False):
    try:
        fmt = "phylip-sequential" if strict else "phylip-relaxed"
        alignment = AlignIO.read(input_file, "fasta")
        AlignIO.write(alignment, output_file, fmt)
        print(
            f"Converted to non-interleaved {'strict' if strict else 'relaxed'} PHYLIP: '{output_file}'"
        )
    except Exception as e:
        print(f"Error during conversion: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert aligned FASTA to non-interleaved PHYLIP format."
    )
    parser.add_argument("input", help="Path to the input aligned FASTA file")
    parser.add_argument("output", help="Path to the output PHYLIP file")
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="Use strict sequential PHYLIP format (10-character max sequence IDs)",
    )
    args = parser.parse_args()
    convert_fasta_to_phylip(args.input, args.output, strict=args.sequential)


if __name__ == "__main__":
    main()
