#!/usr/bin/env python3
import argparse
import sys
from Bio import AlignIO
from collections import Counter
from itertools import combinations
from io import StringIO


def load_alignment(file_path=None):
    """Load a FASTA alignment from a file or stdin."""
    try:
        if file_path:
            alignment = AlignIO.read(file_path, "fasta")
        else:
            # Read from stdin
            fasta_data = sys.stdin.read()
            if not fasta_data.strip():
                raise ValueError("No input received from stdin.")

            alignment = AlignIO.read(StringIO(fasta_data), "fasta")
        return alignment

    except Exception as e:
        print(f"Error loading alignment: {e}")
        exit(1)


def calculate_matching_sites(alignment):
    """Calculate the number of fully matching columns in the alignment."""
    matching_sites = 0
    gaps = 0
    alignment_length = alignment.get_alignment_length()
    for i in range(alignment_length):
        column = alignment[:, i]
        column_counts = Counter(column)
        if len(column_counts) == 1:
            matching_sites += 1
        if '-' in column:
            gaps += 1
    return matching_sites, gaps


def pairwise_identity(seq1, seq2):
    """Calculate pairwise identity between two sequences."""
    matches = sum(a == b for a, b in zip(seq1, seq2) if not (a == "-" and b == "-"))
    length = sum(not (a == "-" and b == "-") for a, b in zip(seq1, seq2))
    return (matches / length) * 100 if length > 0 else 0


def calculate_pairwise_identities(alignment):
    """Calculate pairwise identities between sequences."""
    identities = []
    for seq1, seq2 in combinations(alignment, 2):
        identity = pairwise_identity(seq1.seq, seq2.seq)
        identities.append((seq1.id, seq2.id, identity))
    return identities


def print_alignment_stats(alignment, alignment_name, short_output):
    """Print alignment statistics to stdout."""
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    matching_sites = calculate_matching_sites(alignment)[0]
    gaps = calculate_matching_sites(alignment)[1]
    pairwise_identities = calculate_pairwise_identities(alignment)
    if short_output is True:
        print(
            f"{alignment_name} {matching_sites} {alignment_length} {matching_sites/alignment_length*100:.2f}% {num_sequences}"
        )
        return 0

    print("============================================")
    print("=========== Alignment Statistics ===========")
    print(f"{alignment_name}")
    print(f"Number of sequences: {num_sequences}")
    print(f"Alignment length: {alignment_length}")
    print(f"Number of matching sites: {matching_sites}")
    print(f"Sites containing gaps: {gaps}")
    print(f"Sites mismatching: {alignment_length-matching_sites-gaps}")
    print(f"Percent of matching sites: {matching_sites/alignment_length*100:.2f}%")
    print("\n=== Pairwise Identities ===")
    for seq1_id, seq2_id, identity in pairwise_identities:
        print(f"{seq1_id} ↔ {seq2_id}: {identity:.2f}%")
    print("============================================\n")


def main():
    """Main function to handle arguments and run analysis."""
    parser = argparse.ArgumentParser(
        description="Calculate statistics from a FASTA alignment file or stdin."
    )
    parser.add_argument(
        "alignment_file",
        nargs="?",
        help="Path to the FASTA alignment file (optional, defaults to stdin if not provided).",
    )
    parser.add_argument(
        "-s",
        "--short_output",
        action="store_true",
        help="Prints just the number of sites matching and the total alignment length.",
    )
    args = parser.parse_args()
    alignment = load_alignment(args.alignment_file)
    print_alignment_stats(alignment, args.alignment_file, args.short_output)


if __name__ == "__main__":
    main()
