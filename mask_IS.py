#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import argparse


def replace_sequence_with_N(sequence, pos_start, pos_end):
    length_of_masking = pos_end - pos_start + 1
    return Seq(sequence[: pos_start - 1] + 'N' * length_of_masking + sequence[pos_end:])


def main():
    parser = argparse.ArgumentParser(
        description='Masks a genome given a coordinates file. The coordinates file should have the form: "Start_position,End_position".'
    )
    parser.add_argument('path_to_genome', type=str, help='Path to genome to mask.')
    parser.add_argument(
        'path_to_coordinates', type=str, help='Path to coordinates to mask.'
    )
    parser.add_argument(
        'path_to_masked_genome', type=str, help='Path to output masked genome.'
    )
    args = parser.parse_args()
    # Get list of coords for IS starting and ending positions
    list_of_spans_to_mask = []
    with open(args.path_to_coordinates, "r") as file:
        for line in file:
            # Split the line into columns using ',' as the delimiter
            columns = line.strip().split(',')
            start = int(columns[0])
            end = int(columns[1])
            list_of_spans_to_mask.append((start, end))
    # Read the input FASTA file and create a new FASTA file with replacements
    with open(args.path_to_masked_genome, "w") as output_handle:
        input_genome = SeqIO.read(args.path_to_genome, "fasta")
        for mask_between_positions in list_of_spans_to_mask:
            (mask_start, mask_end) = mask_between_positions
            new_sequence = replace_sequence_with_N(
                str(input_genome.seq), mask_start, mask_end
            )
            input_genome.seq = new_sequence
        SeqIO.write(input_genome, output_handle, "fasta")


if __name__ == "__main__":
    main()
