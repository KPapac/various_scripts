#!/usr/bin/env python3
from BCBio import GFF
from Bio import SeqIO
import sys
import os


def convert_gbk_to_gff(gbk_file, gff_file):
    # Parse all records in the GenBank file
    with open(gbk_file, "r") as input_handle:
        records = list(SeqIO.parse(input_handle, "genbank"))
    if not records:
        raise ValueError("No records found in the GenBank file.")

    # Write all records to the GFF3 file
    with open(gff_file, "w") as output_handle:
        GFF.write(records, output_handle)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gbk_to_gff3.py input.gbk output.gff")
        sys.exit(1)
    gbk_path = sys.argv[1]
    gff_path = sys.argv[2]
    if not os.path.exists(gbk_path):
        print(f"Error: GenBank file '{gbk_path}' not found.")
        sys.exit(1)
    try:
        convert_gbk_to_gff(gbk_path, gff_path)
        print(f"Conversion complete: '{gff_path}' created.")
    except Exception as e:
        print(f"Error during conversion: {e}")
