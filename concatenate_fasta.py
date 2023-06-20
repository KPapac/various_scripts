#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pyfasta import Fasta

# Parser stuff
parser = argparse.ArgumentParser(description='Concatenates the contigs of a draft assembly into one contig.')
parser.add_argument("input_fasta", help='Path to the draft assembly, whose contigs I will concatenate.')
parser.add_argument("header_fasta", help='Header to use for the fasta output.')
parser.add_argument("output_fasta", help='Path at which to write the concatenated sequence.')

args = parser.parse_args()

# Code
f=Fasta(args.input_fasta)

contig_names=sorted(f.keys())
list_of_contigs = [ str(f[contig_header]) for contig_header in contig_names ]
concatenated_sequence = ''.join(list_of_contigs)

with open(args.output_fasta, 'w') as fout:
    fout.write(f">{args.header_fasta}\n{concatenated_sequence}")

