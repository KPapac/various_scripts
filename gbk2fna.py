#!/usr/bin/env python3
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(prog='gbk2fna', description='Converts gbk to fna.')
parser.add_argument('gbk_filename')
parser.add_argument('fna_filename')
args = parser.parse_args()
with open(args.gbk_filename, "r") as gbk, open(args.fna_filename, "w") as fna:
    for seq_record in SeqIO.parse(gbk, "genbank"):
        fna.write(
            ">%s %s\n%s\n" % (seq_record.id, seq_record.description, seq_record.seq)
        )
