#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import GenBank

def summarise_GBK(gbkFile):
    feature_summary=dict()
    with open(gbkFile, 'r') as handle:
        for record in GenBank.parse(handle):
            for feature in record.features:
                if feature.key not in feature_summary.keys():
                    feature_summary[feature.key]=1 
                else: feature_summary[feature.key]+=1
    return f"{gbkFile}: {feature_summary}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse GBK files and optionally write the output to a file.")
    parser.add_argument("input_files", nargs="+", help="Input GBK files to parse")
    parser.add_argument("-o", "--output_file", help="Path to the output file (optional)", default=None)

    args = parser.parse_args()

    summary=''
    for gbk_file in args.input_files:
        summary+=summarise_GBK(gbk_file) + '\n'
    if args.output_file == None: print(summary)
    else: 
        with open(args.output_file, 'w') as f:
            f.write(summary)