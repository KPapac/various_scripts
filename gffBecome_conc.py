#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('gff_file', type=str, help='')
    parser.add_argument('output', type=str, help='')
    args = parser.parse_args()
    return args

def merge_gff(path_to_gff, output_path):
    # Sorting gff
    gff = pd.read_csv(path_to_gff, sep = '\t', comment='#', header=None)
    gff['contig_number'] = gff[0].str.extract(r'(\d+)').astype(int)
    gff = gff.sort_values(['contig_number',3])
    comments=''

    # Getting contig sizes to a dictionary
    contig_sizes={}
    length=0
    with open(path_to_gff) as f:
        for line in f.readlines():
            if '#' not in line: break
            else:
                # Keeping comments
                comments+=line
                match = re.match(r'##sequence-region (\S+) \d+ (\d+)', line)
                if match:
                    contig_name = match.group(1)
                    contig_sizes[contig_name]=length
                    length += int(match.group(2))

    # Appending contig sizes to dataframe
    gff['contig_sizes_to_add']=gff[0].map(contig_sizes)
    # Adding previous coordinates
    gff[3] += gff['contig_sizes_to_add']
    gff[4] += gff['contig_sizes_to_add']
    gff=gff.drop(['contig_number', 'contig_sizes_to_add'], axis=1)
    with open(output_path, 'w') as f:
        f.write(comments)
        gff.to_csv(f, sep='\t', index=False, header=False)

def main():
    # Parsing arguments
    args = parse_arguments()
    input_gff= args.gff_file
    output_path = args.output
    merge_gff(input_gff,output_path)

if __name__ == "__main__":
    main()
