#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse 
from Bio import SeqIO

def get_entry_by_locus_tag(gbk_file, locus_tag):

    with open(gbk_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type=='CDS':
                    if "locus_tag" in feature.qualifiers and feature.qualifiers["locus_tag"][0] == locus_tag:
                        return (locus_tag, feature.qualifiers["product"][0])

    # If it loops through all the GBK and does not find your locus tag:
    return (locus_tag, None)

def main():
    parser = argparse.ArgumentParser(description="Search GenBank file for entries with a specified locus tag.")
    parser.add_argument("gbk_file", help="Path to the GenBank file")
    parser.add_argument("locus_tag", help="Locus tag to search for")

    args = parser.parse_args()

    gbk_file_path = args.gbk_file
    desired_locus_tag = args.locus_tag

    #gbk_file_path = "wMel.gbk" 
    #desired_locus_tag = "WD4002"

    print(get_entry_by_locus_tag(gbk_file_path, desired_locus_tag))

if __name__ == "__main__":
    main()

