#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os.path import basename, splitext
import re
import copy
from Bio import GenBank
from Bio.GenBank.Record import Record as GBK_Record
import csv


def read_data_to_undo_concatenation(input_gbk,tsv_file_name):
    '''
    Reads the tsv made by gbk_m2sEntry.py to acquire data to 
    deconcatenate.
    '''
    deconc_data=[]
    with open('test','r') as handle:
        tsv_data = csv.reader(handle, delimiter="\t")
        for line in tsv_data:
            deconc_data.append(line)
        return deconc_data   

def split_gbk(input_gbk, prefix):
    '''
    Reads the concatenated gbk and extracts from the end, the 
    contig entries. So the concatenated gbk should be exhausted in the end.
    '''
    gbk_records=[]
    # Gets data base on which to split the gbk file
    tsv_data=read_data_to_undo_concatenation(input_gbk, prefix)

    # Reads concatenated gbk file.
    with open(input_gbk) as handle:
        conc_gbk=GenBank.read(handle)

        # Reading the data from last contig to first
        for contig_info in tsv_data[::-1]:
            # Makes a new gbk that will have a record for a contig in it.
            new_record=copy.deepcopy(conc_gbk)
            new_record.locus=contig_info[0]
            contig_length=int(contig_info[1])
            new_record.sequence=conc_gbk.sequence[int(conc_gbk.size)-contig_length:]
            new_record.size=str(contig_length)

            new_record.features=fix_feature_coords(conc_gbk, new_record, contig_length)

            # Removing sequence from conc_gbk that I got, until I exhaust the entry
            conc_gbk.sequence=conc_gbk.sequence[:int(conc_gbk.size)-contig_length]
            conc_gbk.size=str(int(conc_gbk.size)-contig_length)


            # Gather gbk entries for each contig to a list
            gbk_records.append(new_record)
            break
        # Raise error if the concatenate is not consumed.
        # if len(conc_gbk.sequence) != 0 and conc_gbk.size !=0:
        #     raise Exception("Something is wrong, parts of the gbk are not parsed.")
    return gbk_records

def fix_feature_coords(concatenated_gbk, contig_gbk, contig_length):
    '''
    Reads features from conc_entry, returns features with new coordinates 
    for new_entry.
    '''
    # Getting features for the contig
    contig_features=[]

    #Need to add +1, cause python counts from 0
    contig_positions=tuple(range(int(concatenated_gbk.size)-contig_length+1,int(concatenated_gbk.size)+1))

    for feature in concatenated_gbk.features:
        start, end = re.findall(r'\d+', feature.location)
        start, end = int(start), int(end)

        # Fixing coordinates of contig features 
        if start and end in contig_positions:
            new_start=start-min(contig_positions)+1
            new_end=end-min(contig_positions)+1

            if re.match('complement',feature.location):
                fixed_feature=GenBank.Record.Feature(key=feature.key, location=f'complement({new_start}..{new_end})')
            else: 
                fixed_feature=GenBank.Record.Feature(key=feature.key, location=f'{new_start}..{new_end}')
            # Keeps qualifiers from previous entry
            fixed_feature.qualifiers=feature.qualifiers

            contig_features.append(fixed_feature)
    return contig_features

def main(input_gbk):
    # Define basename for gbk file. Used for saving/reading tsv and 
    # writing locus in gbk.
    prefix=splitext(basename(input_gbk))[0]
    read_data_to_undo_concatenation(input_gbk, prefix)
    my_entries=split_gbk(input_gbk, prefix)
    for item in my_entries[::-1]:
        print(item)

main(input_gbk="./conc.gbk")
