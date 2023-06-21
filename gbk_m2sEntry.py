#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from os.path import basename, splitext
import re
from Bio import GenBank
import copy
import argparse

def parse_arguments():
    '''
    Parser of arguments.
    '''
    parser = argparse.ArgumentParser(description='GBK file parser')
    parser.add_argument('gbk_file', type=str, help='Path to the GBK file')
    parser.add_argument('-o', '--output', type=str, help='Output path', default='')
    args = parser.parse_args()
    return args

def make_data_to_undo_concatenation(input_gbk,prefix):
    '''
    Exports locus and size attributes of Bio.GenBank object to a tsv and
    returns a multi gbk file as an Bio.GenBank object.
    '''
    tsv_file_name=prefix+'.tsv'
    with open(input_gbk) as handle:
        with open(tsv_file_name, 'w') as tsv:
            for record in GenBank.parse(handle):
                tsv.write(f'{record.locus}\t{record.size}\n')
 
def concatenate_gbk(input_gbk, prefix):
    '''Returns a new Bio.GenBank object that is the concatenate of the input.'''
    new_record=GenBank.Record.Record()
    # Must set to 0 and not '0' or I get a value error when making the first addition
    new_record.size=0

    with open(input_gbk) as handle:
        for record in GenBank.parse(handle):

            new_record.sequence+=record.sequence
            new_record.features+=fix_coordinates_of_features(record.features, int(new_record.size))
            new_record.size=str(int(new_record.size)+int(record.size))
        
        # Setting the rest of the gbk paramaters that are fixed.
        new_record.locus=record.locus
        new_record.topology=record.topology
        new_record.source=record.source
        new_record.residue_type=record.residue_type
        new_record.organism=record.organism
        new_record.molecule_type=record.molecule_type
        new_record.definition=record.definition
        new_record.comment=record.comment
        new_record.date=record.date
    
    return new_record

def fix_coordinates_of_features(feature_list, genome_size):
    '''
    Takes as input a list of GenBank.Record.Feature objects, makes new 
    coordinates for each object by adding the current genome size, and returns a list of 
    GenBank.Record.Feature objects with the new coordinates.
    '''
    new_feature_list=[]
    for feature in feature_list:
        start, end = re.findall(r'\d+', feature.location)
        new_start=int(start)+int(genome_size)
        new_end=int(end)+int(genome_size)
        if re.match('complement',feature.location):
            fixed_feature=GenBank.Record.Feature(key=feature.key, location=f'complement({new_start}..{new_end})')
        else: 
            fixed_feature=GenBank.Record.Feature(key=feature.key, location=f'{new_start}..{new_end}')
        # Keeps qualifiers from previous entry
        fixed_feature.qualifiers=feature.qualifiers
        new_feature_list.append(fixed_feature)
        
    return new_feature_list

def export_gbk_to_file(object_to_export, output_path):
    '''
    If I provide the -o argument, save gbk to the output_path.
    '''
    with open(output_path, 'w') as outf:
        sys.stdout = outf
        print(object_to_export)

def main():
    # Parsing arguments
    args = parse_arguments()
    input_gbk= args.gbk_file
    output_path = args.output

    # Define basename for gbk file. Used for saving tsv and 
    # writing locus in gbk.
    prefix=splitext(basename(input_gbk))[0] if output_path is '' else splitext(output_path)[0]
    make_data_to_undo_concatenation(input_gbk, prefix)
    conc_gbk=concatenate_gbk(input_gbk, prefix)

    if output_path is not '': export_gbk_to_file(conc_gbk, output_path) 
    else: print(conc_gbk)
        
main()
