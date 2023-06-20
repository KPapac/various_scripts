#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os.path import basename, splitext
import re
from Bio import GenBank
import copy
# from pprint import pprint
#pprint(vars(entry))


def make_data_to_undo_concatenation(input_gbk,tsv_file_name):
    '''
    Exports locus and size attributes of Bio.GenBank object to a tsv and
    returns a multi gbk file as an Bio.GenBank object.
    '''
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

def main(input_gbk):
    # Define basename for gbk file. Used for saving tsv and 
    # writing locus in gbk.
    prefix=splitext(basename(input_gbk))[0]
    make_data_to_undo_concatenation(input_gbk, prefix)
    conc_gbk=concatenate_gbk(input_gbk, prefix)
    # print(type(conc_gbk))
    print(conc_gbk)
    # pprint(vars(conc_gbk))
        
main(input_gbk="./test.gbk")
