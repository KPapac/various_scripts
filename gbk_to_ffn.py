#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Reads through gene features of gbk file, discards pseudogenes. Returns a fasta file, or fasta sequences to std output.')
parser.add_argument("path_to_gbk", help='Path to gbk file.')
parser.add_argument("-o", "--output", help='Path to ffn output file.', default=None)
parser.add_argument("-faa", "--to_protein", help="Outputs amino acid sequences instead of nucleotide.", action='store_true')
args = parser.parse_args()

if args.output!=None and os.path.isfile(args.output): sys.exit(f'Tried to append sequences to existing file {args.output}. Maybe delete this file first?')

def gbk_to_ffn(feature_name, feature_seq, strand, out_path):
    '''Checks if feature is in + or in - strand and prints the
    fasta entry in std output.'''
    if out_path == None:
        if strand == -1: print(f">{feature_name}\n{feature_seq.reverse_complement()}")
        elif strand == 1: print(f">{feature_name}\n{feature_seq}")
        else: raise Exception(f'{feature_name} has no strand orientation. Do not know what to return. Exiting...')
    else:
        with open(out_path, 'a') as outfasta:
            if strand == -1: outfasta.write(f">{feature_name}\n{feature_seq.reverse_complement()}\n")
            elif strand == 1: outfasta.write(f">{feature_name}\n{feature_seq}\n")
            else: raise Exception(f'{feature_name} has no strand orientation. Do not know what to return. Exiting...')

def gbk_to_faa(feature_name, feature_seq, strand, out_path):
    '''Checks if feature is in + or in - strand and prints the
    fasta entry in std output.'''
    if out_path == None:
        if strand == -1:
            aa_seq_rev=feature_seq.reverse_complement().translate(table="Bacterial", cds=True, to_stop=True)
            print(f'>{feature_name}\n{aa_seq_rev}')
        elif strand == 1:
            aa_seq=feature_seq.translate(table="Bacterial", cds=True, to_stop=True)
            print(f'>{feature_name}\n{aa_seq}')
        else: raise Exception(f'{feature_name} has no strand orientation. Do not know what to return. Exiting...')
    
    else:
        with open(out_path, 'a') as outfasta:
            if strand == -1:
                try:
                    aa_seq_rev=feature_seq.reverse_complement().translate(table="Bacterial", cds=True, to_stop=True)
                    outfasta.write(f'>{feature_name}\n{aa_seq_rev}\n')
                except: 
                    #print(f"Could not translate {feature_name}, probably an RNA, skipping...")
                    return -1
            elif strand == 1:
                try:
                    aa_seq=feature_seq.translate(table="Bacterial", cds=True, to_stop=True)
                    outfasta.write(f'>{feature_name}\n{aa_seq}\n')
                except: 
                    #print(f"Could not translate {feature_name}, probably an RNA, skipping...")
                    return -1
            else: raise Exception(f'{feature_name} has no strand orientation. Do not know what to return. Exiting...')


def parse_gbk(path_to_gbk, path_to_ffn=None, faa_val=False):
    '''Reads though gene features of gbk file. Looks if gene feature is pseudogene. Then,
    calls gbk_to_ffn() to output a fasta sequence, either to std output, or to a fasta file.
    '''
    pseudo_genes=0
    partial_genes=0
    skipped_genes = 0
    written_gene_features = 0
    total_gene_features = 0

    records = list(SeqIO.parse(path_to_gbk, "genbank"))
    for record in records:
        for feature in record.features:
            # Get the gene features of the genbank file.
            if feature.type == 'CDS':
                total_gene_features += 1
                # The sequence in the + strand as in the genbank file is:
                feature_name= feature.qualifiers['locus_tag'][0]
                feature_seq = record.seq[feature.location.start:feature.location.end]

                # If pseudogene, then skip.
                if 'pseudo' in feature.qualifiers.keys():
                    pseudo_genes+=1
                    continue
                elif 'pseudogene' in feature.qualifiers.keys():
                    pseudo_genes+=1
                    continue


                # If partial, then skip.
                elif 'partial' in feature.qualifiers.keys():
                    partial_genes+=1
                    continue

                else:
                    try:
                        if faa_val == False: gbk_to_ffn(feature_name, feature_seq, feature.location.strand, path_to_ffn)
                        elif faa_val == True: 
                            gbk_to_faa(feature_name, feature_seq, feature.location.strand, path_to_ffn)
                            if gbk_to_faa == -1: 
                                skipped_genes+=1
                                continue
                        written_gene_features += 1
                    except: f"Failed to get {feature_name}"
                
           # If the feature is not a 'gene', then just continue to the next feature.
            else: continue
    #print(f"{feature_name}\t{total_gene_features}\t{written_gene_features}\t{skipped_genes}")
    with open(f"logs/gbk_to_fasta.log", 'a') as report_f:
        skipped_genes+=(pseudo_genes+partial_genes)
        report_f.write(f"{path_to_gbk}\t{total_gene_features}\t{written_gene_features}\t{pseudo_genes}\t{partial_genes}\t{skipped_genes}\t{written_gene_features+skipped_genes-total_gene_features}\n")
    #print("Genes skipped:", skipped_genes)
    #print("Total genes:", total_gene_features)

if args.output != None: print(f"Writing to {args.output}")
parse_gbk(args.path_to_gbk, args.output, args.to_protein)
