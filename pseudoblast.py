#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

def run_blastx(QUERY="/home/kostas/bases"):
    DB="/space/no_backup/databases/Wolbachia_prot_ncbi_220302.faa"
    cmd = f'blastx -query {QUERY} -db {DB} -evalue 0.01 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe stitle slen"'
    blast_result=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout
    return(blast_result.decode('ascii'))

def save_to_tmp_check_file(QUERY="/home/kostas/bases", FILE="check"):
    DB="/space/no_backup/databases/Wolbachia_prot_ncbi_220302.faa"
    QUERY="/home/kostas/bases"
    cmd = f'blastx -query {QUERY} -db {DB} -evalue 0.01'
    blast_result=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('ascii')
    with open(FILE,'w') as f:
       f.write(blast_result) 

def parse_check_file(CHECK_FILE="check"):
    cmd=r'grep -oP "gi\|.*\|" ~/check | sort | uniq -c | sort -k1'
    parsed_table=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('ascii')
    list_of_candidate_hits=[]
    for line in parsed_table.splitlines():
        (number_of_hits, title) = line.strip().split(sep=" ")
        if int(number_of_hits)>2: 
            list_of_candidate_hits.append((number_of_hits, title))
    return list_of_candidate_hits

def print_results_to_screen(BLAST_RESULT, LIST_OF_HITS):
    if len(LIST_OF_HITS)==0: print("No multiple matches to a subject.") 
    else:
        for hit in LIST_OF_HITS:
            print(hit[1])
            for line in BLAST_RESULT.splitlines():
                blast_columns=line.split(sep="\t")
                if hit[1] in line:
                    msg=f"*{blast_columns[0]}\n\tSubject span: {blast_columns[8]} - {blast_columns[9]}\tFrame:{blast_columns[12]}\t%Identity: {blast_columns[2]} Missmatches: {blast_columns[4]} Gaps: {blast_columns[5]}\tevalue: {blast_columns[10]} Bitscore: {blast_columns[11]}"
                    print(msg)
                    print(f'\tSubject length: {blast_columns[14]}\t{blast_columns[13]}')
            print("-"*len(msg))


def main():
    blast_result=run_blastx()
    save_to_tmp_check_file("/home/kostas/check")
    list_of_possible_hits=parse_check_file("/home/kostas/check")
    print_results_to_screen(blast_result, list_of_possible_hits)
    print("Done.")

if __name__=="__main__":
    main()
