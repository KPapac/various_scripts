#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import argparse

# Parsing arguments
parser = argparse.ArgumentParser(
    prog="backtranslate.py",
    description="Takes a Nucleotide Fasta as input, translates into protein, performs protein alignment with MAFFT and returns the nucleotide alignment.",
    epilog="",
)
parser.add_argument(
    "--genetic_code",
    help="Set the genetic code for protein translation. Default is the Universal code (1). If aligning mitochondrial proteins, or anything else that does not follow the Universal genetic code, then you MUST set this...",
    default=1,
)
parser.add_argument("input_fasta", type=str)
parser.add_argument(
    "-p",
    "--prefix",
    type=str,
    default="out",
    help="Name of output file. Default value: out",
)
parser.add_argument(
    "-T",
    "--threads_mafft",
    type=int,
    default=1,
    help="Threads to use for MAFFT alignment. Default is 1.",
)
args = parser.parse_args()


def translate_dna(dna_seq, code):
    protein_seq = dna_seq.translate(table=code, to_stop=True)
    return protein_seq


def backtranslate_protein(protein_alg, nucleotide_fasta, code):
    prot_char_dict = {}
    # Reading protein alignment
    for record in SeqIO.parse(protein_alg, "fasta"):
        char_array = np.array(list(str(record.seq)), dtype="U1")
        prot_char_dict[record.description] = char_array
    # Reading nucleotides (ungapped, split into codons)
    nucl_char_dict = {}
    for record in SeqIO.parse(nucleotide_fasta, "fasta"):
        seq = record.seq  # important: use .seq
        char_array = np.array(
            ["".join(seq[i: i + 3]) for i in range(0, len(seq), 3)], dtype="U3"
        )
        nucl_char_dict[record.description] = char_array
    for fasta_entry in nucl_char_dict.keys():
        prot = prot_char_dict[fasta_entry]
        nucl = nucl_char_dict[fasta_entry]
        aligned_codons = []
        j = 0  # index into nucl (codon array)
        # Walk along the protein alignment
        for aa in prot:
            if aa == "-":
                # if gap in protein then gap codon
                aligned_codons.append("---")
            else:
                # if amino acid then take next codon if available
                codon = nucl[j]
                j += 1
                aligned_codons.append(codon)
        # Verify: each triplet encodes the correct amino acid
        for i, (codon, aa) in enumerate(zip(aligned_codons, prot)):
            # skip gaps and stop codons
            if aa in ("-", "*") or codon == "---":
                continue

            if Seq(codon).translate(table=code) != aa:
                sys.exit(
                    f"Error: position {i}: codon {codon} does not encode "
                    f"{aa} in {fasta_entry}. Exiting."
                )
        # Join nucleotides to a string
        nucl_char_dict[fasta_entry] = "".join(aligned_codons)
    backtranslated_dna_seq = [
        SeqRecord(Seq(seq), id=head, description="")
        for head, seq in nucl_char_dict.items()
    ]
    return backtranslated_dna_seq


def main():
    # Making a protein FASTA from the nucleotide FASTA
    dna_sequences = list(SeqIO.parse(args.input_fasta, "fasta"))
    protein_sequences = []
    for dna_seq in dna_sequences:
        protein_seq = translate_dna(dna_seq.seq, args.genetic_code)
        protein_record = SeqRecord(
            protein_seq, id=dna_seq.id, description=dna_seq.description
        )
        protein_sequences.append(protein_record)
    SeqIO.write(protein_sequences, f"{args.prefix}_protein_sequences.fasta", "fasta")
    # Running protein alignment
    subprocess.run(
        [
            "mafft-linsi",
            "--thread",
            str(args.threads_mafft),
            f"{args.prefix}_protein_sequences.fasta",
        ],
        stdout=open(f"{args.prefix}_protein_alignment.fasta", "w"),
    )
    # Backtranslating protein alignment
    dna_alignment = backtranslate_protein(
        f"{args.prefix}_protein_alignment.fasta", args.input_fasta, args.genetic_code
    )
    # Write output
    SeqIO.write(dna_alignment, args.prefix + ".mafft", "fasta")
    # Clean intermediate files
    if os.path.exists(f"{args.prefix}_protein_sequences.fasta"):
        os.remove(f"{args.prefix}_protein_sequences.fasta")
    if os.path.exists(f"{args.prefix}_protein_alignment.fasta"):
        os.remove(f"{args.prefix}_protein_alignment.fasta")
    print(
        "If aligning mitochondrial proteins, or anything else that does not follow the Universal genetic code, then you MUST set the --genetic_code argument to have a proper alignment."
    )


if __name__ == "__main__":
    main()
