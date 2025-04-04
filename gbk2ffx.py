#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys


def extract_features(
    gbk_file, output_file=None, faa=False, translation_table=1, locus_tag=None
):
    records = []
    for record in SeqIO.parse(gbk_file, "genbank"):
        feature_counter = 0
        for feature in record.features:
            if feature.type == "gene" or feature.type == "source":
                continue

            if faa:
                if feature.type == "CDS":  # Extract amino acid sequences
                    if "translation" in feature.qualifiers:
                        seq = Seq(feature.qualifiers["translation"][0])
                    else:
                        seq = feature.extract(record.seq).translate(
                            table=translation_table, cds=True
                        )
                else:
                    continue  # Do not print feature if not a CDS

            else:  # Extract nucleotide sequences (FFN)
                seq = feature.extract(record.seq)
            feature_counter += 1
            if locus_tag is None:
                try:
                    fasta_header = feature.qualifiers["locus_tag"][0]
                except:
                    raise Exception(
                        print(
                            f"The feature {feature} may not have a locus tag. Consider supplying the -l argument."
                        )
                    )

            else:
                fasta_header = f"{locus_tag}_{str(feature_counter).zfill(3)}"
            try:
                seq_record = SeqRecord(
                    seq,
                    id=fasta_header,
                    description=f"{feature.qualifiers['product'][0]} from {record.id}",
                )
            except:
                raise Exception(f'Problem with {feature}..')

            records.append(seq_record)
    # Output results
    output_handle = open(output_file, "w") if output_file else sys.stdout
    SeqIO.write(records, output_handle, "fasta")
    if output_file:
        output_handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="Extract FFN or FAA sequences from a GenBank file."
    )
    parser.add_argument("gbk_file", help="Input GenBank (.gbk) file.")
    parser.add_argument("--locus_tag", "-l", help="A locus_tag of your preference.")
    parser.add_argument(
        "-o", "--output", help="Output file path (default: stdout)", default=None
    )
    parser.add_argument(
        "--faa",
        action="store_true",
        help="Extract FAA from CDS features (amino acid sequences), default is FFN.",
    )
    parser.add_argument(
        "--table",
        type=int,
        default=1,
        help="Translation table (default: 1). Relevant for the FAA mode. See tables at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.",
    )
    args = parser.parse_args()
    extract_features(args.gbk_file, args.output, args.faa, args.table, args.locus_tag)


if __name__ == "__main__":
    main()
