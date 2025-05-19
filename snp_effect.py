#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
from Bio import AlignIO
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze variable sites in aligned nucleotide sequences."
    )
    parser.add_argument("fasta_file", help="Aligned FASTA file of nucleotide sequences")
    parser.add_argument(
        "--genetic_code",
        help="Sets the genetic code to use for translating and finding non-synonymous changes.",
        default=1,
    )
    parser.add_argument(
        "reference", help="Reference sequence ID to compare others against"
    )
    parser.add_argument(
        "--skip-cds-check",
        action="store_true",
        help="Skip CDS check and codon-level classification (useful for non-coding RNA alignments)",
    )
    return parser.parse_args()


def check_alignment_validity(records):
    seq_len = len(records[0].seq)
    for record in records:
        if len(record.seq) != seq_len:
            sys.exit("All sequences must be the same length for alignment.")
    if seq_len % 3 != 0:
        sys.exit("CDS alignment length must be a multiple of 3.")


def get_variable_sites(records):
    seq_len = len(records[0].seq)
    variable_nuc_sites = []
    for i in range(seq_len):
        bases = set(record.seq[i] for record in records)
        if len(bases) > 1:
            variable_nuc_sites.append(i)
    return variable_nuc_sites


def classify_variable_sites(records, reference_id, variable_nuc_sites, genetic_code):
    record_dict = {record.id: record for record in records}
    if reference_id not in record_dict:
        sys.exit(f"Reference ID '{reference_id}' not found in the alignment.")
    reference_seq = record_dict[reference_id].seq
    codon_indices = sorted(set(i // 3 for i in variable_nuc_sites))
    synonymous_sites = 0
    nonsynonymous_sites = 0
    nonsyn_table = defaultdict(dict)
    syn_table = defaultdict(dict)
    reference_codons = {}
    for codon_index in codon_indices:
        ref_codon = reference_seq[codon_index * 3: codon_index * 3 + 3]
        ref_aa = str(Seq(ref_codon).translate(table=genetic_code))
        reference_codons[codon_index + 1] = f"{ref_aa} ({ref_codon})"
        codon_has_nonsyn_change = False
        codon_has_syn_change = False
        for record in records:
            if record.id == reference_id:
                continue

            codon = record.seq[codon_index * 3: codon_index * 3 + 3]
            aa = str(Seq(codon).translate(table=genetic_code))
            if codon != ref_codon:
                if aa != ref_aa:
                    codon_has_nonsyn_change = True
                    nonsyn_table[record.id][codon_index + 1] = f"{aa} ({codon})"
                else:
                    codon_has_syn_change = True
                    syn_table[record.id][codon_index + 1] = f". ({codon})"
        if codon_has_nonsyn_change:
            nonsynonymous_sites += 1
        elif codon_has_syn_change:
            synonymous_sites += 1
    return synonymous_sites, nonsynonymous_sites, nonsyn_table, syn_table, reference_codons


def print_summary_noncoding(variable_sites, records, reference_id, total_sites):
    print(f"Variable nucleotide sites: {len(variable_sites)}")
    print(f"Total nucleotide sites in alignment: {total_sites}\n")
    if len(variable_sites) == 0:
        print("No variable sites found.")
        return

    site_indices = [i + 1 for i in variable_sites]  # 1-based
    header = ["Taxon"] + [f"Site{i}" for i in site_indices]
    print("\t".join(header))
    record_dict = {record.id: record for record in records}
    if reference_id not in record_dict:
        sys.exit(f"Reference ID '{reference_id}' not found in the alignment.")
    ref_seq = record_dict[reference_id].seq
    ref_row = [reference_id] + [ref_seq[i] for i in variable_sites]
    print("\t".join(ref_row))
    for taxon in sorted(record_dict.keys()):
        if taxon == reference_id:
            continue

        seq = record_dict[taxon].seq
        row = [taxon]
        for i in variable_sites:
            nt = seq[i]
            row.append("." if nt == ref_seq[i] else nt)
        print("\t".join(row))


def print_summary(
    syn, nonsyn, nonsyn_table, syn_table, reference_codons, reference_id, total_codons
):
    print(f"Synonymous variable sites: {syn}")
    print(f"Non-synonymous variable sites: {nonsyn}\n")
    all_taxa = sorted(set(list(nonsyn_table.keys()) + list(syn_table.keys())))
    if not all_taxa:
        print("No variable codon changes found.")
        print(f"\nTotal codons in alignment: {total_codons}")
        return

    codon_indices = sorted(reference_codons.keys())
    header = ["Taxon"] + [f"Codon{idx}" for idx in codon_indices]
    print("\t".join(header))
    # Reference row
    ref_row = [reference_id] + [reference_codons.get(idx, ".") for idx in codon_indices]
    print("\t".join(ref_row))
    for taxon in all_taxa:
        row = [taxon]
        for idx in codon_indices:
            if idx in nonsyn_table.get(taxon, {}):
                row.append(nonsyn_table[taxon][idx])
            elif idx in syn_table.get(taxon, {}):
                row.append(syn_table[taxon][idx])
            else:
                row.append(".")
        print("\t".join(row))
    print(f"\nTotal codon sites in alignment: {total_codons}")


def main():
    args = parse_args()
    try:
        alignment = AlignIO.read(args.fasta_file, "fasta")
    except Exception as e:
        sys.exit(f"Error reading FASTA file: {e}")
    alignment_length = len(alignment[0].seq)
    variable_sites = get_variable_sites(alignment)
    if args.skip_cds_check:
        print_summary_noncoding(
            variable_sites, alignment, args.reference, alignment_length
        )
    else:
        check_alignment_validity(alignment)
        syn, nonsyn, nonsyn_table, syn_table, ref_codons = classify_variable_sites(
            alignment, args.reference, variable_sites, args.genetic_code
        )
        total_codons = alignment_length // 3
        print_summary(
            syn,
            nonsyn,
            nonsyn_table,
            syn_table,
            ref_codons,
            args.reference,
            total_codons,
        )


if __name__ == "__main__":
    main()
