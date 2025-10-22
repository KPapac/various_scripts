#!/usr/bin/env python3
"""
Translate nucleotide sequences to protein.

- Reads nucleotide sequences from a FASTA (default) or GenBank file.
- Translates using a chosen genetic code (by name or NCBI table ID).
- Optionally keeps stop codons '*' (default: removed).
- Trims trailing bases so the length is a multiple of 3.

Examples
--------
Translate with standard code (remove stops, default):
    python ffn2faa.py input.fasta > proteins.faa

Translate with vertebrate mitochondrial code, keep stops:
    python ffn2faa.py input.fasta -o proteins.faa --code "Vertebrate Mitochondrial" --keep-stops

Read from stdin, write to stdout:
    cat input.fasta | python ffn2faa.py - > proteins.faa
"""

import sys
import argparse
from typing import Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Translate nucleotide sequences to protein using Biopython."
    )
    parser.add_argument(
        "input",
        help="Input nucleotide file (FASTA or GenBank). Use '-' to read from STDIN.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output FASTA file for protein sequences (default: STDOUT).",
    )
    parser.add_argument(
        "--format",
        default="fasta",
        choices=["fasta", "genbank"],
        help="Input file format (default: fasta).",
    )
    parser.add_argument(
        "--code",
        default="Standard",
        help=(
            "Genetic code to use: NCBI table ID (e.g., 1, 2, 5, 11) or name "
            '(e.g., "Standard", "Vertebrate Mitochondrial"). Default: Standard (1).'
        ),
    )
    parser.add_argument(
        "--keep-stops",
        action="store_true",
        help="Keep stop codons '*' in output proteins. Default: remove them.",
    )
    return parser.parse_args()


def resolve_table_id_or_name(table_arg: str) -> Union[int, str]:
    """
    Accept either an integer (as string) for the NCBI table id or a table name.

    Returns the int ID if numeric, otherwise returns the original name string.
    Raises ValueError if the name is not a valid codon table.
    """
    # Try numeric ID first
    try:
        tid = int(table_arg)
        # Validate that the ID exists
        CodonTable.unambiguous_dna_by_id[tid]
        return tid

    except ValueError:
        # Not an int -> treat as name
        pass
    except KeyError:
        raise ValueError(f"Unknown genetic code id: {table_arg!r}")

    # Validate by name
    try:
        CodonTable.unambiguous_dna_by_name[table_arg]
    except KeyError:
        # Try some common normalized variants
        normalized = table_arg.strip().lower()
        for name in CodonTable.unambiguous_dna_by_name:
            if name.lower() == normalized:
                return name

        raise ValueError(f"Unknown genetic code name: {table_arg!r}")

    return table_arg


def translate_record(rec: SeqRecord, table, keep_stops: bool) -> SeqRecord:
    """Translate one nucleotide record to protein, trimming trailing bases."""
    # Remove gaps and whitespace, ensure uppercase
    seq = rec.seq.ungap("-").upper()
    # Trim trailing nucleotides so length is divisible by 3
    if len(seq) % 3 != 0:
        trimmed_len = len(seq) - (len(seq) % 3)
        seq = seq[:trimmed_len]
    # Translate. We *always* include stops in the initial translation,
    # then drop them if keep_stops is False. This avoids premature truncation.
    protein = seq.translate(table=table)
    if not keep_stops:
        # Remove all stop symbols '*'
        protein = protein.replace("*", "")
    desc_extra = f"[translated table={table}]"
    new_desc = (rec.description + " " + desc_extra).strip()
    out_rec = SeqRecord(protein, id=rec.id, name=rec.name, description=new_desc)
    return out_rec


def main():
    args = parse_args()
    try:
        table = resolve_table_id_or_name(args.code)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
    # Input handle
    if args.input == "-":
        in_handle = sys.stdin
    else:
        try:
            in_handle = open(args.input, "r")
        except OSError as e:
            print(f"Error opening input file: {e}", file=sys.stderr)
            sys.exit(2)
    # Output handle
    if args.output == "-" or args.output is None:
        out_handle = sys.stdout
    else:
        try:
            out_handle = open(args.output, "w")
        except OSError as e:
            print(f"Error opening output file: {e}", file=sys.stderr)
            if args.input != "-":
                in_handle.close()
            sys.exit(2)
    try:
        records = SeqIO.parse(in_handle, args.format)
        out_records = (translate_record(r, table, args.keep_stops) for r in records)
        n_written = SeqIO.write(out_records, out_handle, "fasta")
    except Exception as e:
        print(f"Error during translation: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        if args.input != "-":
            in_handle.close()
        if args.output not in ("-", None):
            out_handle.close()
    # Optional: report to stderr
    print(f"Translated {n_written} sequence(s).", file=sys.stderr)


if __name__ == "__main__":
    main()
