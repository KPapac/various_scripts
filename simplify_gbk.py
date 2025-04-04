#!/usr/bin/env python3
import sys


def process_gene_cds():
    output = []
    entry = []
    current_gene = None
    for line in sys.stdin:
        line = line.strip()
        if line.startswith("gene"):
            if current_gene:  # Save previous gene entry
                print(", ".join(current_gene))
            current_gene = [line]  # Start new gene entry
        else:
            if current_gene:
                current_gene.append(line)  # Append metadata to gene entry
    if current_gene:  # Save the last gene entry
        print(", ".join(current_gene))


if __name__ == "__main__":
    process_gene_cds()
