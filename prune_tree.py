#!/usr/bin/env python3
from ete3 import Tree
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Reads a nwk file with multiple trees. Takes each Newick tree and prunes taxa. By default, taxa in the list are kept in the tree. If --exclude is provided, taxa in the list are removed.'
)
# Required argument
parser.add_argument("nwk_in", help='Input Newick file with trees to prune.')
# Create mutually exclusive groups
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "--get_input_taxa",
    action="store_true",
    help="Just prints the taxa from all input trees.",
)
group.add_argument(
    "--prune",
    nargs=2,
    metavar=('list_of_taxa', 'nwk_out'),
    help="Provide a text file with list of taxa, and an output path for the pruned Newick trees.",
)
parser.add_argument(
    "--exclude",
    action="store_true",
    help="If set, taxa in the list are removed instead of kept.",
)
args = parser.parse_args()


def prune_tree(tree_in, list_of_taxa_file):
    with open(list_of_taxa_file, 'r') as f:
        taxa_in_list = {line.strip() for line in f}
    try:
        t = Tree(tree_in, format=2)
    except Exception:
        print(f"Could not open tree from string:\n{tree_in}")
        return None

    taxa_in_tree = {taxon.name for taxon in t.get_descendants() if taxon.name}
    if args.exclude:
        taxa_to_keep = taxa_in_tree - taxa_in_list
    else:
        taxa_to_keep = taxa_in_list
    t.prune(taxa_to_keep, preserve_branch_length=True)
    t.unroot()
    return t.write(format=2)


def print_input_taxa(tree_path):
    input_taxa = set()
    try:
        with open(tree_path, 'r') as f:
            trees = f.readlines()
            for tree_str in trees:
                if tree_str.strip():
                    tree = Tree(tree_str, format=2)
                    for taxon in tree.get_descendants():
                        if taxon.name:
                            input_taxa.add(taxon.name)
    except Exception:
        print("Could not read the nwk file.")
        sys.exit(1)
    for taxon in sorted(input_taxa):
        print(taxon)


def main():
    if args.get_input_taxa:
        print_input_taxa(args.nwk_in)
        sys.exit()
    list_of_taxa, nwk_out = args.prune
    pruned_trees = []
    with open(args.nwk_in, 'r') as trees:
        print("Computing pruned trees...")
        for tree in trees:
            if tree.strip():
                pruned_tree = prune_tree(tree, list_of_taxa)
                if pruned_tree:
                    pruned_trees.append(pruned_tree)
    with open(nwk_out, 'w') as fout:
        print(f"Writing {len(pruned_trees)} trees to {nwk_out}...")
        for tree in pruned_trees:
            fout.write(f"{tree}\n")


if __name__ == "__main__":
    main()
