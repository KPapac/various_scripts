#!/usr/bin/env python3

from ete3 import Tree
import sys
import argparse
import re
from typing import List, Set, Dict


# ---------- helpers ----------
def read_lines(path: str) -> List[str]:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    except Exception:
        print("Could not read the nwk file.", file=sys.stderr)
        sys.exit(1)


def compile_regex(pat: str, ignore_case: bool):
    flags = re.I if ignore_case else 0
    return re.compile(pat, flags)


def leaves_re(tree: Tree, rx: re.Pattern):
    return [lf for lf in tree.get_leaves() if rx.search(lf.name or "")]


def to_ultrametric_by_tip_extension(tree: Tree):
    """Safe fallback: extend only terminal branches to equalize root->leaf distances."""
    leaves = tree.get_leaves()
    if not leaves:
        return
    dists = {lf: tree.get_distance(lf) for lf in leaves}
    Dmax = max(dists.values())
    for lf, d in dists.items():
        lf.dist += Dmax - d


def all_taxa_in_tree(t: Tree) -> Set[str]:
    return {n.name for n in t.get_descendants() if n.name}


def read_tsv_map(tsv_path: str) -> Dict[str, str]:
    """Read 2-column TSV (old_name \t new_name). Ignores blank and '#' lines."""
    mapping: Dict[str, str] = {}
    with open(tsv_path, "r", encoding="utf-8") as fh:
        for i, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    f"Annotate TSV: line {i} doesn't have 2 tab-separated fields."
                )
            old, new = parts[0].strip(), parts[1].strip()
            if not old or not new:
                raise ValueError(f"Annotate TSV: line {i} has empty name(s).")
            mapping[old] = new
    if not mapping:
        raise ValueError("Annotate TSV produced an empty mapping.")
    return mapping


# ---------- core operations ----------
def parse_tree(tree_str: str, fmt: int) -> Tree:
    try:
        return Tree(tree_str, format=fmt)
    except Exception:
        raise ValueError("Failed to parse input tree string.")


def prune_tree_obj(t: Tree, list_of_taxa_file: str, exclude: bool) -> Tree:
    with open(list_of_taxa_file, "r", encoding="utf-8") as f:
        taxa_in_list = {line.strip() for line in f if line.strip()}

    taxa_in_tree = all_taxa_in_tree(t)
    taxa_to_keep = (
        (taxa_in_tree - taxa_in_list) if exclude else (taxa_in_tree & taxa_in_list)
    )

    if not taxa_to_keep:
        raise ValueError("Pruning produced an empty taxa set.")

    t.prune(taxa_to_keep, preserve_branch_length=True)
    return t


def annotate_tree(t: Tree, mapping: Dict[str, str], tree_idx: int):
    """Require full coverage per tree; then rename."""
    leaves = t.get_leaves()
    leaf_names = [lf.name for lf in leaves]
    missing = sorted({name for name in leaf_names if name not in mapping})
    if missing:
        preview = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
        raise ValueError(f"annotate: tree {tree_idx} missing mappings for: {preview}")
    for lf in leaves:
        lf.name = mapping[lf.name]


def root_on_mrca(t: Tree, pat_a: str, pat_b: str, ignore_case: bool, unique: bool):
    rx_a = compile_regex(pat_a, ignore_case)
    rx_b = compile_regex(pat_b, ignore_case)

    A = leaves_re(t, rx_a)
    B = leaves_re(t, rx_b)

    if unique:
        if len(A) != 1 or len(B) != 1:
            raise ValueError(
                f"--unique set but pattern A matched {len(A)} and B matched {len(B)} "
                f"(A={ [lf.name for lf in A] }, B={ [lf.name for lf in B] })"
            )
        mrca = t.get_common_ancestor(A[0], B[0])
    else:
        if not A or not B:
            raise ValueError(
                f"Patterns must each match at least one leaf. " f"A={len(A)} B={len(B)}"
            )
        mrca = t.get_common_ancestor(A + B)

    t.set_outgroup(mrca)


def make_ultrametric(
    t: Tree, do_ultra: bool, strategy: str = None, tree_length: float = None
):
    if not do_ultra:
        return
    used_fallback = False
    try:
        convert = getattr(t, "convert_to_ultrametric", None)
        if convert is None:
            used_fallback = True
        else:
            kwargs = {}
            if strategy:
                kwargs["strategy"] = strategy
            if tree_length is not None:
                kwargs["tree_length"] = float(tree_length)
            convert(**kwargs)
    except Exception:
        used_fallback = True

    if used_fallback:
        to_ultrametric_by_tip_extension(t)


# ---------- CLI ----------
parser = argparse.ArgumentParser(
    description=(
        "Reads a Newick file (one or many trees). "
        "Otherwise, optionally --prune, --annotate, --root, and/or --ultrametric. "
        "If you don't use --prune, you must provide -o/--output."
    )
)

subparsers = parser.add_subparsers(dest="command", required=True)
info_parser = subparsers.add_parser(
    "info",
    help="Returns taxa shared across all trees and taxa not shared across all trees.",
)
# Required positional
info_parser.add_argument(
    "nwk_in", help="Input Newick file (can contain multiple trees)."
)


tool_parser = subparsers.add_parser(
    "tools", help="Prunes, Reroots, Annotates, and/or converts trees to Ultrametric."
)

prune_grp = tool_parser.add_argument_group("Pruning taxa out")
prune_grp.add_argument(
    "--taxa",
    nargs=1,
    help="A list of taxa in a txt file, each line contains one taxon only. By default, only taxa in the list are kept in the trees. Use --exclude to invert (remove only listed taxa).",
)
prune_grp.add_argument(
    "--exclude",
    action="store_true",
)

# annotate
annot_grp = tool_parser.add_argument_group(
    "Annotating with new labels",
    "Renames leaf labels using a 2-column TSV (old_name TAB new_name).",
)
annot_grp.add_argument(
    "--annotate",
    metavar="TSV",
    help="Apply a 2-column TSV mapping to all leaves in each tree (post-prune if pruning). "
    "All leaf names in a tree must exist in the TSV; otherwise the tree errors.",
)

# root
root_grp = tool_parser.add_argument_group(
    "Rerooting Trees",
    "Root on the MRCA of two (regex) taxa sets. Runs after annotate if both are set.",
)
root_grp.add_argument(
    "--root", action="store_true", help="Enable rooting on MRCA of two regex patterns."
)
root_grp.add_argument("--pattern-a", help="Regex for taxon set A (e.g., '^DNEO').")
root_grp.add_argument("--pattern-b", help="Regex for taxon set B (e.g., '^DSTU').")
root_grp.add_argument(
    "--ignore-case",
    action="store_true",
    help="Case-insensitive regex matching for rooting.",
)
root_grp.add_argument(
    "--unique",
    action="store_true",
    help="Require each pattern to match exactly one leaf (else error).",
)

# ultrametric
ultra_grp = tool_parser.add_argument_group(
    "Converting tress to ultrametric",
    'Convert tree to ultrametric. Runs after "root" if both are set.',
)
ultra_grp.add_argument(
    "--ultrametric", action="store_true", help="Enable ultrametric conversion."
)
ultra_grp.add_argument(
    "--strategy",
    default=None,
    help="Strategy for ETE3 convert_to_ultrametric (e.g., 'balanced'). "
    "If unsupported, a safe tip-extension fallback is used.",
)
ultra_grp.add_argument(
    "--tree-length",
    type=float,
    default=None,
    help="Target total tree length for ultrametric conversion (ETE3 only).",
)

# Required positional
tool_parser.add_argument(
    "nwk_in", help="Input Newick file (can contain multiple trees)."
)

tool_parser.add_argument(
    "-f",
    "--format",
    type=int,
    default=2,
    help="ETE3 Newick format code for parsing/writing (default: 2).",
)
tool_parser.add_argument(
    "--skip-errors",
    action="store_true",
    help="On an error for a given tree, log to stderr and continue.",
)
args = parser.parse_args()


def validate_rooting_args():
    if args.root and (not args.pattern_a or not args.pattern_b):
        parser.error("--root requires --pattern-a and --pattern-b.")


# ---------- actions ----------
def print_input_taxa(tree_path: str):
    tree_lines = read_lines(tree_path)
    if not tree_lines:
        print("No trees found in the input file.", file=sys.stderr)
        sys.exit(1)

    tree_taxa_list = []
    for tree_str in tree_lines:
        try:
            tree = Tree(tree_str)
            taxa = all_taxa_in_tree(tree)
            tree_taxa_list.append(taxa)
        except Exception as e:
            print(f"Failed to parse a tree: {e}", file=sys.stderr)

    if not tree_taxa_list:
        print("No parsable trees found.", file=sys.stderr)
        sys.exit(1)

    if len(tree_taxa_list) == 1:
        print("Taxa found in the tree:")
        for taxon in sorted(tree_taxa_list[0]):
            print(taxon)
    else:
        shared_taxa = set.intersection(*tree_taxa_list)
        all_taxa = set.union(*tree_taxa_list)
        print(f"Shared taxa (present in all {len(tree_taxa_list)} trees):")
        for taxon in sorted(shared_taxa):
            print(f"  {taxon}")
        print("\nUnique taxa (appear in some trees but not all):")
        unique_taxa = all_taxa - shared_taxa
        for taxon in sorted(unique_taxa):
            print(f"  {taxon}")


def process_pipeline():
    in_lines = read_lines(args.nwk_in)

    # Inputs for optional steps
    mapping = read_tsv_map(args.annotate) if args.annotate else None
    list_file = args.taxa[0] if args.taxa else None

    out_trees = []
    for idx, tree_str in enumerate(in_lines, start=1):
        try:
            # parse
            t = parse_tree(tree_str, args.format)

            # prune (optional)
            if list_file:
                t = prune_tree_obj(t, list_file, args.exclude)

            # annotate (optional)
            if mapping is not None:
                annotate_tree(t, mapping, idx)

            # root (optional)
            if args.root:
                root_on_mrca(
                    t,
                    pat_a=args.pattern_a,
                    pat_b=args.pattern_b,
                    ignore_case=args.ignore_case,
                    unique=args.unique,
                )

            # ultrametric (optional)
            make_ultrametric(
                t,
                do_ultra=args.ultrametric,
                strategy=args.strategy,
                tree_length=args.tree_length,
            )

            out_trees.append(t.write(format=args.format).strip())

        except Exception as e:
            msg = f"[tree {idx}] error: {e}"
            if args.skip_errors:
                print(msg, file=sys.stderr)
                continue
            else:
                raise

    for newick in out_trees:
        print(newick)
    return 0


def main():
    if args.command == "info":
        print_input_taxa(args.nwk_in)
        return
    process_pipeline()


if __name__ == "__main__":
    main()
