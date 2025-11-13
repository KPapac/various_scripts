#!/usr/bin/env python3
import sys
import re
import argparse
from typing import Iterable, TextIO, Set, List
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter FASTA sequences by name. Reads from one or more files or stdin\n"
            "and prints matching records to stdout.\n\n"
            "Names to search for can be provided as positional arguments, via\n"
            "-f/--list-of-names, or both. All names from both sources are combined\n"
            "into a single list. How these names are interpreted depends on the\n"
            "selected matching mode:\n"
            "  * default (no --partial, no --regex): exact match against record ID\n"
            "    or full FASTA header (description)\n"
            "  * --partial: names are used as plain substrings of ID or header\n"
            "  * --regex: names are interpreted as regular expressions and matched\n"
            "    with re.search() against ID and header\n\n"
            "Options --partial and --regex are mutually exclusive."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Optional positional names (can be empty if -f/--list-of-names is used)
    parser.add_argument(
        "names",
        nargs="*",
        help=(
            "Sequence names or patterns to select. By default, these are used for\n"
            "exact matches against the FASTA record ID or full header. With\n"
            "--partial, they are treated as substrings. With --regex, they are\n"
            "interpreted as regular expressions.\n"
            "You may also provide names via -f/--list-of-names; the two sources\n"
            "are combined."
        ),
    )
    # File with names (one per line)
    parser.add_argument(
        "-f",
        "--list-of-names",
        metavar="FILE",
        help=(
            "File containing one name or pattern per line to search for. Empty\n"
            "lines and lines starting with '#' are ignored. Names from this file\n"
            "are combined with any positional NAMES. They are interpreted\n"
            "according to the selected matching mode (exact / --partial / --regex)."
        ),
    )
    # Optional: input FASTA files (if none given, read from stdin)
    parser.add_argument(
        "-i",
        "--input",
        metavar="FASTA",
        nargs="*",
        default=[],
        help=(
            "Input FASTA file(s). Use '-' or omit this option to read from stdin.\n"
            "Multiple files (and/or '-') may be provided."
        ),
    )
    # Matching mode flags (mutually exclusive)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-p",
        "--partial",
        action="store_true",
        help=(
            "Enable partial matching: a record is selected if any provided name\n"
            "appears as a substring of the record ID or the full header.\n"
            "Mutually exclusive with --regex."
        ),
    )
    group.add_argument(
        "--regex",
        action="store_true",
        help=(
            "Interpret all provided names (both positional NAMES and entries from\n"
            "-f/--list-of-names) as regular expressions. Each pattern is compiled\n"
            "and matched with re.search() against the record ID and full header.\n"
            "Mutually exclusive with --partial."
        ),
    )
    return parser.parse_args()


def load_names_from_file(path: str) -> List[str]:
    """Load names from a flat file, ignoring empty lines and comments."""
    names: List[str] = []
    try:
        with open(path, "r") as fh:
            for line in fh:
                name = line.strip()
                if not name:
                    continue  # ignore empty lines

                if name.startswith("#"):
                    continue  # ignore comment lines

                names.append(name)
    except OSError as e:
        print(f"Error opening name file '{path}': {e}", file=sys.stderr)
        sys.exit(1)
    return names


def records_from_handles(handles: Iterable[TextIO]):
    """Yield all FASTA records from a list of open file handles."""
    for handle in handles:
        for record in SeqIO.parse(handle, "fasta"):
            yield record


def main() -> None:
    args = parse_args()
    # Collect names from command line and optional file
    cli_names: List[str] = list(args.names)
    file_names: List[str] = []
    if args.list_of_names:
        file_names = load_names_from_file(args.list_of_names)
    all_names_list: List[str] = cli_names + file_names
    if not all_names_list:
        print(
            "Error: You must provide at least one name either as positional NAMES "
            "or via -f/--list-of-names.",
            file=sys.stderr,
        )
        sys.exit(1)
    # Prepare matching structures depending on mode
    if args.regex:
        try:
            regex_patterns = [re.compile(p) for p in all_names_list]
        except re.error as e:
            print(f"Invalid regular expression: {e}", file=sys.stderr)
            sys.exit(1)
        names_set: Set[str] = set()  # unused in regex mode
        names_list: List[str] = []  # unused in regex mode
    else:
        names_set: Set[str] = set(all_names_list)  # for exact match
        names_list: List[str] = all_names_list  # for partial matching
    # Decide input sources
    handles: List[TextIO] = []
    if not args.input:
        # No input files provided → read from stdin
        handles = [sys.stdin]
    else:
        # If "-" is among the inputs, read stdin as one of the sources
        for path in args.input:
            if path == "-":
                handles.append(sys.stdin)
            else:
                try:
                    fh = open(path, "r")
                except OSError as e:
                    print(f"Error opening {path}: {e}", file=sys.stderr)
                    continue

                handles.append(fh)
    # Process records and write matches to stdout
    matches_found = 0
    for record in records_from_handles(handles):
        header_id = record.id
        header_desc = record.description
        if args.regex:
            # Regex-based match on ID or full header
            if any(
                p.search(header_id) or p.search(header_desc) for p in regex_patterns
            ):
                SeqIO.write(record, sys.stdout, "fasta")
                matches_found += 1
        elif args.partial:
            # Substring match on ID or full header
            if any(name in header_id or name in header_desc for name in names_list):
                SeqIO.write(record, sys.stdout, "fasta")
                matches_found += 1
        else:
            # Exact match on ID or full header
            if header_id in names_set or header_desc in names_set:
                SeqIO.write(record, sys.stdout, "fasta")
                matches_found += 1
    # Close any file handles we opened (not stdin)
    for h in handles:
        if h is not sys.stdin:
            h.close()
    if matches_found == 0:
        print("No matching sequences found.", file=sys.stderr)


if __name__ == "__main__":
    main()
