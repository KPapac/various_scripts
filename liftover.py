#!/usr/bin/env python3
#"""
#Lift a contig-level GFF3 onto scaffold coordinates using an AGP file.
#Also emit 'gap' features for N-rows in the AGP.
#
#Usage:
#  python lift_gff_to_scaffold.py --agp example_scaffolds.agp --gff example_annotations.gff3 --fasta example_scaffolds.fa --out scaffold_annotations.gff3
#
#Notes:
#- Supports AGP v2.1 W (component) and N (gap) rows.
#- Assumes GFF3 features are on contig seqids that match AGP component_id values.
#- Preserves feature rows and attributes; reprojects coordinates and flips strand if the component orientation is '-'.
#- Adds one 'gap' feature per AGP N-row with attributes: ID=gapX;Gap_type=...;Length=...
#"""
import argparse
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(
        description="Lift contig GFF to scaffold using AGP and add gap features."
    )
    ap.add_argument("--agp", required=True, help="AGP file path")
    ap.add_argument("--gff", required=True, help="Contig-level GFF3 file")
    ap.add_argument(
        "--fasta",
        required=False,
        help="(Optional) Scaffold FASTA to get lengths for sequence-region pragmas",
    )
    ap.add_argument("--out", required=True, help="Output GFF3 path")
    return ap.parse_args()


def read_fasta_lengths(path):
    if not path:
        return {}

    lengths = {}
    try:
        with open(path) as f:
            name = None
            L = 0
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if name is not None:
                        lengths[name] = L
                    name = line[1:].split()[0]
                    L = 0
                else:
                    L += len(line)
            if name is not None:
                lengths[name] = L
    except FileNotFoundError:
        return {}

    return lengths


def parse_agp(agp_path):
    """Return:
       - components: dict[scaffold] -> list of dict rows for 'W' and 'N' parts in ascending object_beg
       - contig_to_spans: dict[contig] -> list of (scaffold, S, E, cb, ce, ori)
    """
    components = defaultdict(list)
    contig_to_spans = defaultdict(list)
    with open(agp_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                # try spaces
                cols = line.rstrip("\n").split()
                if len(cols) < 9:
                    raise ValueError(f"AGP line has <9 columns: {line}")

            obj, obeg, oend, part, ctype = cols[0], int(cols[1]), int(cols[2]), cols[
                3
            ], cols[
                4
            ]
            row = {
                "object": obj,
                "object_beg": obeg,
                "object_end": oend,
                "part_number": part,
                "component_type": ctype,
            }
            if ctype.upper() == "W":
                row.update(
                    {
                        "component_id": cols[5],
                        "component_beg": int(cols[6]),
                        "component_end": int(cols[7]),
                        "orientation": cols[8],
                    }
                )
                contig_to_spans[row["component_id"]].append(
                    (
                        obj,
                        obeg,
                        oend,
                        row["component_beg"],
                        row["component_end"],
                        row["orientation"],
                    )
                )
            elif ctype.upper() == "U":
                # U gap: 6=gap_length, 7=gap_type, 8=linkage or 'yes/no', 9=linkage_evidence (may be omitted in some flavors)
                row.update(
                    {
                        "gap_length": int(cols[5]),
                        "gap_type": cols[6],
                        "linkage": cols[7] if len(cols) > 7 else ".",
                        "linkage_evidence": cols[8] if len(cols) > 8 else ".",
                    }
                )
            else:
                # other components not supported in this minimal example
                pass
            components[obj].append(row)
    # Ensure sorted by object_beg
    for k in components:
        components[k].sort(key= lambda r: r["object_beg"])
    return components, contig_to_spans


def parse_gff(gff_path):
    headers = []
    records = []
    with open(gff_path) as g:
        for line in g:
            if not line.strip():
                continue

            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                # be lenient: skip bad rows
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            records.append(
                {
                    "seqid": seqid,
                    "source": source,
                    "type": ftype,
                    "start": start,
                    "end": end,
                    "score": score,
                    "strand": strand,
                    "phase": phase,
                    "attributes": attrs,
                }
            )
    return headers, records


def flip_strand(s):
    return {"+": "-", "-": "+", ".": "."}.get(s, s)


def lift_record(rec, span):
    """Project a contig-based GFF record rec onto scaffold using AGP span:
       span = (scaffold, S, E, cb, ce, ori)
       Returns a new record dict with updated seqid/start/end/strand.
    """
    scaffold, S, E, cb, ce, ori = span
    start, end = rec["start"], rec["end"]
    if ori == "+":
        new_start = S + (start - cb)
        new_end = S + (end - cb)
        new_strand = rec["strand"]
    else:
        # reverse
        new_start = S + (ce - end)
        new_end = S + (ce - start)
        new_strand = flip_strand(rec["strand"])
    if new_start > new_end:
        new_start, new_end = new_end, new_start
    lifted = rec.copy()
    lifted["seqid"] = scaffold
    lifted["start"] = new_start
    lifted["end"] = new_end
    lifted["strand"] = new_strand
    return lifted


def write_gff(out_path, headers, seqregion, gap_features, lifted_records):
    with open(out_path, "w") as out:
        # gff version header
        wrote_version = any(h.startswith("##gff-version") for h in headers)
        if not wrote_version:
            out.write("##gff-version 3\n")
        for h in headers:
            out.write(h.rstrip("\n") + "\n")
        # sequence regions for scaffolds
        for seqid, L in seqregion.items():
            out.write(f"##sequence-region {seqid} 1 {L}\n")
        # gaps
        for gf in gap_features:
            out.write(
                "\t".join(
                    [
                        gf["seqid"],
                        gf["source"],
                        gf["type"],
                        str(gf["start"]),
                        str(gf["end"]),
                        gf["score"],
                        gf["strand"],
                        gf["phase"],
                        gf["attributes"],
                    ]
                ) +
                "\n"
            )
        # lifted
        for r in lifted_records:
            out.write(
                "\t".join(
                    [
                        r["seqid"],
                        r["source"],
                        r["type"],
                        str(r["start"]),
                        str(r["end"]),
                        r["score"],
                        r["strand"],
                        r["phase"],
                        r["attributes"],
                    ]
                ) +
                "\n"
            )


def main():
    args = parse_args()
    components, contig_spans = parse_agp(args.agp)
    headers, gff_records = parse_gff(args.gff)
    scaffold_lengths = read_fasta_lengths(args.fasta)
    # Build gap features from AGP N-rows
    gap_features = []
    gap_id_counter = 1
    for scaffold, rows in components.items():
        for row in rows:
            if row["component_type"].upper() == "U":
                gf = {
                    "seqid": scaffold,
                    "source": "ragtag",
                    "type": "sequence_gap",
                    "start": row["object_beg"],
                    "end": row["object_end"],
                    "score": ".",
                    "strand": ".",
                    "phase": ".",
                    "attributes": f"ID=gap_{scaffold}_{gap_id_counter};gap_type={row.get('gap_type','scaffold')};estimated_length=unknown;length={row.get('gap_length','.')};linkage={row.get('linkage','.')};linkage_evidence={row.get('linkage_evidence','.')}",
                }
                gap_features.append(gf)
                gap_id_counter += 1
    # Lift records: only those whose seqid exists in contig_spans
    lifted_records = []
    for rec in gff_records:
        contig = rec["seqid"]
        if contig not in contig_spans:
            # pass through unmodified if not a contig in AGP (optional)
            continue

        for span in contig_spans[contig]:
            lifted_records.append(lift_record(rec, span))
    # If no FASTA, derive scaffold length from max object_end
    if not scaffold_lengths:
        for scaffold, rows in components.items():
            max_end = max(r["object_end"] for r in rows)
            scaffold_lengths[scaffold] = max_end
    write_gff(args.out, headers, scaffold_lengths, gap_features, lifted_records)


if __name__ == "__main__":
    sys.exit(main())
