#!/usr/bin/env python3
#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def write_nexus_sequential(
    alignment, output_file, molecule_type="DNA", add_sets_block=False
):
    ntax = len(alignment)
    nchar = alignment.get_alignment_length()
    with open(output_file, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN DATA;\n")
        f.write(f"    DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
        f.write(f"    FORMAT DATATYPE={molecule_type} MISSING=? GAP=- INTERLEAVE=NO;\n")
        f.write("    MATRIX\n")
        for record in alignment:
            f.write(f"{record.id.ljust(15)} {str(record.seq)}\n")
        f.write("    ;\nEND;\n")
        if add_sets_block:
            f.write("\nBEGIN SETS;\n")
            f.write("    [ Define your partitions below. Example: ]\n")
            f.write("    [ CHARSET gene1 = 1-300; ]\n")
            f.write("    [ CHARSET gene2 = 301-600; ]\n")
            f.write("END;\n")


def fasta_to_nexus(input_file, output_file, molecule_type, add_sets_block, interleave):
    try:
        # Read the alignment from FASTA
        alignment = AlignIO.read(input_file, "fasta")
        # Ensure all sequences are the same length (aligned)
        seq_lengths = {len(rec.seq) for rec in alignment}
        if len(seq_lengths) > 1:
            raise ValueError("All sequences must be the same length (aligned)")

        if interleave:
            # Use BioPython's native NEXUS writer (interleaved)
            AlignIO.convert(input_file, "fasta", output_file, "nexus", molecule_type)
            if add_sets_block:
                with open(output_file, "a") as f:
                    f.write("\nBEGIN SETS;\n")
                    f.write("    [ Define your partitions below. Example: ]\n")
                    f.write("    [ CHARSET gene1 = 1-300; ]\n")
                    f.write("    [ CHARSET gene2 = 301-600; ]\n")
                    f.write("END;\n")
        else:
            # Use custom writer for sequential format
            write_nexus_sequential(
                alignment, output_file, molecule_type, add_sets_block
            )
        print(
            f"✅ Converted '{input_file}' to NEXUS format ({'interleaved' if interleave else 'sequential'} mode)"
        )
    except Exception as e:
        print(f"❌ Error: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a FASTA alignment to NEXUS format."
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output NEXUS file")
    parser.add_argument(
        "--type",
        default="DNA",
        help="Molecule type: DNA, RNA, or protein (default: DNA)",
    )
    parser.add_argument(
        "--add-sets",
        action="store_true",
        help="Include an empty BEGIN SETS block with partition comments",
    )
    parser.add_argument(
        "--interleave",
        action="store_true",
        help="Write sequences in interleaved format (default is sequential)",
    )
    args = parser.parse_args()
    fasta_to_nexus(
        args.input, args.output, args.type.upper(), args.add_sets, args.interleave
    )


if __name__ == "__main__":
    main()
