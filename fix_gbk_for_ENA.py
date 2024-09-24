import sys
from Bio.GenBank.Record import Record
from Bio import SeqIO
from Bio.Seq import translate
from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
from io import StringIO
import argparse

# from pprint import pprint
parser = argparse.ArgumentParser(
    description='Fixes potential errors for submission of Chromosomes/Scaffolds to ENA. Writes to an EMBL file.'
)
parser.add_argument('gbk_to_fix', help='Path to gbk file to fix for ENA')
parser.add_argument('-f', '--topology', choices=['linear', 'circular'])
parser.add_argument('embl_out', help='Path to EMBL output file')
args = parser.parse_args()
modified_records = []
for record in SeqIO.parse(args.gbk_to_fix, "genbank"):
    if args.topology == 'linear':
        record.annotations["topology"] = "linear"
    if args.topology == 'circular':
        record.annotations["topology"] = "circular"
    for index, feature in enumerate(record.features):
        # Fixes partial coordinates by adding "<" or ">"
        if 'partial' in feature.qualifiers:
            sequence_expected = feature.extract(record.seq).translate(
                to_stop=True, table="Bacterial"
            )
            if len(
                record.features
            ) <= 3:  # If the contig is small and I have only one gene that is partial, then check manually what to do
                print(
                    f"Please check manually:",
                    feature.qualifiers['locus_tag'],
                    feature.location,
                )
                continue

            elif index == 1 or index == 2:
                feature.location = FeatureLocation(
                    BeforePosition(feature.location.start),
                    feature.location.end,
                    strand=feature.location.strand,
                )
            elif index == len(record.features) - 1 or index == len(record.features) - 2:
                feature.location = FeatureLocation(
                    feature.location.start,
                    AfterPosition(feature.location.end),
                    strand=feature.location.strand,
                )
            else:
                pass
        #        try:
        #            print(
        #                feature.qualifiers['locus_tag'],
        #                feature.location,
        #                f"Total length {len(record.seq)}",
        #            )
        #        except:
        #           pass
        if 'old_locus_tag' in feature.qualifiers:
            del feature.qualifiers['old_locus_tag']
        if 'protein_id' in feature.qualifiers:
            del feature.qualifiers['protein_id']
        if 'blastp_file' in feature.qualifiers:
            del feature.qualifiers['blastp_file']
        if 'blastx_file' in feature.qualifiers:
            del feature.qualifiers['blastx_file']
        if 'ID' in feature.qualifiers:
            del feature.qualifiers['ID']
        if 'Parent' in feature.qualifiers:
            del feature.qualifiers['Parent']
        if 'Name' in feature.qualifiers:
            feature.qualifiers['gene'] = feature.qualifiers.pop('Name')
        if 'pseudogene' in feature.qualifiers:
            if feature.qualifiers['pseudogene'] == ["Unknown"]:
                feature.qualifiers['pseudogene'] = ["unknown"]
        if 'score' in feature.qualifiers:
            del feature.qualifiers['score']
        if 'color' in feature.qualifiers:
            del feature.qualifiers['color']
        if 'eC_number' in feature.qualifiers:
            feature.qualifiers['EC_number'] = feature.qualifiers.pop('eC_number')
        if 'Note' in feature.qualifiers:
            feature.qualifiers['note'] = feature.qualifiers.pop('Note')
        if feature.type == "CDS" and 'partial' in feature.qualifiers:  # Pass this if/elif/else if partial, since translation will be incomplete and thus not a proper CDS.
            sequence_expected = str(
                feature.extract(record.seq).translate(to_stop=True, table="Bacterial")
            )
            if 'translation' not in feature.qualifiers:
                feature.qualifiers['translation'] = [sequence_expected]
            if feature.qualifiers['translation'] != [sequence_expected]:
                feature.qualifiers['translation'] = [sequence_expected]
            del feature.qualifiers['partial']
        elif feature.type == "CDS" and not 'partial' in feature.qualifiers:
            sequence_expected = str(
                feature.extract(record.seq).translate(
                    to_stop=True, table="Bacterial", cds=True
                )
            )
            if 'translation' not in feature.qualifiers:
                feature.qualifiers['translation'] = [sequence_expected]
            if feature.qualifiers['translation'] != [sequence_expected]:
                feature.qualifiers['translation'] = [sequence_expected]
            else:
                pass
        if feature.type == "misc_RNA":
            if 'accession' in feature.qualifiers:
                feature.qualifiers[
                    'db_xref'
                ] = f"RFAM:{feature.qualifiers.pop('accession')}"
    # Store the modified record
    modified_records.append(record)
# Write the modified records to an in-memory string (as GenBank format)
fixed_gbk = StringIO()
SeqIO.write(modified_records, fixed_gbk, "genbank")
fixed_gbk.seek(0)  # Reset the StringIO object for reading
SeqIO.convert(fixed_gbk, 'genbank', args.embl_out, 'embl')
