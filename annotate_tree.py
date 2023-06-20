#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('tree')
args = parser.parse_args()

with open('/space/no_backup/Kostas/wolbachia_clustering/data/meta/sample_dictionary') as inp1:
    metadata = [item.split(',') for item in inp1.read().splitlines()]

#with open('iqtree_partition_file.contree') as inp2:
with open(args.tree) as inp2:
    mytree=inp2.read()

for entry in metadata:
    mytree = re.sub(entry[0], f"{entry[1]}_{entry[0]}", mytree)

print(mytree)
# with open('annotated.nex', 'w') as out:
#     out.write(mytree)
