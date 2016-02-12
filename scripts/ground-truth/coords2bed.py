#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--coords', required=True)
parser.add_argument('-b', '--bed', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

out = open(args.bed, 'w')

with open(args.coords) as f:
  f.readline()
  f.readline()
  for line in f:
    fields = line.split('|')
    start_ref, end_ref = fields[0].strip().split()
    start_ctg, end_ctg = fields[1].strip().split()
    len_ref, len_ctg = fields[2].strip().split()
    idty = float(fields[3])
    ref_name, ctg_name = fields[4].strip().split()
    ref_name = ref_name.split('_')[0]

    if idty < 95.00: continue

    out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ref_name, start_ref, end_ref, ctg_name, start_ctg, end_ctg, len_ref, len_ctg))
