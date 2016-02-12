#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bed', required=True)
parser.add_argument('-c', '--csv', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

intervals = dict()
with open(args.bed) as f:
  for line in f:
    fields = line.strip().split()
    ctg = fields[3]
    if ctg not in intervals: intervals[ctg] = list()
    # ivl_str = '%s-%s-%s' % (fields[0], fields[1], fields[2])
    ivl = fields[:3]
    intervals[ctg].append((fields[0], int(fields[1]), int(fields[2])))

with open(args.csv, 'w') as out:
  for ctg, ivls in intervals.iteritems():
    # print ivls
    ivls = [i for i in ivls if i[2]-i[1] > 1000]
    sorted_ivls = sorted(ivls, key=lambda x: x[2]-x[1], reverse=True)
    ivl_strs = ['%s-%d-%d' % (i[0], i[1], i[2]) for i in sorted_ivls]
    ivl_str = ':'.join(ivl_strs[:3])
    out.write('%s,%s\n' % (ctg, ivl_str))
