#!/usr/bin/python

import argparse
import pysam
import re

# -----------------------------------------------------------------------------
# arguments

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-c', '--containment', required=True)
parser.add_argument('-t', '--threshold', type=int, default=75)
parser.add_argument('-s', '--shift', type=int, default=0)
# parser.add_argument('-m', '--map', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# load map

# print 'Loading well map...'
#
# well_map = dict()
# n = 0
# with open(args.map) as mapfile:
# 	mapfile.readline()
# 	for line in mapfile:
# 		_, r, w = line.split()
# 		well_map[r] = int(w)
# 		n += 1
# 		if n % 1000000 == 0:
# 			print n
#
# print 'Well map loaded.'			

# -----------------------------------------------------------------------------
# bam -> containment

# open alignment bam file
samfile = pysam.Samfile(args.bam, "rb")

# print some stats
print "Mapped reads:", samfile.mapped
print "Unmpped reads:", samfile.unmapped

# # load contig names and lengths
# contigs = {name : length for name, length 
# 		   in zip(samfile.references, samfile.lengths)}

# parse file
contig_interval_map = dict()
contig_well_suport_counts = dict()

n = 0

for read in samfile:
  if not read.positions:
    continue

  if read.mapping_quality < 30:
    continue

  name = read.qname
  well = int(name.split('_')[0][4:]) + args.shift
  # t = ':'.join(name.split(':')[3:])
  # if t not in well_map: continue
  # well = well_map[t]
  contig = samfile.getrname(read.tid)
  start, end = min(read.positions), max(read.positions)

  if contig not in contig_interval_map:
    contig_interval_map[contig] = dict()
    contig_well_suport_counts[contig] = dict()
  if well not in contig_interval_map[contig]:
    contig_interval_map[contig][well] = [[start], [end]]
    contig_well_suport_counts[contig][well] = 0

  start_positions, end_positions = contig_interval_map[contig][well]

  if start < start_positions[-1]:
    start_positions.append(start)
    start_positions.sort()
    start_positions = start_positions[:10]
  if end > end_positions[0]:
    end_positions.append(end)
    end_positions.sort()
    end_positions = end_positions[1:]
    # interval[1] = end

  contig_interval_map[contig][well] = start_positions, end_positions
  contig_well_suport_counts[contig][well] += 1

  n += 1
  if n % 1000000 == 0:
    print 'reads processed: ', n, samfile.mapped

# write containment
with open(args.containment, 'w') as c:
  for contig, interval_map in contig_interval_map.iteritems():
    for well, interval in interval_map.iteritems():
      if contig_well_suport_counts[contig][well] > args.threshold:
        start, end = interval[0][-1], interval[1][0]
        c.write('W\t%s\t%d\t%d\t%d\n' % (contig, well, start, end))

for contig, well_support_map in contig_well_suport_counts.iteritems():
  print contig
  for well, support in well_support_map.iteritems():
    print '\t', well, ':', support
