#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import itertools
import pickle
import pysam

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-f', '--fasta', required=True)
parser.add_argument('-e', '--edges', required=True)
parser.add_argument('-l', '--rlen', type=int, default=101)
parser.add_argument('-s', '--insert-size', type=int, default=200)

args = parser.parse_args()

# ----------------------------------------------------------------------------

bamfile = pysam.Samfile(args.bam, 'rb')
fasta = pysam.FastaFile(args.fasta)

ctg_lengths = {ctg: length for ctg, length in zip(fasta.references, fasta.lengths)}

# ----------------------------------------------------------------------------

# load full connections
R1 = dict()
R2 = dict()
n_ref = len(bamfile.references)
d = args.insert_size / 2
d_cutoff = 2*d + args.rlen
for i, ref in enumerate(bamfile.references):
  if i % 1000 == 0: print '%d/%d' % (i, n_ref)
  ref_len = ctg_lengths[ref]
  if ref_len > d_cutoff:
    reads = itertools.chain(
      bamfile.fetch(ref, 0, d),
      bamfile.fetch(ref, ref_len-d, ref_len)
    )
  else:
    # continue
    reads = bamfile.fetch(ref)
  for read in reads:
    if not read.is_paired: continue
    if read.mapq < 10: continue
    if read.reference_id != read.next_reference_id:
      name1 = bamfile.getrname(read.reference_id)
      name2 = bamfile.getrname(read.next_reference_id)
      len1 = ctg_lengths[name1]
      len2 = ctg_lengths[name2]
      start1 = read.reference_start
      start2 = read.next_reference_start

      # matches = sum([l for (o,l) in read.cigartuples if o == 0])
      # if matches < 20: continue
      tdict = R1 if read.is_read1 else R2

      if read.qname not in tdict:
        tdict[read.qname] = list()

      strand = 'R' if read.is_reverse else 'S'

      if ref_len <= d_cutoff:
        tdict[read.qname].append((name1, 'H', strand, start1))
        tdict[read.qname].append((name1, 'T', strand, start1))
      else:
        if start1 < d:
          tdict[read.qname].append((name1, 'H', strand, start1))
        elif start1 >= len1 - d:
          tdict[read.qname].append((name1, 'T', strand, start1))
        else:
          continue
          # print start1, len1, len1 - start1
          # exit(1)

simple_edge_counts = dict()
edge_dists = dict()
n_strange = 0
n_missing = 0
n_valid = 0
for read1, links in R1.iteritems():
  fi = read1.split('.')
  read2 = '.'.join(fi[:2]) + '.2'

  if read2 not in R2: 
    n_missing += 1
    continue

  for ctg1, conn1, strand1, start1 in links:
    for ctg2, conn2, strand2, start2 in R2[read2]:
      # note: we consider certain types of edges invalid
      if (conn1 == conn2 and strand1 != strand2) or \
         (conn1 != conn2 and strand1 == strand2):
        n_strange += 1
        continue

      simple_edge = frozenset([ctg1,ctg2])
      if simple_edge not in simple_edge_counts:
        simple_edge_counts[simple_edge] = 0
      simple_edge_counts[simple_edge] += 1

      strand = 'S' if strand1 == strand2 else 'R'

      # determine distance
      # we may choose keep d > 10 to simplify some downstream analysis
      if conn1 == 'H':
        o1 = start1 + args.rlen
      elif conn1 == 'T':
        o1 = ctg_lengths[ctg1] - start1
      if conn2 == 'H':
        o2 = start2
      elif conn2 == 'T':
        o2 = ctg_lengths[ctg2] - start2
      # d = max(10, args.insert_size - o1 - o2)
      d = args.insert_size - o1 - o2
      if o1 < 25 or o2 < 25: continue
      if o1 + o2 < args.insert_size + 6*25: continue

      link1 = (ctg1, conn1)
      link2 = (ctg2, conn2)
      edge_l = [link1, link2, strand]
      edge = frozenset(edge_l)
      if edge not in edge_dists:
        edge_dists[edge] = list()
      edge_dists[edge].append((edge_l, d))
      n_valid += 1

print 'Strange reads (skipped): %d' % n_strange
print 'Reads whose mate was not good: %d' % n_missing
print 'Read paired that were used: %d' % n_valid
print 'Connections established: %d' % len(edge_dists)

# this is for validation later on
pickle.dump(simple_edge_counts, open('pe_accurate_edge_counts.pkl', 'wb'))

# write down the edges:
with open(args.edges, 'w') as out:
  for dists in edge_dists.values():
    link1, link2, strand = dists[0][0]
    ctg1, conn1 = link1
    ctg2, conn2 = link2

    count = len(dists)
    avg_dist = sum([float(d[1]) for d in dists]) / count

    # write result
    out.write('{ctg1}\t{ctg2}\t{conn1}\t{conn2}\t{strand}\t{count}\t{dist}\n'.format(
      ctg1=ctg1,ctg2=ctg2,conn1=conn1,conn2=conn2,strand=strand,count=count,dist=avg_dist)
    )

