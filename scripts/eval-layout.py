#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-l', '--layout', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

def parse_intervals(ivls, ctg_len):
  # find top chromosome
  chroms = set([ivl[0] for ivl in ivls])
  chrom_len = {chrom :
              sum([ivl[2]-ivl[2] for ivl in ivls if ivl[0] == chrom])
              for chrom in chroms}
  chrom = sorted(list(chroms), key=lambda x:chrom_len[x],
                reverse=True)[0]

  # find ivls for this chrom
  ivls = [ivl for ivl in ivls if ivl[0] == chrom]
  ivls = sorted(ivls, key=lambda x:x[2]-x[1])

  # keep merging until we exceed contig length
  ilen = lambda x: x[2]-x[1]
  ivl = ivls.pop()
  while ivls:
    ivl0 = ivls.pop()
    merged_ivl = ( ivl[0], min(ivl[1],ivl0[1]),
                           max(ivl[2],ivl0[2]) )
    if ilen(merged_ivl) > ctg_len:
      break
    else:
      ivl = merged_ivl

  return ivl

def parse_ivl_str(ivls_str):
  ifields = ivls_str.split(',')
  ivls = list()
  for istr in ifields:
    f1, f2 = istr.split(':')
    s, e = f2.split('-')
    ivls.append( (int(f1), int(s), int(e)) )
  return ivls

def eval_ivls(ivls):
  # find largest ivl:
  ilen = lambda x: x[2]-x[1]
  ix, max_i = sorted(enumerate(ivls),
                    key=lambda x: ilen(x[1]))[-1]

  # take it as having the correct orientation
  # and look at other's contigs position w.r.t. to it
  n_correct1, n_wrong1 = eval_ivl_seq(ivls[ix:])
  n_correct2, n_wrong2 = eval_ivl_seq(ivls[:ix+1])

  n_correct = n_correct1 + n_correct2 - ilen(max_i)
  n_wrong = n_wrong1 + n_wrong2

  return n_correct, n_wrong

def comes_after(ivl1, ivl2):
  if ivl1[0] != ivl2[0]: return False
  if ivl1[1] < ivl2[1] and ivl2[1] > ivl1[2] - 5000:
    return True
  else:
    return False

def eval_ivl_seq(ivls):
  ilen = lambda x: x[2]-x[1]
  n_correct, n_wrong = 0, 0
  last_correct_ivl = ivls[0]
  for i, ivl in enumerate(ivls):
    if i == 0:
      n_correct += ilen(ivl)
      last_correct_ivl = ivls[0]
      continue
    if comes_after(last_correct_ivl, ivl):
      last_correct_ivl = ivl
      n_correct += ilen(ivl)
    else:
      n_wrong += ilen(ivl)

  return n_correct, n_wrong

bp_total = 0
bp_verifiable = 0
bp_correct = 0
ctg_intervals = dict()
with open(args.layout) as f:
  for line in f:
    fields = line.strip().split()
    ctg = int(fields[0])
    ctg_intervals[ctg] = list()
    for cstr in fields[1:]:
      ctg_fi = cstr.split(';')
      clen = int(ctg_fi[2])
      istr = ctg_fi[1]
      bp_total += clen

      if istr:
        ivls = parse_ivl_str(istr)
        ivl = parse_intervals(ivls, clen)
        ctg_intervals[ctg].append(ivl)

        n_correct1, n_wrong1 = eval_ivls(ctg_intervals[ctg])
        n_correct2, n_wrong2 = eval_ivls(ctg_intervals[ctg][::-1])
        if n_correct1 > n_correct2:
          n_correct, n_wrong = n_correct1, n_wrong1
        else:
          n_correct, n_wrong = n_correct2, n_wrong2
        bp_correct += n_correct
        bp_verifiable += n_correct + n_wrong

print bp_correct, bp_verifiable, bp_total




