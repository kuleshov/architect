#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from libkuleshov.stats import n50

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
  n_correct1, n_wrong1, profile1 = eval_ivl_seq(ivls[ix:],dir='fwd')
  n_correct2, n_wrong2, profile2 = eval_ivl_seq(ivls[:ix+1],dir='bwd')

  n_correct = n_correct1 + n_correct2 - ilen(max_i)
  n_wrong = n_wrong1 + n_wrong2
  print profile1, profile2
  profile = profile2 + profile1[1:]

  return n_correct, n_wrong, profile

def comes_after(ivl1, ivl2):
  if ivl1[0] != ivl2[0]: return False
  if ivl1[1] < ivl2[1] and ivl2[1] > ivl1[2] - 5000:
    return True
  else:
    return False

def comes_before(ivl1, ivl2):
  if ivl1[0] != ivl2[0]: return False
  if ivl2[2] < ivl1[2] and ivl2[1] < ivl1[2] + 5000:
    return True
  else:
    return False

def eval_ivl_seq(ivls, dir):
  ilen = lambda x: x[2]-x[1]
  n_correct, n_wrong = 0, 0
  profile = []
  if dir == 'fwd':
    print 'F',
    last_correct_ivl = ivls[0]
    for i, ivl in enumerate(ivls):
      print ivl,
      if i == 0:
        n_correct += ilen(ivl)
        last_correct_ivl = ivl
        profile.append('K')
        continue
      if comes_after(last_correct_ivl, ivl):
        last_correct_ivl = ivl
        n_correct += ilen(ivl)
        profile.append('K')
        print 'K',
      else:
        n_wrong += ilen(ivl)
        profile.append('E')
        print 'E',
    print n_correct, n_wrong
  elif dir == 'bwd':
    print 'B',
    for i, ivl in enumerate(reversed(ivls)):
      print ivl,
      if i == 0:
        n_correct += ilen(ivl)
        last_correct_ivl = ivl
        profile = ['K']
        continue
      if comes_before(last_correct_ivl, ivl):
        last_correct_ivl = ivl
        n_correct += ilen(ivl)
        profile = ['K'] + profile
        print 'K',
      else:
        n_wrong += ilen(ivl)
        profile = ['E'] + profile
        print 'E',
    print n_correct, n_wrong

  return n_correct, n_wrong, profile

bp_total = 0
bp_verifiable = 0
bp_correct = 0
profiles = list()
with open(args.layout) as f:
  for line in f:
    fields = line.strip().split()
    cluster_id = int(fields[0])
    cluster_ivls = list()
    cluster_len = 0

    all_lengths = list()
    all_ids = list()
    verifiable_ids = list()

    for cstr in fields[1:]:
      ctg_fi = cstr.split(';')
      cid = int(ctg_fi[0])
      clen = int(ctg_fi[2])
      istr = ctg_fi[1]
      bp_total += clen
      cluster_len += clen

      all_ids.append(cid)
      all_lengths.append(clen)

      if istr:
        verifiable_ids.append(cid)
        ivls = parse_ivl_str(istr)
        ivl = parse_intervals(ivls, clen)
        cluster_ivls.append(ivl)

    n_correct, n_wrong = 0,0
    if cluster_ivls:
      print cluster_ivls
      n_correct1, n_wrong1, profile1 = eval_ivls(cluster_ivls)
      n_correct2, n_wrong2, profile2 = eval_ivls(cluster_ivls[::-1])
      print 'p1', profile1
      print 'p2', profile2
      if n_correct1 > n_correct2:
        n_correct, n_wrong, profile = n_correct1, n_wrong1, profile1
      else:
        n_correct, n_wrong, profile = n_correct2, n_wrong2, profile2[::-1]
    bp_correct += n_correct
    bp_verifiable += n_correct + n_wrong
    profile_dict = dict(zip(verifiable_ids, profile))
    full_profile = [(cid, clen, profile_dict.get(cid,'N'))
                   for cid, clen in zip(all_ids, all_lengths)]
    profiles.append(full_profile)
    print cluster_len, n_correct, n_wrong

print bp_correct, bp_verifiable, bp_total

# n50 calculations
bp_total = 0
uncorrected_lengths = \
  [sum([clen for (cid, clen, cprof) in profile]) for profile in profiles]
lengths = list()
for profile in profiles:
  bp_total += sum([clen for (cid, clen, cprof) in profile])
  curr_len = 0
  for _, clen, cprof in profile:
    if cprof in ('K', 'N'):
      curr_len += clen
    elif cprof == 'E':
      if curr_len:
        lengths.append(curr_len)
      curr_len = 0
      lengths.append(clen)

  if curr_len:
    lengths.append(curr_len)

print bp_total, n50(uncorrected_lengths), n50(lengths)





