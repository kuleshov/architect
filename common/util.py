import random

# ----------------------------------------------------------------------------
# stats

def n50(l):
  total_length = sum(l)
  half_length = total_length / 2
  length = 0

  l.sort(reverse=True)

  while length < half_length:
    x = l.pop()
    length += x

  return x

def nx(L,x):
  assert 0.0 <= x <= 1.0
  assert L
  l=list(L)
  total_length = float(sum(l))
  x_length = total_length * (1.0-x)
  length = 0

  l.sort(reverse=True)

  while length < x_length:
    t = l.pop()
    length += t

  return t


# ----------------------------------------------------------------------------
# dna manipulation

def complement(char):
  if char == 'A':
    return 'T'
  elif char == 'T':
    return 'A'
  elif char == 'C':
    return 'G'
  elif char == 'G':
    return 'C'
  elif char == 'N':
    return 'N'
  else:
    raise Exception("ERROR: Invalid DNA character: '%s'" % char)

def complement_string(string):
  return ''.join([complement(x) for x in string])

def reverse_complement(string):
  return complement_string(string)[::-1]

# ----------------------------------------------------------------------------
# misc

def peek(S):
  x = S.pop()
  S.add(x)
  return x

def weighted_choice(choices):
  total = sum(w for c, w in choices)
  r = random.uniform(0, total)
  upto = 0
  for c, w in choices:
    if upto + w >= r:
      return c
    upto += w
  assert False, "Shouldn't get here"

def reverse_string(s):
  return s[::-1]