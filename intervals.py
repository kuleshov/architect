###############################################################################
## BASIC INTERNVAL FUNCTIONS

def merge_intervals(I):
	merged_I = list()
	sorted_I = sorted(I)
	current_i = sorted_I[0]

	for i in sorted_I[1:]:
		if overlaps_sorted(current_i, i):
			current_i = _merge(current_i, i)
		else:
			merged_I.append(current_i)
			current_i = i

	merged_I.append(current_i)

	return merged_I

def overlaps_sorted(a, b):
	# assumes a is before b
	return a[0] == b[0] and b[1] <= a[2]

def overlap(I1, I2, tol=0):
	for ivl1 in I1:
		for ivl2 in I2:
			if overlaps(ivl1, ivl2, tol):
				return True
	return False

def overlaps(a, b, tol=0):
	return a[0] == b[0] and \
	(
		(a[1] - tol <= b[1] <= a[2] + tol) or 
		(a[1] - tol <= b[2] <= a[2] + tol)
	)

def union(ivl1, ivl2):
	assert ivl1[0] == ivl2[0]
	start = min(ivl1[1], ivl2[1])
	end = max(ivl1[2], ivl2[2])
	return (ivl1[0], start, end)

def shift(ivl, shift):
	return ivl[0], ivl[1]+shift, ivl[2]+shift

def _merge(a, b):
  return (a[0], a[1], b[2])

###############################################################################
## PARSE INTERVALS ASSOCIATED WITH A VERTEX

def parse_interval(w, ctg):
	fields = ctg.split('_')
	chrom, coords = fields[1].split(':')
	start, end = (int(x) for x in coords.split('-'))

	# correction for how we generated the reads. bed format is 1-based
	start -= 1
	end -= 1

	internal_start = int(fields[2])
	return int(chrom), start + internal_start, start + internal_start + 199


def get_true_intervals(w, ctgs):
	I = list()
	for ctg in ctgs:
		I.append(parse_interval(w, ctg))

	if I:
		return merge_intervals(I)

def get_head_intervals(w, ctgs):
	I = list()
	head_ctgs = [ctg for ctg in  ctgs if w.metadata['contig_starts'][ctg] != -1 and \
										 w.metadata['contig_starts'][ctg] < 200 ]
	for ctg in head_ctgs:
		I.append(parse_interval(w, ctg))

	if I:
		return merge_intervals(I)
	else:
		return list()

def get_tail_intervals(w, ctgs):
	I = list()
	w_len = len(w.seq)
	tail_ctgs = [ctg for ctg in ctgs if w.metadata['contig_ends'][ctg] != -1 and \
		  							    w.metadata['contig_ends'][ctg] > w_len - 200 ]
	for ctg in tail_ctgs:
		I.append((int(chrom), start + internal_start, start + internal_start + 199))

	if I:
		return merge_intervals(I)
	else:
		return list()

def print_true_coordinates(w, ctg):
	print parse_interval(w, ctg)

def print_true_intervals(w, ctgs):
	I = get_true_intervals(w, ctgs)
	merged_I = merge_intervals(I)
	for i in merged_I:
		print '\t', i
