from pmpp import spurious_connection, get_wells_by_edge
from intervals import get_true_intervals, overlaps
from visualize import print_connection, print_repeat

# -----------------------------------------------------------------------------
# visualization of spurious edges

def examine_connections(g, conservative='very'):
	num_spurious = 0
	wrong_calls = 0
	for e in g.edges:
		spurious = spurious_connection(e, conservative)
		valid = valid_connection(e)

		if spurious_connection(e, conservative):
			num_spurious += 1
			if valid:
				# false positive
				print '>>> WRONG CALL (FALSE POSITIVE):'
				print_connection(e)
				wrong_calls += 1

		# # false negative
		# if not valid and not spurious:
		# 	print '>>> NOT DETECTED (FALSE NEGATIVE):
		# 	print_connection(e)

	print 'Called %d spurious edges' % num_spurious
	print 'Made %d wrong calls' % wrong_calls

def valid_connection(e):
	"""Determines whether connection is valid.

	Note that may not be triggered in some cases when it should.
	For example, when we have A - B - C; A is a repat, B, C connect to C
	but C has no direct connection to A.
	"""
	v1, v2 = e.v1, e.v2
	I_v1 = get_true_intervals(v1, v1.metadata['contigs'])
	I_v2 = get_true_intervals(v2, v2.metadata['contigs'])

	wrong_call = False
	for i1 in I_v1:
		for i2 in I_v2:
			if overlaps(i1, i2):
				wrong_call = True
				break

	return wrong_call

# -----------------------------------------------------------------------------
# repeat visualization 

def examine_repeats(g):
	MAX_V = 50
	for count, v in enumerate(g.vertices):
		# if count >= MAX_V: break
		T = v.tail_edges
		H = v.head_edges

		if len(T) == 0 or len(H) == 0: continue
		if len(T) + len(H) <= 2: continue

		wells_by_edge = get_wells_by_edge(v, v.edges)
		print_repeat(v, wells_by_edge)