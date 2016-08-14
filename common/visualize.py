"""Helper functions for analyzing the assembly graph. 

These may be helpful to debug or examine assemblies.
"""

import os
import sys

import graph.intervals as intervals
from graph.intervals import get_head_intervals, get_tail_intervals
from common.util import reverse_complement

# ----------------------------------------------------------------------------
# print vertex, edge information

def print_vertex(v):
	print 'Vertex: %d' % v.id
	print 'True intervals:', v.intervals
	print 'Wells:', ', '.join(['%d:%d-%d' % (w, v.well_interval(w)[0], v.well_interval(w)[1]) for w in v.wells])
	T_wells = set(v.tail_wells)
	H_wells = set(v.head_wells)
	print 'Head edges:'
	for e in v.head_edges:
		w = e.other_vertex(v)
		w_wells = w.head_wells if e.connection[w] == 'H' else w.tail_wells
		wells = ', '.join(sorted(['%d:%d-%d' % (well, v.well_interval(well)[0], v.well_interval(well)[1]) 
														for well in w_wells if well in H_wells]))
		print e.id, len(w_wells), wells
	print 'Tail edges:'
	for e in v.tail_edges:
		w = e.other_vertex(v)
		w_wells = w.head_wells if e.connection[w] == 'H' else w.tail_wells
		wells = ', '.join(sorted(['%d:%d-%d' % (well, v.well_interval(well)[0], v.well_interval(well)[1]) 
														 for well in w_wells if well in T_wells]))
		print e.id, len(w_wells), wells
	print

def print_connection(e):
	v1, v2 = e.v1, e.v2

	if e.connection[v1] == 'H':
		v1_wells = v1.head_wells
	elif e.connection[v1] == 'T':
		v1_wells = v1.tail_wells
	if e.connection[v2] == 'H':
		v2_wells = v2.head_wells
	elif e.connection[v2] == 'T':
		v2_wells = v2.tail_wells

	if e.is_overlap_edge:
		print 'Edge %d: v1 overlap: %d %d %s, v2 overlap: %d %d %s, ori: %d' \
		% (e.id, e.ovl_start[v1], e.ovl_end[v1], e.connection[v1],
				  e.ovl_start[v2], e.ovl_end[v2], e.connection[v2],
				  e.orientation)
	elif e.is_scaffold_edge:
		print 'Edge %d: %s%s, ori: %d' % (e.id, e.connection[v1], e.connection[v2], e.orientation)
	print 'V1: id: %d, %d bp' % (v1.id, len(v1))
	print 'V1 wells:', ','.join(sorted([str(w) for w in v1_wells]))
	print 'V2: id: %d, %d bp' % (v2.id, len(v2))
	print 'V2 wells:', ','.join(sorted([str(w) for w in v2_wells]))

def print_repeat(v, wells_by_edge):
	print 'Vertex %d: %d contigs, %d bp' % (v.id, len(v.metadata['contigs']), len(v))
	print_true_intervals(v, v.metadata['contigs'])
	print 'V_wells:', ','.join([str(well) for well in v.get_wells()])

	T_wells = {well for e in v.tail_edges for well in wells_by_edge[e]}
	H_wells = {well for e in v.head_edges for well in wells_by_edge[e]}

	for e in v.tail_edges:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in wells_by_edge[e] if well in H_wells])
		print '%d: %d contigs\t %d bp\t wells %s' % (w.id, len(w.metadata['contigs']), len(w), wells)
		print_true_intervals(w, w.metadata['contigs']) # keep only edge contigs?

	print '-------------------------------'

	for e in v.head_edges:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in wells_by_edge[e] if well in T_wells])
		print '%d: %d contigs\t %d bp\t wells %s' % (w.id, len(w.metadata['contigs']), len(w), wells)
		print_true_intervals(w, w.metadata['contigs']) # keep only edge contigs?

	print

# ----------------------------------------------------------------------------
# visualize graph sub-components

def visualize_well_correctness(g):
	for v1 in g.vertices:
		if len(v1) < 5000: continue
		for v2 in g.vertices:
			if len(v2) < 5000: continue
			if v1.id >= v2.id: continue
			for conn1 in ('H', 'T'):
				for conn2 in ('H', 'T'):
					common_wells, v1_wells, v2_wells = _get_wells_between_v(v1, v2, conn1, conn2)
					if not v1_wells or not v2_wells: continue
					frac_common1 = float(len(common_wells)) / float(len(v1_wells))
					frac_common2 = float(len(common_wells)) / float(len(v2_wells))
					if len(common_wells) >= 4 and max(frac_common1, frac_common2) > 0.33 \
																		and min(frac_common1, frac_common2) > 0.33 :
						print
						print '='* 80
						print v1.id, v2.id, '%s%s' % (conn1, conn2), intervals.overlap(v1.intervals, v2.intervals)
						print v1.id, len(v1), v1.intervals
						print v1_wells
						print v2.id, len(v2), v2.intervals
						print v2_wells
						print common_wells

def _get_wells_between_v(v1, v2, conn1, conn2):
  if conn1 == 'H':
    v1_wells = v1.head_wells
  elif conn1 == 'T':
    v1_wells = v1.tail_wells
  if conn2 == 'H':
    v2_wells = v2.head_wells
  elif conn2 == 'T':
    v2_wells = v2.tail_wells

  return v1_wells & v2_wells, v1_wells, v2_wells

# ----------------------------------------------------------------------------
# create dot files

def to_graphviz_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			dot.write('%d [label = "%d"]\n' % (v.id, len(v.seq)))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [label= "%s%s %d"]\n' \
				% (v1.id, v2.id, e.connection[v1], e.connection[v2], e.orientation))

		dot.write('}\n')

def to_abyss_explorer_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		idto_num = dict()
		for v in g.vertices:
			dot.write('"%d+" [l=%d C=1]\n' % (v.id, len(v.seq)))
			dot.write('"%d-" [l=%d C=1]\n' % (v.id, len(v.seq)))

		for e in g.edges:
			v1, v2, ori = e.v1, e.v2, e.v2_orientation
			if ori == 0:
				dot.write('"%d+" -> "%d+"\n' % (v1.id, v2.id))
			elif ori == 1:
				dot.write('"%d+" -> "%d-"\n' % (v1.id, v2.id))
			else:
				raise Exception("ERROR: Invalid orientation")

		dot.write('}\n')
