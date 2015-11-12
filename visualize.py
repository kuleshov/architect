import os
import sys

import intervals
from intervals import get_head_intervals, get_tail_intervals
from libkuleshov.dna import reverse_complement

# ----------------------------------------------------------------------------
# print vertex, edge information

def print_vertex(v):
	print '================================='
	print 'Vertex: %d' % v.id
	print 'True intervals:', v.intervals
	print 'wells:', ', '.join(['%d:%d-%d' % (w, v.well_interval(w)[0], v.well_interval(w)[1]) for w in v.wells])
	T_wells = set(v.tail_wells)
	H_wells = set(v.head_wells)
	print 'H:'
	for e in v.head_edges:
		w = e.other_vertex(v)
		w_wells = w.head_wells if e.connection[w] == 'H' else w.tail_wells
		wells = ', '.join(sorted(['%d:%d-%d' % (well, v.well_interval(well)[0], v.well_interval(well)[1]) 
														for well in w_wells if well in H_wells]))
		print e.id, len(w_wells), wells
	print 'T:'
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
	# v1_wells = get_wells(v1)
	# v2_wells = get_wells(v2)

	print '================================'
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
	print '================================'
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
						# print_vertex(v1)
						# print_vertex(v2)

def visualize_assembly(v, I, genome):
	print v.seq
	for i, ctg in enumerate(v.metadata['contigs']):
		print ctg, v.metadata['contig_starts'][ctg], v.metadata['contig_ends'][ctg]
		print genome[I[i][0]][I[i][1]:I[i][2]+1]

def visualize_connection(e):
	v1, v2 = e.v1, e.v2
	print e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation
	print v1.seq
	print v2.seq
	print v1.seq[e.ovl_start[v1]:e.ovl_end[v1]+1]
	print v2.seq[e.ovl_start[v2]:e.ovl_end[v2]+1]


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
				exit("ERROR: Invalid orientation")

		dot.write('}\n')

def to_graphviz_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			dot.write('%d [label = "%d"]\n' % (v.id, len(v.seq)))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id, v2.id, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')

def to_graphviz_dot_with_intervals(g, dot_file=sys.stdout):
	# genome = load_genome()
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I = get_true_intervals(v, v.metadata['contigs'])
			color = "black"
			if len(I) > 1:
				color = "red"
			if v.id == 3764187:
				color = "green"

			dot.write('%d [label = "%d:%d:%s", color="%s"]\n' % (v.id, v.id, len(v), ','.join([str(i) for i in I]), color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			# if v1.id == 3706567 and v2.id == 3712515:
			# 	visualize_connection(e)
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %d %s%s %d"]\n' % (v1.id, v2.id, e.id, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')

def to_graphviz_dot_with_connections(g, resolver, dot_file=sys.stdout):
	# figure out which color to use for which edge
	COLORS = ['red', 'blue', 'green', 'cyan', 'yellow', 'magenta', 'purple', 'pink', 'chocolate']
	color_map = dict()
	for v in g.vertices:
		i = 0
		pairs_to_resolve = resolver(v, g)
		print v.id, len(pairs_to_resolve)
		for e_h, e_t in pairs_to_resolve:
			if v.id == 3754668:
				print e_h.id, e_t.id
			if e_h not in color_map:
				color_map[e_h] = list()
			if e_t not in color_map:
				color_map[e_t] = list()

			color_map[e_h].append(COLORS[i])
			color_map[e_t].append(COLORS[i])
			i += 1

	# print the graph
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I = get_true_intervals(v, v.metadata['contigs'])
			if len(I) > 1:
				color = "red"
			else:
				color = "black"

			dot.write('%d [label = "%d:%s", color="%s"]\n' % (v.id, v.id, ','.join([str(i) for i in I]), color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [color="%s" label= "%d %d %d %d %d %s%s %d"]\n' % (v1.id, v2.id, ':'.join(color_map.get(e, ['black'])), e.id, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')

def to_graphviz_dot_with_markup(g, V_sets, E_sets, dot_file=sys.stdout):
	COLORS = ['blue', 'green', 'cyan', 'yellow', 'magenta', 'purple', 'pink', 'chocolate']
	v_color_map = dict()
	for i, V in enumerate(V_sets):
		for v in V:
			v_color_map[v] = COLORS[i]

	e_color_map = dict()
	for i, E in enumerate(E_sets):
		for e in E:
			e_color_map[e] = COLORS[i]

	# print the graph
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I = get_true_intervals(v, v.metadata['contigs'])
			if len(I) > 1:
				color = "red"
			else:
				color = "black"

			dot.write('%d [label = "%d:%d:%s", color="%s"]\n' % (v.id, v.id, len(v), ','.join([str(i) for i in I]), v_color_map.get(v,color)))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [color="%s" label= "%d %d %d %d %d %s%s %d"]\n' % (v1.id, v2.id, e_color_map.get(e, 'black'), e.id, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')

def to_graphviz_dot_with_double_intervals(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I_head = get_head_intervals(v, v.metadata['contig_starts'].keys())
			I_tail = get_tail_intervals(v, v.metadata['contig_ends'].keys())
			if len(I_head) > 1 or len(I_tail) > 1:
				color = "red"
			else:
				color = "black"

			if I_head: 
				head_label = ','.join([str(i) for i in I_head])
			else:
				head_label = ''

			if I_tail: 
				tail_label = ','.join([str(i) for i in I_tail])
			else:
				tail_label = ''

			label = head_label + '\n' + tail_label + '\n' + str(len(v))
			dot.write('%d [label = "%s", color="%s"]\n' % (v.id, label, color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			# if v1.id == 3706567 and v2.id == 3712515:
			# 	visualize_connection(e)
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id, v2.id, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')