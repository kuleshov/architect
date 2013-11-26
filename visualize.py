import os
import sys

from intervals import print_true_intervals, get_true_intervals, get_head_intervals, get_tail_intervals
from libkuleshov.dna import reverse_complement

###############################################################################
## PRINT INFORMATION FOR A VERTEX, EDGE

# this probably was moved here and never tested: get_tail_wells is not imported
def print_vertex(x):
	print '================================='
	print 'Vertex: %d' % v.id_
	print [(e_tail.id_, e_head.id_) for (e_tail, e_head) in pairs_to_resolve]
	print ','.join([str(w) for w in v_wells])
	print_true_intervals(v, v.metadata['contigs'])
	print
	print 'BEFORE:'
	print 'H:'
	for e in H:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in w.get_tail_wells() if well in T_wells])
		print e.id_, '\t', wells
		print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_ends'][ctg] > len(w) - 1000])
	print 'T:'
	for e in T:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in w.get_head_wells() if well in H_wells])
		print e.id_, '\t', wells
		print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_starts'][ctg] < 1000])
		# for ctg in w.metadata['contigs'][:15]:
		# 	if w.metadata['contig_ends'][ctg] < 1000:
		# 		print_true_coordinates(w, ctg)
	print

def print_connection(e):
	v1, v2 = e.v1, e.v2

	if e.connection[v1] == 'H':
		v1_wells = v1.get_head_wells()
	elif e.connection[v1] == 'T':
		v1_wells = v1.get_tail_wells()
	if e.connection[v2] == 'H':
		v2_wells = v2.get_head_wells()
	elif e.connection[v2] == 'T':
		v2_wells = v2.get_tail_wells()
	# v1_wells = get_wells(v1)
	# v2_wells = get_wells(v2)

	I_v1 = get_true_intervals(v1, v1.metadata['contigs'])
	I_v2 = get_true_intervals(v2, v2.metadata['contigs'])

	print '================================'
	print 'Edge %d: v1 overlap: %d %d %s, v2 overlap: %d %d %s, ori: %d' \
	% (e.id_, e.ovl_start[v1], e.ovl_end[v1], e.connection[v1],
			  e.ovl_start[v2], e.ovl_end[v2], e.connection[v2],
			  e.v2_orientation)
	print 'V1: id: %d ,%d contigs, %d bp' % (v1.id_, len(v1.metadata['contigs']), len(v1))
	print 'V1 wells:', ','.join([str(w) for w in v1_wells])
	print 'I1', I_v1
	print 'V2: id: %d ,%d contigs, %d bp' % (v2.id_, len(v2.metadata['contigs']), len(v2))
	print 'V2 wells:', ','.join([str(w) for w in v2_wells])
	print 'I2', I_v2

def print_repeat(v, wells_by_edge):
	print '================================'
	print 'Vertex %d: %d contigs, %d bp' % (v.id_, len(v.metadata['contigs']), len(v))
	print 'V_wells:', ','.join([str(well) for well in v.get_wells()])

	T_wells = {well for e in v.tail_edges for well in wells_by_edge[e]}
	H_wells = {well for e in v.head_edges for well in wells_by_edge[e]}

	for e in v.tail_edges:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in wells_by_edge[e] if well in H_wells])
		print '%d: %d contigs\t %d bp\t wells %s' % (w.id_, len(w.metadata['contigs']), len(w), wells)
		print_true_intervals(w, w.metadata['contigs']) # keep only edge contigs?

	print '-------------------------------'

	for e in v.head_edges:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in wells_by_edge[e] if well in T_wells])
		print '%d: %d contigs\t %d bp\t wells %s' % (w.id_, len(w.metadata['contigs']), len(w), wells)
		print_true_intervals(w, w.metadata['contigs']) # keep only edge contigs?

	print

###############################################################################
## VISUALIZE SUB-COMPONENTS OF THE GRAPH

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

###############################################################################
## SAVE TO DOT

def to_abyss_explorer_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		id_to_num = dict()
		for v in g.vertices:
			dot.write('"%d+" [l=%d C=1]\n' % (v.id_, len(v.seq)))
			dot.write('"%d-" [l=%d C=1]\n' % (v.id_, len(v.seq)))

		for e in g.edges:
			v1, v2, ori = e.v1, e.v2, e.v2_orientation
			if ori == 0:
				dot.write('"%d+" -> "%d+"\n' % (v1.id_, v2.id_))
			elif ori == 1:
				dot.write('"%d+" -> "%d-"\n' % (v1.id_, v2.id_))
			else:
				exit("ERROR: Invalid orientation")

		dot.write('}\n')

def to_graphviz_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			dot.write('%d [label = "%d"]\n' % (v.id_, len(v.seq)))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id_, v2.id_, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')

def to_graphviz_dot_with_intervals(g, dot_file=sys.stdout):
	# genome = load_genome()
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I = get_true_intervals(v, v.metadata['contigs'])
			# if v.id_== 3599553:
			# 	visualize_assembly(v, I, genome)
			if len(I) > 1:
				color = "red"
				# examine_repeats(I, genome)
				# if is_misassembly(v, I, genome):
					# print_containment(v, I, genome)
			else:
				color = "black"

			dot.write('%d [label = "%d:%s", color="%s"]\n' % (v.id_, v.id_, ','.join([str(i) for i in I]), color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			# if v1.id_ == 3706567 and v2.id_ == 3712515:
			# 	visualize_connection(e)
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %d %s%s %d"]\n' % (v1.id_, v2.id_, e.id_, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

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
			dot.write('%d [label = "%s", color="%s"]\n' % (v.id_, label, color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			# if v1.id_ == 3706567 and v2.id_ == 3712515:
			# 	visualize_connection(e)
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id_, v2.id_, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')