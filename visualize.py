from intervals import get_intervals

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
		wells = ','.join([str(well) for well in get_tail_wells(w) if well in T_wells])
		print e.id_, '\t', wells
		print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_ends'][ctg] > len(w) - 1000])
	print 'T:'
	for e in T:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in get_head_wells(w) if well in H_wells])
		print e.id_, '\t', wells
		print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_starts'][ctg] < 1000])
		# for ctg in w.metadata['contigs'][:15]:
		# 	if w.metadata['contig_ends'][ctg] < 1000:
		# 		print_true_coordinates(w, ctg)
	print

def get_true_intervals(w, ctgs):
	I = list()
	for ctg in ctgs:
		fields = ctg.split('_')
		chrom, coords = fields[1].split(':')
		start, end = (int(x) for x in coords.split('-'))
		internal_start = int(fields[2])

		I.append((int(chrom), start + internal_start, start + internal_start + 200))

	return get_intervals(I)

def print_true_coordinates(w, ctg):
	fields = ctg.split('_')
	chrom, coords = fields[1].split(':')
	start, end = (int(x) for x in coords.split('-'))
	internal_start = int(fields[2])

	print chrom, start + internal_start, start + internal_start + 200

def print_true_intervals(w, ctgs):
	I = list()
	for ctg in ctgs:
		fields = ctg.split('_')
		chrom, coords = fields[1].split(':')
		start, end = (int(x) for x in coords.split('-'))
		internal_start = int(fields[2])

		I.append((int(chrom), start + internal_start, start + internal_start + 200))

	merged_I = get_intervals(I)
	for i in merged_I:
		print '\t', i

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

def to_graphviz_dot_with_intervals(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I = get_true_intervals(v, v.metadata['contigs'])
			if len(I) > 1:
				color = "red"
			else:
				color = "black"

			dot.write('%d [label = "%s", color="%s"]\n' % (v.id_, ','.join([str(i) for i in I]), color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id_, v2.id_, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')