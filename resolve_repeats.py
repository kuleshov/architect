def examine_repeats(g, overlap_store):
	for v in g.vertices:
		T = v.tail_edges

		print '----------------------------'
		print 'Vertex %d: %d contigs, %d bp' % (v.id_, len(v.metadata['contigs']), len(v))

		for e in T:
			w = e.other_vertex(v)
			print '%d: %d contigs\t %d bp\t %d ovl' % (w.id_, len(w.metadata['contigs']), len(w), e.ovl_end[w] - e.ovl_start[w]),
			print '\t%s%s \t %d v->w \t %d w->v' % (e.connection[v], e.connection[w], count_overlaping_fragments(v, w, overlap_store), count_overlaping_fragments(w, v, overlap_store))

		print

def count_overlaping_fragments(v, w, overlap_store):
	""" Counts number of contigs that overlap both v and w. """

	# Count number of contigs of v that overlap with a contig of w
	num_overlaping_contigs = 0
	WC = set(w.metadata['contigs'])
	for c in v.metadata['contigs']:
		if overlap_store[c] & WC:
			num_overlaping_contigs += 1

	# Count number of contigs of w that overlap with a contig of v
	# VC = set(v.metadata['contigs'])
	# for c in w.metadata['contigs']:
	# 	if overlap_store[c] & VC:
	# 		num_overlaping_contigs += 1

	return num_overlaping_contigs

def resolve_short_repeats(g):
	for v in g.vertices:
		e = overlap_decider(v, v.tail_edges)
		if e:
			for f in (v.tail_edges - set([e])):
				g.remove_edge(f)

		e = overlap_decider(v, v.head_edges)
		if e:
			for f in (v.head_edges - set([e])):
				g.remove_edge(f)

def overlap_decider(v, E):
	if len(E) < 2:
		return None

	edge_overlaps = list()

	for e in E:
		overlap = e.ovl_end[v] - e.ovl_start[v]
		edge_overlaps.append((e, overlap))

	edge_overlaps = sorted(edge_overlaps, key=lambda x: x[1], reverse=True)
	assert edge_overlaps[0][1] - edge_overlaps[1][1] >= 0

	if edge_overlaps[0][1] - edge_overlaps[1][1] > 50:
		return edge_overlaps[0][0]
	else:
		return None

def load_overlap_store(asqg_file):
	overlap_store = dict()
	with open(asqg_file) as asqg:
		for line in asqg:
			if not line.startswith('ED'): continue
		
			fields = line.strip().split()

			assert fields[0] == 'ED'

			c1, c2 = fields[1], fields[2]

			if c1 not in overlap_store:
				overlap_store[c1] = set()
			if c2 not in overlap_store:
				overlap_store[c2] = set()

			overlap_store[c1].add(c2)
			overlap_store[c2].add(c1)
	
	return overlap_store