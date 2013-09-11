from string_graph import OverlapVertex
from visualize import print_true_intervals

def examine_repeats(g):
	MAX_V = 50
	for count, v in enumerate(g.vertices):
		if count >= MAX_V: break
		T = v.tail_edges
		H = v.head_edges

		if len(T) <= 1 or len(H) <= 1: continue

		print '================================'
		print 'Vertex %d: %d contigs, %d bp' % (v.id_, len(v.metadata['contigs']), len(v))

		T_wells = {well for e in T for well in get_head_wells(e.other_vertex(v))}
		H_wells = {well for e in H for well in get_tail_wells(e.other_vertex(v))}

		for e in T:
			w = e.other_vertex(v)
			wells = ','.join([str(well) for well in get_head_wells(w) if well in H_wells])
			print '%d: %d contigs\t %d bp\t wells %s' % (w.id_, len(w.metadata['contigs']), len(w), wells)
			for ctg in w.metadata['contigs'][:15]:
				if w.metadata['contig_starts'][ctg] < 100:
					print_true_coordinates(w, ctg)
					# fields = ctg.split('_')
					# print fields[0], fields[1], w.metadata['contig_starts'][ctg], w.metadata['contig_ends'][ctg]
			print

		print '-------------------------------'

		for e in H:
			w = e.other_vertex(v)
			wells = ','.join([str(well) for well in get_tail_wells(w) if well in T_wells])
			print '%d: %d contigs\t %d bp\t wells %s' % (w.id_, len(w.metadata['contigs']), len(w), wells)
			for ctg in w.metadata['contigs'][:15]:
				if w.metadata['contig_ends'][ctg] > len(w) - 100:
					print_true_coordinates(w, ctg)
					# fields = ctg.split('_')
					# print fields[0], fields[1], w.metadata['contig_starts'][ctg], w.metadata['contig_ends'][ctg]
			print

		print

def resolve_repeats(g):
	for v in g.vertices:
		if len(v.head_edges) > 1 or len(v.tail_edges) > 1:
			try_to_resolve(v, g)

def try_to_resolve(v, g):
	H = {e for e in v.head_edges if e.v1 != e.v2}
	T = {e for e in v.tail_edges if e.v1 != e.v2}

	# for e in (H | T):
	# 	assert v == e.v1 or v == e.v2

	v_wells = get_wells(v)
	T_wells = {well for e in T for well in get_head_wells(e.other_vertex(v))}
	H_wells = {well for e in H for well in get_tail_wells(e.other_vertex(v))}
	T_wells_map = {e: get_head_wells(e.other_vertex(v)) for e in T}
	H_wells_map = {e: get_tail_wells(e.other_vertex(v)) for e in H}

	T_support_map = get_unique_support(T_wells_map)
	H_support_map = get_unique_support(H_wells_map)

	# print [(well, e.id_) for well, e in H_support_map.iteritems() if e is not None]
	# print [(well, e.id_) for well, e in T_support_map.iteritems() if e is not None]

	# e_head <-> e_tail if exists well w
	# s.t. w occurs only in e_head in H
	# and only in e_tail in T
	H_to_T = {e: set() for e in H}
	T_to_H = {e: set() for e in T}
	for well in v_wells:
		# if has_support(w, v):
		e_head = H_support_map.get(well, None)
		e_tail = T_support_map.get(well, None)
		if e_head and e_tail:
			# print well, e_head.id_, e_tail.id_
			H_to_T[e_head].add(e_tail)
			T_to_H[e_tail].add(e_head)
	
	# only keep pairs that e_head <-> e_tail
	# in which both e have degree 1
	pairs_to_resolve = set()
	for e_head in H:
		# print e_head.id_, [f.id_ for f in H_to_T[e_head]]
		if len(H_to_T[e_head]) == 1:
			e_tail = list(H_to_T[e_head])[0]
			e_reverse = list(T_to_H[e_tail])[0]
			if len(T_to_H[e_tail]) == 1 and e_head == e_reverse:
				pairs_to_resolve.add((e_tail,e_head))
	
	if pairs_to_resolve:
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
			print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_ends'][ctg] > len(w) - 100])
		print 'T:'
		for e in T:
			w = e.other_vertex(v)
			wells = ','.join([str(well) for well in get_head_wells(w) if well in H_wells])
			print e.id_, '\t', wells
			print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_starts'][ctg] < 100])
			# for ctg in w.metadata['contigs'][:15]:
			# 	if w.metadata['contig_ends'][ctg] < 100:
			# 		print_true_coordinates(w, ctg)
		print

	for e_tail, e_head in pairs_to_resolve:
		if len(v.head_edges) > 1 or len(v.tail_edges) > 1:
			v_new = resolve(v, e_tail, e_head, g)
			print 'AFTER:'
			print 'H:'
			for e in v.head_edges:
				w = e.other_vertex(v)
				wells = ','.join([str(well) for well in get_tail_wells(w) if well in T_wells])
				print e.id_, '\t', wells
			print 'T:'
			for e in v.tail_edges:
				w = e.other_vertex(v)
				wells = ','.join([str(well) for well in get_head_wells(w) if well in H_wells])
				print e.id_, '\t', wells
			print '+'
			w = e_head.other_vertex(v_new)
			wells = ','.join([str(well) for well in get_tail_wells(w) if well in T_wells])
			print e_head.id_, '\t', wells
			w = e_tail.other_vertex(v_new)
			wells = ','.join([str(well) for well in get_head_wells(w) if well in H_wells])
			print e_tail.id_, '\t', wells


def resolve(v, e_tail, e_head, g):
	v_new = OverlapVertex(g.vertex_id_generator.get_id(), str(v.seq))
	v_new.head_edges.add(e_head)
	v_new.tail_edges.add(e_tail)
	g.add_vertex(v_new)

	print "resolving", v.id_, v_new.id_, e_tail.id_, e_head.id_

	v_new.metadata = v.metadata.copy()

	e_head.replace(v, v_new)
	e_tail.replace(v, v_new)

	v.head_edges.remove(e_head)
	v.tail_edges.remove(e_tail)

	return v_new

def get_unique_support(edges_to_wells):
	""" 
	Returns well -> edge dict D where D[w] = e means w occurs only in e. 
	"""
	unique_wells_to_edges = dict()
	for e, W in edges_to_wells.iteritems():
		for w in W:
			if w not in unique_wells_to_edges:
				unique_wells_to_edges[w] = e
			else:
				unique_wells_to_edges[w] = None

	return unique_wells_to_edges

def has_support(w, v):
	# head_ctgs = [ctg in v.metadata['contigs'] if v.metadata['contig_starts'][ctg] < 2000]
	# tail_ctgs = [ctg in v.metadata['contigs'] if v.metadata['contig_starts'][ctg] > len(v) - 2000]
	# wells = {get_well(ctg) for ctg in head_ctgs} | {get_well(ctg) for ctg in tail_ctgs}
	
	well_counts = dict()
	for ctg in v.metadata['contigs']:
		well = get_well(ctg)
		if well not in well_counts:
			well_counts[well] = 0
		well_counts += 1

	return well_counts[w] >= 1

def get_wells(v):
	return {get_well(ctg) for ctg in v.metadata['contigs']}

def get_head_wells(v):
	return {get_well(ctg) for ctg in v.metadata['contigs']
			if v.metadata['contig_starts'][ctg] < 100}

def get_tail_wells(v):
	len_v = len(v)
	return {get_well(ctg) for ctg in v.metadata['contigs']
			if v.metadata['contig_ends'][ctg] > len_v - 100}

def get_well(ctg):
	assert ctg.startswith('well')
	fields = ctg.split('_')
	return fields[0][4:]

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