from string_graph import OverlapVertex
from visualize import print_true_intervals, get_true_intervals
from intervals import overlaps

###############################################################################
## REMOVAL OF SPURIOUS EDGES

def print_connection(e):
	v1, v2 = e.v1, e.v2

	if e.connection[v1] == 'H':
		v1_wells = get_head_wells(v1)
	elif e.connection[v1] == 'T':
		v1_wells = get_tail_wells(v1)
	if e.connection[v2] == 'H':
		v2_wells = get_head_wells(v2)
	elif e.connection[v2] == 'T':
		v2_wells = get_tail_wells(v2)
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

def examine_connections(g, conservative=True):
	num_spurious = 0
	wrong_calls = 0
	for e in g.edges:
		v1, v2 = e.v1, e.v2

		I_v1 = get_true_intervals(v1, v1.metadata['contigs'])
		I_v2 = get_true_intervals(v2, v2.metadata['contigs'])

		if spurious_connection(e, conservative):
			num_spurious += 1

			wrong_call = False
			for i1 in I_v1:
				for i2 in I_v2:
					if overlaps(i1, i2):
						wrong_call = True
						break
			if wrong_call:
				print '>>> WRONG CALL:'
				print_connection(e)
				wrong_calls += 1

		# if true spurious connection, print:
		spurious = True
		for i1 in I_v1:
			for i2 in I_v2:
				if overlaps(i1, i2):
					spurious = False

		if spurious:
			if spurious_connection(e, conservative):
				print '>>> CALLED:'
			print_connection(e)

	print 'Called %d spurious edges' % num_spurious
	print 'Made %d wrong calls' % wrong_calls

def spurious_connection(e, conservative=True):
	v1, v2 = e.v1, e.v2

	# temporary hack:
	if (e.connection[v1] == 'H' and len(v1.head_edges) == 1):
		if e.connection[v2] == 'T' and len(v2.tail_edges) == 1 \
		or e.connection[v2] == 'H' and len(v2.head_edges) == 1:
			# print e.id_, 'path'
			return False
	elif (e.connection[v1] == 'T' and len(v1.tail_edges) == 1):
		if e.connection[v2] == 'T' and len(v2.tail_edges) == 1 \
		or e.connection[v2] == 'H' and len(v2.head_edges) == 1:
			# print e.id_, 'path'
			return False

	if e.connection[v1] == 'H':
		v1_wells = get_head_wells(v1)
	elif e.connection[v1] == 'T':
		v1_wells = get_tail_wells(v1)
	if e.connection[v2] == 'H':
		v2_wells = get_head_wells(v2)
	elif e.connection[v2] == 'T':
		v2_wells = get_tail_wells(v2)
	# v1_wells = get_wells(v1)
	# v2_wells = get_wells(v2)

	common_wells = v1_wells & v2_wells
	v1_fraction = len(common_wells) / float(len(v1_wells))
	v2_fraction = len(common_wells) / float(len(v2_wells))

	# len(v1_wells) + len(v2_wells) < 3 ---> 114 calls, 10 wrong, to FN
	# avg len: 30k (from 12k), 11 conn. comp.
	# below ---> 96 calls, 3 wrong
	# avg len: 26k, 9 conn. comp.
	# IDEA: do this in several steps, b/c larger contigs -> more wells

	if conservative:
		if v1_fraction > 0.25 or v2_fraction > 0.25 \
		or len(v1_wells) < 3 or len(v2_wells) < 3:
			return False
		else:
			return True
	else:
		print e.id_, v1_fraction, v2_fraction, min(len(v1_wells), len(v2_wells))

		if min(len(v1_wells), len(v2_wells)) == 1:
			if min(v1_fraction, v2_fraction) >= 0.25 \
			or len(v1_wells) + len(v2_wells) <= 3:				
				return False
			else:
				return True
		else:
			if v1_fraction >= 0.25 or v2_fraction >= 0.25 \
			or len(v1_wells) + len(v2_wells) <= 3:				
				return False
			else:
				return True

def delete_spurious_edges(g, conservative=True):
	for e in g.edges:
		if spurious_connection(e, conservative):
			g.remove_edge(e)


###############################################################################
## REPEAT RESOLUTION

def examine_repeats(g):
	MAX_V = 50
	for count, v in enumerate(g.vertices):
		# if count >= MAX_V: break
		T = v.tail_edges
		H = v.head_edges

		if len(T) == 0 or len(H) == 0: continue
		if len(T) + len(H) <= 2: continue

		print_repeat(v)

def print_repeat(v):
	wells_by_edge = get_wells_by_edge(v, v.edges)

	print '================================'
	print 'Vertex %d: %d contigs, %d bp' % (v.id_, len(v.metadata['contigs']), len(v))
	print 'V_wells:', ','.join([str(well) for well in get_wells(v)])

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

def resolve_repeats(g, wells='edges'):
	for v in g.vertices:
		if len(v.head_edges) > 1 or len(v.tail_edges) > 1:
			try_to_resolve(v, g, wells=wells)

def get_wells_by_edge(v, E):
	wells_by_edge = {e: set() for e in E}
	for e in E:
		w = e.other_vertex(v)
		if e.connection[w] == 'H':
			wells_by_edge[e] = get_head_wells(w)
		elif e.connection[w] == 'T':
			wells_by_edge[e] = get_tail_wells(w)
		else:
			exit("Error.")

	return wells_by_edge

def try_to_resolve(v, g, wells='edges'):
	H = {e for e in v.head_edges if e.v1 != e.v2}
	T = {e for e in v.tail_edges if e.v1 != e.v2}

	# place the wells of neighbouring vertices into a dict, such that
	# H_wells_map[e] = wells(e.other_vertex(v))
	v_wells = get_wells(v)
	H_wells_map = get_wells_by_edge(v, H)
	T_wells_map = get_wells_by_edge(v, T)
	common_wells = [well for e in H for well in H_wells_map[e]] + \
				   [well for e in T for well in T_wells_map[e]]

	# calculate which wells occur uniquely on edges on either side
	# i.e. H_support_map[w] = e iff e is the only edge in H that contains w
	T_support_map = get_unique_support(T_wells_map)
	H_support_map = get_unique_support(H_wells_map)

	# compute weights for the edges in the graph H <-> T
	edge_pair_weights = {(e_head, e_tail): 0 for e_head in H for e_tail in T}
	# for well in v_wells:
	for well in common_wells:
		# if has_support(w, v):
		e_head = H_support_map.get(well, None)
		e_tail = T_support_map.get(well, None)
		if e_head and e_tail:
			edge_pair_weights[(e_head, e_tail)] += 1
	
	# connect a e_head -> e_tail pair if it is supported by >= 3 wells than
	# the second best option
	pairs_to_resolve = set()
	for edge_pair, weight in edge_pair_weights.iteritems():
		if weight == 0: continue
		e_head, e_tail = edge_pair
		alternative_pairs = {(e_h, e_t) for e_h in H for e_t in T if e_h == e_head or e_t == e_tail}
		sorted_pairs = sorted(alternative_pairs, key=lambda x: edge_pair_weights[x], reverse=True)

		if sorted_pairs[0] != edge_pair: continue
		if len(sorted_pairs) == 1 \
		or edge_pair_weights[sorted_pairs[0]] - edge_pair_weights[sorted_pairs[1]] >= 3:
			pairs_to_resolve.add((e_tail, e_head))

	print '>>> RESOLVING:'
	print_repeat(v)
	print {(e_h.id_, e_t.id_, edge_pair_weights[(e_h, e_t)]) for e_h in H for e_t in T}

	H_wells = {well for e in H for well in H_wells_map[e]}
	T_wells = {well for e in T for well in T_wells_map[e]}

	for e_tail, e_head in pairs_to_resolve:
		if len(v.head_edges) > 1 or len(v.tail_edges) > 1:
			v_new = resolve(v, e_tail, e_head, g)
			print 'AFTER:'
			print 'H:'
			for e in v.head_edges:
				if e not in H: continue
				wells = ','.join([str(well) for well in H_wells_map[e] if well in T_wells])
				print e.id_, '\t', wells
			print 'T:'
			for e in v.tail_edges:
				if e not in T: continue
				wells = ','.join([str(well) for well in T_wells_map[e] if well in H_wells])
				print e.id_, '\t', wells
			print '+'
			wells = ','.join([str(well) for well in H_wells_map[e_head] if well in T_wells])
			print e_head.id_, '\t', wells
			wells = ','.join([str(well) for well in T_wells_map[e_tail] if well in H_wells])
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
	return {get_well(ctg) for ctg in v.metadata['contig_starts']
			if v.metadata['contig_starts'][ctg] < 500}

def get_tail_wells(v):
	len_v = len(v)
	return {get_well(ctg) for ctg in v.metadata['contig_ends']
			if v.metadata['contig_ends'][ctg] > len_v - 500}

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