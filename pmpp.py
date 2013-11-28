from string_graph import OverlapVertex
from intervals import print_true_intervals, get_true_intervals
from visualize import print_repeat

# -----------------------------------------------------------------------------
# spurious edge removal

def delete_spurious_edges(g, conservative='very'):
	for e in g.edges:
		if spurious_connection(e, conservative):
			g.remove_edge(e)

def spurious_connection(e, conservative='very'):
	v1, v2 = e.v1, e.v2

	# temporary hack:
	# don't call spurious connections that can be contracted
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
		v1_wells = v1.get_head_wells()
	elif e.connection[v1] == 'T':
		v1_wells = v1.get_tail_wells()
	if e.connection[v2] == 'H':
		v2_wells = v2.get_head_wells()
	elif e.connection[v2] == 'T':
		v2_wells = v2.get_tail_wells()
	# v1_wells = get_wells(v1)
	# v2_wells = get_wells(v2)

	common_wells = v1_wells & v2_wells
	v1_fraction = len(common_wells) / float(len(v1_wells))
	v2_fraction = len(common_wells) / float(len(v2_wells))

	if conservative == 'very':
		# make sure this is an "extra" edge for each vertex

		unbalanced_v1 = True
		unbalanced_v2 = True
		if e.connection[v1] == 'H':
			if len(v1.head_edges) <= len(v1.tail_edges): unbalanced_v1 = False
		elif e.connection[v1] == 'T':
			if len(v1.tail_edges) <= len(v1.head_edges): unbalanced_v1 = False

		if e.connection[v2] == 'H':
			if len(v2.head_edges) <= len(v2.tail_edges): unbalanced_v2 = False
		elif e.connection[v2] == 'T':
			if len(v2.tail_edges) <= len(v2.head_edges): unbalanced_v2 = False

		# if (not unbalanced_v2) and (not unbalanced_v1):
		# 	return False

		if v1_fraction > 0.0 or v2_fraction > 0.0 \
		or len(v1_wells) < 4 or len(v2_wells) < 4:
			return False
		else:
			return True

	elif conservative == 'yes':
		if v1_fraction > 0.25 or v2_fraction > 0.25 \
		or len(v1_wells) < 3 or len(v2_wells) < 3:
			return False
		else:
			return True

	elif conservative == 'no':
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

def get_neighborhood(g, v, e):
	"""Returns all nodes that overlap with v and have a path through e."""
	
	# get other vertex
	w = e.other_vertex(v)

	# get part of w that overlaps with v:
	ovl = e.ovl_start[w], e.ovl_end[w]

	# get nodes that overlap with w using BFS
	visited = set([v, w])
	neighborhood = set([w])
	N_w, visited = visit_neighborhood(g, w, visited, ovl)

	neighborhood.update(N_w)
	return neighborhood

def visit_neighborhood(g, v, visited, ovl):
	"""Returns nodes that overlap with v at ovl.

	Performs a BFS starting at v.
	"""
	N = set()
	to_visit = list()
	# check if the neighbors of v intersect it at the desired overlap
	for e in v.edges:
		n = e.other_vertex(v)

		# skip visited vertices
		if n in visited:
			continue
		else:
			visited.add(n)

		# intersect ovl with overlap of n in v
		nv_ovl = e.ovl_start[v], e.ovl_end[v]

		if nv_ovl[0] <= ovl[0] <= nv_ovl[1]:
			common_ovl = (ovl[0], min(nv_ovl[1], ovl[1]))
		elif nv_ovl[0] <= ovl[1] <= nv_ovl[1]:
			common_ovl = (max(nv_ovl[0], ovl[0]), ovl[1])
		else:
			continue

		# if we are still here, the is a common overlap

		# add to set of neighbors:
		N.add(n)

		# determine the coordinates of common_ovl in n:
		start_shift = common_ovl[0] - nv_ovl[0]
		end_shift = nv_ovl[1] - common_ovl[1]
		shifted_ovl_in_n = e.ovl_start[n] + start_shift, e.ovl_end[n] - end_shift

		# mark this node to be visited later:
		to_visit.append((n, shifted_ovl_in_n))

	# now, visit all neighbors recusively:
	for n, n_ovl in to_visit:
		N_n, visited = visit_neighborhood(g, n, visited, n_ovl)
		N.update(N_n)
			
	# return all neighbors found and all nodes visited
	return N, visited

# -----------------------------------------------------------------------------
# repeat resolution

def resolve_repeats(g, V, wells='edges'):
	for v in V:
		if len(v.head_edges) > 1 or len(v.tail_edges) > 1:
			# pairs_to_resolve = try_to_resolve(v, g, wells=wells)
			# resolve_pairs(g, v, pairs_to_resolve)
			H_to_T = try_to_resolve_new(v, g)
			resolve_pairs_new(g, v, pairs_to_resolve)

def try_to_resolve_new(v, g):
	"""Function to match wells to each other."""
	# H = {e.other_vertex(v) for e in v.head_edges if e.v1 != e.v2}
	# T = {e.other_vertex(v) for e in v.tail_edges if e.v1 != e.v2}
	H = {e for e in v.head_edges if e.v1 != e.v2}
	T = {e for e in v.tail_edges if e.v1 != e.v2}

	# place the wells of neighbouring vertices into a dict, such that
	# H_wells_map = {v: v.get_wells() for v in H}
	# T_wells_map = {v: v.get_wells() for v in H}
	H_wells_map = get_wells_by_edge(v, H)
	T_wells_map = get_wells_by_edge(v, T)

	pairs_to_resolve = set()
	for e_h in H_wells_map:
		for e_t in T_wells_map:
			# if they have >=3 wells in common, connect them:
			if len(H_wells_map[e_h] & T_wells_map[e_t]) >= 4:
				pairs_to_resolve.add((e_h, e_t))

	return pairs_to_resolve

def try_to_resolve(v, g, wells='edges'):
	H = {e for e in v.head_edges if e.v1 != e.v2}
	T = {e for e in v.tail_edges if e.v1 != e.v2}

	# place the wells of neighbouring vertices into a dict, such that
	# H_wells_map[e] = wells(e.other_vertex(v))
	v_wells = v.get_wells()
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
		# if has_support(well, v):
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
		or (edge_pair_weights[sorted_pairs[0]] - edge_pair_weights[sorted_pairs[1]] >= 3):
			# and edge_pair_weights[sorted_pairs[1]] <= 3):
			pairs_to_resolve.add((e_tail, e_head))

	print '>>> RESOLVING:'
	wells_by_edge = get_wells_by_edge(v, v.edges)
	print_repeat(v, wells_by_edge)
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

# -----------------------------------------------------------------------------
# well-related helpers

def get_wells_by_edge(v, E):
	wells_by_edge = {e: set() for e in E}
	for e in E:
		w = e.other_vertex(v)
		if e.connection[w] == 'H':
			wells_by_edge[e] = w.get_head_wells()
		elif e.connection[w] == 'T':
			wells_by_edge[e] = w.get_tail_wells()
		else:
			exit("Error.")

	return wells_by_edge

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

	def get_well(ctg):
		assert ctg.startswith('well')
		fields = ctg.split('_')
		return fields[0][4:]
	
	well_counts = dict()
	for ctg in v.metadata['contigs']:
		well = get_well(ctg)
		if well not in well_counts:
			well_counts[well] = 0
		well_counts[well] += 1

	return w in well_counts