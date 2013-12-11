from collections import deque

def compute_traversals(g):
	# compute neighborhoods
	compute_neighborhoods(g)

	# compute valid edge pairs
	validate_edge_pairs(g)

	# compute non-repeats (currently, longest vertices)
	sorted_V = sorted(g.vertices, reverse=True)
	v_start = sorted_V[0]

	# traverse components from non-repeats
	E = traverse(v_start)

	return E

def traverse(start_v):
	# build initial neighborhood around start_v
	# i.e. decide which edges you can follow
	start_E = set()
	for e in start_v.edges:
		if e.connection[start_v] == 'H':
			W = start_v.get_head_wells()
		elif e.connection[start_v] == 'T':
			W = start_v.get_tail_wells()

		N = start_v.metadata['neighborhood'][e]
		if validate_vertex(N, W):
			start_E.add(e)

	# follow the initial edges if well information is okay
	# perform a BFS to decide where to follow
	chosen = set(start_E)
	to_visit = deque(start_E)

	# while there are edges whose neighbors still need to be explored
	while to_visit:
		e = to_visit.popleft()
		E_v1, E_v2 = e.v1.metadata['connections'][e], e.v2.metadata['connections'][e]

		for E in (E_v1, E_v2):
			for e in E:
				if e not in chosen:
					to_visit.append(e)
				chosen.add(e)

	return chosen

def validate_edge_pairs(g):
	"""Determines which pairs around a vertex can be joined.

	Currently results stored in vertex metadata, but this should
	be done properly, once the algorithm is worked out.
	"""

	for v in g.vertices:
		connections = {e:set() for e in v.edges}
		for e in v.edges:
			u = e.other_vertex(v)
			if e.connection[u] == 'H':
				W = v.get_head_wells()
			elif e.connection[u] == 'T':
				W = v.get_tail_wells()

			for f in v.edges:
				if e == f: continue
				N = v.metadata['neighborhood'][f]
				if validate_vertex(N, W):
					connections[e].add(f)

		v.metadata['connections'] = connections

def get_wells_by_edge(v, E):
	"""Returns a map e -> W of wells present at e.other_vertex(v)

	Helper function.
	"""
	
	wells_by_edge = {e: set() for e in E}
	for e in E:
		w = e.other_vertex(v)
		if e.connection[w] == 'H':
			wells_by_edge[e] = w.get_head_wells()
		elif e.connection[w] == 'T':
			wells_by_edge[e] = w.get_tail_wells()
		else:
			raise Exception()

	return wells_by_edge

def validate_vertex(N, W):
	"""Determines whether a neighbor n in N shares wells in W.

	Used when initializing traversal to determine from which edges
	to start walk.
	"""

	for n, C in N:
		if C == 'H':
			W_n = n.get_head_wells()
		elif C == 'T':
			W_n = n.get_tail_wells()
		else:
			raise Exception()

		if len(W_n & W) >= 4:
			return True

	return False

# -----------------------------------------------------------------------------
# neighborhood calculations

def compute_neighborhoods(g):
	"""Computes extended neighborhood of each vertex.

	Determines the set of vertices W_e that overlap with each
	vertex v and that can be reached through edge e.

	Note that we shouldn't need this if we don't throw away
	transitive edges in SGA.
	"""

	for v in g.vertices:
		v.metadata['neighborhood'] = dict()
		for e in v.edges:
			W_e = get_neighborhood(g, v, e)
			v.metadata['neighborhood'][e] = W_e

def get_neighborhood(g, v, e):
	"""Returns all nodes that overlap with v and have a path through e."""
	
	# get other vertex
	w = e.other_vertex(v)

	# get part of w that overlaps with v:
	ovl = e.ovl_start[w], e.ovl_end[w]

	# get nodes that overlap with w using BFS
	visited = set([v, w])
	neighborhood = set([(w, e.connection[w])])
	N_w, visited = visit_neighborhood(g, w, visited, ovl)

	neighborhood.update(N_w)
	return neighborhood

def visit_neighborhood(g, v, visited, ovl):
	"""Returns nodes that overlap with v at ovl.

	Performs a BFS starting at v. Output set N contains tuples
	(n, C), where C is either 'H' or 'T' and indicates if the overlap
	with v is in the head or tail of n.
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
		N.add((n, e.connection[n]))

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