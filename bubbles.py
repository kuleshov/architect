from collections import deque
from contraction import contract_edges

def simplify_graph(g):
	for i in xrange(5):
		print i
		edges_to_contract = pop_triangles(g)
		contract_edges(g, edges_to_contract)

def check_bubbles(g, forward_pointers):
	v0 = g.vertices_by_id[3764304]
	e0 = g.edges_by_id[188]

	start_v, end_v = detect_bubble_from_pointers(g, v0, e0, forward_pointers, 500)

	print start_v.id_
	print end_v.id_

# -----------------------------------------------------------------------------
# Pop "transitive" bubbles

def get_triangle(v):
	# bubbles = list()
	for n1 in v.suffix_neighbors:
		for n2 in v.suffix_neighbors:
			if n1 == n2: continue
			if len(n2.edges) == 2:
				if v in n2.prefix_neighbors and n1 in n2.suffix_neighbors \
				or v in n2.suffix_neighbors and n1 in n2.prefix_neighbors:
					# bubbles.append((v, n1, n2))
					return v, n1, n2

	for n1 in v.prefix_neighbors:
		for n2 in v.prefix_neighbors:
			if n1 == n2: continue
			if len(n2.edges) == 2:
				if v in n2.suffix_neighbors and n1 in n2.prefix_neighbors \
				or v in n2.prefix_neighbors and n1 in n2.suffix_neighbors:
					# bubbles.append((v, n1, n2))
					return v, n1, n2

	return None

def pop_triangles(g):
	edges_to_contract = set()
	for v in g.vertices:
		T = get_triangle(v)
		if T:
			x, y, t = T
		
			# always prefer the longest path
			straight_length = path_length((x, y))
			transitive_length = path_length((x, t, y))
			assert transitive_length > straight_length

			# only try to resolve short tandem repeats:
			if transitive_length - straight_length < 200:
			
				# delete the straight path edge
				e_straight = x.edge_to_vertex(y)
				e_trans1 = x.edge_to_vertex(t)
				e_trans2= t.edge_to_vertex(y)

				print 'removing edge %d' % e_straight.id_
				g.remove_edge(e_straight)
				edges_to_contract.add(e_trans1)
				edges_to_contract.add(e_trans2)


# -----------------------------------------------------------------------------
# BFS stuff

def detect_bubble_from_pointers(g, v0, e0, forward_pointers, max_distance):
	# do a BFS from v0 until you reach max. distance
	# and collect all nodes at the boundary
	V, B = get_border_vertices(g, v0, e0, forward_pointers, max_distance)

	print 'V'
	for v in V:
		print v.id_
	print 'B'
	for b in B:
		print b.id_

	# take on vertex
	t = B.pop()

	print 't'
	print t.id_

	# get edges on path to t
	edges_on_path, vertices_on_path = path_via_bfs(g, v0, e0, forward_pointers, t)
	print 'P'
	for p in vertices_on_path:
		print p.id_ 

	# get all nodes via BFS from v0 without tacking edges on this path
	skipped_vertices = set(vertices_on_path)
	# B REDEFINED???
	_, B_new = get_border_vertices(g, v0, e0, forward_pointers, max_distance, 
							   skipped_vertices)

	print 'B_new'
	for b in B_new:
		print b.id_

	# if one of these nodes touches boundary, we cannot resolve this bubble
	resolvable = True
	if B_new & B:
		resolvable = False

	print 'extreme vertex calculation'

	# otherwise, the extreme vertex is the last on the path
	# that intersects with N
	extreme_vertex = None
	if resolvable:
		for p in reversed(vertices_on_path):
			print p.id_
			if p in B:
				extreme_vertex = p

	print '---'

	return v0, extreme_vertex

def path_via_bfs(g, v0, e0, forward_pointers, t):
	visited = set([v0])
	to_expand = deque([(v0, e0)])
	edge_back_pointers = dict()
	vertex_back_pointers = dict()

	while to_expand:
		v, e = to_expand.popleft()
		W = [(f.other_vertex(v), f) for f in forward_pointers[e]]

		for w, f in W:
			if w in visited:
				continue
			elif w == t:
				edge_back_pointers[f] = e
				vertex_back_pointers[w] = v

				e_curr = f
				edges_on_path = list()
				while e_curr in edge_back_pointers:
					edges_on_path.append(e_curr)
					e_curr = edge_back_pointers[e_curr]

				v_curr = w
				vertices_on_path = list()
				while v_curr in vertex_back_pointers:
					vertices_on_path.append(v_curr)
					v_curr = vertex_back_pointers[v_curr]

				return edges_on_path, vertices_on_path
			else:
				edge_back_pointers[f] = e
				vertex_back_pointers[w] = v
				to_expand.append((w, f))

	return None, None

def get_border_vertices(g, v0, e0, forward_pointers, max_distance, 
						skipped_vertices=set()):
	# input is like this:
	# ----> v0 ---->
	# (e0)
	visited = set([v0])
	# we remember vertex to expand, ~shortest distance from v0,
	# and the edge through which we're coming in
	to_expand = deque([(v0, 0, e0)])
	border_vertices = set()

	while to_expand:
		v, dist, e = to_expand.popleft()
		W = [(f.other_vertex(v), f) for f in forward_pointers[e]]
		
		if not W:
			border_vertices.add(v)
		else:
			WW = [(w, distance(v, w, f), f) for w, f in W 
				  if w not in visited and w not in skipped_vertices]

		if not WW: border_vertices.add(v)

		for w, d, f in sorted(WW, key=lambda x: x[1]):
			if w in visited: 
				continue
			if w in skipped_vertices: 
				continue
			elif d > max_distance:
				border_vertices.add(w)
				continue
			else:
				visited.add(w)
				to_expand.append((w, d, f))

	return visited, border_vertices

def distance(v1, v2, e):
	# we always compute distances from the edge of v0
	# to the start of the read v1
	# so this function gives you the distances between
	# the edges of v1, v2 in the given direction
	if e.connection[v2] == 'H':
		return len(v2) - e.ovl_end[v2]
	elif e.connection[v1] == 'T':
		return len(v2) - e.ovl_start[v2]


# ----------------------------------------------------------------------------
# helpers

def path_length(v_list):
	v0 = v_list[0]
	L = len(v0)
	for v in v_list[1:]:
		assert v in v0.neighbors
		e = v0.edge_to_vertex(v)
		ovl_len = e.ovl_end[v] - e.ovl_start[v] + 1
		assert ovl_len > 0
		L += (len(v) - ovl_len)

	return L
