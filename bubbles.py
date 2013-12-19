from contraction import contract_edges

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

def simplify_graph(g):
	for i in xrange(5):
		print i
		edges_to_contract = pop_triangles(g)
		contract_edges(g, edges_to_contract)

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
