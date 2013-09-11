def pop_triangles(g):
	# pop triangles
	for v in g.vertices:
		for x in v.prefix_neighbors:
			for y in v.prefix_neighbors:
				if x == y: continue
				e = get_neighboring_edge(x,y)
				if e:
					g.remove_edge(e)

		for x in v.suffix_neighbors:
			for y in v.prefix_neighbors:
				if x == y: continue
				e = get_neighboring_edge(x,y)
				if e:
					g.remove_edge(e)


def get_neighboring_edge(v,w):
	for e in v.head_edges:
		if w == e.v1 or w == e.v2:
			return e
	for e in v.tail_edges:
		if w == e.v1 or w == e.v2:
			return e

	return None

def detect_cycle(g, start_v):
	stack = list()
	stack.append(start_v)
	visited = set()
	back_pointers = dict()
	prev_v = None

	print "cycle detection..."

	while (stack):
		v = stack.pop()
		print '~', v.id_, start_v in stack

		if v not in visited:
			back_pointers[v] = prev_v
		
		N = v.get_neighbors()
		N.discard(prev_v)
		
		if start_v in N: recover_cycle(v, back_pointers)

		prev_v = v
		stack.extend(N)

	print len(visited)

def recover_cycle(v, back_pointers):
	while v:
		print '>', v.id_
		v = back_pointers[v]