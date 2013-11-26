from libkuleshov.debug import keyboard

###############################################################################
## GRAPH
## Implementation: Store dict of id->vertex. Store adjacency info in vertex
## and edge objects.

class Graph(object):
	def __init__(self):
		self.vertices_by_id = dict()
		self.edges_by_id = dict()
		self.metadata = dict()

		self.vertex_id_generator = IdGenerator()
		self.edge_id_generator = IdGenerator()

	@property
	def vertices(self):
		return self.vertices_by_id.values()

	@property
	def edges(self):
		return self.edges_by_id.values()

	def get_edge(self, id_):
		return self.edges_by_id[id_]

	def add_vertex(self, v):
		self.vertices_by_id[v.id_] = v

	def remove_vertex(self, v):
		for e in v.head_edges:
			g.remove(e)
		for e in v.tail_ddges:
			g.remove(e)

		self.vertices_by_id.pop(v.id_)

	def add_edge(self, e):
		self.edges_by_id[e.id_] = e

	def remove_edge(self, e):
		v1, v2 = e.v1, e.v2
		
		self.edges_by_id.pop(e.id_)
		v1.disconnect_edge(e)
		v2.disconnect_edge(e)

	def remove_vertex(self, v):
		self.vertices_by_id.pop(v.id_)

	def count_connected_components(self):
		to_visit = set(self.vertices)
		num_components = 0

		while(to_visit):
			edge_nodes = set([to_visit.pop()])
			while(edge_nodes):
				v = edge_nodes.pop()
				to_visit.discard(v)
				edge_nodes |= (v.neighbors & to_visit)
			num_components += 1

		return num_components

class Vertex(object):
	def __init__(self, id_):
		self.id_=id_
		self.metadata = dict()
		
	def __eq__(self,v):
		return self.id_ == v.id_

	def __hash__(self):
		return hash(self.id_)

	def __eq__(self, v):
		return self.id_ == v.id_

	def __ne__(self, v):
		return self.id_ != v.id_

	## get_neighbors abstract method

	## disconnect edge_abstract method

	## method for dealing with wells

	def get_wells(self):
		return {Vertex.get_well(ctg) for ctg in self.metadata['contigs']}

	def get_head_wells(self):
		return {Vertex.get_well(ctg) for ctg in self.metadata['contig_starts']
				if self.metadata['contig_starts'][ctg] < 500}

	def get_tail_wells(self):
		len_v = len(self)
		return {Vertex.get_well(ctg) for ctg in self.metadata['contig_ends']
				if self.metadata['contig_ends'][ctg] > len_v - 500}

	@staticmethod
	def get_well(ctg):
		assert ctg.startswith('well')
		fields = ctg.split('_')
		return fields[0][4:]

class Edge(object):
	def __init__(self, id_, v1, v2):
		self.id_ = id_
		self.v1 = v1
		self.v2 = v2
		self.metadata = dict()
		
	def __eq__(self, e):
		return self.id_ == e.id_

	def other_vertex(self, v):
		if v == self.v1:
			return self.v2
		elif v == self.v2:
			return self.v1
		else:
			exit("ERROR: Vertex not found")

	def __hash__(self):
		return hash(self.id_)

class IdGenerator:
	def __init__(self, start=0):
		self.counter = start

	def get_id(self):
		i = self.counter
		self.counter += 1
		return i

	def set_counter(self, n):
		self.counter = n

