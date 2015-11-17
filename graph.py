import sys
import intervals
import networkx as nx
from libkuleshov.debug import keyboard

# ----------------------------------------------------------------------------
# graph-related classes. internally, stores dict of id->vertex. 
# adjacency info is in vertex and edge objects.

class Graph(object):
	def __init__(self):
		self.vertices_by_id = dict()
		self.edges_by_id = dict()
		self.metadata = dict()

		self.vertex_id_generator = IdGenerator()
		self.edge_id_generator = IdGenerator()
		self._graph = nx.Graph()

	@property
	def vertices(self):
		return self.vertices_by_id.values()

	@property
	def edges(self):
		return self.edges_by_id.values()

	def vertex_from_id(self, vid):
		return self.vertices_by_id[vid]

	def get_edge(self, id_):
		return self.edges_by_id[id_]

	def add_vertex(self, v):
		self.vertices_by_id[v.id] = v
		vg1, vg2 = (v.id, 'H'), (v.id, 'T')
		self._graph.add_nodes_from([vg1, vg2])
		self._graph.add_edge(vg1, vg2)

	def remove_vertex(self, v):
		E_head = [e for e in v.head_edges]
		E_tail = [e for e in v.tail_edges]
		for e in E_head:
			self.remove_edge(e)
		for e in E_tail:
			self.remove_edge(e)

		self.remove_vertex_from_index(v)

	def remove_vertex_from_index(self, v):
		self.vertices_by_id.pop(v.id)
		vg1, vg2 = (v.id, 'H'), (v.id, 'T')
		self._graph.remove_node(vg1)
		self._graph.remove_node(vg2)

	def add_edge(self, e):
		self.edges_by_id[e.id] = e
		v1, v2 = e.v1, e.v2
		assert v1 != v2
		vg1 = (v1.id, e.connection[v1])
		vg2 = (v2.id, e.connection[v2])
		self._graph.add_edge(vg1, vg2)

	def remove_edge(self, e):
		v1, v2 = e.v1, e.v2

		vg1 = (v1.id, e.connection[v1])
		vg2 = (v2.id, e.connection[v2])
		self._graph.remove_edge(vg1, vg2)
		
		self.edges_by_id.pop(e.id)
		v1.disconnect_edge(e)
		v2.disconnect_edge(e)

	def reconnect_edge(self, e, v, w):
		pass

	def count_connected_components(self):
		return nx.number_connected_components(self._graph)
	
	@property
	def nxgraph(self):
	  return self._graph

# ----------------------------------------------------------------------------

class Vertex(object):
	def __init__(self, id_):
		self._my_id=id_
		self._metadata = {'wells': dict(), 'intervals': list(), 
											'contigs': list()}
		
	def __eq__(self,v):
		return self._my_id == v._my_id

	def __hash__(self):
		return hash(self.id)

	def __eq__(self, v):
		return self.id == v.id

	def __ne__(self, v):
		return self.id != v.id

	@property
	def id(self):
		return self._my_id

	#TODO: remove this after refactoring code that used metadata directly
	@property
	def metadata(self):
	  return self._metadata
	
	## methods for dealing with wells

	def add_well(self, w, start, end):
		ivl = (0, start, end)
		if w not in self._metadata['wells']:
			self._metadata['wells'][w] = ivl
		else:
			ivl0 = self._metadata['wells'][w]
			self._metadata['wells'][w] = intervals.union(ivl0, ivl)

	def shift_well(self, w, shift):
		self._metadata['wells'][w][1] += shift
		self._metadata['wells'][w][2] += shift

	def well_interval(self, w):
		if w in self._metadata['wells']:
			return self._metadata['wells'][w][1], self._metadata['wells'][w][2]
		else:
			return None

	def set_well_interval(self, w, s, e):
		assert w in self._metadata['wells']
		self._metadata['wells'][w] = (0, s, e)

	@property
	def wells(self):
	  return self._metadata['wells'].keys()

	@property
	def head_wells(self, d=500):
		return {w for w, i in self.metadata['wells'].iteritems() if i[1] < d}

	@property
	def tail_wells(self, d=500):
		len_v = len(self)
		return { w for w, i in self.metadata['wells'].iteritems()
						 if i[2] > len_v - d }

	@staticmethod
	def get_well(ctg):
		assert ctg.startswith('well')
		fields = ctg.split('_')
		return fields[0][4:]

	## methods for dealing with intervals

	def add_interval(self, ivl):
		self._metadata['intervals'].append(ivl)
		self._merge_intervals()

	def shift_intervals(self, shift):
		ivls = [(c, s+shift, e+shift) for (c,s,e) in self._metadata['intervals']]
		self._metadata['intervals'] = ivls

	@property
	def intervals(self):
		return self._metadata['intervals']

	def print_true_intervals(self):
		I = self.get_true_intervals()
		merged_I = intervals.merge_intervals(I)
		for i in merged_I:
			print '\t', i

	def _merge_intervals(self):
		I = self._metadata['intervals']
		self._metadata['intervals'] = intervals.merge_intervals(I)

	## methods for dealing with internal contigs

	@property
	def contigs(self):
	  return self._metadata['contigs']

	def initialize_contigs(self):
		self._metadata['contigs'] = [(self.id, self.intervals, len(self), 'S')]

	def set_contigs_from_vertices(self, v1, v2):
		self._metadata['contigs'] = v1.contigs + v2.contigs

	def flip_contigs(self):
		def _flip_ori(ctg):
			new_ori = 'R' if ctg[3] == 'S' else 'S'
			return (ctg[0], ctg[1], ctg[2], new_ori)

		self._metadata['contigs'] \
			= list(reversed([_flip_ori(ctg) for ctg in self.contigs]))

# ----------------------------------------------------------------------------

class Edge(object):
	def __init__(self, id_, v1, v2):
		super(Edge, self).__init__()
		self._id = id_
		self._v1 = v1
		self._v2 = v2
		self.metadata = dict()
		
	def __eq__(self, e):
		return self._id == e.id

	def other_vertex(self, v):
		if v == self.v1:
			return self.v2
		elif v == self.v2:
			return self.v1
		else:
			raise Exception("ERROR: Vertex %d not found in v.other_vertex" % v.id)

	def __hash__(self):
		return hash(self._id)

	@property
	def id(self):
	  return self._id

	@property
	def v1(self):
	  return self._v1

	@v1.setter
	def v1(self, w):
	  self._v1 = w

	@property
	def v2(self):
	  return self._v2

	@v2.setter
	def v2(self, w):
	  self._v2 = w

class IdGenerator:
	def __init__(self, start=0):
		self.counter = start

	def get_id(self):
		i = self.counter
		self.counter += 1
		return i

	def set_counter(self, n):
		self.counter = n
