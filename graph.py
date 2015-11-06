import sys
import intervals
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
		print >> sys.stderr, 'xx', v.id_
		for e in v.head_edges:
			self.remove(e)
		for e in v.tail_edges:
			self.remove(e)

		self.vertices_by_id.pop(v.id_)

	def remove_vertex_from_index(self, v):
		self.vertices_by_id.pop(v.id_)

	def add_edge(self, e):
		self.edges_by_id[e.id_] = e

	def remove_edge(self, e):
		v1, v2 = e.v1, e.v2
		
		self.edges_by_id.pop(e.id_)
		v1.disconnect_edge(e)
		v2.disconnect_edge(e)

		# if v in v1.neighbors:
		# 	e2 = '...', v1.edge_to_vertex(v)
		# 	print v1.id_, e2.id_, e.id_
		# if v in v2.neighbors:
		# 	e2 = v2.edge_to_vertex(v)
		# 	print '...', v2.id_, e2.id_, e.id_
		# assert v not in v1.neighbors
		# assert v not in v2.neighbors

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

	## methods for dealing with wells

	def get_wells(self):
		return {Vertex.get_well(ctg) for ctg in self.metadata['contigs']}

	def get_head_wells(self, d=500):
		return {Vertex.get_well(ctg) for ctg in self.metadata['contig_starts']
				if self.metadata['contig_starts'][ctg] < d}

	def get_tail_wells(self, d=500):
		len_v = len(self)
		return {Vertex.get_well(ctg) for ctg in self.metadata['contig_ends']
				if self.metadata['contig_ends'][ctg] > len_v - d}

	def get_head_history(self):
		return {ctg: dict(well=Vertex.get_well(ctg), pos=self.metadata['contig_starts'][ctg])
				for ctg in self.metadata['contig_starts']
				if self.metadata['contig_starts'][ctg] < 4000}

	def get_tail_history(self):
		len_v = len(self)
		return {ctg: dict(well=Vertex.get_well(ctg), pos=len_v-self.metadata['contig_ends'][ctg]-1)
				for ctg in self.metadata['contig_ends']
				if self.metadata['contig_ends'][ctg] > len_v - 4000}

	@staticmethod
	def get_well(ctg):
		assert ctg.startswith('well')
		fields = ctg.split('_')
		return fields[0][4:]

	## methods for dealing with intervals

	def get_true_intervals(self):
		I = list()
		for ctg in self.metadata['contigs']:
			I.append(Vertex.parse_interval(ctg))

		if I:
			return intervals.merge_intervals(I)

	def get_head_intervals(self):
		I = list()
		head_ctgs = [ctg for ctg in self.metadata['contigs'] 
					 if w.metadata['contig_starts'][ctg] != -1 and \
						w.metadata['contig_starts'][ctg] < 200 ]
		for ctg in head_ctgs:
			I.append(Vertex.parse_interval(w, ctg))

		if I:
			return intervals.merge_intervals(I)
		else:
			return list()

	def get_tail_intervals(w, ctgs):
		I = list()
		w_len = len(w.seq)
		tail_ctgs = [ctg for ctg in self.metadata['contigs']
					 if w.metadata['contig_ends'][ctg] != -1 and \
			  			w.metadata['contig_ends'][ctg] > w_len - 200 ]
		for ctg in tail_ctgs:
			I.append((int(chrom), start + internal_start, start + internal_start + 199))

		if I:
			return intervals.merge_intervals(I)
		else:
			return list()

	def print_true_intervals(self):
		I = self.get_true_intervals()
		merged_I = intervals.merge_intervals(I)
		for i in merged_I:
			print '\t', i

	@staticmethod
	def parse_interval(ctg):
		fields = ctg.split('_')
		chrom, coords = fields[1].split(':')
		start, end = (int(x) for x in coords.split('-'))

		# correction for how we generated the reads. bed format is 1-based
		start -= 1
		end -= 1

		internal_start = int(fields[2])
		return int(chrom), start + internal_start, start + internal_start + 199

class Edge(object):
	def __init__(self, id_, v1, v2):
		self.id_ = id_
		self._v1 = v1
		self._v2 = v2
		self.metadata = dict()
		
	def __eq__(self, e):
		return self.id_ == e.id_

	def other_vertex(self, v):
		if v == self.v1:
			return self.v2
		elif v == self.v2:
			return self.v1
		else:
			raise Exception("ERROR: Vertex %d not found in v.other_vertex" % v.id_)

	def __hash__(self):
		return hash(self.id_)

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

