import sys

from bisect import bisect_left, bisect_right
from graph import Vertex, Edge, Graph

from libkuleshov.misc import reverse_string
from libkuleshov.dna import reverse_complement
from libkuleshov.debug import keyboard

class OverlapVertex(Vertex):
	def __init__(self, id_, seq):
		super(OverlapVertex, self).__init__(id_)

		self.seq = seq

		self.head_edges = set() # edges that connect to the head of the segment
		self.tail_edges = set() # edges that connect to the tail of the segment

	def __len__(self):
		return len(self.seq)

	@property
	def edges(self):
		return self.head_edges | self.tail_edges

	def edge_to_vertex(self, v, direction='none'):
		if direction == 'head':
			E = self.head_edges
		elif direction == 'tail':
			E = self.tail_edges
		else:
			E = self.edges

		for e in E:
			if e.other_vertex(self) == v:
				return e

		raise Exception("Error: Vertex %d not found from vertex %d" % (v.id_, self.id_))

	@property
	def prefix_neighbors(self):
		N = set()
		for e in self.head_edges:
			N.add(e.v1)
			N.add(e.v2)

		if N:
			N.remove(self)

		return N

	@property
	def suffix_neighbors(self):
		N = set()
		for e in self.tail_edges:
			N.add(e.v1)
			N.add(e.v2)

		if N:
			N.remove(self)

		return N

	@property
	def neighbors(self):
		return self.prefix_neighbors | self.suffix_neighbors

	def disconnect_edge(self, e):
		print >> sys.stderr, 'deleting', self.id_, e.id_
		self.head_edges.discard(e)
		self.tail_edges.discard(e)
	
class OverlapEdge(Edge):
	def __init__(self, id_, v1, v2, 
				 v1_ovl_start, v1_ovl_end, v1_len, 
				 v2_ovl_start, v2_ovl_end, v2_len, 
				 v2_orientation):

		# figure out the type of connection (L->R), (R->L)
		if v1_ovl_start != 0:
			connection = 'LR'
		elif v1_ovl_start == 0:
			connection = 'RL'
		else:
			exit()

		if connection == 'LR':
			super(OverlapEdge, self).__init__(id_, v1, v2)
		elif connection == 'RL':
			super(OverlapEdge, self).__init__(id_, v2, v1)

		self.ovl_start = {v1: v1_ovl_start, v2: v2_ovl_start}
		self.ovl_end = {v1: v1_ovl_end, v2: v2_ovl_end}

		assert v1_len == len(v1)
		assert v2_len == len(v2)
		length = {v1: v1_len, v2: v2_len}
		
		self.v2_orientation = v2_orientation

		self.connection = dict()
		for v in (v1, v2):
			if self.ovl_start[v] == 0:
				self.connection[v] = 'H'
			elif self.ovl_end[v] == length[v]-1:
				self.connection[v] = 'T'
			else:
				exit("ERROR: Reads don't overlap at endpoints")

	def shift_overlap(self, v, ovl_shift):
		self.ovl_start[v] += ovl_shift
		self.ovl_end[v] += ovl_shift

	def flip_connection(self, v):
		v_len = len(v)

		if self.ovl_start[v] == 0:
			self.ovl_start[v] = v_len - self.ovl_end[v] - 1
			self.ovl_end[v] = v_len-1
		elif self.ovl_end[v] == v_len-1:
			self.ovl_end[v] = v_len - self.ovl_start[v] - 1
			self.ovl_start[v] = 0
		else:
			print self.ovl_start, self.ovl_end
			print v_len
			exit("ERROR: Invalid overlap at edge")

		if self.connection[v] == 'H':
			self.connection[v] = 'T'
		elif self.connection[v] == 'T':
			self.connection[v] = 'H'
		else:
			exit("ERROR")

	def replace(self, v, w):
		self.ovl_start[w] = self.ovl_start[v]
		self.ovl_end[w] = self.ovl_end[v]
		self.connection[w] = self.connection[v]
		
		if self.v1 == v:
			self.v1 = w
		elif self.v2 == v:
			self.v2 = w

		self.ovl_start.pop(v)
		self.ovl_end.pop(v)
		self.connection.pop(v)

	def flip(self):
		""" Flips edge from H->T to T->H and vice versa. """

		self.v1, self.v2 = self.v2, self.v1

class Sequence(object):
	def __init__(self, seq):
		self.seq = seq

	def __len__(self):
		return len(self.seq)

	def __eq__(self, s):
		if self.seq == s.seq or self.seq == reverse_complement(s.seq):
			return True
		else:
			return False

	def __ne__(self, s):
		return (not (self == s))

	def __add__(self, x):
		return Sequence(self.seq + x.seq)

	def __getitem__(self, x):
		return Sequence(self.seq[x])

	def __str__(self):
		return self.seq.__str__()

	def reverse(self):
		self.seq = reverse_string(self.seq)

	def get_reverse(self):
		return Sequence(reverse_string(self.seq))

###############################################################################
## GRAPH TESTING METHODS

def no_diedge(v):
	N = set()
	for e in v.head_edges:
		n1, n2 = e.v1, e.v2

		if n1 != v and n1 in N:
			exit("ERROR: Diedge found")
		if n2 != v and n2 in N:
			exit("ERROR: Diedge found")

		N.add(n1)
		N.add(n2)

	return True

###############################################################################
## BREAKING LONGS CONTIGS INTO PIECES

def test_break_vertex():
	g = Graph()
	v0 = OverlapVertex(g.vertex_id_generator.get_id(), 'XABCD')
	big_v = OverlapVertex(g.vertex_id_generator.get_id(), 'ABCDEFGHIJ')
	big_v.metadata['contigs'] = list()
	big_v.metadata['contig_starts'] = dict()
	big_v.metadata['contig_ends'] = dict()

	big_v.metadata['contigs'].append('c1')
	big_v.metadata['contig_starts']['c1'] = 0
	big_v.metadata['contig_ends']['c1'] = 3

	big_v.metadata['contigs'].append('c2')
	big_v.metadata['contig_starts']['c2'] = 2
	big_v.metadata['contig_ends']['c2'] = 5

	big_v.metadata['contigs'].append('c3')
	big_v.metadata['contig_starts']['c3'] = 4
	big_v.metadata['contig_ends']['c3'] = 7

	big_v.metadata['contigs'].append('c4')
	big_v.metadata['contig_starts']['c4'] = 6
	big_v.metadata['contig_ends']['c4'] = 9

	big_v.metadata['contigs'].append('c5')
	big_v.metadata['contig_starts']['c5'] = -1
	big_v.metadata['contig_ends']['c5'] = -1

	v1 = OverlapVertex(g.vertex_id_generator.get_id(), 'GHIJY')

	e0 = OverlapEdge(g.edge_id_generator.get_id(), v0, big_v,
					 1, 4, 5,
					 0, 3, 10,
					 0)

	e1 = OverlapEdge(g.edge_id_generator.get_id(), big_v, v1,
					 6, 9, 10,
					 0, 3, 5,
					 0)

	v0.tail_edges.add(e0)
	big_v.head_edges.add(e0)
	big_v.tail_edges.add(e1)
	v1.head_edges.add(e1)

	g.add_vertex(v0)
	g.add_vertex(big_v)
	g.add_vertex(v1)
	g.add_edge(e0)
	g.add_edge(e1)

	break_contigs(g)

	print 'VERTICES:'
	for v in g.vertices:
		print v.id_, v.seq
		if v.id_ not in (0, 2	):
			for ctg in v.metadata['contigs']:
				print ctg, 
				if ctg != 'c5':
					print v.metadata['contig_starts'][ctg], v.metadata['contig_ends'][ctg],
				print
	print 

	print 'EDGES:'
	for e in g.edges:
		print e.id_, e.v1.id_, e.ovl_start[e.v1], e.ovl_end[e.v1], e.v2.id_, e.ovl_start[e.v2], e.ovl_end[e.v2]

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def break_vertex(v, g):
	current_v = v
	current_len = len(v)
	# take pieces off from the front:
	while current_len > 40000:
		current_contigs = current_v.metadata['contigs']
		later_contig_starts = sorted([v for v in current_v.metadata['contig_starts'].values() 
									  if v != -1 and v > 2500])
		later_contig_ends = sorted([v for v in current_v.metadata['contig_ends'].values() 
							  if v != -1 and v > 2500])

		# print contig_starts
		# print contig_ends

		# determine positions at which to perform the breaks
		# positions will be considered inclusive.
		# if min(later_contig_ends) < 20000:
		# 	print 'WHOAH'
		# 	new_end = find_le(later_contig_ends, 20000)		# bactrack until first contig end
		# else:
		# 	new_end = 20000

		# if max(later_contig_starts) > 10000:
		# 	exit('FAIL')
		# 	new_start = find_ge(later_contig_starts, 10000)	# go forward until first contig start
		# else:
		# 	new_start = 10000

		new_end = 20000
		new_start = 10000

		# create new vertices:
		small_v = OverlapVertex(g.vertex_id_generator.get_id(), current_v.seq[:new_end+1])
		small_v.metadata['contigs'] = list()
		small_v.metadata['contig_starts'] = dict()
		small_v.metadata['contig_ends'] = dict()

		big_v = OverlapVertex(g.vertex_id_generator.get_id(), current_v.seq[new_start:])
		big_v.metadata['contigs'] = list()
		big_v.metadata['contig_starts'] = dict()
		big_v.metadata['contig_ends'] = dict()

		# populate inner contigs
		for ctg in current_contigs:
			if ctg in current_v.metadata['contig_starts'] \
			and current_v.metadata['contig_starts'][ctg] != -1:
				start = current_v.metadata['contig_starts'][ctg]
				end = current_v.metadata['contig_ends'][ctg]

				if end <= new_end:
					small_v.metadata['contigs'].append(ctg)
					small_v.metadata['contig_starts'][ctg] = start
					small_v.metadata['contig_ends'][ctg] = end

				if new_start <= start:
					big_v.metadata['contigs'].append(ctg)
					big_v.metadata['contig_starts'][ctg] = (start - new_start)
					big_v.metadata['contig_ends'][ctg] = (end - new_start)
				
			else:
				small_v.metadata['contigs'].append(ctg)
				big_v.metadata['contigs'].append(ctg)

		# create new edge
		new_edge = OverlapEdge(g.edge_id_generator.get_id(), small_v, big_v,
							   new_start, len(small_v)-1, len(small_v),
							   0, new_end - new_start, len(big_v),
							   0)

		# attach new edge
		small_v.tail_edges.add(new_edge)
		big_v.head_edges.add(new_edge)

		# reconnect earlier edges:
		for e in current_v.head_edges:
			e.replace(current_v, small_v)
			small_v.head_edges.add(e)

		for e in current_v.tail_edges:
			e.replace(current_v, big_v)
			e.shift_overlap(big_v, -new_start)
			big_v.tail_edges.add(e)

		# remove old vertex:
		g.remove_vertex(current_v)

		# insert vertices and edge into graph
		g.add_vertex(big_v)
		g.add_vertex(small_v)
		g.add_edge(new_edge)

		# move pointer to the remained of the piece
		current_v = big_v
		current_len = len(big_v)

def break_contigs(g):
	V = list(g.vertices)
	
	for v in V:
		if len(v) > 6:
			break_vertex(v, g)