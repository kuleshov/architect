from graph import Vertex, Edge

from libkuleshov.misc import reverse_string
from libkuleshov.dna import reverse_complement
from libkuleshov.debug import keyboard

## call it string graph!

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
