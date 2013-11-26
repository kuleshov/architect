# testneighborhood.py
import sys
import unittest
sys.path.append('..')

from pmpp import get_neighborhood
from graph_load import load_from_sga_asqg, save_graph, load_graph

class TestNeighborhood(unittest.TestCase):
	def setUp(self):
		self.g = load_graph('./data/2.pruned.asqg', 
	        	   	        './data/2.pruned.containment')

	def tearDown(self):
		pass

	def testNb(self):
		v = self.g.vertices_by_id[3764187]
		e = self.g.edges_by_id[47]
		N = get_neighborhood(self.g, v, e)
		for n in N:
			print '>>>', n.id_

if __name__ == '__main__':
	unittest.main()