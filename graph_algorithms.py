import itertools
import networkx as nx

from scaffolder import make_wellscaff_edges #, contract_edges
from contraction import contract_edges

# ----------------------------------------------------------------------------

def scaffold_via_wells_mst(g):
  # delete all existing edges from the graph
  E = g.edges
  for e in E:
    g.remove_edge(e)

  # create new edges whenever vertices have similar well profiles
  n_edges = make_wellscaff_edges(g)
  print '%d scaffold edges created via wells.' % n_edges
  
  # construct well-based scaffold graph in networkx format
  nxg = g.nxgraph
  # nxg = _construct_graph(g)

  # find the maxinum spanning forest
  msf = nx.minimum_spanning_tree(nxg)

  # remove edges not in forest
  E = [e for e in g.edges]
  n_removed = 0
  for e in E:
    e_nx = ( (e.v1.id, e.connection[e.v1]), (e.v2.id, e.connection[e.v2]) )
    if not msf.has_edge(*e_nx):
      g.remove_edge(e)
      n_removed += 1

  print '%d edges not in MST removed.' % n_removed

  for v in g.vertices:
    v.initialize_contigs()
  
  # contract edges
  n_contracted = contract_edges(g, store_layout=True)
  print '%d edges contracted.' % n_contracted

  # recompute the maxinum spanning forest
  msf = nx.minimum_spanning_tree(g.nxgraph)

  # for each tree in forest:
  n_removed = 0
  trunk = list()
  for mst in nx.connected_component_subgraphs(msf):
    # add to mst trunk
    trunk.extend(_mst_trunk(mst))

  # remove edges not in trunk:
  E = [e for e in g.edges]
  trunk = set(trunk)
  n_removed = 0
  for e in E:
    e_nx = ( (e.v1.id, e.connection[e.v1]), (e.v2.id, e.connection[e.v2]) )
    if not e_nx in trunk:
      g.remove_edge(e)
      n_removed += 1

  print '%d edges not in trunk removed.' % n_removed

  # contract one last time
  n_contracted = contract_edges(g, store_layout=True)
  print '%d edges contracted.' % n_contracted


# ----------------------------------------------------------------------------
# helpers

def _mst_trunk(mst):
  # compute shortest path distances between nodes
  all_pairs_shortest_dists = nx.shortest_path_length(mst)

  # determine the pair that corresponds to the longest distance
  all_pairs = itertools.product(mst.nodes_iter(), mst.nodes_iter())
  s, t = max(all_pairs, key=lambda p: all_pairs_shortest_dists[p[0]][p[1]])

  # return the path between s and t
  return nx.shortest_path(mst,s,t)

def _construct_graph(g):
  nxg = g.nxgraph.copy()

  # clear all edges except the internal ones
  nxg.remove_edges_from(nxg.edges)
  for v in g.vertices:
    nxg.add_edge( (v.id,'H'), (v.id,'T') )

  edge_generator = (
    ( (v1, conn1), (v2, conn2), -_edge_weight(v1, v2, conn1, conn2) )
    for v1 in g.vertices for v2 in g.vertices 
    for conn1 in ('H', 'T') for conn2 in ('H', 'T')
    if (v1.id > v2.id) and len(v1) > 5000 and len(v2) < 5000
  )

  nxg.add_weighted_edges_from(edge_generator)

  return nxg