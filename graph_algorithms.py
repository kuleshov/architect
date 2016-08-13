"""Graph manipulation algorithms

Currently, this module contains an alternative way to prune the scaffold graph
using Maximum Spanning Trees, in a way that is similar to the FragScaff 
algorithm.

Our main approach seems to currently work better, so we currently don't use 
this code.
"""

import itertools
import networkx as nx

from scaffolder import make_wellscaff_edges #, contract_edges
from contraction import contract_edges

# ----------------------------------------------------------------------------

def full_scaffold_via_wells_mst(g):
  # delete all existing edges from the graph
  E = g.edges
  for e in E:
    g.remove_edge(e)

  # create new edges whenever vertices have similar well profiles
  n_edges = make_wellscaff_edges(g, min_common=4, min_thr=0.1)
  print '%d scaffold edges created via wells.' % n_edges

  scaffold_via_wells_mst(g)

def scaffold_via_wells_mst(g):
  # initialize internal contig labels (used for downstream qc)
  for v in g.vertices:
      v.initialize_contigs()
  
  # construct well-based scaffold graph in networkx format
  nxg = g.nxgraph
  # nxg = _construct_graph(g)

  # weigh edges according to how many wells they are sharing:
  _reweigh_edges(nxg, g, type_='wells')

  # find the maxinum spanning forest
  msf = nx.minimum_spanning_tree(nxg)

  # keep simplifying the graph until the msf has no branching nodes:
  n_iter = 1
  while _has_branches(msf) and n_iter <= 10:
    print 'MSF simplificaiton iteration %d' % n_iter

    # print '...', max(msf.degree(weight=None).values())
    # print '...', sorted(msf.degree(weight=None).iteritems(), key=lambda x: x[1], reverse=True)[:10]
    # vg = sorted(msf.degree(weight=None).iteritems(), key=lambda x: x[1], reverse=True)[0][0]
    # v = g.vertex_from_id(vg[0])
    # N = [n.id for n in g.vertices if v in n.neighbors]
    # print ',,,', N
    # print msf.neighbors(v)


    # remove edges of g not selected in forest MSF
    E = [e for e in g.edges]
    n_removed = 0
    for e in E:
      e_nx = ( (e.v1.id, e.connection[e.v1]), (e.v2.id, e.connection[e.v2]) )
      if not msf.has_edge(*e_nx):
        g.remove_edge(e)
        n_removed += 1

    print '%d edges not in MST removed.' % n_removed
    
    # contract edges
    n_contracted = contract_edges(g, store_layout=True)
    print '%d edges contracted.' % n_contracted

    # now we are going to compute the trunk

    # get the networkx graph again
    nxg = g.nxgraph
    _reweigh_edges(nxg, g, type_='wells') # FIXME: do this once

    # recompute the maxinum spanning forest
    msf = nx.minimum_spanning_tree(g.nxgraph)

    # for each tree in forest:
    trunk = list()
    for mst in nx.connected_component_subgraphs(msf):
      # add to mst trunk
      if len(mst) >= 4:
        trunk.extend(_mst_trunk(mst, g))

    # remove edges not in trunk:
    E = [e for e in g.edges]
    print trunk
    trunk_v = set([v[0] for v in trunk])
    n_removed = 0
    for e in E:
      v1_id, v2_id = e.v1.id, e.v2.id
      if v1_id not in trunk_v or v2_id not in trunk_v:
        g.remove_edge(e)
        n_removed += 1

    if n_iter >= 4: keyboard()

    print '%d edges not in trunk removed.' % n_removed

    # contract one last time
    n_contracted = contract_edges(g, store_layout=True)
    print '%d edges contracted.' % n_contracted

    # construct well-based scaffold graph in networkx format
    nxg = g.nxgraph
    # nxg = _construct_graph(g)

    # weigh edges according to how many wells they are sharing:
    _reweigh_edges(nxg, g, type_='wells')

    # find the maxinum spanning forest
    msf = nx.minimum_spanning_tree(nxg)

    n_iter += 1


# ----------------------------------------------------------------------------
# helpers

def _mst_trunk(mst, g):
  # weigh edges according to their distance
  _reweigh_edges(mst, g, type_='lengths')

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

def _has_branches(nxg):
  if max(nxg.degree(weight=None).values()) > 2:
    return True
  else:
    return False

def _reweigh_edges(nxg, g, type_):
  for u, v, d in nxg.edges(data=True):
    if u[0] == v[0]:
      if type_ == 'wells':
        d['weight'] = -_transfer_fn(1.)-1000
      elif type_ == 'lengths':
        d['weight'] = len(g.vertex_from_id(u[0]))
    else:
      if type_ == 'wells':
        v1 = g.vertex_from_id(u[0])
        v2 = g.vertex_from_id(v[0])
        d['weight'] = _edge_weight(v1, v2, u[1], v[1])
      elif type_ == 'lengths':
        d['weight'] = 0.

def _edge_weight(v1, v2, conn1, conn2):
  if conn1 == 'H':
    v1_wells = v1.head_wells
  elif conn1 == 'T':
    v1_wells = v1.tail_wells
  if conn2 == 'H':
    v2_wells = v2.head_wells
  elif conn2 == 'T':
    v2_wells = v2.tail_wells

  common_wells = v1_wells & v2_wells

  if not v1_wells or not v2_wells: 
    return 0.

  all_wells = v1_wells | v2_wells
  frac_common = float(len(common_wells)) / float(len(all_wells))

  return _transfer_fn(frac_common)

def _transfer_fn(fraction):
  return 1./(1.01-fraction) - 1.