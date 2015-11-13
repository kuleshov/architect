import networkx as nx

import intervals
from string_graph import ScaffoldEdge
from contraction import contract_edges

from visualize import print_vertex

# ----------------------------------------------------------------------------
# classical scaffolding

def prune_scaffold_edges(g):
  sorted_v = sorted(g.vertices, key=lambda x: len(x.seq), reverse=True)
  n_pruned = 0
  for v in sorted_v:
    # look at the head side
    if len(v.head_edges) >= 2 and all(e.is_scaffold_edge for e in v.head_edges):
      sorted_e = sorted(v.head_edges, key=lambda x: x.support, reverse=True)
      e1, e2 = sorted_e[0], sorted_e[1]
      abs_support_delta = e1.support - e2.support
      rel_support_delta = float(e2.support) / float(e1.support)
      if (abs_support_delta >= 3) and rel_support_delta < 0.7:
        n_pruned += len(v.head_edges) - 1
        _select_edge(g, v.head_edges, e1)

    # look at the tail side
    if len(v.tail_edges) >= 2 and all(e.is_scaffold_edge for e in v.tail_edges):
      sorted_e = sorted(v.tail_edges, key=lambda x: x.support, reverse=True)
      e1, e2 = sorted_e[0], sorted_e[1]
      abs_support_delta = e1.support - e2.support
      rel_support_delta = float(e2.support) / float(e1.support)
      if (abs_support_delta >= 3) and rel_support_delta < 0.7:
        n_pruned += len(v.tail_edges) - 1
        _select_edge(g, v.tail_edges, e1)

  return n_pruned

def prune_scaffold_edges_via_wells(g):
  sorted_v = sorted(g.vertices, key=lambda x: len(x.seq), reverse=True)
  n_pruned = 0
  for v in sorted_v:
    # look at the head side
    if len(v.head_edges) >= 2 and all(e.is_scaffold_edge for e in v.head_edges):
      supported_edges = [e for e in v.head_edges if len(_get_common_wells(e)) >= 4]
      if len(supported_edges) == 1:
        n_pruned += len(v.head_edges) - 1
        _select_edge(g, v.head_edges, supported_edges[0])  

    # look at the tail side
    if len(v.tail_edges) >= 2 and all(e.is_scaffold_edge for e in v.tail_edges):
      supported_edges = [e for e in v.tail_edges if len(_get_common_wells(e)) >= 4]
      if len(supported_edges) == 1:
        n_pruned += len(v.tail_edges) - 1
        _select_edge(g, v.tail_edges, supported_edges[0])

  return n_pruned

def cut_tips(g, d=500):
  n_cut = 0
  for v in g.vertices:
    if len(v.seq) < d:
      continue

    head_edge_num = len(v.head_edges)
    tail_edge_num = len(v.tail_edges)

    if min(head_edge_num, tail_edge_num) == 0 and max(head_edge_num, tail_edge_num) == 1:
      g.remove_vertex(v)
      n_cut += 1

  return n_cut

# ----------------------------------------------------------------------------
# well-based scaffolding

def scaffold_via_wells(g):
  # delete all existing edges from the graph
  E = g.edges
  for e in E:
    g.remove_edge(e)

  # create new edges whenever vertices have similar well profiles
  n_edges = _make_wellscaff_edges(g)
  print '%d scaffold edges created via wells.' % n_edges

  # _inspect_new_edges(g)

  for v in g.vertices:
    v.initialize_contigs()

  # contract edges
  n_contracted = contract_edges(g, store_layout=True)
  print '%d edges contracted.' % n_contracted

def _make_wellscaff_edges(g):
  n_edges = 0
  for v1 in g.vertices:
    if len(v1) < 5000: continue
    for v2 in g.vertices:
      if len(v2) < 5000: continue
      if v1.id >= v2.id: continue
      # uncomment this to only connect "true edges"
      # if not intervals.overlap(v1.intervals, v2.intervals): continue
      for conn1 in ('H', 'T'):
        for conn2 in ('H', 'T'):
          common_wells, v1_wells, v2_wells = _get_wells_between_v(v1, v2, conn1, conn2)
          if not v1_wells or not v2_wells: continue
          frac_common1 = float(len(common_wells)) / float(len(v1_wells))
          frac_common2 = float(len(common_wells)) / float(len(v2_wells))
          if len(common_wells) >= 4 and max(frac_common1, frac_common2) > 0.33 \
                                    and min(frac_common1, frac_common2) > 0.33 :
            ori = 0 if conn1 != conn2 else 1
            j = g.edge_id_generator.get_id()
            e = ScaffoldEdge(j, v1, v2, conn1, conn2, ori, 3000)
            for v in (v1, v2):
              if e.connection[v] == 'H':
                v.head_edges.add(e)
              elif e.connection[v] == 'T':
                v.tail_edges.add(e)
              else:
                exit('ERROR: Invalid edge connection!')
            g.add_edge(e)
            n_edges += 1

  return n_edges

# ----------------------------------------------------------------------------

def _inspect_new_edges(g):
  for e in g.edges:
    v1, v2 = e.v1, e.v2

    contractable = False
    if e.connection[v1] == 'H' and len(v1.head_edges) == 1:
      if e.connection[v2] == 'H' and len(v2.head_edges) == 1:
        contractable = True
      elif e.connection[v2] == 'T' and len(v2.tail_edges) == 1:
        contractable = True
    elif e.connection[v1] == 'T' and len(v1.tail_edges) == 1:
      if e.connection[v2] == 'H' and len(v2.head_edges) == 1:
        contractable = True
      elif e.connection[v2] == 'T' and len(v2.tail_edges) == 1:
        contractable = True

    print
    print e.id, v1.id, v2.id, e.connection[v1], e.connection[v2], contractable
    print_vertex(v1)
    print_vertex(v2)

def examine_scaffold_ambiguities(g):
  sorted_v = sorted(g.vertices, key=lambda x: len(x.seq), reverse=True)
  n_pruned = 0
  for v in sorted_v:
    # look at the head side
    if len(v.head_edges) >= 2 and all(e.is_scaffold_edge for e in v.head_edges):
      sorted_e = sorted(v.head_edges, key=lambda x: x.support, reverse=True)
      e1, e2 = sorted_e[0], sorted_e[1]
      abs_support_delta = e1.support - e2.support
      rel_support_delta = float(e2.support) / float(e1.support)
      w1, w2 = e1.other_vertex(v), e2.other_vertex(v)

      if e1.connection[v] == 'H':
        v_wells1 = v.head_wells
      elif e1.connection[v] == 'T':
        v_wells1 = v.tail_wells
      if e2.connection[v] == 'H':
        v_wells2 = v.head_wells
      elif e2.connection[v] == 'T':
        v_wells2 = v.tail_wells

      if e1.connection[w1] == 'H':
        w1_wells = w1.head_wells
        w1_connectivity = len(w1.head_edges)
      elif e1.connection[w1] == 'T':
        w1_wells = w1.tail_wells
        w1_connectivity = len(w1.tail_edges)
      if e2.connection[w2] == 'H':
        w2_wells = w2.head_wells
      elif e2.connection[w2] == 'T':
        w2_wells = w2.tail_wells

      print v.id, 'H', len(v), v.intervals
      print 'v wells:', ','.join([str(w) for w in v_wells1])
      print 'e1/w1:', e1.id, w1.id, len(w1), w1.intervals 
      print 'w1 wells:', ','.join([str(w) for w in w1_wells])
      print 'e2/w2:', e2.id, w2.id, len(w2), w2.intervals 
      print 'w2 wells:', ','.join([str(w) for w in w2_wells])
      print e1.support, e2.support, abs_support_delta, rel_support_delta, w1_connectivity
      print

    # look at the tail side
    if len(v.tail_edges) >= 2 and all(e.is_scaffold_edge for e in v.tail_edges):
      sorted_e = sorted(v.tail_edges, key=lambda x: x.support, reverse=True)
      e1, e2 = sorted_e[0], sorted_e[1]
      abs_support_delta = e1.support - e2.support
      rel_support_delta = float(e2.support) / float(e1.support)
      w1, w2 = e1.other_vertex(v), e2.other_vertex(v)
      
      if e1.connection[v] == 'H':
        v_wells1 = v.head_wells
      elif e1.connection[v] == 'T':
        v_wells1 = v.tail_wells
      if e2.connection[v] == 'H':
        v_wells2 = v.head_wells
      elif e2.connection[v] == 'T':
        v_wells2 = v.tail_wells

      if e1.connection[w1] == 'H':
        w1_wells = w1.head_wells
        w1_connectivity = len(w1.head_edges)
      elif e1.connection[w1] == 'T':
        w1_wells = w1.tail_wells
        w1_connectivity = len(w1.tail_edges)
      if e2.connection[w2] == 'H':
        w2_wells = w2.head_wells
      elif e2.connection[w2] == 'T':
        w2_wells = w2.tail_wells

      print v.id, 'H', len(v), v.intervals
      print 'v wells:', ','.join([str(w) for w in v_wells1])
      print 'e1/w1:', e1.id, w1.id, len(w1), w1.intervals 
      print 'w1 wells:', ','.join([str(w) for w in w1_wells])
      print 'e2/w2:', e2.id, w2.id, len(w2), w2.intervals 
      print 'w2 wells:', ','.join([str(w) for w in w2_wells])
      print e1.support, e2.support, abs_support_delta, rel_support_delta, w1_connectivity
      print

# ----------------------------------------------------------------------------
# helpers

def _select_edge(g, E, e_selected):
  to_remove = [e for e in E]
  for e in to_remove:
    if e != e_selected:
      assert e.v1 != e.v2
      g.remove_edge(e)

def _get_wells_between_v(v1, v2, conn1, conn2):
  if conn1 == 'H':
    v1_wells = v1.head_wells
  elif conn1 == 'T':
    v1_wells = v1.tail_wells
  if conn2 == 'H':
    v2_wells = v2.head_wells
  elif conn2 == 'T':
    v2_wells = v2.tail_wells

  return v1_wells & v2_wells, v1_wells, v2_wells

def _get_common_wells(e):
  v1, v2 = e.v1, e.v2

  if e.connection[v1] == 'H':
    v1_wells = v1.head_wells
  elif e.connection[v1] == 'T':
    v1_wells = v1.tail_wells
  if e.connection[v2] == 'H':
    v2_wells = v2.head_wells
  elif e.connection[v2] == 'T':
    v2_wells = v2.tail_wells

  return v1_wells & v2_wells