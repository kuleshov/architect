import sys
import logging

def print_stats(g):
  stats_str1 = "Graph stats: %d vertices, %d edges, %d connected components" \
              % (len(g.vertices), len(g.edges), g.count_connected_components())
  stats_str2 = "Graph stats: N50: %d; Max achievable N50: %d; Total length: %d" \
              % (_graph_n50(g), g.idealized_n50(), _graph_tot(g))
           
  logging.info(stats_str1)
  logging.info(stats_str2)

  # stats.write("Vertices: %d\n" % len(g.vertices))
  # stats.write("Edges: %d\n" % len(g.edges))
  # stats.write("Connected components: %d\n" % g.count_connected_components())
  # stats.write("N50: %d\n" % _graph_n50(g))
  # stats.write("Idealized N50: %d\n" % g.idealized_n50())
  # stats.write("Total/Average contig length: %d\t%d\n" % _graph_avg(g))

# ----------------------------------------------------------------------------
# helpers

def _graph_n50(g):
  return _n50([len(v) for v in g.vertices])

def _graph_tot(g):
  L = [len(v) for v  in g.vertices]
  return sum(L)

def _graph_avg(g):
  L = [len(v) for v  in g.vertices]
  return float(sum(L)) / len(L), sum(L)

def _n50(l, total_length=None):
  if not total_length:
    total_length = sum(l)
  
  half_length = total_length / 2
  length = 0

  l.sort(reverse=True)

  while length < half_length:
    x = l.pop()
    length += x

  return x