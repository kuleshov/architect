import sys

def print_stats(g, stats_file=sys.stdout):
  if stats_file == sys.stdout:
    stats = stats_file
  else:
    stats = open(stats_file, 'w')

  stats.write("Vertices: %d\n" % len(g.vertices))
  stats.write("Edges: %d\n" % len(g.edges))
  stats.write("Connected components: %d\n" % g.count_connected_components())
  stats.write("N50: %d\n" % graph_n50(g))
  stats.write("Idealized N50: %d\n" % g.idealized_n50())
  stats.write("Average contig length: %d\n" % graph_avg(g))

  if stats != sys.stdout:
    stats.close()

def graph_n50(g):
	return _n50([len(v) for v in g.vertices])

def graph_avg(g):
	L = [len(v) for v  in g.vertices]
	return float(sum(L)) / len(L)

# ----------------------------------------------------------------------------
# helpers

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