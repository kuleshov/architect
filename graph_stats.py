import sys

from libkuleshov.stats import n50

def print_stats(g, stats_file=sys.stdout):
	with open(stats_file, 'w') as stats:
		stats.write("Vertices: %d" % len(g.vertices))
		stats.write("Edges: %d" % len(g.edges))
		stats.write("Connected components: %d" % g.count_connected_components())
		stats.write("N50: %d" % graph_n50(g))
		stats.write("Average contig length: %d" % graph_avg(g))

def graph_n50(g):
	return n50([len(v) for v in g.vertices])

def graph_avg(g):
	L = [len(v) for v  in g.vertices]
	return float(sum(L)) / len(L)