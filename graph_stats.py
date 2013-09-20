from libkuleshov.stats import n50

def print_stats(g):
	print "Vertices:", len(g.vertices)
	print "Edges:", len(g.edges)
	print "Connected components:", g.count_connected_components()
	print "N50:", graph_n50(g)
	print "Average contig length:", graph_avg(g)

def graph_n50(g):
	return n50([len(v) for v in g.vertices])

def graph_avg(g):
	L = [len(v) for v  in g.vertices]
	return float(sum(L)) / len(L)