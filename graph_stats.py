from libkuleshov.stats import n50

def graph_n50(g):
	return n50([len(v) for v in g.vertices])

def graph_avg(g):
	L = [len(v) for v  in g.vertices]
	return float(sum(L)) / len(L)