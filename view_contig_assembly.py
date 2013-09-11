from graph_load import load_from_sga_asqg, to_graphviz_dot, save_graph, load_graph
from graph_stats import graph_n50, graph_avg
from intervals import get_intervals

##############################################################################			
## LOAD

g = load_graph('graph.asqg', 'graph.containment')

##############################################################################			

def get_true_intervals(w, ctgs):
	I = list()
	for ctg in ctgs:
		fields = ctg.split('_')
		chrom, coords = fields[1].split(':')
		start, end = (int(x) for x in coords.split('-'))
		internal_start = int(fields[2])

		I.append((int(chrom), start + internal_start, start + internal_start + 200))

	return get_intervals(I)

##############################################################################			

I = list()
for v in g.vertices:
	v_intervals = get_true_intervals(v, v.metadata['contigs'])
	for i in v_intervals:
		I.append((i, v.id_))

for i in sorted(I):
	print i[0], '\t', i[1]