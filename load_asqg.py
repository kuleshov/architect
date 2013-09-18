#!/usr/bin/env python
import argparse
from graph_load import load_from_sga_asqg, save_graph, load_graph
from visualize import to_graphviz_dot, to_graphviz_dot_with_intervals, \
					  examine_misassemblies
from graph_stats import graph_n50, graph_avg

from contraction import contract_edges
from pmpp import resolve_repeats, examine_repeats, delete_spurious_edges, get_well, examine_connections
from remove_transitive import pop_triangles

from libkuleshov.debug import keyboard

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--asqg')
parser.add_argument('--all_asqg')

args = parser.parse_args()

def find_read(g, id_):
	for v in g.vertices:
		if id_ in v.metadata['contigs']:
			for ctg in v.metadata['contigs']:
				print v.metadata['contig_starts'][ctg], get_well(ctg)
			return v


##############################################################################			
## STATS

def print_stats(g):
	print "Vertices:", len(g.vertices)
	print "Edges:", len(g.edges)
	print "Connected components:", g.count_connected_components()
	print "N50:", graph_n50(g)
	print "Average contig length:", graph_avg(g)

##############################################################################			
## LOAD

g = load_graph('graph.asqg', 'graph.containment')
print_stats(g)

# g = load_from_sga_asqg(args.asqg)
# print_stats(g)

# ##############################################################################			
# ## CONTRACT PATHS

contract_edges(g)	
print_stats(g)

# ##############################################################################			
# ## DELETE SPURIOUS EDGES

delete_spurious_edges(g)
# # examine_connections(g)

# contract_edges(g)
# print_stats(g)

# # examine_misassemblies(g)

# save_graph(g, 'graph.asqg', 'graph.containment')
# # exit()


##############################################################################			
## RESOLVE REPEATS

for i in xrange(4):
	resolve_repeats(g, wells='edges')

# # # contract_edges(g)
examine_repeats(g)
# print_stats(g)

# save_graph(g, 'graph.asqg', 'graph.containment')

# examine_connections(g, conservative=False)
# delete_spurious_edges(g, conservative=False)

to_graphviz_dot_with_intervals(g, 'out.dot')

exit()

##############################################################################			
## REMOVE TRANSITIVE EDGES:

pop_triangles(g)
print_stats(g)
examine_misassemblies(g)

save_graph(g, 'graph.asqg', 'graph.containment')

to_graphviz_dot_with_intervals(g, 'out.dot')
exit()

###############################################################################
## RESOLVE REPEATS

print "Resolving repeats..."

for i in xrange(100):
	resolve_repeats(g)
	contract_edges(g)

print_stats(g)
examine_repeats(g)
print_stats(g)

###############################################################################
## REMOVE SMALL VERTICES

# for v in g.vertices:
# 	if len(v) < 500:
# 		g.remove_vertex(v)

# print_stats(g)

# ##############################################################################			
# ## COMPUTE LENGTH OF WRONG AND CORRECT OVERLAPS

# L_bad = list()
# L_good = list()
# for e in g.edges:
# 	v1, v2 = e.v1, e.v2
# 	if e.connection[v1] == e.connection[v2]:
# 		# erroneous edge!
# 		L_bad.append(abs(e.ovl_end[v1] - e.ovl_start[v1]))
# 	else:
# 		# could be good edge:
# 		L_good.append(abs(e.ovl_end[v1] - e.ovl_start[v1]))

# print float(sum(L_good)) / len(L_good)
# # print float(sum(L_bad)) / len(L_bad)


# to_graphviz_dot(g, 'out.dot')
to_graphviz_dot_with_intervals(g, 'out.dot')
# save_graph(g, 'graph.asqg', 'graph.containment')

###############################################################################
## LOAD OVERLAP STORE

# overlap_store = load_overlap_store(args.all_asqg)

# examine_repeats(g, overlap_store)