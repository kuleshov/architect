#!/usr/bin/env python
import argparse

from graph_load import load_from_sga_asqg, save_graph, load_graph, save_graph_to_fasta
from visualize import to_graphviz_dot, to_graphviz_dot_with_intervals, \
					  examine_misassemblies, to_graphviz_dot_with_double_intervals
from graph_stats import print_stats
from string_graph import test_break_vertex, break_contigs

from contraction import contract_edges
from pmpp import resolve_repeats, examine_repeats, delete_spurious_edges, get_well, examine_connections
from remove_transitive import pop_triangles
from resolve_repeats import resolve_short_repeats

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
## LOAD

CHECKPOINTDIR = 'checkpoints/working/'
# CHECKPOINTDIR = 'checkpoints/130916/130916'

# g = load_from_sga_asqg(args.asqg)

g = load_graph('../rerun_sga/reads.pre_simplify', 
			   'graph.containment')

# g = load_graph(CHECKPOINTDIR + 'graph.04.repruned.asqg', 
			  # CHECKPOINTDIR + 'graph.04.repruned.containment')

print_stats(g)
to_graphviz_dot_with_intervals(g, 'out.dot')
exit()

# contract_edges(g)
# break_contigs(g)
# to_graphviz_dot_with_double_intervals(g, 'out.dot')
# save_graph_to_fasta(g, 'graph.fasta', 'graph.containment')
# exit()
# contract_edges(g)
# examine_misassemblies(g)

print_stats(g)

save_graph(g, CHECKPOINTDIR + 'graph.00.loaded.asqg', 
			  CHECKPOINTDIR + 'graph.00.loaded.containment')

exit()

# ##############################################################################			
# ## CONTRACT PATHS

# contract_edges(g)	
# print_stats(g)

# save_graph(g, CHECKPOINTDIR + 'graph.01.contracted.asqg', 
# 			  CHECKPOINTDIR + 'graph.01.contracted.containment')

# ##############################################################################			
# ## DELETE SPURIOUS EDGES

# # for i in xrange(4):
# # 	examine_connections(g, conservative='very')
# # 	delete_spurious_edges(g, conservative='very')
# # 	contract_edges(g)

# for i in xrange(2):
# 	examine_connections(g, conservative='yes')
# 	delete_spurious_edges(g, conservative='yes')

# # contract_edges(g)
# print_stats(g)
# to_graphviz_dot_with_intervals(g, 'out.dot')

# # save_graph(g, CHECKPOINTDIR + 'graph.04.repruned.asqg', 
# # 			  CHECKPOINTDIR + 'graph.04.repruned.containment')

# # examine_misassemblies(g)

##############################################################################			
## RESOLVE REPEATS

# for i in xrange(4):
# 	resolve_repeats(g, wells='edges')

resolve_short_repeats(g)
# resolve_repeats(g)
# examine_repeats(g)
print_stats(g)

to_graphviz_dot_with_intervals(g, 'out.dot')

save_graph(g, CHECKPOINTDIR + 'graph.05.resolved.asqg', 
			  CHECKPOINTDIR + 'graph.05.resolved.containment')

exit()

# save_graph(g, 'graph.asqg', 'graph.containment')

# # examine_connections(g, conservative=False)
# delete_spurious_edges(g, conservative=False)

# save_graph(g, CHECKPOINTDIR + 'graph.04.repruned.asqg', 
# 			  CHECKPOINTDIR + 'graph.04.repruned.containment')

contract_edges(g)
to_graphviz_dot_with_intervals(g, 'out.dot')
print_stats(g)

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

