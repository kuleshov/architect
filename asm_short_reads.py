#!/usr/bin/env python
import argparse
from graph_load import load_from_asqg, to_graphviz_dot
from graph_stats import graph_n50, graph_avg

from contraction import *
from pmpp import resolve_repeats, examine_repeats
from remove_transitive import pop_triangles

from libkuleshov.debug import keyboard

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--asqg')
parser.add_argument('--all_asqg')

args = parser.parse_args()

##############################################################################			
## STATS

def print_stats(g):
	print "Vertices:", len(g.vertices)
	print "Edges:", len(g.edges)
	print "Connected components:", g.count_connected_components()
	print "N50:", graph_n50(g)
	print "Average contig length:", graph_avg(g)

g = load_from_asqg(args.asqg)
print_stats(g)

for v in g.vertices:
	assert len(v.metadata['contigs']) > 0

##############################################################################			
## CONTRACT PATHS

print "Contracting paths..."

candidate_edges = set(g.edges)

while candidate_edges:
	e = candidate_edges.pop()
	# if e.v1.id_ == '6074-1_1044-1': keyboard()
	# E = g.edges.copy()
	if can_be_contracted(e): 
		contract_edge(g,e)

for v in g.vertices:
	assert len(v.metadata['contigs']) > 0
	
print_stats(g)

##############################################################################			
## REMOVE TRANSITIVE EDGES:

print "Removing triangles..."

pop_triangles(g)

candidate_edges = set(g.edges)

while candidate_edges:
	e = candidate_edges.pop()
	# if e.v1.id_ == '6074-1_1044-1': keyboard()
	# E = g.edges.copy()
	if can_be_contracted(e): 
		contract_edge(g,e)
	
print_stats(g)

###############################################################################
## RESOLVE REPEATS

print "Resolving repeats..."

resolve_repeats(g)

candidate_edges = set(g.edges)

while candidate_edges:
	e = candidate_edges.pop()
	# if e.v1.id_ == '6074-1_1044-1': keyboard()
	# E = g.edges.copy()
	if can_be_contracted(e): 
		contract_edge(g,e)

print_stats(g)

examine_repeats(g)

# ###############################################################################
# ## REMOVE SMALL VERTICES

# for v in g.vertices:
# 	if len(v) < 8000:
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


to_graphviz_dot(g, 'out.dot')

###############################################################################
## LOAD OVERLAP STORE

# overlap_store = load_overlap_store(args.all_asqg)

# examine_repeats(g, overlap_store)