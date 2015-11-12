#!/usr/bin/env python
import sys
import argparse
import logging

# get graph from SGA's asqg; then load/store from internal asqg-like format
from graph_load import load_from_sga_asqg, load_from_asqg, \
											 load_from_fasta_tsv, \
											 save_graph, save_fasta, save_layout, \
											 unpickle_graph, pickle_graph

# save file in dot format
import visualize
from visualize import to_graphviz_dot, to_graphviz_dot_with_intervals, \
         			  to_graphviz_dot_with_double_intervals, \
         			  to_graphviz_dot_with_connections, \
         			  to_graphviz_dot_with_markup

# visualize parts of assembly graph
from visualize import print_connection, print_repeat

# examine misassemblies
from misassemblies import examine_misassemblies

# collect statistics about the graph
from graph_stats import print_stats

# contract edges that are similar
from contraction import contract_edges, remove_loops, remove_parallel_edges

# graph simplification operation that involve counting wells (PMMP)
from pmpp import resolve_repeats, delete_spurious_edges, get_wells_by_edge, \
				 try_to_resolve_new

# traverse graph according to well information
from traversal import compute_traversals

# method to verify graph correctness
from verificator import examine_repeats, examine_connections

# pop bubbles
from bubbles import pop_triangles, simplify_graph, check_bubbles, \
					detect_bubble_from_pointers

# resolve short repeats in the graph by using lengths
from resolve_repeats import resolve_short_repeats

from scaffolder import prune_scaffold_edges, cut_tips, \
											 prune_scaffold_edges_via_wells, \
											 scaffold_via_wells, \
											 examine_scaffold_ambiguities

from libkuleshov.debug import keyboard
from intervals import print_true_intervals

# ----------------------------------------------------------------------------

def main():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(title='Commands')

	## LOAD

	load_parser = subparsers.add_parser('load')
	load_parser.set_defaults(func=load)

	load_parser.add_argument('--asqg', required=True)
	load_parser.add_argument('--out', required=True)
	load_parser.add_argument('--log')

	## CONTRACT

	contract_parser = subparsers.add_parser('contract')
	contract_parser.set_defaults(func=contract)

	contract_parser.add_argument('--inp', required=True)
	contract_parser.add_argument('--out', required=True)
	contract_parser.add_argument('--dot')
	contract_parser.add_argument('--stats')
	contract_parser.add_argument('--log')
	contract_parser.add_argument('--masm', action='store_true')

	## DELETE SPURIOUS

	spurious_parser = subparsers.add_parser('prune')
	spurious_parser.set_defaults(func=spurious)

	spurious_parser.add_argument('--inp', required=True)
	spurious_parser.add_argument('--out', required=True)
	spurious_parser.add_argument('--dot')
	spurious_parser.add_argument('--stats')
	spurious_parser.add_argument('--log')
	spurious_parser.add_argument('--masm', action='store_true')

	## RESOLVE SHORT REPEATS

	short_repeats_parser = subparsers.add_parser('resolve-short-repeats')
	short_repeats_parser.set_defaults(func=short_repeats)

	short_repeats_parser.add_argument('--inp', required=True)
	short_repeats_parser.add_argument('--out', required=True)
	short_repeats_parser.add_argument('--dot')
	short_repeats_parser.add_argument('--stats')
	short_repeats_parser.add_argument('--log')
	short_repeats_parser.add_argument('--masm', action='store_true')

	## RESOLVE REPEATS BY WELLS

	repeats_parser = subparsers.add_parser('resolve')
	repeats_parser.set_defaults(func=repeats)

	repeats_parser.add_argument('--inp', required=True)
	repeats_parser.add_argument('--out', required=True)
	repeats_parser.add_argument('--dot')
	repeats_parser.add_argument('--stats')
	repeats_parser.add_argument('--log')
	repeats_parser.add_argument('--masm', action='store_true')

	## TRAVERSE REGIONS CONNECTED BY WELLS

	traverse_parser = subparsers.add_parser('traverse')
	traverse_parser.set_defaults(func=traverse)

	traverse_parser.add_argument('--inp', required=True)
	traverse_parser.add_argument('--out', required=True)
	traverse_parser.add_argument('--dot')
	traverse_parser.add_argument('--stats')
	traverse_parser.add_argument('--log')
	traverse_parser.add_argument('--masm', action='store_true')

	## POP BUBBLES

	bubbles_parser = subparsers.add_parser('pop-bubbles')
	bubbles_parser.set_defaults(func=bubbles)

	bubbles_parser.add_argument('--inp', required=True)
	bubbles_parser.add_argument('--out', required=True)
	bubbles_parser.add_argument('--dot')
	bubbles_parser.add_argument('--stats')
	bubbles_parser.add_argument('--log')
	bubbles_parser.add_argument('--masm', action='store_true')

	## SCAFFOLDER

	scaffold_parser = subparsers.add_parser('scaffold')
	scaffold_parser.set_defaults(func=scaffold)

	scaffold_parser.add_argument('--fasta', required=True)
	scaffold_parser.add_argument('--edges', required=True)
	scaffold_parser.add_argument('--containment')
	scaffold_parser.add_argument('--log')

	wellscaffold_parser = subparsers.add_parser('wellwellscaffold')
	wellscaffold_parser.set_defaults(func=wellscaffold)

	wellscaffold_parser.add_argument('--pkl', required=True)
	wellscaffold_parser.add_argument('--log')

	## VIEW STATISTICS

	view_parser = subparsers.add_parser('view')
	view_parser.set_defaults(func=view)

	view_parser.add_argument('--fasta')
	view_parser.add_argument('--edges')
	view_parser.add_argument('--containment')
	view_parser.add_argument('--inp')
	view_parser.add_argument('--edge', type=int)
	view_parser.add_argument('--vertex', type=int)
	view_parser.add_argument('--dot')
	view_parser.add_argument('--log')

	args = parser.parse_args()

	if args.log:
		logging.basicConfig(filename=args.log, 
							filemode='a',
							level=logging.INFO)

	logging.info('ok')

	args.func(args)

# ----------------------------------------------------------------------------
# entry point functions

def load(args):
	g = load_from_sga_asqg(args.asqg)

	save_graph(g, args.out + '.asqg', 
        		  args.out + '.containment')

def contract(args):
	g = load_from_asqg(args.inp + '.asqg', args.inp + '.containment')

	contract_edges(g)
	save_optional_output(g, args)

	save_graph(g, args.out + '.asqg', 
        		  args.out + '.containment')

def spurious(args):
	g = load_from_asqg(args.inp + '.asqg', 
        	       args.inp + '.containment')

	for i in xrange(4):
	  examine_connections(g, conservative='very')
	  delete_spurious_edges(g, conservative='very')
	  contract_edges(g)

	print 'SECOND RUN'

	examine_connections(g, conservative='yes')
	delete_spurious_edges(g, conservative='yes')
	contract_edges(g)

	save_optional_output(g, args)
	examine_misassemblies(g)

	save_graph(g, args.out + '.asqg', 
        		  args.out + '.containment')

def short_repeats(args):
	g = load_from_asqg(args.inp + '.asqg', 
        	       args.inp + '.containment')

	resolve_short_repeats(g)
	save_optional_output(g, args)

	save_graph(g, args.out + '.asqg', 
        		  args.out + '.containment')

def repeats(args):
	g = load_from_asqg(args.inp + '.asqg', 
        	       args.inp + '.containment')

	V_odd = [v for v in g.vertices if len(v.edges) % 2 == 1]

	resolve_repeats(g, V_odd)
	contract_edges(g)
	# resolve_repeats(g, V_odd)

	V_even = [v for v in g.vertices if len(v.edges) % 2 == 0]

	# for v in V_even:
	# 	if v.id_ == 3762621:
	# 		print_true_intervals(v, [ctg for ctg in v.metadata['contigs']])
	# 		# keyboard()

	# resolve_repeats(g, V_even)
	to_graphviz_dot_with_connections(g, try_to_resolve_new, args.dot)
	examine_repeats(g)
	# save_optional_output(g, args)

	save_graph(g, args.out + '.asqg', 
        		  args.out + '.containment')

def traverse(args):
	g = load_from_asqg(args.inp + '.asqg', 
        	       args.inp + '.containment')


	E, forward_pointers = compute_traversals(g)
	check_bubbles(g, forward_pointers)
	to_graphviz_dot_with_markup(g, [[]], [E], args.dot)


def bubbles(args):
	g = load_from_asqg(args.inp + '.asqg', 
        	       args.inp + '.containment')

	E, forward_pointers = compute_traversals(g)

	# keyboard()

	v0 = g.vertices_by_id[3764088]
	e0 = g.edges_by_id[106]

	v0, v1 = detect_bubble_from_pointers(g, v0, e0, forward_pointers, 1000)

	print v0.id_, v1.id_

	to_graphviz_dot_with_markup(g, [[]], [E], args.dot)

	save_graph(g, args.out + '.asqg', 
        		  args.out + '.containment')

def view(args):
	# g = load_from_asqg(args.inp + '.asqg', 
  #        	       args.inp + '.containment')
	g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
	visualize.visualize_well_correctness(g)

	if args.edge:
		print_connection(g.get_edge(args.edge))
	elif args.dot:
		to_graphviz_dot_with_connections(g, try_to_resolve_new, args.dot)

def scaffold(args):
	g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
	print_stats(g)
	contract_edges(g, store_layout=True)
	print_stats(g)
	save_fasta(g, 'contracted.fasta')
	# n_cut = cut_tips(g)
	# print '%d tips were cut.' % n_cut
	# examine_scaffold_ambiguities(g)
	n_pruned = prune_scaffold_edges(g)
	print '%d edges pruned.' % n_pruned
	n_pruned = prune_scaffold_edges_via_wells(g)
	print '%d edges pruned via wells.' % n_pruned
	contract_edges(g)
	print_stats(g)
	save_fasta(g, 'pruned.fasta')
	scaffold_via_wells(g)
	print_stats(g)
	save_fasta(g, 'pmmp.fasta')
	save_layout(g, 'pmmp.layout')
	# pickle_graph(g, 'scaffolded.pkl')

def wellscaffold(args):
	g = unpickle_graph('scaffolded.pkl')
	scaffold_via_wells(g)
	print_stats(g)
	save_fasta(g, 'pmmp.fasta')

# ----------------------------------------------------------------------------
# helpers

def save_optional_output(g, args):
	if args.dot:
		to_graphviz_dot_with_intervals(g, args.dot)
	if args.stats:
		print_stats(g, args.stats)
	if args.masm:
		examine_misassemblies(g)


if __name__ == '__main__':
	main()