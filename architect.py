#!/usr/bin/env python
import sys
import argparse
import logging
import pickle

# save/load graphs
from graph_load import load_from_sga_asqg, load_from_asqg, \
											 load_from_fasta_tsv, \
											 save_graph, save_fasta, save_layout, save_bandage_gfa, \
											 save_to_fasta_tsv, \
											 unpickle_graph, pickle_graph

# save file in dot format
import visualize
from visualize import to_graphviz_dot, to_graphviz_dot_with_intervals, \
         			  to_graphviz_dot_with_double_intervals, \
         			  to_graphviz_dot_with_connections, \
         			  to_graphviz_dot_with_markup

# visualize parts of assembly graph
from visualize import print_connection, print_repeat

# collect statistics about the graph
from graph_stats import print_stats

# contract edges that are similar
from contraction import contract_edges, remove_loops, remove_parallel_edges

# method to verify graph correctness
from verificator import examine_repeats, examine_connections

from scaffolder import prune_scaffold_edges, cut_tips, \
											 prune_scaffold_edges_via_wells, \
											 scaffold_via_wells, \
											 examine_scaffold_ambiguities, \
											 make_wellscaff_edges

from graph_algorithms import scaffold_via_wells_mst

# ----------------------------------------------------------------------------

def main():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(title='Commands')

	## SCAFFOLDER

	scaffold_parser = subparsers.add_parser('scaffold',
											help='Order contigs/scaffold using wells')
	scaffold_parser.set_defaults(func=scaffold)

	scaffold_parser.add_argument('--fasta', required=True,
		help='Input scaffolds/contigs')
	scaffold_parser.add_argument('--edges', required=True,
		help='Known paired-end or overlap connections')
	scaffold_parser.add_argument('--containment', required=True,
		help='Container hits and various meta-data')
	scaffold_parser.add_argument('--out', required=True,
		help='Prefix for the ouput files')
	scaffold_parser.add_argument('--min-ctg-len', type=int, default=0,
		help='Discard contigs smaller than this length (def: 0)')
	scaffold_parser.add_argument('--log',
		help='Save stdout to log file')

	## VIEW STATISTICS

	view_parser = subparsers.add_parser('view',
											help='View statistics about the graph')
	view_parser.set_defaults(func=view)

	view_parser.add_argument('--fasta', required=True,
		help='Input scaffolds/contigs')
	view_parser.add_argument('--edges', required=True,
		help='Known paired-end or overlap connections')
	view_parser.add_argument('--containment', required=True,
		help='Container hits and various meta-data')
	view_parser.add_argument('--edge', type=int,
		help='View neighborhood around particular edge id')
	view_parser.add_argument('--vertex', type=int,
		help='View neighborhood around particular vertex id')
	view_parser.add_argument('--dot',
		help='Dotfile name for visualization')
	view_parser.add_argument('--log',
		help='Save stdout to log file')

	args = parser.parse_args()

	if args.log:
		logging.basicConfig(filename=args.log, 
							filemode='a',
							level=logging.INFO)

	args.func(args)

# ----------------------------------------------------------------------------
# entry point functions

def view(args):
	# g = load_from_asqg(args.inp + '.asqg', 
  #        	       args.inp + '.containment')
	g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
	save_to_fasta_tsv(g, 'test.fasta', 'test.tsv', 'test.containment')
	g = load_from_fasta_tsv('test.fasta', 'test.tsv', 'test.containment')

	visualize.visualize_well_correctness(g)

	if args.edge:
		print_connection(g.get_edge(args.edge))
	elif args.dot:
		to_graphviz_dot_with_connections(g, try_to_resolve_new, args.dot)

def scaffold(args):
	if False:
		g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
		print_stats(g)

		# delete small vertices
		if args.min_ctg_len:
			n_removed = 0
			for v in g.vertices:
				if len(v.seq) < args.min_ctg_len:
					g.remove_vertex(v)
					n_removed += 1
			print 'Removed %d vertices' % n_removed
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
		# save_fasta(g, 'pruned.fasta')
		
		# delete all existing edges from the graph
		E = g.edges
		for e in E:
			g.remove_edge(e)

		# create new edges whenever vertices have similar well profiles
		# n_edges = make_wellscaff_edges(g)
		n_edges = make_wellscaff_edges(g, min_common=3, min_thr=0.2)
		print '%d scaffold edges created via wells.' % n_edges

		save_to_fasta_tsv(g, 'wellscaff.fasta', 'wellscaff.tsv', 'wellscaff.containment')

	g = load_from_fasta_tsv('wellscaff.fasta', 'wellscaff.tsv', 'wellscaff.containment',
													min_supp=1)
	print_stats(g)
	# scaffold_via_wells(g)
	scaffold_via_wells(g, min_common=3, min_thr=0.2)
	# scaffold_via_wells_mst(g)
	print_stats(g)
	save_fasta(g, '%s.fasta' % args.out)
	save_layout(g, '%s.layout' % args.out)
	# pickle_graph(g, 'scaffolded.pkl')
	# save_bandage_gfa(g, '%s.gfa' % args.out)


# ----------------------------------------------------------------------------
# helpers

def save_optional_output(g, args):
	if args.dot:
		to_graphviz_dot_with_intervals(g, args.dot)
	if args.stats:
		print_stats(g, args.stats)


if __name__ == '__main__':
	main()