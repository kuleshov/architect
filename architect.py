#!/usr/bin/env python
import sys
import argparse
import logging
import pickle

# save/load graphs
from graph_load import load_from_fasta_tsv, \
                       save_fasta, save_layout, save_bandage_gfa, \
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

# scaffolding subroutines
from scaffolder import prune_scaffold_edges, cut_tips, \
                       prune_scaffold_edges_via_wells, \
                       prune_via_wells, \
                       make_wellscaff_edges

# ----------------------------------------------------------------------------

def main():
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(title='Commands')

  # scaffolder

  scaffold_parser = subparsers.add_parser('scaffold',
   help='Scaffold a set of input contigs')
  scaffold_parser.set_defaults(func=scaffold)

  scaffold_parser.add_argument('--fasta', required=True,
    help='Input scaffolds/contigs')
  scaffold_parser.add_argument('--edges', required=False,
    help='Known paired-end or overlap connections')
  scaffold_parser.add_argument('--containment', required=True,
    help='Container hits and various meta-data')
  scaffold_parser.add_argument('--out', required=True,
    help='Prefix for the ouput files')
  scaffold_parser.add_argument('--min-ctg-len', type=int, default=0,
    help='Discard contigs smaller than this length (def: 0)')
  scaffold_parser.add_argument('--cut-tip-len', type=int, default=0,
    help='Cut tips smaller than this length')
  scaffold_parser.add_argument('--log',
    help='Save stdout to log file')

  # viewer

  view_parser = subparsers.add_parser('view',
   help='Display information about the scaffold graph')
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

  FORMAT = '[architect] %(message)s'
  logging.basicConfig(filename=args.log, 
   filemode='a',
   level=logging.INFO,
   format=FORMAT)

  args.func(args)

# ----------------------------------------------------------------------------
# entry point functions

def view(args):
  g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
  save_to_fasta_tsv(g, 'test.fasta', 'test.tsv', 'test.containment')
  g = load_from_fasta_tsv('test.fasta', 'test.tsv', 'test.containment')

  visualize.visualize_well_correctness(g)

  if args.edge:
    print_connection(g.get_edge(args.edge))
  elif args.dot:
    to_graphviz_dot_with_connections(g, try_to_resolve_new, args.dot)

def scaffold(args):
  logging.info('Creating the scaffold graph')
  g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
  logging.info('Statistics of the loaded graph:')
  print_stats(g)

  # delete small vertices
  if args.min_ctg_len:
    logging.info('Removing vertices smaller than %d bp' % args.min_ctg_len)
    n_removed = 0
    for v in g.vertices:
      if len(v.seq) < args.min_ctg_len:
        g.remove_vertex(v)
        n_removed += 1
        logging.info('Removed %d vertices' % n_removed)
        print_stats(g)

  # prune scaffold edges
  if g.edges:
    logging.info('Simplifying the graph using paired-end reads')
    logging.info('Contracting unambigous paths')
    contract_edges(g, store_layout=True)
    logging.info('Statistics after contraction:')
    print_stats(g)
    save_fasta(g, 'contracted.fasta')
    
    if args.cut_tip_len:
      n_cut = cut_tips(g, d=args.cut_tip_len)
      logging.info('Cut %d tips shorter than %d bp' % (n_cut, args.cut_tip_len))
    # examine_scaffold_ambiguities(g)
    
    logging.info('Pruning edges with low support')
    n_pruned1 = prune_scaffold_edges(g)
    n_pruned2 = prune_scaffold_edges_via_wells(g)
    logging.info('%d edges pruned' % (n_pruned1 + n_pruned2))

    logging.info('Contracting unambigous paths')
    n_contracted = contract_edges(g)
    logging.info('Statistics after contraction:')
    print_stats(g)

  # delete all existing edges from the graph
  E = g.edges
  for e in E:
    g.remove_edge(e)

  # create new edges whenever vertices have similar well profiles
  logging.info('Creating edges from read clouds')
  n_edges = make_wellscaff_edges(g, min_common=3, min_thr=0.2)
  logging.info('%d scaffold edges from read clouds' % n_edges)

  # save_to_fasta_tsv(g, 'wellscaff.fasta', 'wellscaff.tsv', 'wellscaff.containment')
  # g = load_from_fasta_tsv('wellscaff.fasta', 'wellscaff.tsv', 'wellscaff.containment',
  #                         min_supp=1)

  logging.info('Pruning edges with low support')
  n_pruned = prune_via_wells(g, min_common=3, min_thr=0.2)
  logging.info('%d edges pruned' % n_pruned)

  logging.info('Contracting unambigous paths')
  n_contracted = contract_edges(g, store_layout=True)
  logging.info('Statistics after contraction:')
  print_stats(g)

  logging.info('Saving scaffolding results')
  save_fasta(g, '%s.fasta' % args.out)
  save_layout(g, '%s.layout' % args.out)

# ----------------------------------------------------------------------------

if __name__ == '__main__':
  main()
