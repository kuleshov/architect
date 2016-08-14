#!/usr/bin/env python
import sys
import argparse
import logging
import pickle

# save/load graphs
from graph.load import load_from_fasta_tsv, \
                       save_fasta, save_ordering, save_bandage_gfa, \
                       save_to_fasta_tsv, \
                       save_bandage_gfa

# save file in dot format
from common import visualize
from common.visualize import to_graphviz_dot, to_abyss_explorer_dot

# visualize parts of assembly graph
from common.visualize import print_connection, print_vertex

# collect statistics about the graph
from graph.stats import print_stats

# contract edges that are similar
from algorithms.contraction import contract_edges, remove_loops, \
                                   remove_parallel_edges

# scaffolding subroutines
from algorithms.scaffolder import prune_scaffold_edges, cut_tips, \
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
    help='Known paired-end or read cloud connections')
  scaffold_parser.add_argument('--containment', required=True,
    help='Container hits and various meta-data')
  scaffold_parser.add_argument('--out', required=True,
    help='Prefix for the ouput files')
  scaffold_parser.add_argument('--min-ctg-len', type=int, default=0,
    help='Discard contigs smaller than this length (def: 0)')
  scaffold_parser.add_argument('--cut-tip-len', type=int, default=0,
    help='Cut tips smaller than this length')
  scaffold_parser.add_argument('--pe-abs-thr', type=int, default=3,
    help='Threshold for absolute support when pruning paired-end edges')
  scaffold_parser.add_argument('--pe-rel-thr', type=float, default=0.7,
    help='Threshold for relative support when pruning paired-end edges')
  scaffold_parser.add_argument('--pe-rc-rel-thr', type=float, default=0.5,
    help='Threshold for relative support for read-cloud / paired-end pruning')
  scaffold_parser.add_argument('--rc-abs-thr', type=int, default=3,
    help='Minimum support for create read-cloud based edge')
  scaffold_parser.add_argument('--rc-rel-edge-thr', type=float, default=0.2,
    help='Threshold for relative support when creating read-cloud based edges')
  scaffold_parser.add_argument('--rc-rel-prun-thr', type=float, default=0.2,
    help='Threshold for relative support when pruning read-cloud based edges')
  scaffold_parser.add_argument('--log',
    help='Save stdout to log file')

  # viewer

  view_parser = subparsers.add_parser('view',
   help='Display information about the scaffold graph')
  view_parser.set_defaults(func=view)

  view_parser.add_argument('--fasta', required=True,
    help='Input scaffolds/contigs')
  view_parser.add_argument('--edges', required=True,
    help='Known paired-end or read cloud connections')
  view_parser.add_argument('--containment', required=True,
    help='Container hits and various meta-data')
  view_parser.add_argument('--edge', type=int,
    help='View neighborhood around particular edge id')
  view_parser.add_argument('--vertex', type=int,
    help='View neighborhood around particular vertex id')
  view_parser.add_argument('--check-correctness', action='store_true',
    help='Check if edges are correct using true intervals')
  view_parser.add_argument('--dot',
    help='Dotfile for visualization in Abyss Explorer')
  view_parser.add_argument('--gfa',
    help='GFA file for visualization in Bandage')
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
  # save_to_fasta_tsv(g, 'test.fasta', 'test.tsv', 'test.containment')
  # g = load_from_fasta_tsv('test.fasta', 'test.tsv', 'test.containment')

  if args.check_correctness:
    visualize_correctness(g)

  if args.vertex:
    print_vertex(g.vertex_from_id(args.vertex))
  if args.edge:
    print_connection(g.get_edge(args.edge))
  
  if args.dot:
    to_graphviz_dot(g, args.dot)
  if args.gfa:
    save_bandage_gfa(g, args.gfa)

def scaffold(args):
  logging.info('Creating the scaffold graph')
  g = load_from_fasta_tsv(args.fasta, args.edges, args.containment)
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
    contract_edges(g, store_ordering=True)
    print_stats(g)
    save_fasta(g, 'contracted.fasta')
    
    if args.cut_tip_len:
      n_cut = cut_tips(g, d=args.cut_tip_len)
      logging.info('Cut %d tips shorter than %d bp' \
                    % (n_cut, args.cut_tip_len))
    
    logging.info('Pruning edges with low support')
    n_pruned1 = prune_scaffold_edges(g, abs_support_thr=args.pe_abs_thr, 
                                        rel_support_thr=args.pe_rel_thr)
    n_pruned2 = prune_scaffold_edges_via_wells(g, thr=args.pe_rc_rel_thr)
    logging.info('%d edges pruned' % (n_pruned1 + n_pruned2))

    logging.info('Contracting unambigous paths')
    n_contracted = contract_edges(g)
    print_stats(g)

  # delete all existing edges from the graph
  E = g.edges
  for e in E:
    g.remove_edge(e)

  # create new edges whenever vertices have similar well profiles
  logging.info('Creating edges from read clouds')
  n_edges = make_wellscaff_edges(g, min_common=args.rc_abs_thr, 
                                    min_thr=args.rc_rel_edge_thr)
  logging.info('%d scaffold edges from read clouds' % n_edges)

  logging.info('Auto-saving graph with prefix %s.wellscaff' % args.out)
  save_to_fasta_tsv(g, '%s.wellscaff.fasta' % args.out, 
                       '%s.wellscaff.tsv' % args.out, 
                       '%s.wellscaff.containment' % args.out)
  g = load_from_fasta_tsv('%s.wellscaff.fasta' % args.out, 
                          '%s.wellscaff.tsv' % args.out, 
                          '%s.wellscaff.containment' % args.out, 
                           min_supp=1)

  logging.info('Pruning edges with low support')
  n_pruned = prune_via_wells(g, min_common=args.rc_abs_thr, 
                                min_thr=args.rc_rel_prun_thr)
  logging.info('%d edges pruned' % n_pruned)

  logging.info('Contracting unambigous paths')
  n_contracted = contract_edges(g, store_ordering=True)
  print_stats(g)

  logging.info('Saving scaffolding results')
  save_fasta(g, '%s.fasta' % args.out)
  save_ordering(g, '%s.ordering' % args.out)

# ----------------------------------------------------------------------------

if __name__ == '__main__':
  main()
