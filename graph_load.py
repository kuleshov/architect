import pickle
import pysam

from graph import Graph
from string_graph import OverlapVertex, OverlapEdge, ScaffoldEdge, no_diedge
from libkuleshov.dna import reverse_complement
from libkuleshov.debug import keyboard

# ----------------------------------------------------------------------------
# constants

ASQG_VERTEX_REC = 'VT'
ASQG_EDGE_REC = 'ED'
ASQG_HEAD_REC = 'HT'

TSV_LEFT = 'L'
TSV_RIGHT = 'R'
TSV_SAME = 'S'
TSV_REVERSE = 'R'

CTMT_WELL_REC = 'W'
CTMT_IVL_REC = 'R'

# ----------------------------------------------------------------------------
# main functions

def load_from_sga_asqg(asqg_path,):
	g = _load_from_sga_asqg(asqg_path)
	return g

def load_from_asqg(asqg_path, containment_path=None):
	g, vertices_by_contig = _load_from_asqg(asqg_path)
	if containment_path:
		_load_containment(g, containment_path, vertices_by_contig)
	return g

def load_from_fasta_tsv(fasta_path, tsv_path, containment_path=None):
	g, vertices_by_contig = _load_from_fasta(fasta_path)
	_load_edges_from_tsv(g, tsv_path, vertices_by_contig)
	if containment_path:
		_load_containment(g, containment_path, vertices_by_contig)
	return g

def unpickle_graph(g, pickle_path):
	with open(pickle_path, 'rb') as f:
		g = pickle.load(f)

	return g

def save_graph(g, asqg_path, containment_path):
	_write_asqg(g, asqg_path)
	_write_containment(g, containment_path)

def save_graph_to_fasta(g, fasta_path, containment_path):
	_write_fasta(g, fasta_path)
	_write_containment(g, containment_path)

def save_fasta(g, fasta_file):
	with open(fasta_file, 'w') as fasta:
		for v in g.vertices:
			fasta.write('>' + str(v.id) + '\n')
			fasta.write(v.seq + '\n')

def pickle_graph(g, pickle_path):
	with open(pickle_path, 'wb') as f:
		pickle.dump(g, f)

# ----------------------------------------------------------------------------
# loading the graph

def _load_from_sga_asqg(asqg_path):
	# Assumptions:
	#	- vertices are listed first
	#	- no contained edges
	#	- reads overlap at endpoints only

	g = Graph()
	vertices_by_contig = dict()

	with open(asqg_path) as asqg:
		for line in asqg:
			fields = line.strip().split()

			## VERTEX
			if fields[0] == ASQG_VERTEX_REC:
				i = g.vertex_id_generator.get_id()
				v = OverlapVertex(i, fields[2])
				v.metadata['contigs'] = [fields[1]]
				v.metadata['contig_starts'] = {fields[1]: 0}
				v.metadata['contig_ends'] = {fields[1]: len(v)-1}
				vertices_by_contig[fields[1]] = v
				g.add_vertex(v)

			## EDGE
			elif fields[0] == ASQG_EDGE_REC:
				v1 = vertices_by_contig[fields[1]]
				v2 = vertices_by_contig[fields[2]]
				v1_ovl_start 	= int(fields[3])
				v1_ovl_end 		= int(fields[4])
				v1_len 			= int(fields[5])
				v2_ovl_start 	= int(fields[6])
				v2_ovl_end 		= int(fields[7])
				v2_len 			= int(fields[8])
				v2_orientation 	= int(fields[9])

				# basic sanity checking:
				# if int(fields[10]) != 0: print "WARNING: Non-perfect overlap found"
				if v1_ovl_start == 0 and v1_ovl_end == v1_len-1: 
					print 'WARNING: Contained read found! Skipping.'
					continue
					# exit("ERROR: Contained read found")
				if v2_ovl_start == 0 and v2_ovl_end == v2_len-1: 
					print 'WARNING: Contained read found! Skipping.'
					continue
					# exit("ERROR: Contained read found")

				# do the reads actually ovelap?
				if v2_orientation == 0:
					assert (v1.seq[v1_ovl_start:v1_ovl_end+1] == v2.seq[v2_ovl_start:v2_ovl_end+1])
				elif v2_orientation == 1:
					assert v1.seq[v1_ovl_start:v1_ovl_end+1] == \
						reverse_complement(v2.seq[v2_ovl_start:v2_ovl_end+1])
					# keyboard()
				else:
					exit("ERROR: Invalid orientation")

				j = g.edge_id_generator.get_id()

				e = OverlapEdge(j, v1, v2, v1_ovl_start, v1_ovl_end, v1_len, v2_ovl_start, v2_ovl_end, v2_len, v2_orientation)

				g.add_edge(e)

				for v in (v1, v2):
					if e.connection[v] == 'H':
						v.head_edges.add(e)
					elif e.connection[v] == 'T':
						v.tail_edges.add(e)
					else:
						exit('ERROR: Invalid edge connection!')

			## HEADER
			elif fields[0] == ASQG_HEAD_REC:
				g.metadata['asqg_header'] = line.strip()

		for v in g.vertices:
			assert no_diedge(v)

	return g, vertices_by_contig

def _load_from_asqg(asqg_path):
	g = Graph()
	vertices_by_contig = dict()

	with open(asqg_path) as asqg:
		for line in asqg:
			fields = line.strip().split()

			## VERTEX
			if fields[0] == ASQG_VERTEX_REC:
				i = int(fields[1])
				v = OverlapVertex(i, fields[2])
				vertices_by_contig[fields[1]] = v
				g.add_vertex(v)

			## EDGE
			elif fields[0] == 'ED':
				v1 = vertices_by_contig[fields[1]]
				v2 = vertices_by_contig[fields[2]]
				v1_ovl_start 	= int(fields[3])
				v1_ovl_end 		= int(fields[4])
				v1_len 			= int(fields[5])
				v2_ovl_start 	= int(fields[6])
				v2_ovl_end 		= int(fields[7])
				v2_len 			= int(fields[8])
				v2_orientation 	= int(fields[9])

				# basic sanity checking:
				# if int(fields[10]) != 0: print "WARNING: Non-perfect overlap found"
				if v1_ovl_start == 0 and v1_ovl_end == v1_len-1: 
					print 'WARNING: Contained read found! Skipping.'
					continue
					# exit("ERROR: Contained read found")
				if v2_ovl_start == 0 and v2_ovl_end == v2_len-1: 
					print 'WARNING: Contained read found! Skipping.'
					continue
					# exit("ERROR: Contained read found")

				# do the reads actually ovelap?
				if v2_orientation == 0:
					assert (v1.seq[v1_ovl_start:v1_ovl_end+1] == v2.seq[v2_ovl_start:v2_ovl_end+1])
				elif v2_orientation == 1:
					assert v1.seq[v1_ovl_start:v1_ovl_end+1] == \
						reverse_complement(v2.seq[v2_ovl_start:v2_ovl_end+1])
					# keyboard()
				else:
					exit("ERROR: Invalid orientation")

				if v1_ovl_start == 0:
					v1_connection = 'H'
				elif v1_ovl_end == v1_len-1:
					v1_connection = 'T'
				else:
					exit("ERROR: Reads don't overlap at endpoints")

				if v2_ovl_start == 0:
					v2_connection = 'H'
				elif v2_ovl_end == v2_len-1:
					v2_connection = 'T'
				else:
					exit("ERROR: Reads don't overlap at endpoints")

				# # loops are useless:
				# if v1 == v2:
				# 	print "WARNING: Discarded a loop"
				# 	continue

				# some connections don't make sense

				# HH and TT overlaps only make sense if the reads have *opposite* orientation
				if v1_connection == v2_connection and v2_connection == 0:
					exit("WARNING: Nonsense overlap found")
					continue

				# HT and TH overlaps only make sense if the reads have *the same* orientation
				if v1_connection != v2_connection and v2_connection == 1:
					exit("WARNING: Nonsense overlap found")
					continue

				j = g.edge_id_generator.get_id()

				e = OverlapEdge(j, v1, v2, v1_ovl_start, v1_ovl_end, v1_len, v2_ovl_start, v2_ovl_end, v2_len, v2_orientation)

				g.add_edge(e)

				for v in (v1, v2):
					if e.connection[v] == 'H':
						v.head_edges.add(e)
					elif e.connection[v] == 'T':
						v.tail_edges.add(e)
					else:
						exit('ERROR: Invalid edge connection!')

			## HEADER
			elif fields[0] == ASQG_HEAD_REC:
				g.metadata['asqg_header'] = line.strip()

		for v in g.vertices:
			assert no_diedge(v)

	max_v_id = max(g.vertices_by_id.keys())
	# max_e_id = max(g.edges_by_id.keys())

	g.vertex_id_generator.set_counter(max_v_id+1)
	# g.edge_id_generator.set_counter(max_e_id+1)

	return g, vertices_by_contig

def _load_from_fasta(fasta_path):
	g = Graph()
	vertices_by_contig = dict()

	fasta = pysam.FastaFile(fasta_path)
	n = len(fasta.references)
	for i, ctg in enumerate(fasta.references):
		# if i % 1000 == 0: print '%d/%d' % (i, n)
		id_ = g.vertex_id_generator.get_id()
		seq = fasta.fetch(ctg)
		v = OverlapVertex(id_, seq)
		assert ctg not in vertices_by_contig
		vertices_by_contig[ctg] = v
		g.add_vertex(v)

	print 'Contigs loaded.'

	return g, vertices_by_contig

def _load_edges_from_tsv(g, tsv_path, vertices_by_contig=None):
	if not vertices_by_contig:
		vertices_by_contig = {v.id : v for v in g.vertices}

	with open(tsv_path) as tsv:
		for line in tsv:
			ctg1, ctg2, c1, c2, o, spt, d = line.strip().split()
			if int(spt) < 3:
				continue

			v1 = vertices_by_contig[ctg1]
			v2 = vertices_by_contig[ctg2]

			if c1 == TSV_LEFT:
				conn1 = 'H'
			elif c1 == TSV_RIGHT:
				conn1 = 'T'
			else:
				raise ValueError('Invalid connection value in .tsv')

			if c2 == TSV_LEFT:
				conn2 = 'H'
			elif c2 == TSV_RIGHT:
				conn2 = 'T'
			else:
				raise ValueError('Invalid connection value in .tsv')

			if o == TSV_SAME:
				ori = 0
			elif o == TSV_REVERSE:
				ori = 1
			else:
				raise ValueError('Invalid orientation value in .tsv')

			j = g.edge_id_generator.get_id()
			# FIXME: user proper distance
			e = ScaffoldEdge(j, v1, v2, conn1, conn2, ori, 25)
			# e = ScaffoldEdge(j, v1, v2, conn1, conn2, ori, int(d))
			e.support = int(spt)

			g.add_edge(e)

			for v in (v1, v2):
				if e.connection[v] == 'H':
					v.head_edges.add(e)
				elif e.connection[v] == 'T':
					v.tail_edges.add(e)
				else:
					exit('ERROR: Invalid edge connection!')

	print 'Edge connections loaded.'

	return g

def _load_containment(g, containment_file, vertices_by_contig=None):
	if not vertices_by_contig:
		vertices_by_contig = {v.id : v for v in g.vertices}

	with open(containment_file) as in_:
		for line in in_:
			fields = line.split()

			v = vertices_by_contig.get(fields[1], None)
			if not v:
				print 'WARNING: Vertex not found:', fields[1]
				continue
			
			if fields[0] == CTMT_WELL_REC:
				well, start, end = int(fields[2]), int(fields[3]), int(fields[4])
				v.add_well(well, start, end)

			elif fields[0] == CTMT_IVL_REC:
				ivl = (int(fields[2]), int(fields[3]), int(fields[4]))
				v.add_interval(ivl)

			else:
				print 'WARNING: Invalid record type found:', fields[0]

# ----------------------------------------------------------------------------
# saving the graph

def _write_asqg(g, asqg_file):
	with open(asqg_file, 'w') as asqg:
		for v in g.vertices:
			asqg.write('VT\t%d\t%s\n' % (v.id, v.seq))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			asqg.write('ED\t{v1}\t{v2}\t{v1os}\t{v1oe}\t{v1l}\t{v2os}\t{v2oe}\t{v2l}\t{ori}\n'.format
				(v1=v1.id, v2=v2.id, 
				 v1os=e.ovl_start[v1], v1oe=e.ovl_end[v1], v1l=len(v1), 
				 v2os=e.ovl_start[v2], v2oe=e.ovl_end[v2], v2l=len(v2), ori=e.v2_orientation))

def _write_containment(g, containment_file):
	with open(containment_file, 'w') as out:
		for v in g.vertices:
			for ctg in v.metadata['contigs']:
				out.write('%d\t%s\t%d\t%d\n' % (v.id, ctg, v.metadata['contig_starts'].get(ctg,-1), 
																				v.metadata['contig_ends'].get(ctg,-1)))