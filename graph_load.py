from graph import Graph
from string_graph import OverlapVertex, OverlapEdge, no_diedge
from libkuleshov.dna import reverse_complement
from libkuleshov.debug import keyboard

def load_from_sga_asqg(asqg_path):
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
			if fields[0] == 'VT':
				i = g.vertex_id_generator.get_id()
				v = OverlapVertex(i, fields[2])
				v.metadata['contigs'] = [fields[1]]
				v.metadata['contig_starts'] = {fields[1]: 0}
				v.metadata['contig_ends'] = {fields[1]: len(v)-1}
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

				# if v1_ovl_start == 0:
				# 	v1_connection = 'H'
				# elif v1_ovl_end == v1_len-1:
				# 	v1_connection = 'T'
				# else:
				# 	exit("ERROR: Reads don't overlap at endpoints")

				# if v2_ovl_start == 0:
				# 	v2_connection = 'H'
				# elif v2_ovl_end == v2_len-1:
				# 	v2_connection = 'T'
				# else:
				# 	exit("ERROR: Reads don't overlap at endpoints")

				# some connections don't make sense

				# HH and TT overlaps only make sense if the reads have *opposite* orientation
				# if v1_connection == v2_connection and v2_connection == 0:
				# 	exit("WARNING: Nonsense overlap found")
				# 	continue

				# # HT and TH overlaps only make sense if the reads have *the same* orientation
				# if v1_connection != v2_connection and v2_connection == 1:
				# 	exit("WARNING: Nonsense overlap found")
				# 	continue

				j = g.edge_id_generator.get_id()

				# change H->T connections so they're T->H (i.e. left to right)
				# if v2_connection == 'T' and v1_connection == 'H':
				# 	# tail(V2) -> head(V1)
				# 	# so make direction be v2 -> v1
				# 	e = OverlapEdge(j, v2, v1, v2_ovl_start, v2_ovl_end, v2_len, v1_ovl_start, v1_ovl_end, v1_len, v2_orientation)
				# elif connection == 'LR'
				# 	# don't change orientation
				# 	e = OverlapEdge(j, v1, v2, v1_ovl_start, v1_ovl_end, v1_len, v2_ovl_start, v2_ovl_end, v2_len, v2_orientation)
				# else:
				# 	exit()

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
			elif fields[0] == 'HT':
				g.metadata['asqg_header'] = line.strip()

		for v in g.vertices:
			assert no_diedge(v)

	return g

def load_graph(asqg_path, containment_file):
	g = load_asqg(asqg_path)
	load_containment(g, containment_file)
	return g

def load_asqg(asqg_path):
	g = Graph()
	vertices_by_contig = dict()

	with open(asqg_path) as asqg:
		for line in asqg:
			fields = line.strip().split()

			## VERTEX
			if fields[0] == 'VT':
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
			elif fields[0] == 'HT':
				g.metadata['asqg_header'] = line.strip()

		for v in g.vertices:
			assert no_diedge(v)

	max_v_id = max(g.vertices_by_id.keys())
	# max_e_id = max(g.edges_by_id.keys())

	g.vertex_id_generator.set_counter(max_v_id+1)
	# g.edge_id_generator.set_counter(max_e_id+1)

	return g

def load_containment(g, containment_file):
	with open(containment_file) as in_:
		for line in in_:
			fields = line.split()
			v = g.vertices_by_id[int(fields[0])]
			ctg = fields[1]

			if 'contigs' not in v.metadata:
				v.metadata['contigs'] = list()
				v.metadata['contig_starts'] = dict()
				v.metadata['contig_ends'] = dict()
			v.metadata['contigs'].append(ctg)

			if fields[2] != '-1':
				v.metadata['contig_starts'][ctg] = int(fields[2])
				v.metadata['contig_ends'][ctg] = int(fields[3])

def save_graph(g, asqg_file, containment_file):
	write_asqg(g, asqg_file)
	write_containment(g, containment_file)

def write_asqg(g, asqg_file):
	with open(asqg_file, 'w') as asqg:
		for v in g.vertices:
			asqg.write('VT\t%d\t%s\n' % (v.id_, v.seq))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			asqg.write('ED\t{v1}\t{v2}\t{v1os}\t{v1oe}\t{v1l}\t{v2os}\t{v2oe}\t{v2l}\t{ori}\n'.format
				(v1=v1.id_, v2=v2.id_, 
				 v1os=e.ovl_start[v1], v1oe=e.ovl_end[v1], v1l=len(v1), 
				 v2os=e.ovl_start[v2], v2oe=e.ovl_end[v2], v2l=len(v2), ori=e.v2_orientation))

def write_containment(g, containment_file):
	with open(containment_file, 'w') as out:
		for v in g.vertices:
			for ctg in v.metadata['contigs']:
				out.write('%d\t%s\t%d\t%d\n' % (v.id_, ctg, v.metadata['contig_starts'].get(ctg,-1), v.metadata['contig_ends'].get(ctg,-1)))

def to_stabq(g, stabq_file):
	with open(stabq_file) as stabq:
		for v in vertices:
			stabq.write('>%s\n' % str(v.id_))
			stabq.write('%s\n' % v.seq)