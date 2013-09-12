import os

from intervals import get_intervals
from libkuleshov.dna import reverse_complement

def print_vertex(x):
	print '================================='
	print 'Vertex: %d' % v.id_
	print [(e_tail.id_, e_head.id_) for (e_tail, e_head) in pairs_to_resolve]
	print ','.join([str(w) for w in v_wells])
	print_true_intervals(v, v.metadata['contigs'])
	print
	print 'BEFORE:'
	print 'H:'
	for e in H:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in get_tail_wells(w) if well in T_wells])
		print e.id_, '\t', wells
		print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_ends'][ctg] > len(w) - 1000])
	print 'T:'
	for e in T:
		w = e.other_vertex(v)
		wells = ','.join([str(well) for well in get_head_wells(w) if well in H_wells])
		print e.id_, '\t', wells
		print_true_intervals(w, [ctg for ctg in w.metadata['contigs'] if w.metadata['contig_starts'][ctg] < 1000])
		# for ctg in w.metadata['contigs'][:15]:
		# 	if w.metadata['contig_ends'][ctg] < 1000:
		# 		print_true_coordinates(w, ctg)
	print

###############################################################################
## DETERMINE THE TRUE REFERENCE SEQUENCES THAT ARE CONTAINED IN A CONTIG

def get_true_intervals(w, ctgs):
	I = list()
	for ctg in ctgs:
		fields = ctg.split('_')
		chrom, coords = fields[1].split(':')
		start, end = (int(x) for x in coords.split('-'))

		# correction for how we generated the reads. bed format is 1-based
		start -= 1
		end -= 1

		internal_start = int(fields[2])

		I.append((int(chrom), start + internal_start, start + internal_start + 199))

	return get_intervals(I)

def print_true_coordinates(w, ctg):
	fields = ctg.split('_')
	chrom, coords = fields[1].split(':')
	start, end = (int(x) for x in coords.split('-'))

	# correction for how we generated the reads. bed format is 1-based
	start -= 1
	end -= 1

	internal_start = int(fields[2])

	print chrom, start + internal_start, start + internal_start + 199

def print_true_intervals(w, ctgs):
	I = list()
	for ctg in ctgs:
		fields = ctg.split('_')
		chrom, coords = fields[1].split(':')
		start, end = (int(x) for x in coords.split('-'))
		
		# correction for how we generated the reads. bed format is 1-based
		start -= 1
		end -= 1

		internal_start = int(fields[2])
		I.append((int(chrom), start + internal_start, start + internal_start + 199))

	merged_I = get_intervals(I)
	for i in merged_I:
		print '\t', i

###############################################################################
## EXAMINE POTENTIAL MISSASEMBLIES

def load_genome():
	GENOME_PATH="/home/kuleshov/metagenomica/simulate/rs/simulate_short_reads/renamed.ref"
	FASTX_PATH="/home/kuleshov/lib/fastx"

	genome = dict()
	with os.popen('cat %s | %s/fasta_formatter' % (GENOME_PATH, FASTX_PATH)) as genome_file:
		current_ctg = ""
		for line in genome_file:
			line.strip()
			if line.startswith('>'):
				current_ctg = int(line[1:])
			else:
				assert current_ctg and current_ctg not in genome
				genome[current_ctg] = line

	return genome

def common_overlap(i1, i2):
	""" Find the common overlap of i1, i2. """
	# NOTE: Doesn't work b/c we don't consider rev. compl.

	# i2 -> i1?
	for i1_start_in_i2 in xrange(len(i2)):
		# l = len(i2) - i1_start_in_i2
		l = len(i2[i1_start_in_i2:])
		if i2[i1_start_in_i2:] == i1[:l]:
			return i1[:l]

	# i1 -> i2?
	for i2_start_in_i1 in xrange(len(i1)):
		# l = len(i1) - i2_start_in_i1
		l = len(i1[i2_start_in_i1:])
		if i1[i2_start_in_i1:] == i2[:l]:
			return i2[:l]

def get_containments(v, true_intervals, genome):
	"""
	Returns a list of that descripbes how true_intervals are contained 
	in the assembled sequence of v.
	A containment is a tuple (start_pos, end_pos, orientation).
	"""
	containments = list()
	assembly = v.seq
	for i in true_intervals:
		# print i, i[2]-i[1]+1,
		true_region = genome[i[0]][i[1]:i[2]+1]

		# try same strand:
		for j in xrange(len(assembly)):
			if true_region == assembly[j:j+len(true_region)]:
				containments.append((j, j+len(true_region)-1, 0))
				break

		# try opposite strands:
		true_region_opp = reverse_complement(true_region)
		for j in xrange(len(assembly)):
			if true_region_opp == assembly[j:j+len(true_region_opp)]:
				containments.append((j, j+len(true_region_opp)-1, 1))
				break

		# containments.append((-1, -1, -1))

	return containments

def print_containment(v, true_intervals, genome):
	# if len(v.seq) > 250: return
	print '-----------------------------------------------------------'
	print 'ASSEMBLY:', v.id_, len(v.seq)
	# print v.seq
	print 'VERIFYING INTERVALS:', ', '.join([str(i) for i in true_intervals])

	containments = get_containments(v, true_intervals, genome)
	for c, intv in zip(containments, true_intervals):
		print intv, intv[2]-intv[1]+1, c

def examine_repeats(intervals, genome):
	print '-----------------------------------------------------------'
	print 'VERIFYING INTERVALS:', ', '.join([str(i) for i in intervals])
	for i in intervals:
		print i, i[2] - i[1] + 1

		# if i[2] - i[1] < 100:
		# 	print genome[i[0]][i[1]:i[2]+1]
		# else:
		# 	print

	if len(intervals) == 2:
		a = intervals[0]
		b = intervals[1]
		# need to consider rev. compl. in this function:
		o = common_overlap(genome[a[0]][a[1]:a[2]+1], genome[b[0]][b[1]:b[2]+1])
		if o:
			print o
		else:
			print "overlap not found"

def is_misassembly(v, true_intervals, genome):
	"""
	If a true genomic region does not touch the 200bp-long edges of a contig,
	we consider that it is a misassembly.
	"""
	L = len(v.seq)
	containments = get_containments(v, true_intervals, genome)
	for c, intv in zip(containments, true_intervals):
		if c[0] > 200 or c[1] < L - 200:
			return True

	return False

###############################################################################
## VISUALIZE SUB-COMPONENTS OF THE GRAPH

def visualize_assembly(v, I, genome):
	print v.seq
	for i, ctg in enumerate(v.metadata['contigs']):
		print ctg, v.metadata['contig_starts'][ctg], v.metadata['contig_ends'][ctg]
		print genome[I[i][0]][I[i][1]:I[i][2]+1]

def visualize_connection(e):
	v1, v2 = e.v1, e.v2
	print e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation
	print v1.seq
	print v2.seq
	print v1.seq[e.ovl_start[v1]:e.ovl_end[v1]+1]
	print v2.seq[e.ovl_start[v2]:e.ovl_end[v2]+1]

###############################################################################
## SAVE TO DOT

def to_abyss_explorer_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		id_to_num = dict()
		for v in g.vertices:
			dot.write('"%d+" [l=%d C=1]\n' % (v.id_, len(v.seq)))
			dot.write('"%d-" [l=%d C=1]\n' % (v.id_, len(v.seq)))

		for e in g.edges:
			v1, v2, ori = e.v1, e.v2, e.v2_orientation
			if ori == 0:
				dot.write('"%d+" -> "%d+"\n' % (v1.id_, v2.id_))
			elif ori == 1:
				dot.write('"%d+" -> "%d-"\n' % (v1.id_, v2.id_))
			else:
				exit("ERROR: Invalid orientation")

		dot.write('}\n')

def to_graphviz_dot(g, dot_file):
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			dot.write('%d [label = "%d"]\n' % (v.id_, len(v.seq)))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id_, v2.id_, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')

def to_graphviz_dot_with_intervals(g, dot_file):
	genome = load_genome()
	with open(dot_file, 'w') as dot:
		dot.write('digraph adj {\n')
		for v in g.vertices:
			I = get_true_intervals(v, v.metadata['contigs'])
			# if v.id_== 3599553:
			# 	visualize_assembly(v, I, genome)
			if len(I) > 1:
				color = "red"
				# examine_repeats(I, genome)
				# if is_misassembly(v, I, genome):
					# print_containment(v, I, genome)
			else:
				color = "black"

			dot.write('%d [label = "%s", color="%s"]\n' % (v.id_, ','.join([str(i) for i in I]), color))
		for e in g.edges:
			v1, v2 = e.v1, e.v2
			# if v1.id_ == 3706567 and v2.id_ == 3712515:
			# 	visualize_connection(e)
			dot.write('"%d" -> "%d" [label= "%d %d %d %d %s%s %d"]\n' % (v1.id_, v2.id_, e.ovl_start[v1], e.ovl_end[v1], e.ovl_start[v2], e.ovl_end[v2], e.connection[v1], e.connection[v2], e.v2_orientation))

		dot.write('}\n')