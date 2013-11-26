import os

from intervals import get_true_intervals
from libkuleshov.dna import reverse_complement

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

def examine_misassemblies(g):
	genome = load_genome()
	for v in g.vertices:
		I = get_true_intervals(v, v.metadata['contigs'])
		# if v.id_== 3599553:
		# 	visualize_assembly(v, I, genome)
		if len(I) > 1:
			# examine_repeats(I, genome)
			if is_misassembly(v, I, genome):
				print_containment(v, I, genome)