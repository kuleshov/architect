def get_intervals(I):
	merged_I = list()
	sorted_I = sorted(I)
	current_i = sorted_I[0]

	for i in sorted_I[1:]:
		if overlaps_sorted(current_i, i):
			current_i = merge(current_i, i)
		else:
			merged_I.append(current_i)
			current_i = i

	merged_I.append(current_i)

	return merged_I

def overlaps_sorted(a, b):
	# assumes a is before b
	return a[0] == b[0] and b[1] <= a[2]

def merge(a, b):
	return (a[0], a[1], b[2])