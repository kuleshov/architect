COORDS:=/path/to/quast_results/latest/contigs_reports/nucmer_output/scaffolds.coords.filtered
BEDTOOLS:=/path/to/bedtools
FASTA:=/path/to/Drosophila_melanogaster.BDGP6.dna.toplevel.fa

scaffolds.bed: $(COORDS)
	python ./coords2bed.py -c $< -b $@.tmp
	mv $@.tmp $@

scaffolds.sorted.bed: scaffolds.bed
	$(BEDTOOLS) sort -i $< > $@.tmp
	mv $@.tmp $@

scaffolds.merged.d1.bed: scaffolds.sorted.bed
	$(BEDTOOLS) merge -i $< -d 1 > $@.tmp
	mv $@.tmp $@

scaffolds.merged.d200.bed: scaffolds.sorted.bed
	$(BEDTOOLS) merge -i $< -d 200 > $@.tmp
	mv $@.tmp $@

scaffolds.csv: scaffolds.bed
	python ./bed2csv.py -b $< -c $@.tmp
	mv $@.tmp $@

scaffolds.containment: scaffolds.bed
	python ./bed2containment.py -b $< -c $@.tmp -r $(FASTA)
	mv $@.tmp $@
