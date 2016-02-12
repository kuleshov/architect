ARCH:=/path/to/architect/architect.py
ARCHDIR:=/path/to/architect
FASTA:=/path/to/scaffolds.fasta
TSV:=/path/to/besst.edges.tsv
CTM:=/path/to/drosophila.containment


wells: $(CTM)
	$(ARCH) scaffold \
	  --fasta $(FASTA) \
	  --edges $(TSV) \
	  --containment $< \
	  --out drosophila

view: $(CTM)
	$(ARCH) view \
	  --fasta $(FASTA) \
	  --edges $(TSV) \
	  --containment $<

check-layout: drosophila.layout
	python $(ARCHDIR)/scripts/eval-layout.py -l $<

clean:
	-rm *.fasta *.containment *.fasta.fai *.layout *.tsv *.gfa
