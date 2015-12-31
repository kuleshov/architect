BAM=/home/kuleshov/scaffolder/drosophila/align/sub0.25-to-spades-shotgun/reads.rg.bam
BAMPARSE=./dro.s025.E500.o5000.bamParse
LINKFILE=./dro.s025.r1.fragScaff.r1.links.txt
FASTA=/home/kuleshov/scaffolder/drosophila/asm/spades-shotgun/out/scaffolds.fasta

THREADS=5

for j in 1.00 1.25 1.50 1.75 2 2.5 3; do
  echo $j
  perl ~/bin/asm/fragScaff.pl \
    -B $BAMPARSE \
    -K $LINKFILE \
    -O dro.s025.j$j \
    -j $j \
    -F $FASTA \
    -P F \
    -G R \
    -E 500 \
    -b 0 \
    -t $THREADS \
    -v
    # -C 50
    #-A \
done
