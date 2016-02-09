Architect
=========

Architect is a genomic scaffolder aimed at synthetic long read cloud sequencing technologies
such as Illumina Tru-Seq synthetic long reads or the 10X GemCode platform.

## Requirements

Architect is implement in Python and requires

* `pysam >= 0.82`
* `networkx >= 1.10`

## Input data

Architect takes as input:
* Genomic contigs in `fasta` format assembled using a standard (short-read) assembler.
* A mapping of read clouds to contigs in `bam` format.
* Optionally, an aligment of paired-end reads to the contigs.

## Usage

Architect is run as follows.

```
usage: architect.py scaffold [-h] --fasta FASTA --edges EDGES
                             [--containment CONTAINMENT] [--log LOG]

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA
  --edges EDGES
  --containment CONTAINMENT
  --log LOG
```

A `containment` file encodes container hits in the genome. The `edges` file encodes paired-end read information. The `fasta` file contains the pre-assembled contigs.

The input files are generated from `bam` aligments using their corresponding scripts in the `/bam` folder.

## Results

We used Architect to assemble the genomes of *D. melanogaster* and *C. elegans* as well as two gut metgenomic datasets.
Architect took as input standard short read assemblies augmented with raw short reads cloud based on the Illumina TSLRs technology that were subsampled to various depths.
We found that the scaffolder produced up to 5x improvements in contig contiguity without increasing the misassembly rate, and using between 4-20x less sequencing data.