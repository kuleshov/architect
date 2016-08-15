Architect
=========

Architect is a genomic scaffolder aimed at synthetic long read and read-cloud sequencing technologies
such as Illumina Tru-Seq or the 10X platform. 

We describe Architect in detail in the following paper.

```
V. Kuleshov et al., Genome assembly from synthetic long read clouds, Bioinformatics (2016) 32 (12): i216-i224.
```

## Requirements

Architect is implemented in Python and requires

* `pysam >= 0.82`
* `networkx >= 1.10`

## Installation

To install Architect, clone this repo and add it to your `PYTHONPATH`.

```
git clone https://github.com/kuleshov/architect.git
cd architect
export PYTHONPATH=$PYTHONPATH:`pwd`
```

You may now run the program as `python /path/to/architect.py`.

## Input data

Architect takes as input:
* Genomic contigs in `fasta` format assembled using a standard (short-read) assembler.
* A mapping of read clouds to contigs in `bam` format.
* Optionally, an alignment of paired-end reads to the contigs.

## Usage

Architect is run as follows.

```
usage: architect.py scaffold [-h] --fasta FASTA [--edges EDGES] --containment
                             CONTAINMENT --out OUT [--min-ctg-len MIN_CTG_LEN]
                             [--cut-tip-len CUT_TIP_LEN]
                             [--pe-abs-thr PE_ABS_THR]
                             [--pe-rel-thr PE_REL_THR]
                             [--pe-rc-rel-thr PE_RC_REL_THR]
                             [--rc-abs-thr RC_ABS_THR]
                             [--rc-rel-edge-thr RC_REL_EDGE_THR]
                             [--rc-rel-prun-thr RC_REL_PRUN_THR] [--log LOG]

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         Input scaffolds/contigs
  --edges EDGES         Known paired-end or read cloud connections
  --containment CONTAINMENT
                        Container hits and various meta-data
  --out OUT             Prefix for the ouput files
  --min-ctg-len MIN_CTG_LEN
                        Discard contigs smaller than this length (def: 0)
  --cut-tip-len CUT_TIP_LEN
                        Cut tips smaller than this length
  --pe-abs-thr PE_ABS_THR
                        Threshold for absolute support when pruning paired-end
                        edges
  --pe-rel-thr PE_REL_THR
                        Threshold for relative support when pruning paired-end
                        edges
  --pe-rc-rel-thr PE_RC_REL_THR
                        Threshold for relative support for read-cloud /
                        paired-end pruning
  --rc-abs-thr RC_ABS_THR
                        Minimum support for create read-cloud based edge
  --rc-rel-edge-thr RC_REL_EDGE_THR
                        Threshold for relative support when creating read-
                        cloud based edges
  --rc-rel-prun-thr RC_REL_PRUN_THR
                        Threshold for relative support when pruning read-cloud
                        based edges
  --log LOG             Save stdout to log file
```

A `containment` file encodes container hits in the genome. The `edges` file encodes paired-end read information. The `fasta` file contains the pre-assembled contigs.

The input files are generated from `bam` alignments using scripts in the `/bam` folder.

The user may tune Architect via various parameters described above. At the moment, we have set their defaults to values that were found to work well with Illumina TruSeq Synthetic Long Read datasets.

## Results

We used Architect to assemble the genomes of *D. melanogaster* and *C. elegans* as well as two gut metagenomic datasets.
Architect took as input standard short read assemblies augmented with raw short reads cloud based on the Illumina TSLRs technology that were subsampled to various depths.
We found that the scaffolder produced up to 5x improvements in contig contiguity without increasing the misassembly rate, and using between 4-20x less sequencing data.

## Documentation

For more information on running Architect, have a look at the [wiki](https://github.com/kuleshov/architect/wiki). You may find there:

* General [usage information](https://github.com/kuleshov/architect/wiki/Usage-Information)
* An [example](https://github.com/kuleshov/architect/wiki/D.-melanogaster-Example) of how to use Architect to scaffold the genome of D. melanogaster
* Information on the [file formats](https://github.com/kuleshov/architect/wiki/File-Format-Information) used by Architect