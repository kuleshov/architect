Architect
=========

**NOTE**: We will release the full version of Architect after ISMB. We suggest waiting until then before using the program. Until then, this repo contains a demo that can reproduce results from the paper.

***

Architect is a genomic scaffolder aimed at synthetic long read and read-cloud sequencing technologies
such as Illumina Tru-Seq synthetic long reads or the 10X GemCode platform.

## Requirements

Architect is implemented in Python and requires

* `pysam >= 0.82`
* `networkx >= 1.10`

## Input data

Architect takes as input:
* Genomic contigs in `fasta` format assembled using a standard (short-read) assembler.
* A mapping of read clouds to contigs in `bam` format.
* Optionally, an alignment of paired-end reads to the contigs.

## Usage

Architect is run as follows.

```
usage: architect.py scaffold [-h] --fasta FASTA --edges EDGES --containment
                             CONTAINMENT --out OUT [--min-ctg-len MIN_CTG_LEN]
                             [--log LOG]

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         Input scaffolds/contigs
  --edges EDGES         Known paired-end or overlap connections
  --containment CONTAINMENT
                        Container hits and various meta-data
  --out OUT             Prefix for the ouput files
  --min-ctg-len MIN_CTG_LEN
                        Discard contigs smaller than this length (def: 0)
  --log LOG             Save stdout to log file
```

A `containment` file encodes container hits in the genome. The `edges` file encodes paired-end read information. The `fasta` file contains the pre-assembled contigs.

The input files are generated from `bam` aligments using their corresponding scripts in the `/bam` folder.

## Results

We used Architect to assemble the genomes of *D. melanogaster* and *C. elegans* as well as two gut metagenomic datasets.
Architect took as input standard short read assemblies augmented with raw short reads cloud based on the Illumina TSLRs technology that were subsampled to various depths.
We found that the scaffolder produced up to 5x improvements in contig contiguity without increasing the misassembly rate, and using between 4-20x less sequencing data.

### Scaffolding the genome of D. melanogaster

Let's now see how well Architect can assemble the genome of D. melanogaster.
We will take as input a set of contigs assembled from regular short reads 
using SPAdes. In addition, we will use as input pre-computed links
between contigs that were derived from an alignment of paired-end reads,
as well as "container hits" obtained from aligning a TLSR read cloud library, 
in which the short reads have been subsampled to 25%.

More concretely, we will run Architect on 3 files:

* `scaffolds.fasta`: scaffolds obtained via SPAdes from short reads
* `edges.tsv`: paired-end links obtained from the same reads
* `drosophila.containment`: container hits obtained from aligning a subsampled TSLR read cloud library

These files can be downloaded from `https://stanford.box.com/s/10crf36gvf998ylhjo59pr7xswwb0uga`.

### Running Architect

We may run the scaffolder on this data as follows:

```
python architect.py scaffold \
    --fasta scaffolds.fasta \
    --edges edges.tsv \
    --containment drosophila.containment \
    --out drosophila
```

This will produce two files:

* `drosophila.layout`: which will contain the scaffold orderings.
* `drosophila.fasta`: scaffolds determined by the ordering; the number of N's between contigs is unreliable (may be off by +/-2kbp) and may confuse aligners

The format of the layout file is:
```
<ordering id>   <ctg_id1;containment_string1;orientation1>    <ctg_id2;containment_string2;orientation2>
```
A containment string, if available, indicates the regions of the (known) reference genome to which a contig maps.
It can be used for validation; its format is `ctg:start-end`.
Finally, `orientation` is either `R` or `S`, indicating whether the contig came from the reverse or the forward strand, respectively.

### Verifying the accuracy of the method

We can use the information in the `layout` file to verify our accuracy.
This is done using an evaluation script:
```
python scripts/eval-layout.py -l drosophila.layout
```

This will print statistics for each contig. The most important line
is the last one; it prints out five numbers:
```
<num_orderings> <num_misorderings> <bp_verifiable> <n50> <na50>
```
Here, the misorderings are the equivalent of a misassembly (in the QUAST sense).
The NA50 is also an extension of the QUAST definition to orderings.

In our case, this result should read:
```
57567 40 122423762 262797 252190
```

The alignment information in `drosophila.layout` is
taken from `drosophila.containment`. This file contains a mapping of
most contigs to the known Drosophila reference. This information
may be used for internald diagnostics, as well as for evaluating the final layout.

### Scripts for generating the input files

The input files `edges.tsv` and `drosophila.containment` are generated 
via scripts in the `/scripts` subfolder.

The paired-end link file `edges-tsv` is generated using `pe-connections.py`
from an alignment of paired-end reads to the input `scaffolds.fasta`.
We used the same paired-end reads as for de-novo assembly to generate 
`edges.tsv`.

The file `drosophila.containment` is actually a concatenation of two files:

* A map of container hits; these correspond to entries with `W` in the first column.
* An alignment of contigs to the reference; these are lines with `R` in the first column.

The container hits are obtained using the script `bam_to_containment.py`
from an alignment of TSLR read clouds to `scaffolds.fasta`.

The true intervals are converted directly from the true alignments
identified by the QUAST tool. The scripts for generating them can
be found in the `scripts/ground-truth` subfolder.
