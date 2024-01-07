# plume - pseudoaligner for metagenomics

### WARNING - ALPHA MODE

This software is in alpha. Perhaps even pre-alpha. It works, but it is poorly documented and subject to change. 

### What does plume do? 

Plume is a fast, lightweight pseudoaligner for **calculating coverages/abundances** for metagenomic reads against **genomes/contigs**. Plume is similar to [kallisto](https://github.com/pachterlab/kallisto) or [salmon](https://github.com/COMBINE-lab/salmon) in spirit, but instead of quantifying transcripts from RNA-seq, it quantifies/calculates coverage for metagenomics. 

The inputs/outputs of plume are similar to [CoverM](https://github.com/wwood/CoverM), but plume is a new algorithm written from scratch. 

#### Performance - speed

* CoverM (with minimap2) took about **10 minutes (with 4 threads)** to quantify a metagenome against 180 genomes.
* Plume took **2 minutes to index, but 20 seconds (with 4 threads)** to quantify the same metagenome.

I estimate that plume ranges from 1-5x faster (one-to-one coverage calculation) to 10x faster (all-to-all) than minimap2 for coverage calculation; it will depend on how much mapping vs indexing you need. The performance will also change as I optimize/tweak.

#### Performance - accuracy

The coverage is pretty concordant with CoverM (using minimap2) on genomes/contigs. I quickly tested plume for metagenomic binning coverage; the results were worse than BWA but passable. 

## Quick start

### Requirements
1. Make sure [rust](https://www.rust-lang.org/tools/install) is installed.
2. Standard toolchain (make, GCC, and make)

```sh
git clone https://github.com/bluenote-1577/plume
cd plume
cargo install --path .
plume -h
```

### Indexing

Plume requires contigs/genomes/reads to be indexed before usage. 

```sh
# indexing contigs inside a genome/assembly. outputs 2 indices
plume index -c contigs1.fa contigs2.fa

# indexing a collection of genomes. outputs 1 index
plume index -g genome1.fa genome2.fa

# indexing paired-end reads
plume index -1 A_1.fq B_1.fq -2 A_2.fq B_2.fq

# or single-end reads
plume index -s reads1.fq reads2.fq
```

### Quantification/Coverage calculation

```sh
# coverage
plume quant contigs.fa.plref A_1.fq.paired.plreads > output.tsv

# relative (sequence) abundance
plume quant -r *.plref *.plreads > abund.tsv
```

## Output format

The output is a TSV file that looks like the following:

```
Reference_file	Genome_or_contigs	sample1.fastq.gz.plreads	sample2.fastq.gz.plreads
temp.fa.plref	GG662017.1	0.0000	0.0000
temp.fa.plref	GG662016.1	0.0000	0.3170
temp.fa.plref	GG662015.1	0.0000	0.0000
```
