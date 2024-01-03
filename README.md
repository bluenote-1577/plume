# plume - pseudoaligner for metagenomics

### WARNING - ALPHA MODE

This software is in alpha. Perhaps even pre-alpha. It works, but it is poorly documented and subject to change. 

### What does plume do? 

Plume is a fast, lightweight pseudoaligner for **calculating coverages/abundances** for metagenomic reads against **genomes/contigs**. 

The inputs/outputs of plume are similar to [CoverM](https://github.com/wwood/CoverM), but plume is a new algorithm written from scratch. 

### Why plume?

Plume is similar to [kallisto](https://github.com/pachterlab/kallisto) or [salmon](https://github.com/COMBINE-lab/salmon), but it offers a few engineering and methodological innovations that make it suitable for metagenomics:

1. Faster all-to-all coverage calculation 
2. Long-read aware pseudoalignment 
3. Optimized parameters for metagenomics

#### Performance

* CoverM (with minimap2) took about **10 minutes (with 4 threads)** to quantify a metagenome against 180 genomes.
* Plume took **2 minutes to index, but 20 seconds (with 4 threads)** to quantify the same metagenome.

I estimate that plume ranges from 1-5x faster (one-to-one coverage calculation) to 20-30x faster (all-to-all) than minimap2 for coverage calculation; it will depend on how much mapping vs indexing you need. 

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
plume quant contigs*.plref reads*.plreads -t 20 > all-to-all.tsv

# relative (sequence) abundance
plume quant -r *.plref *.plreads > abund.tsv
```

## Output format

The output is a TSV file that looks like the following:

```
Genome_or_contigs	sample1.fq.plreads  smaple2.fq.plreads
temp/GCA_900542375.1_genomic.fna	0.0486	0.0182
temp/GCF_003865035.1_genomic.fna	1.4225	0.0167
temp/GCA_000020605.1_genomic.fna	0.3205	0.0111
temp/GCF_000313565.1_genomic.fna	0.0153	0.0258
```

## Citation

Forthcoming (eventually) 
