use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "Plume is a lightweight pseudoaligner for metagenomics, allowing for coverage/abundance calculation of reads against genomes or contigs.", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Index reads and genomes/contigs
    Index(IndexArgs),
    /// Calculate coverage/relative abundance
    Quant(QuantArgs),
}


#[derive(Args, Default)]
pub struct IndexArgs {
    //#[clap(multiple=true, help = "")]
    //pub files: Vec<String>,
        //#[clap(short,long, default_value = "", help_heading = "OUTPUT", help = "If specified, output indices will have this string as an attached prefix.")]
    //pub prefix: String,
    
    //INPUTS----------------------
    #[clap(short,long="contigs", multiple=true, help_heading = "INPUT", help = "Contigs to index (fasta/gzip/bz2). Each group of contigs is indexed.")]
    pub contigs: Vec<String>,
    #[clap(short,long="genomes", help_heading = "INPUT", multiple=true, help = "Treat each file (fasta/gzip/bz2) as a genome and combine into a single index. File will be named <--genome-out>.plref")]
    pub genomes: Vec<String>,
    #[clap(short='G',long="genome-out", default_value = "indexed_genomes", help_heading = "INPUT", help = "Index name for --genomes option")]
    pub genome_out: String,
    #[clap(short='1',long="first-pair", multiple=true, help_heading = "INPUT", help = "First pairs in paired end reads e.g. S1_1.fq S2_1.fq")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pair", multiple=true, help_heading = "INPUT", help = "Second pairs in paried end reads e.g. S1_2.fq S2_2.fq")]
    pub second_pair: Vec<String>,
    #[clap(short='s',long="single-end", multiple=true, help_heading = "INPUT", help = "Single-end reads (e.g. long-reads)")]
    pub single_reads: Vec<String>,

    //ALGORITHM----------------------
    #[clap(short, default_value_t = 19,help_heading = "ALGORITHM", help ="k-mer size")]
    pub k: usize,
    #[clap(short, default_value_t = 15, help_heading = "ALGORITHM", help = "Window size for minimizers. Density is 2/(w+1)")]
    pub w: usize,
    #[clap(long, hidden=true, help_heading="ALGORITHM", help = "Use FracMinHash seeds")]
    pub fmh: bool,
    #[clap(long, help_heading="ALGORITHM", help = "Use open syncmer seeds (default)")]
    pub window_sync: bool,
    #[clap(long, help_heading="ALGORITHM", help = "Use minimizer seeds")]
    pub minimizer: bool,

    
    //MISC
    #[clap(short,long, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace-level output for debugging")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug-level output for debugging")]
    pub debug: bool,
    //#[clap(short,long="list", multiple=true, help_heading = "INPUT", help = "Input files in a list instead of in the command line.")]
    //pub list: Option<String>,

}

#[derive(Args)]
pub struct QuantArgs {
    #[clap(multiple=true, help = "Pre-sketched query or sample files. Raw fastq/fasta also allowed but presketching is recommended; see sylph sketch for more info")]
    pub files: Vec<String>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,

    #[clap(long="debug", help = "Debug output for debugging")]
    pub debug: bool,

    #[clap(long="trace", help = "Trace output for debugging")]
    pub trace: bool,
    #[clap(short,long="relative-abundance", help_heading = "OUTPUT", help = "Output relative abundances instead of coverages")]
    pub relative_abundance: bool,

    #[clap(long="no-supp", help_heading = "OUTPUT", help = "Don't use supplementary pseudomappings for short-reads")]
    pub no_supp: bool,

}
