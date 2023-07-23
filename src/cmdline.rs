use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "TODO", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// TODO
    Index(IndexArgs),
    /// TODO
    Quant(QuantArgs),
}


#[derive(Args, Default)]
pub struct IndexArgs {
    #[clap(multiple=true, help = "")]
    pub files: Vec<String>,
    #[clap(short,long,default_value = "reference_index", help_heading = "OUTPUT", help = "")]
    pub reference_prefix: String,
    #[clap(short,long="individual-records", help_heading = "INPUT", help = "Use individual records (e.g. contigs) as queries instead")]
    pub individual: bool,
    #[clap(short, default_value_t = 21,help_heading = "ALGORITHM", help ="")]
    pub k: usize,
    #[clap(short, default_value_t = 40, help_heading = "ALGORITHM", help = "")]
    pub w: usize,
    #[clap(short,long, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output for debugging")]
    pub trace: bool,
    #[clap(short='1',long="first-pair", multiple=true, help_heading = "INPUT", help = "First pairs in paired end reads e.g. S1_1.fq S2_1.fq")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pair", multiple=true, help_heading = "INPUT", help = "Second pairs in paried end reads e.g. S1_2.fq S2_2.fq")]
    pub second_pair: Vec<String>,
    #[clap(long)]
    pub fmh: bool,
    #[clap(long)]
    pub window_sync: bool,
}

#[derive(Args)]
pub struct QuantArgs {
    #[clap(multiple=true, help = "Pre-sketched query or sample files. Raw fastq/fasta also allowed but presketching is recommended; see sylph sketch for more info")]
    pub files: Vec<String>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output for debugging")]
    pub trace: bool,
}
