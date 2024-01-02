use clap::Parser;
use plume::cmdline::*;
use plume::index;
use plume::quant;

//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though. 
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;


fn main() {
    let cli = Cli::parse();
    match cli.mode {
        Mode::Index(index_args) => index::index(index_args),
        Mode::Quant(quant_args) => quant::quant(quant_args),
    }
}
