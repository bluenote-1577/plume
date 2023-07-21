use crate::cmdline::*;
use std::time::Instant;
use crate::types::*;
use fxhash::FxHashSet;
use log::*;
use rayon::prelude::*;
use needletail::parse_fastx_file;
use crate::seeding::*;

pub fn index(args: IndexArgs){

    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else {
        level = log::LevelFilter::Info;
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    if args.files.is_empty()
        && args.first_pair.is_empty()
        && args.second_pair.is_empty()
    {
        error!("No input sequences found. Exiting.");
        std::process::exit(1);
    }

    let mut read_files = vec![];
    let mut ref_files = vec![];
    for file in args.files{
        if is_fastq(&file){
            read_files.push(file);
        }
        else if is_fasta(&file){
            ref_files.push(file);
        }
        else{
            warn!("{} does not have a gz, bz2 fasta/fastq extensions -- skipping.", &file);
        }
    }
    for ref_file in ref_files{
        sketch_genome_individual(args.w, args.k, &ref_file);
    }
}

fn sketch_genome_individual(
    w: usize,
    k: usize,
    ref_file: &str,
) -> Vec<RefSequenceIndex> {
    let reader = parse_fastx_file(&ref_file);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        return vec![];
    } else {
        let mut reader = reader.unwrap();
        let mut return_vec = vec![];
        while let Some(record) = reader.next() {
            let mut return_genome_sketch = RefSequenceIndex::default();
            return_genome_sketch.w = w;
            return_genome_sketch.k = k;
            return_genome_sketch.file_name = ref_file.to_string();
            if record.is_ok() {
                let mut kmer_vec = vec![];
                let record = record.expect(&format!("Invalid record for file {}", ref_file));
                let contig_name = String::from_utf8_lossy(record.id()).to_string();
                return_genome_sketch.first_contig_name = contig_name;
                let seq = record.seq();

                let start_t = Instant::now();
                minimizer_seeds(&seq, w, k, &mut kmer_vec);
                println!("{:?}",Instant::now() - start_t);
//                let start_t = Instant::now();
//                fmh_seeds(&seq, &mut kmer_vec, w/2, k);
//                println!("{:?}",Instant::now() - start_t);

                let mut kmer_set = FxHashSet::default();
                let mut duplicate_set = FxHashSet::default();
                let mut new_vec = Vec::with_capacity(kmer_vec.len());
                for km in kmer_vec.iter() {
                    if !kmer_set.contains(&km) {
                        kmer_set.insert(km);
                    } else {
                        duplicate_set.insert(km);
                    }
                }
                return_genome_sketch.genome_kmers = new_vec;
                return_vec.push(return_genome_sketch);
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return vec![];
            }
        }
        return return_vec;
    }
}

pub fn is_fastq(file: &str) -> bool {
    if file.ends_with(".fq")
        || file.ends_with(".fnq")
        || file.ends_with(".fastq")
        || file.ends_with(".fq.gz")
        || file.ends_with(".fnq.gz")
        || file.ends_with(".fastq.gz")
        || file.ends_with(".fastq.bz2")
        || file.ends_with(".fq.bz2")
        || file.ends_with(".fnq.bz2")

    {
        return true;
    } else {
        return false;
    }
}

pub fn is_fasta(file: &str) -> bool {
    if file.ends_with(".fa")
        || file.ends_with(".fna")
        || file.ends_with(".fasta")
        || file.ends_with(".fa.gz")
        || file.ends_with(".fna.gz")
        || file.ends_with(".fasta.gz")
        || file.ends_with(".fa.bz2")
        || file.ends_with(".fna.bz2")
        || file.ends_with(".fasta.bz2")

    {
        return true;
    } else {
        return false;
    }
}
