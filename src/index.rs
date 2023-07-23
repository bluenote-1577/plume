use crate::cmdline::*;
use std::time::Instant;
use std::fs::File;
use std::io::BufWriter;
use crate::types::*;
use fxhash::*;
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

    let start_t_initial = Instant::now();
    let mut read_files = vec![];
    let mut ref_files = vec![];
    for file in args.files.iter(){
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
    ref_files.into_par_iter().for_each(|ref_file|{
        let genome_sketch = sketch_genome_individual(args.w, args.k, &ref_file, &args);
        if let Some(ind) = genome_sketch{
            encode_and_dump(&args, ind, &ref_file);
        }
    });

    read_files.into_par_iter().for_each(|read_file|{ 
        let read_index = sketch_reads(args.w, args.k, &read_file, &args);
        if let Some(ind) = read_index{
            dump_read_indices(&args, ind, &read_file);
        }
    });

    log::info!("Indexing time {:?}", Instant::now() - start_t_initial);
}

fn encode_and_dump(args: &IndexArgs, index: Vec<RefSequenceIndex>, ref_file: &str){
    if index.is_empty(){
        warn!("{} has no valid records.", ref_file);
        return;
    }
    let mut table_ind = FxHashMap::default();
    let mut contig_names = vec![];
    let mut lens = vec![];
    let mut file_name = String::new();
    let mut w = 0;
    let mut k = 0;
    for (i,ref_ind) in index.into_iter().enumerate(){
        for kmer in ref_ind.genome_kmers{
            let ref_corresp = table_ind.entry(kmer).or_insert(vec![]);
            ref_corresp.push(i as u32);
        }
        contig_names.push(ref_ind.first_contig_name);
        lens.push(ref_ind.len);
        file_name = ref_ind.file_name;
        w = ref_ind.w;
        k = ref_ind.k;
    }
    let vec_table : Vec<(Kmer,Vec<u32>)> = table_ind.into_iter().collect();
    let enc = ShadeRefIndexEncode{ vec_table, file_name, contig_names, w, k, lens};

    let dump_str = format!("{}.shref", ref_file); 
    let mut read_sk_file = BufWriter::new(
        File::create(&dump_str)
            .expect(&format!("{} path not valid, exiting.", dump_str)),
    );
    bincode::serialize_into(&mut read_sk_file, &enc).unwrap();

}


fn dump_read_indices(args: &IndexArgs, index: ReadsIndex, read_file: &str){
    let dump_str = format!("{}.shread", read_file); 
    let mut read_sk_file = BufWriter::new(
        File::create(&dump_str)
            .expect(&format!("{} path not valid, exiting.", dump_str)),
    );
    bincode::serialize_into(&mut read_sk_file, &index).unwrap();
}

fn sketch_reads(
    w: usize, 
    k: usize,
    read_file: &str,
    args: &IndexArgs,
) -> Option<ReadsIndex>{
    let reader = parse_fastx_file(&read_file);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", read_file);
        return None;
    } else {
        let mut reader = reader.unwrap();
        let mut return_read_index = ReadsIndex::default();
        return_read_index.w = w;
        return_read_index.k = k;
        return_read_index.file_name = read_file.to_string();
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let mut kmer_vec = vec![];
                let record = record.expect(&format!("Invalid record for file {}", read_file));
                let seq = record.seq();
                let seq_len = seq.len();

                if args.window_sync{
                    window_sync_seeds(&seq, w, k, &mut kmer_vec);
                }
                else if args.fmh{
                    fmh_seeds(&seq, w, k, &mut kmer_vec);
                }
                else{
                    minimizer_seeds(&seq, w, k, &mut kmer_vec);
                }
                
                return_read_index.read_kmers.push((kmer_vec,seq_len as u32));
            } else {
                warn!("File {} is not a valid fasta/fastq file", read_file);
                return None
            }
        }
        return Some(return_read_index);
    }
}


fn sketch_genome_individual(
    w: usize,
    k: usize,
    ref_file: &str,
    args: &IndexArgs
) -> Option<Vec<RefSequenceIndex>> {
    let reader = parse_fastx_file(&ref_file);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        return None;
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
                return_genome_sketch.len = seq.len();

                if args.window_sync{
                    window_sync_seeds(&seq, w, k, &mut kmer_vec);
                }
                else if args.fmh{
                    fmh_seeds(&seq, w, k, &mut kmer_vec);
                }
                else{
                    minimizer_seeds(&seq, w, k, &mut kmer_vec);
                }

                let mut kmer_set = FxHashSet::default();
                let mut new_vec = Vec::with_capacity(kmer_vec.len());
                for km in kmer_vec.iter() {
                    if !kmer_set.contains(&km) {
                        new_vec.push(*km);
                        kmer_set.insert(km);
                    } 
                }
                return_genome_sketch.genome_kmers = new_vec;
                return_vec.push(return_genome_sketch);
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return None;
            }
        }
        return Some(return_vec);
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
