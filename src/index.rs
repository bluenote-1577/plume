use crate::cmdline::*;
use crate::constants::*;
use crate::io::*;
use crate::seeding::*;
use crate::types::*;
use fxhash::*;
use log::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use smallvec::*;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::sync::Mutex;
use std::time::Instant;

pub fn index(args: IndexArgs) {
    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
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

    check_valid_inputs(&args);

    let start_t_initial = Instant::now();

    //Do paired reads
    info!("Indexing reference/genome files...");

    let genome_sketch = sketch_genome_all(args.w, args.k, &args.genomes, &args);
    if let Some(ind) = genome_sketch {
        encode_and_dump(&args, ind, "genomes_combined");
    }

    args.contigs.par_iter().for_each(|ref_file| {
        let genome_sketch = sketch_contig_assembly(args.w, args.k, ref_file, &args);
        if let Some(ind) = genome_sketch {
            encode_and_dump(&args, ind, &ref_file);
        }
    });

    info!("Indexing read files...");
    args.single_reads.par_iter().for_each(|read_file| {
        sketch_reads_stream(args.w, args.k, read_file, &args);
    });

    (0..args.first_pair.len()).into_par_iter().for_each(|i| {
        let read_file1 = &args.first_pair[i];
        let read_file2 = &args.second_pair[i];
        sketch_paired_reads_stream(args.w, args.k, &read_file1, &read_file2, &args);
    });

    log::info!("Indexing time {:?}", Instant::now() - start_t_initial);
}

fn check_valid_inputs(args: &IndexArgs) {
    if args.contigs.is_empty()
        && args.first_pair.is_empty()
        && args.second_pair.is_empty()
        && args.genomes.is_empty()
        && args.single_reads.is_empty()
    {
        error!("No input sequences found. Exiting; use -h/--help for help.");
        std::process::exit(1);
    }

    if args.first_pair.len() != args.second_pair.len() {
        error!("Different number of files in -1 option compared to -2 option. Not all paired-end reads are paired. Exiting.");
        std::process::exit(1);
    }
}

fn encode_and_dump(_args: &IndexArgs, index: Vec<RefSequenceIndex>, ref_file: &str) {
    if index.is_empty() {
        warn!("{} has no valid records.", ref_file);
        return;
    }
    let mut table_ind = FxHashMap::default();
    let mut contig_names = vec![];
    let mut lens = vec![];
    let mut file_name = String::new();
    let mut w = 0;
    let mut k = 0;
    let mut num_kmers_in_contigs = vec![];
    for (i, ref_ind) in index.into_iter().enumerate() {
        num_kmers_in_contigs.push(ref_ind.genome_kmers.len());
        for kmer in ref_ind.genome_kmers {
            let ref_corresp = table_ind.entry(kmer).or_insert(smallvec![]);
            ref_corresp.push(i as u32);
        }
        contig_names.push(ref_ind.first_contig_name);
        lens.push(ref_ind.len);
        file_name = ref_ind.file_name;
        w = ref_ind.w;
        k = ref_ind.k;
    }
    let vec_table: Vec<(Kmer, SmallVec<[u32; 1]>)> = table_ind.into_iter().collect();
    let enc = ShadeRefIndexEncode {
        vec_table,
        file_name,
        contig_names,
        num_kmers_in_contigs,
        w,
        k,
        lens,
    };

    let dump_str = format!("{}.{}", ref_file, REF_IND_SUFFIX);
    let mut read_sk_file = BufWriter::new(
        File::create(&dump_str).expect(&format!("{} path not valid, exiting.", dump_str)),
    );
    bincode::serialize_into(&mut read_sk_file, &enc).unwrap();
    log::info!("Finished indexing references at {}", &dump_str);
}

//fn dump_read_indices(_args: &IndexArgs, index: ReadsIndex, read_file: &str) {
//    let dump_str = format!("{}.{}", read_file, READ_IND_SUFFIX);
//    let mut read_sk_file = BufWriter::new(
//        File::create(&dump_str).expect(&format!("{} path not valid, exiting.", dump_str)),
//    );
//    bincode::serialize_into(&mut read_sk_file, &index).unwrap();
//    log::info!("Finished indexing reads at {}", &dump_str);
//}


fn sketch_paired_reads_stream(
    w: usize,
    k: usize,
    read_file1: &str,
    read_file2: &str,
    args: &IndexArgs,
) {
    let reader1 = parse_fastx_file(&read_file1);
    let reader2 = parse_fastx_file(&read_file2);
    if !reader1.is_ok() || !reader2.is_ok() {
        warn!(
            "{} or {} is not a valid fasta/fastq file; skipping.",
            read_file1, read_file2
        );
        return;
    } else {
        let mut reader1 = reader1.unwrap();
        let mut reader2 = reader2.unwrap();
        let mut return_read_index = ReadsIndex::default();
        return_read_index.w = w;
        return_read_index.k = k;
        return_read_index.file_name = read_file1.to_string();

        //Start serializing
        let dump_str = format!("{}.paired.{}", read_file1, READ_IND_SUFFIX);
        let opt = start_serialize(&return_read_index, &dump_str);
        if opt.is_none() {
            return;
        }
        let mut writer = opt.unwrap();
        let mut read_kmers = vec![];
        let mut num_bases = 0;

        let mut counter = 0;
        while true{
            let mut pass1 = false;
            if let Some(record) = reader1.next() {
                if record.is_ok() {
                    pass1 = true;
                    let mut kmer_vec = vec![];
                    let record = record.expect(&format!("Invalid record for file {}", read_file1));
                    let seq = record.seq();
                    let seq_len = seq.len();

                    if args.window_sync {
                        window_sync_seeds(&seq, w, k, &mut kmer_vec);
                    } else if args.fmh {
                        fmh_seeds(&seq, w, k, &mut kmer_vec);
                    } else if args.minimizer {
                        minimizer_seeds(&seq, w, k, &mut kmer_vec);
                    } else{
                        window_sync_seeds(&seq, w, k, &mut kmer_vec);
                    }

                    read_kmers
                        .push((kmer_vec, seq_len as u32));
                } else {
                    warn!("File {} is not a valid fasta/fastq file", read_file1);
                    return;
                }
            }
            let mut pass2 = false;
            if let Some(record) = reader2.next() {
                if record.is_ok() {
                    pass2 = true;
                    let mut kmer_vec = vec![];
                    let record = record.expect(&format!("Invalid record for file {}", read_file1));
                    let seq = record.seq();
                    let seq_len = seq.len();

                    if args.window_sync {
                        window_sync_seeds(&seq, w, k, &mut kmer_vec);
                    } else if args.fmh {
                        fmh_seeds(&seq, w, k, &mut kmer_vec);
                    } else if args.minimizer{
                        minimizer_seeds(&seq, w, k, &mut kmer_vec);
                    } else{
                        window_sync_seeds(&seq, w, k, &mut kmer_vec);
                    }

                    if counter >= read_kmers.len() {
                        warn!("File {}/{} have different number of records and is not a valid fastq pair. The index {} may be invalid.", read_file1, read_file2, dump_str);
                        return;
                    }
                    read_kmers[counter].0.extend(kmer_vec);
                    read_kmers[counter].1 += seq_len as u32;
                    num_bases += read_kmers[counter].1;
                    counter += 1;
                } else {
                    warn!("File {} is not a valid fasta/fastq file", read_file1);
                    return;
                }
            }
            else if counter != read_kmers.len(){
                warn!("File {}/{} have different number of records and is not a valid fastq pair. The index {} may be invalid.", read_file1, read_file2, dump_str);
                return;

            }
            if num_bases as usize > BATCH_SIZE_BASE {
                let dump_kmers = std::mem::take(&mut read_kmers);
                serialize_kmers(dump_kmers, &mut writer);
                counter = 0;
                num_bases = 0;
            }
            if !pass1 && !pass2{
                break;
            }
        }
        if !read_kmers.is_empty() {
            serialize_kmers(read_kmers, &mut writer);
        }
        log::info!("Finished indexing reads at {}", &dump_str);
    }
}

fn sketch_reads_stream(w: usize, k: usize, read_file: &str, args: &IndexArgs) {
    let reader = parse_fastx_file(&read_file);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", read_file);
        return;
    } else {
        let mut reader = reader.unwrap();
        let mut return_read_index = ReadsIndex::default();
        return_read_index.w = w;
        return_read_index.k = k;
        return_read_index.file_name = read_file.to_string();

        //Start serializing
        let dump_str = format!("{}.{}", read_file, READ_IND_SUFFIX);
        let opt = start_serialize(&return_read_index, &dump_str);
        if opt.is_none() {
            return;
        }
        let mut writer = opt.unwrap();
        let mut read_kmers = vec![];
        let mut num_bases = 0;

        //Streaming deserialization
        while let Some(record) = reader.next() {
            if record.is_ok() {
                
                let mut kmer_vec = vec![];
                let record = record.expect(&format!("Invalid record for file {}", read_file));
                let seq = record.seq();
                let seq_len = seq.len();
                num_bases += seq_len;

                if args.window_sync {
                    window_sync_seeds(&seq, w, k, &mut kmer_vec);
                } else if args.fmh {
                    fmh_seeds(&seq, w, k, &mut kmer_vec);
                } else if args.minimizer{
                    minimizer_seeds(&seq, w, k, &mut kmer_vec);
                } else{
                    window_sync_seeds(&seq, w, k, &mut kmer_vec);
                }

                read_kmers.push((kmer_vec, seq_len as u32));
            } else {
                warn!("File {} is not a valid fasta/fastq file", read_file);
                return;
            }
            if num_bases as usize> BATCH_SIZE_BASE {
                let dump_kmers = std::mem::take(&mut read_kmers);
                serialize_kmers(dump_kmers, &mut writer);
                num_bases = 0;
            }
        }

        if !read_kmers.is_empty() {
            serialize_kmers(read_kmers, &mut writer);
        }
        log::info!("Finished indexing reads at {}", &dump_str);
    }

}

fn sketch_contig_assembly(
    w: usize,
    k: usize,
    ref_file: &str,
    args: &IndexArgs,
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

                if args.window_sync {
                    window_sync_seeds(&seq, w, k, &mut kmer_vec);
                } else if args.fmh {
                    fmh_seeds(&seq, w, k, &mut kmer_vec);
                } else if args.minimizer {
                    minimizer_seeds(&seq, w, k, &mut kmer_vec);
                } else{
                    window_sync_seeds(&seq, w, k, &mut kmer_vec);
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

fn sketch_genome_all(
    w: usize,
    k: usize,
    ref_files: &Vec<String>,
    args: &IndexArgs,
) -> Option<Vec<RefSequenceIndex>> {
    if ref_files.is_empty() {
        return None;
    }
    let return_vec = Mutex::new(vec![]);
    ref_files.par_iter().for_each(|ref_file| {
        let mut return_genome_sketch = RefSequenceIndex::default();
        return_genome_sketch.w = w;
        return_genome_sketch.k = k;
        let reader = parse_fastx_file(&ref_file);
        return_genome_sketch.first_contig_name = ref_file.to_string();
        return_genome_sketch.file_name = ref_file.to_string();
        let mut kmer_set = FxHashSet::default();
        let mut kmer_vec = vec![];
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        } else {
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record = record.expect(&format!("Invalid record for file {}", ref_file));
                    let seq = record.seq();
                    return_genome_sketch.len += seq.len();

                    if args.window_sync {
                        window_sync_seeds(&seq, w, k, &mut kmer_vec);
                    } else if args.fmh {
                        fmh_seeds(&seq, w, k, &mut kmer_vec);
                    } else if args.minimizer{
                        minimizer_seeds(&seq, w, k, &mut kmer_vec);
                    } else{
                        window_sync_seeds(&seq, w, k, &mut kmer_vec);
                    }
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                }
            }
        }
        let mut new_vec = Vec::with_capacity(kmer_vec.len());
        for km in kmer_vec.iter() {
            if !kmer_set.contains(&km) {
                new_vec.push(*km);
                kmer_set.insert(km);
            }
        }
        return_genome_sketch.genome_kmers = new_vec;

        return_vec.lock().unwrap().push(return_genome_sketch);
    });
    return Some(
        return_vec
            .into_inner()
            .expect("Could not unlock mutex while multiprocessing sketching"),
    );
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
