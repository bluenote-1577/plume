use crate::constants::*;
use crate::io::*;
use crate::{cmdline::*, types::*};
use fxhash::*;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::sync::mpsc;
use std::sync::mpsc::TryRecvError;
use std::sync::{Arc, Mutex};
use std::thread;

pub fn quant(args: QuantArgs) {
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

    if args.files.is_empty() {
        error!("No input sequences found. Exiting; use -h/--help for help.");
        std::process::exit(1);
    }

    let mut ref_index_files = vec![];
    let mut read_index_files = vec![];

    for file in args.files.iter() {
        if file.ends_with(REF_IND_SUFFIX) {
            ref_index_files.push(file);
        } else if file.ends_with(READ_IND_SUFFIX) {
            read_index_files.push(file);
        } else {
            warn!("{} file extension is not a plume index.", &file);
        }
    }

    //Parallelize this later? TODO. Not sure
    //the correct way to load things into memory.
    //Let's just load things one by one for now...

    for ref_index_file in ref_index_files.iter(){
        log::info!("Loading index for {} ...", ref_index_file);
        let ref_index = load_and_decode(ref_index_file);
        let mut result_mat = vec![vec![0.; read_index_files.len()]; ref_index.contig_names.len()];

        log::info!("Pseudoaligning reads ...");
        for (i, read_index_file) in read_index_files.iter().enumerate() {
            let mappings = streaming_pseudoalign(&args, read_index_file, &ref_index);

            //let mappings = pseudoalign(&ref_index, &read_index, i, j);
            let abund_result = em_alg(&mappings);
            if args.relative_abundance {
                for k in 0..abund_result.len() {
                    let cov = abund_result[k] * 100. * mappings.number_of_bases_mapped as f64
                        / mappings.number_of_bases_total as f64;
                    result_mat[k][i] = cov;
                }
            } else {
                for k in 0..abund_result.len() {
                    let modifier;
                    if ref_index.lens[k] > 400{
                        modifier = 50.;
                    }
                    else{
                        modifier = 0.;
                    }
                    let cov = abund_result[k] * mappings.number_of_bases_mapped as f64
                        / (ref_index.lens[k] as f64 - modifier);
                    result_mat[k][i] = cov;
                }
            }
        }

        log::info!("Quantification complete for {}", ref_index_file);
        print!("Sequence_names");
        for read_index_file in read_index_files.iter() {
            print!("\t{}", read_index_file);
        }
        println!("");
        for i in 0..result_mat.len() {
            print!("{}", ref_index.contig_names[i].split(' ').collect::<Vec<_>>()[0]);
            for j in 0..result_mat[i].len() {
                print!("\t{:.4}", result_mat[i][j]);
            }
            println!("");
        }
    }
}

fn pseudoalign_reads(
    args: &QuantArgs,
    read_kmers: Vec<(Vec<Kmer>, u32)>,
    ref_index: &ShadeRefIndex,
    equiv_class_matrix: Arc<Mutex<FxHashMap<Vec<u32>, f64>>>,
    number_of_bases_mapped: Arc<Mutex<usize>>,
    number_of_bases_total: Arc<Mutex<usize>>,
) {
    let inv_index = &ref_index.inv_table;
    let test_supp = !args.no_supp;
    read_kmers.par_iter().for_each(|(kmers, len)| {
        {
            {
                *number_of_bases_total.lock().unwrap() += *len as usize;
            }
            let long = *len as usize > LEN_SUPP_CUTOFF;
            let mut ref_hit_map = FxHashMap::default();
            let mut ref_hit_ranges = FxHashMap::default();
            for (i, kmer) in kmers.iter().enumerate() {
                if inv_index.contains_key(kmer) {
                    let hit_refs = &inv_index[kmer];
                    for reference in hit_refs {
                        let count = ref_hit_map.entry(*reference).or_insert(0);
                        *count += 1;
                        let range = ref_hit_ranges.entry(*reference).or_insert((i, i));
                        range.1 = i;
                    }
                }
            }
            let mut vec_hit_map: Vec<(u32, usize)> = ref_hit_map.into_iter().collect();
            vec_hit_map.sort_unstable_by(|x, y| y.1.cmp(&x.1));
            //Take best equivalent mappings right now
            let mut last_count = 0;
            let mut equiv_class = vec![];
            let mut supp_alns = vec![];
            let mut mapped = true;
            let kmer_pass;
            if long {
                kmer_pass = 5;
            } else {
                let k = ref_index.k as f64;
                let c = (ref_index.w + 1) as f64 / 2.;
                let eff_mult = f64::min(1., k / c);
                let eff_ss = *len as f64 / k * eff_mult;
                let expected =
                    (0.95_f64).powi(ref_index.k.try_into().unwrap()) * *len as f64 / k as f64 * 1.5;
                kmer_pass = usize::max((expected / 2.0) as usize, 2);
//                if hash(&kmers[0]) % 10000 == 0 {
//                    log::trace!("{} - len {}, num_kmer {}", kmer_pass, len, kmers.len());
//                }
            }
            if vec_hit_map.len() == 0 || vec_hit_map[0].1 < kmer_pass {
                mapped = false;
            }
            if mapped {
                if long && test_supp{
                    let mut overlap_vecs: Vec<(usize, usize)> = vec![];
                    let max_count = vec_hit_map[0].1;
                    for (reference, count) in vec_hit_map.iter() {
                        let mut overlap = false;
                        let range = ref_hit_ranges[reference];
                        for old_range in overlap_vecs.iter() {
                            if range.0 > old_range.1 || range.1 < old_range.0 {
                            } else {
                                overlap = true;
                                break;
                            }
                        }
                        if last_count == 0 {
                            equiv_class.push(*reference);
                            overlap_vecs.push(ref_hit_ranges[reference]);
                        } else if *count == max_count {
                            equiv_class.push(*reference);
                            overlap_vecs.push(ref_hit_ranges[reference]);
                        } else if !overlap && count >= &kmer_pass{
                            supp_alns.push(*reference);
                            overlap_vecs.push(ref_hit_ranges[reference]);
                        }
                        last_count = *count;
                    }
//                    if supp_alns.len() > 0{
//                        dbg!(&supp_alns, &overlap_vecs, &equiv_class);
//                    }
                } else {
                    for (reference, count) in vec_hit_map.iter() {
                        if last_count != 0 && *count < last_count {
                            break;
                        } else {
                            equiv_class.push(*reference)
                        }
                        last_count = *count;
                    }
                }
            }
            if mapped {
                let mut lock = equiv_class_matrix.lock().unwrap();
                if long && test_supp{
                    let range = ref_hit_ranges[&equiv_class[0]];
                    let class_val = lock.entry(equiv_class).or_insert(0.);
                    *class_val += *len as f64 * (range.1 - range.0 + 1) as f64 / kmers.len() as f64;

                    for supp_aln in supp_alns {
                        let range = ref_hit_ranges[&supp_aln];
                        let class_val = lock.entry(vec![supp_aln]).or_insert(0.);
                        *class_val +=
                            *len as f64 * (range.1 - range.0 + 1) as f64 / kmers.len() as f64;
                    }
                } else {
                    let class_val = lock.entry(equiv_class).or_insert(0.);
                    *class_val += *len as f64;
                }
                if mapped {
                    *number_of_bases_mapped.lock().unwrap() += *len as usize;
                }
            }
            //*class_val += 1.;
        }
    });
}

fn streaming_pseudoalign(args: &QuantArgs, read_index_file: &str, ref_index: &ShadeRefIndex) -> ReadRefMappings {
    let equiv_class_matrix = Arc::new(Mutex::new(FxHashMap::default()));
    let number_of_bases_mapped = Arc::new(Mutex::new(0));
    let number_of_bases_total = Arc::new(Mutex::new(0));

    let (sender, receiver) = mpsc::channel();
    let receiver = Arc::new(Mutex::new(receiver));

    // Producer thread
    let sender_clone = sender.clone();
    let read_index_file = read_index_file.to_owned();

    let producer_handle = thread::spawn(move || {
        let opt = start_deserialize(read_index_file);
        if let Some((mut reader, _read_index)) = opt {
            let mut reads = read_batch(&mut reader, 100_000);
            while !reads.is_empty() {
                match sender_clone.send(Some(reads)) {
                    Ok(_) => reads = read_batch(&mut reader, 100_00),
                    Err(_) => break, // Receiver has dropped, exit loop
                }
            }
            let _ = sender_clone.send(None);
        }
    });

    // Consumer thread
    let e = equiv_class_matrix.clone();
    let f = number_of_bases_mapped.clone();
    let g = number_of_bases_total.clone();
    thread::scope(|s| {
        s.spawn(move || {
            loop {
                let recv_result = receiver.lock().expect("Receiver lock poisoned");
                match recv_result.try_recv() {
                    Ok(Some(kmers)) => {
                        pseudoalign_reads(&args, kmers, &ref_index, e.clone(), f.clone(), g.clone());
                        log::trace!("Batch aligned");
                    }
                    Ok(None) => {
                        break
                    }
                    Err(TryRecvError::Empty) => {
                        thread::yield_now();
                    }
                    Err(_) => break, // Producer has finished
                }
            }
        });
    });

    producer_handle.join().unwrap();


    let equiv_class_matrix = Arc::try_unwrap(equiv_class_matrix)
        .expect("Arc unwrap failed")
        .into_inner()
        .expect("Mutex poisoned");
    let number_of_bases_mapped = Arc::try_unwrap(number_of_bases_mapped)
        .expect("Arc unwrap failed")
        .into_inner()
        .expect("Mutex poisoned");
    let number_of_bases_total = Arc::try_unwrap(number_of_bases_total)
        .expect("Arc unwrap failed")
        .into_inner()
        .expect("Mutex poisoned");

    ReadRefMappings {
        number_of_refs: ref_index.contig_names.len(),
        equiv_class_matrix,
        read_index: 0,
        ref_index: 0,
        number_of_bases_mapped,
        number_of_bases_total,
    }
}

fn load_and_decode(ref_index_file: &str) -> ShadeRefIndex {
    let buf_reader = BufReader::with_capacity(
        10_000_000,
        File::open(ref_index_file).expect(&format!("{} was not a valid index", ref_index_file)),
    );
    let read_index: ShadeRefIndexEncode = bincode::deserialize_from(buf_reader).expect(&format!(
        "reference index {} is not valid.",
        &ref_index_file
    ));
    let inv_table = read_index.vec_table.into_iter().collect();
    return ShadeRefIndex {
        inv_table,
        file_name: read_index.file_name,
        contig_names: read_index.contig_names,
        num_kmers_in_contigs: read_index.num_kmers_in_contigs,
        w: read_index.w,
        k: read_index.k,
        lens: read_index.lens,
    };
}

//fn pseudoalign(
//    ref_index: &ShadeRefIndex,
//    read_index: &ReadsIndex,
//    i: usize,
//    j: usize,
//) -> ReadRefMappings {
//    let equiv_class_matrix = Mutex::new(FxHashMap::default());
//    let inv_index = &ref_index.inv_table;
//    let number_of_bases_mapped = Mutex::new(0);
//    let number_of_bases_total = Mutex::new(0);
//    read_index.read_kmers.par_iter().for_each(|(kmers, len)| {
//        {
//            *number_of_bases_total.lock().unwrap() += *len as usize;
//        }
//        let long = *len as usize > LEN_SUPP_CUTOFF;
//        let mut ref_hit_map = FxHashMap::default();
//        let mut ref_hit_ranges = FxHashMap::default();
//        for (i, kmer) in kmers.iter().enumerate() {
//            if inv_index.contains_key(kmer) {
//                let hit_refs = &inv_index[kmer];
//                for reference in hit_refs {
//                    let count = ref_hit_map.entry(*reference).or_insert(0);
//                    *count += 1;
//                    let range = ref_hit_ranges.entry(*reference).or_insert((i, i));
//                    range.1 = i;
//                }
//            }
//        }
//        let mut vec_hit_map: Vec<(u32, usize)> = ref_hit_map.into_iter().collect();
//        vec_hit_map.sort_unstable_by(|x, y| y.1.cmp(&x.1));
//        //Take best equivalent mappings right now
//        let mut last_count = 0;
//        let mut equiv_class = vec![];
//        let mut supp_alns = vec![];
//        let mut mapped = true;
//        let kmer_pass;
//        if long {
//            kmer_pass = kmers.len() / 10;
//        } else {
//            let k = ref_index.k as f64;
//            let c = (ref_index.w + 1) as f64 / 2.;
//            let eff_mult = f64::min(1., k / c);
//            let eff_ss = *len as f64 / k * eff_mult;
//            let expected =
//                (0.95_f64).powi(ref_index.k.try_into().unwrap()) * *len as f64 / k as f64 * 1.5;
//            kmer_pass = usize::max((expected / 2.0) as usize, 2);
////            if hash(&kmers[0]) % 10000 == 0 {
////                log::trace!("{} - len {}, num_kmer {}", kmer_pass, len, kmers.len());
////            }
//        }
//        if vec_hit_map.len() == 0 || vec_hit_map[0].1 < kmer_pass {
//            mapped = false;
//        }
//        if mapped {
//            if long {
//                let mut overlap_vecs: Vec<(usize, usize)> = vec![];
//                let max_count = vec_hit_map[0].1;
//                for (reference, count) in vec_hit_map.iter() {
//                    let mut overlap = false;
//                    let range = ref_hit_ranges[reference];
//                    for old_range in overlap_vecs.iter() {
//                        if range.0 > old_range.1 || range.1 < old_range.0 {
//                        } else {
//                            overlap = true;
//                            break;
//                        }
//                    }
//                    if last_count == 0 {
//                        equiv_class.push(*reference);
//                        overlap_vecs.push(ref_hit_ranges[reference]);
//                    } else if *count == max_count {
//                        equiv_class.push(*reference);
//                        overlap_vecs.push(ref_hit_ranges[reference]);
//                    } else if !overlap {
//                        supp_alns.push(*reference);
//                        overlap_vecs.push(ref_hit_ranges[reference]);
//                    }
//                    last_count = *count;
//                }
//                //dbg!(&supp_alns, &overlap_vecs);
//            } else {
//                for (reference, count) in vec_hit_map.iter() {
//                    if last_count != 0 && *count < last_count {
//                        break;
//                    } else {
//                        equiv_class.push(*reference)
//                    }
//                    last_count = *count;
//                }
//            }
//        }
//        if mapped {
//            let mut lock = equiv_class_matrix.lock().unwrap();
//            if long {
//                let range = ref_hit_ranges[&equiv_class[0]];
//                let class_val = lock.entry(equiv_class).or_insert(0.);
//                *class_val += *len as f64 * (range.1 - range.0 + 1) as f64 / kmers.len() as f64;
//
//                for supp_aln in supp_alns {
//                    let range = ref_hit_ranges[&supp_aln];
//                    let class_val = lock.entry(vec![supp_aln]).or_insert(0.);
//                    *class_val += *len as f64 * (range.1 - range.0 + 1) as f64 / kmers.len() as f64;
//                }
//            } else {
//                let class_val = lock.entry(equiv_class).or_insert(0.);
//                *class_val += *len as f64;
//            }
//            if mapped {
//                *number_of_bases_mapped.lock().unwrap() += *len as usize;
//            }
//        }
//        //*class_val += 1.;
//    });
//    let equiv_class_matrix = equiv_class_matrix.into_inner().unwrap();
//    let number_of_bases_mapped = number_of_bases_mapped.into_inner().unwrap();
//    let number_of_bases_total = number_of_bases_total.into_inner().unwrap();
//    let ret_mappings = ReadRefMappings {
//        number_of_refs: ref_index.contig_names.len(),
//        equiv_class_matrix,
//        read_index: i,
//        ref_index: j,
//        number_of_bases_mapped,
//        number_of_bases_total,
//    };
//    log::trace!("{:?}", &ret_mappings);
//
//    return ret_mappings;
//}

fn em_alg(result: &ReadRefMappings) -> Vec<f64> {
    let mut ns = vec![1. / result.number_of_refs as f64; result.number_of_refs];
    let mut next_ns = vec![0.; result.number_of_refs];

    for _ in 0..100{
        for (class, probs) in result.equiv_class_matrix.iter() {
            let mut normal_factor = 0.;
            for index in class.iter() {
                normal_factor += ns[*index as usize];
            }
            for index in class.iter() {
                next_ns[*index as usize] += probs * ns[*index as usize] / normal_factor;
            }
        }
        let total_sum: f64 = next_ns.iter().sum();
        ns = next_ns.into_iter().map(|x| x / total_sum).collect();
        next_ns = vec![0.; result.number_of_refs];
    }

    return ns;
}
