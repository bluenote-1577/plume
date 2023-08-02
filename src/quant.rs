use crate::{cmdline::*, types::*};
use crate::constants::*;
use fxhash::*;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;
use log::*;
use rayon::prelude::*;

pub fn quant(args: QuantArgs){

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
    {
        error!("No input sequences found. Exiting; use -h/--help for help.");
        std::process::exit(1);
    }


    let mut ref_index_files = vec![];
    let mut read_index_files = vec![];

    for file in args.files.iter() {
        if file.ends_with(".shref") {
            ref_index_files.push(file);
        } else if file.ends_with(".shread") {
            read_index_files.push(file);
        } else {
            warn!(
                "{} file extension is not a shade sketch.",
                &file
            );
        }
    }

    
    //Parallelize this later? TODO. Not sure 
    //the correct way to load things into memory. 
    //Let's just load things one by one for now...

    for (j,ref_index_file) in ref_index_files.iter().enumerate(){

        let ref_index = load_and_decode(ref_index_file);
        let mut result_mat = vec![vec![0.;read_index_files.len()]; ref_index.contig_names.len()];

        for (i,read_index_file) in read_index_files.iter().enumerate(){
            let buf_reader = BufReader::with_capacity(10_000_000, File::open(read_index_file).expect(&format!("{} was not a valid index", read_index_file)));
        let read_index: ReadsIndex = bincode::deserialize_from(buf_reader)
            .expect(&format!(
                "read index {} is not valid.",
                &read_index_file
            ));

            let mappings = pseudoalign(&ref_index, &read_index, i, j);
            let abund_result = em_alg(&mappings);
            if args.relative_abundance{
                for k in 0..abund_result.len(){
                    let cov = abund_result[k] * 100. * mappings.number_of_bases_mapped as f64 / mappings.number_of_bases_total as f64;
                    result_mat[k][i] = cov;
                }
            }
            else{
                for k in 0..abund_result.len(){
                    let cov = abund_result[k] * mappings.number_of_bases_mapped as f64 / ref_index.lens[k] as f64;
                    result_mat[k][i] = cov;
                }
            }
        }
        for i in 0..result_mat.len(){
            print!("{}", ref_index.contig_names[i]);
            for j in 0..result_mat[i].len(){
                print!("\t{:.4}", result_mat[i][j]);
            }
            println!("");
        }
    }
}

fn load_and_decode(ref_index_file: &str) -> ShadeRefIndex{
    let buf_reader = BufReader::with_capacity(10_000_000, File::open(ref_index_file).expect(&format!("{} was not a valid index", ref_index_file)));
    let read_index: ShadeRefIndexEncode = bincode::deserialize_from(buf_reader)
        .expect(&format!(
            "reference index {} is not valid.",
            &ref_index_file
        ));
    let inv_table = read_index.vec_table.into_iter().collect();
    return ShadeRefIndex { inv_table, file_name: read_index.file_name, contig_names: read_index.contig_names, w: read_index.w, k: read_index.k, lens: read_index.lens}

}

fn pseudoalign(ref_index: &ShadeRefIndex, read_index: &ReadsIndex, i: usize, j: usize) -> ReadRefMappings{
    let equiv_class_matrix = Mutex::new(FxHashMap::default());
    let inv_index = &ref_index.inv_table;
    let number_of_bases_mapped = Mutex::new(0);
    let number_of_bases_total = Mutex::new(0);
    read_index.read_kmers.par_iter().for_each(|(kmers, len)|{
        {
            *number_of_bases_total.lock().unwrap() += *len as usize;
        }
        let long = *len as usize > LEN_SUPP_CUTOFF;
        let mut ref_hit_map = FxHashMap::default();
        let mut ref_hit_ranges = FxHashMap::default();
        for (i,kmer) in kmers.iter().enumerate(){
            if inv_index.contains_key(kmer){
                let hit_refs = &inv_index[kmer];
                for reference in hit_refs{
                    let count = ref_hit_map.entry(*reference).or_insert(0);
                    *count += 1;
                    if long{
                        let range = ref_hit_ranges.entry(*reference).or_insert((i,i));
                        range.1 = i;
                    }
                }
            }
        }
        let mut vec_hit_map: Vec<(u32,usize)> = ref_hit_map.into_iter().collect();
        vec_hit_map.sort_unstable_by(|x,y| y.1.cmp(&x.1));
        //Take best equivalent mappings right now
        let mut last_count = 0;
        let mut equiv_class = vec![];
        let mut supp_alns = vec![];
        let mut mapped = true;
        if vec_hit_map.len() == 0 ||  vec_hit_map[0].1 < 2{
            mapped = false;
        }
        if mapped{
            if long{
                let mut overlap_vecs: Vec<(usize,usize)> = vec![];
                let max_count = vec_hit_map[0].1;
                for (reference, count) in vec_hit_map.iter(){
                    let mut overlap = false;
                    let range = ref_hit_ranges[reference];
                    for old_range in overlap_vecs.iter(){
                        if range.0 > old_range.1 || range.1 < old_range.0{
                        }
                        else{
                            overlap = true;
                            break;
                        }
                    }
                    if last_count == 0{
                        equiv_class.push(*reference);
                        overlap_vecs.push(ref_hit_ranges[reference]);
                    }
                    else if *count == max_count {
                        equiv_class.push(*reference);
                        overlap_vecs.push(ref_hit_ranges[reference]);
                    }
                    else if !overlap{
                        supp_alns.push(*reference);
                        overlap_vecs.push(ref_hit_ranges[reference]);
                    }
                    last_count = *count;
                }
                //dbg!(&supp_alns, &overlap_vecs);
            }
            else{
                for (reference, count) in vec_hit_map.iter(){
                    if last_count != 0 && *count < last_count{
                        break
                    }
                    else{
                        equiv_class.push(*reference)
                    }
                    last_count = *count;
                }
            }
        }
        if mapped{
            let mut lock = equiv_class_matrix.lock().unwrap();
            if long{
                let range = ref_hit_ranges[&equiv_class[0]];
                let class_val = lock.entry(equiv_class).or_insert(0.);
                *class_val += *len as f64 * (range.1 - range.0 + 1) as f64 / kmers.len() as f64;

                for supp_aln in supp_alns{
                    let range = ref_hit_ranges[&supp_aln];
                    let class_val = lock.entry(vec![supp_aln]).or_insert(0.);
                    *class_val += *len as f64 * (range.1 - range.0 + 1) as f64 / kmers.len() as f64;
                }
            }
            else{
                let class_val = lock.entry(equiv_class).or_insert(0.);
                *class_val += *len as f64;
            }
            if mapped{
                *number_of_bases_mapped.lock().unwrap() += *len as usize;
            }
        }
        //*class_val += 1.;
});
    let equiv_class_matrix = equiv_class_matrix.into_inner().unwrap();
    let number_of_bases_mapped = number_of_bases_mapped.into_inner().unwrap();
    let number_of_bases_total = number_of_bases_total.into_inner().unwrap();
    let ret_mappings = ReadRefMappings {number_of_refs: ref_index.contig_names.len(), equiv_class_matrix, read_index: i, ref_index: j, number_of_bases_mapped, number_of_bases_total};
    log::trace!("{:?}", &ret_mappings);

    return ret_mappings
}

fn em_alg(result: &ReadRefMappings) -> Vec<f64>{

    let mut ns = vec![1./result.number_of_refs as f64 ;result.number_of_refs];
    let mut next_ns = vec![0.; result.number_of_refs];

    for _ in 0..100{
        for (class, probs) in result.equiv_class_matrix.iter(){
            let mut normal_factor = 0.;
            for index in class.iter(){
                normal_factor += ns[*index as usize];
            }
            for index in class.iter(){
                next_ns[*index as usize] += probs * ns[*index as usize] / normal_factor;
            }
        }
        let total_sum :f64 = next_ns.iter().sum();
        ns = next_ns.into_iter().map(|x| x/total_sum).collect();
        next_ns = vec![0.; result.number_of_refs];
    }

    return ns;

}
