use crate::{cmdline::*, types::*};
use fxhash::*;
use std::fs::File;
use std::io::BufReader;
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
        error!("No input sequences found. Exiting.");
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

    for (i,read_index_file) in read_index_files.iter().enumerate(){
        let buf_reader = BufReader::with_capacity(10_000_000, File::open(read_index_file).expect(&format!("{} was not a valid index", read_index_file)));
        let read_index: ReadsIndex = bincode::deserialize_from(buf_reader)
            .expect(&format!(
                "read index {} is not valid.",
                &read_index_file
            ));

        for (j,ref_index_file) in ref_index_files.iter().enumerate(){
            let ref_index = load_and_decode(ref_index_file);
            let mappings = pseudoalign(&ref_index, &read_index, i, j);
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
    let mut equiv_class_matrix = FxHashMap::default();
    let inv_index = &ref_index.inv_table;
    for (kmers, len) in read_index.read_kmers.iter(){
        let mut ref_hit_map = FxHashMap::default();
        for kmer in kmers{
            if inv_index.contains_key(kmer){
                let hit_refs = &inv_index[kmer];
                for reference in hit_refs{
                    let count = ref_hit_map.entry(*reference).or_insert(0);
                    *count += 1;
                }
            }
        }
        let mut vec_hit_map: Vec<(u32,usize)> = ref_hit_map.into_iter().collect();
        vec_hit_map.sort_unstable_by(|x,y| y.1.cmp(&x.1));
        //Take best equivalent mappings right now
        let mut last_count = 0;
        let mut equiv_class = vec![];
        for (reference, count) in vec_hit_map.iter(){
            if last_count != 0 && *count != last_count{
                break
            }
            else{
                equiv_class.push(*reference)
            }
            last_count = *count;
        }
        let class_val = equiv_class_matrix.entry(equiv_class).or_insert(0.);
        *class_val += 1.;

    }
    dbg!(&equiv_class_matrix);
    return ReadRefMappings { equiv_class_matrix, read_index: i, ref_index: j}
}
