//Various DNA lookup tables and hashing methods are taken from miniprot by Heng Li. Attached below is their license:
//The MIT License

// **** miniprot LICENSE ***
//Copyright (c) 2022-     Dana-Farber Cancer Institute
//
//Permission is hereby granted, free of charge, to any person obtaining
//a copy of this software and associated documentation files (the
//"Software"), to deal in the Software without restriction, including
//without limitation the rights to use, copy, modify, merge, publish,
//distribute, sublicense, and/or sell copies of the Software, and to
//permit persons to whom the Software is furnished to do so, subject to
//the following conditions:
//
//The above copyright notice and this permission notice shall be
//included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//******************************
//
//
use serde::{Deserialize, Serialize};
use fxhash::*;
pub type Kmer = u64;
pub const  DNA_TO_AA: [u8; 64] =
            *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

pub const BYTE_TO_SEQ: [u32; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];


#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct ShadeRefIndex{
    pub inv_table: FxHashMap<Kmer,Vec<u32>>,
    pub file_name: String,
    pub contig_names: Vec<String>,
    pub w: usize,
    pub k :usize,
    pub lens: Vec<usize>
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct ShadeRefIndexEncode{
    pub vec_table: Vec<(Kmer,Vec<u32>)>,
    pub file_name: String,
    pub contig_names: Vec<String>,
    pub w: usize,
    pub k :usize,
    pub lens: Vec<usize>
}


#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct RefSequenceIndex{
    pub genome_kmers: Vec<Kmer>,
    pub file_name: String,
    pub first_contig_name: String,
    pub w: usize,
    pub k: usize,
    pub len: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct ReadsIndex{
    pub read_kmers: Vec<(Vec<Kmer>,u32)>,
    pub file_name: String,
    pub w: usize,
    pub k: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct ReadRefMappings{
    pub equiv_class_matrix: FxHashMap<Vec<u32>, f64>,
    pub read_index: usize,
    pub ref_index: usize,
}
