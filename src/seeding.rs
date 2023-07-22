use crate::types::*;

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

pub fn decode(byte: u64) -> u8 {
    if byte == 0 {
        return b'A';
    } else if byte == 1 {
        return b'C';
    } else if byte == 2 {
        return b'G';
    } else if byte == 3 {
        return b'T';
    } else {
        panic!("decoding failed")
    }
}
pub fn print_string(kmer: u64, k: usize) {
    let mut bytes = vec![];
    let mask = 3;
    for i in 0..k {
        let val = kmer >> 2 * i;
        let val = val & mask;
        bytes.push(decode(val));
    }
    dbg!(std::str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
}
#[inline]
fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn fmh_seeds(string: &[u8],  w: usize, k: usize,kmer_vec: &mut Vec<u64>) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }
    let c = (w+1)/2;

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold_marker = u64::MAX / (c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);

        if hash_marker < threshold_marker {
            kmer_vec.push(hash_marker as u64);
        }
    }
}

pub fn window_sync_seeds(string: &[u8], minimizer_w: usize, k: usize, kmer_vec: &mut Vec<u64>) {
    let c = (minimizer_w + 1)/2;
    let s = 11;
    let mut w = c;
    if w%2 == 0{
        w += 1;
    }
    if string.len() < k {
        return;
    }

    type MarkerBits = u64;
    let marker_k = k;

    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let mut rolling_kmer_f_smer: MarkerBits = 0;
    let mut rolling_kmer_r_smer: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));

    let s_reverse_shift_dist = 2 * (s - 1);
    let s_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * s);
    let s_rev_mask = !(3 << (2 * s - 2));

    let len = string.len();
    let mut ring_buffer = vec![0u64; w];
    let mut minimum_pos = usize::MAX;
    let mut min_val = u64::MAX;

    //substring length is c + s - 1
    //index of last kmer is c + s - k;
    //if w = 2; offset is 1 = c - 1
    //first k bases are 
    for i in 0..len{
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;

        rolling_kmer_f_smer <<= 2;
        rolling_kmer_f_smer |= nuc_f;
        rolling_kmer_f_smer &= s_mask;
        rolling_kmer_r_smer >>= 2;
        rolling_kmer_r_smer &= s_rev_mask;
        rolling_kmer_r_smer |= nuc_r << s_reverse_shift_dist;

        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;

        if i >= w + s - k - 1{
            let old_nuc_r = 3 - BYTE_TO_SEQ[string[(i + k + 1 - s - w)] as usize] as u64;
            rolling_kmer_r_marker >>= 2;
            rolling_kmer_r_marker &= marker_rev_mask;
            rolling_kmer_r_marker |= old_nuc_r << marker_reverse_shift_dist;
        }

        if i >= s-1 {
            let buffer_ind = (i - (s - 1)) % w;
            let canonical_smer = rolling_kmer_f_smer < rolling_kmer_r_smer;
            let canonical_kmer_smer = if canonical_smer {
                rolling_kmer_f_smer
            } else {
                rolling_kmer_r_smer
            };

            let hash_smer = mm_hash64(canonical_kmer_smer);
            ring_buffer[buffer_ind] = hash_smer;
            if hash_smer < min_val{
                min_val = hash_smer;
                minimum_pos = buffer_ind;
            }
            else if buffer_ind == minimum_pos{
                minimum_pos = position_min(&ring_buffer).unwrap();
                min_val = ring_buffer[minimum_pos];
            }
            if i >= s + w - 2{
                let middle_sync_pos = (i - (s - 1) + c/2 + 1) % w;
                if minimum_pos == middle_sync_pos{
                    let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
                    let canonical_kmer_marker = if canonical_marker {
                        rolling_kmer_f_marker
                    } else {
                        rolling_kmer_r_marker
                    };

                    kmer_vec.push(mm_hash64(canonical_kmer_marker));
                }
            }
        }
    }
}

pub fn minimizer_seeds(string: &[u8], w: usize, k: usize, kmer_vec: &mut Vec<u64>) {
    if string.len() < k {
        return;
    }

    type MarkerBits = u64;
    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    let mut ring_buffer = vec![0u64; w];
    let mut minimum_pos = usize::MAX;
    let mut min_val = u64::MAX;

    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    let mut start_mini_calc = false;
    for i in marker_k - 1..len {
        let buffer_ind = (i - (marker_k - 1)) % w;
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);
        let mut new_min = false;

        if hash_marker < min_val {
            min_val = hash_marker;
            minimum_pos = buffer_ind;
            new_min = true;
        }

        ring_buffer[buffer_ind] = hash_marker;

        if start_mini_calc {
            if new_min == true{
                kmer_vec.push(min_val);
            }
            else{
                if minimum_pos == buffer_ind{
                    minimum_pos = position_min(&ring_buffer).unwrap();
                    min_val = ring_buffer[minimum_pos];
                    kmer_vec.push(min_val);
                }
            }
        } else {
            //k = 16. i = 15 is when the rolling k-mer gets populated, ring_buffer length = 1. If ring buffer
            //length (w) = 5, then i = 19 is when the buffer is filled. so i = k + w - 2 is the
            //starting point
            if i < k + w - 2 {
                continue;
            } else {
                start_mini_calc = true;
                kmer_vec.push(min_val);
            }
        }
    }
}

pub fn fmh_seeds_positions(
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold_marker = u64::MAX / (c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);

        if hash_marker < threshold_marker {
            kmer_vec.push((contig_number, i, hash_marker as u64));
        }
    }
}
