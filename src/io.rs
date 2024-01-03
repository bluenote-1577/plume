use crate::constants::*;
use crate::types::*;
use log::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

pub fn start_serialize(reads_index: &ReadsIndex, filename: &str) -> Option<BufWriter<File>> {
    //TODO
    //    if Path::new(filename).exists(){
    //        log::warn!("Index for {} already exists. Skipping", filename);
    //        return None;
    //    }
    // Open a file in write mode.
    let file = File::create(filename).ok()?;
    let mut writer = BufWriter::new(file);
    // Serialize `file_name`.
    // Write the length of the string and then the string itself.
    writer
        .write_all(&(reads_index.file_name.len() as u64).to_le_bytes())
        .ok()?;
    writer.write_all(reads_index.file_name.as_bytes()).ok()?;

    // Serialize `w` and `k`.
    writer.write_all(&reads_index.w.to_le_bytes()).ok()?;
    writer.write_all(&reads_index.k.to_le_bytes()).ok()?;

    return Some(writer);
}

pub fn serialize_kmers(reads: Vec<(Vec<Kmer>, u32)>, writer: &mut BufWriter<File>) {
    // Serialize `read_kmers`.
    // Write the length of the vector.
    for (kmers, count) in &reads {
        // Write the length of the `kmers` vector and then its contents.
        writer
            .write_all(&(kmers.len() as u64).to_le_bytes())
            .expect("Kmer serializing failure");
        for kmer in kmers {
            writer
                .write_all(&kmer.to_le_bytes())
                .expect("Kmer serializing failure");
        }
        // Write the `count`.
        writer
            .write_all(&count.to_le_bytes())
            .expect("Kmer serializing failure");
    }

    // Make sure all data is flushed to the file.
    writer.flush().expect("Flush failure");
}

pub fn start_deserialize(filename: String) -> Option<(BufReader<File>, ReadsIndex)> {
    // Open the file in read mode.
    let file = File::open(&filename).ok()?;
    let mut reader = BufReader::new(file);
    let mut buffer = [0u8; 8]; // Helper buffer for reading u64 values.
                               // Deserialize `file_name`.
    reader.read_exact(&mut buffer).ok()?;
    let file_name_len = u64::from_le_bytes(buffer) as usize;
    let mut file_name_bytes = vec![0u8; file_name_len];
    reader.read_exact(&mut file_name_bytes).ok()?;
    let file_name = String::from_utf8(file_name_bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .ok()?;

    // Deserialize `w` and `k`.
    reader.read_exact(&mut buffer).ok()?;
    let w = u64::from_le_bytes(buffer) as usize;

    reader.read_exact(&mut buffer).ok()?;
    let k = u64::from_le_bytes(buffer) as usize;

    return Some((
        reader,
        ReadsIndex {
            file_name,
            w,
            k,
        },
    ));
}

pub fn read_batch(reader: &mut BufReader<File>, batch_size: usize) -> Vec<(Vec<Kmer>, u32)> {
    let mut buffer = [0u8; 8]; // Helper buffer for reading u64 values.
                               // Deserialize `read_kmers`.
    let mut read_kmers: Vec<(Vec<Kmer>, u32)> = vec![];
    for _ in 0..batch_size {
        // Read the length of the `kmers` vector.
        if reader.read_exact(&mut buffer).is_err() {
            break;
        }
        let kmers_len = u64::from_le_bytes(buffer) as usize;

        let mut kmers = Vec::with_capacity(kmers_len);
        for _ in 0..kmers_len {
            let mut kmer_buffer = [0u8; 8];
            reader
                .read_exact(&mut kmer_buffer)
                .expect("Deserializing reads failed");
            kmers.push(Kmer::from_le_bytes(kmer_buffer));
        }

        // Read the `length`.
        let mut buffer = [0u8; 4]; // Helper buffer for reading u64 values.
        reader
            .read_exact(&mut buffer)
            .expect("Deserializing reads failed");
        let length = u32::from_le_bytes(buffer);

        read_kmers.push((kmers, length));
    }
    return read_kmers;
}
