/**
 * file: fastq.rs
 * desc: Write simulated reads to a FASTQ file.
 */
use std::fs;
use std::io::Write;

use crate::genome;
use crate::simulate::SimulatedRead;
use crate::util;

pub fn write_to_fastq(
    uuid: &genome::UUID,
    reads: &Vec<SimulatedRead>,
    output: &str,
    header_format: &str,
    append: bool,
) -> Result<(), String> {
    // Set up file for writing
    let mut file = match fs::OpenOptions::new()
        .write(true)
        .append(append)
        .create(true)
        .open(output)
    {
        Ok(f) => f,
        Err(e) => return Err(format!("{}", e)),
    };

    for read in reads.iter() {
        // Use the provided header format and interpolate where necessary
        let fwd_header = header_format
            .replace("{:genome_id:}", &uuid.to_string())
            .replace("{:read_id:}", &read.id.to_string())
            .replace("{:pair:}", "1");
        let rev_header = header_format
            .replace("{:genome_id:}", &uuid.to_string())
            .replace("{:read_id:}", &read.id.to_string())
            .replace("{:pair:}", "2");

        // Again optionize all these b/c we don't care if the writes fail for the msot part
        // write the forward read
        file.write_all(fwd_header.as_bytes()).ok();
        file.write_all("\n".as_bytes()).ok();
        file.write_all(&read.forward.sequence).ok();
        file.write_all("\n".as_bytes()).ok();
        file.write_all("+".as_bytes()).ok();
        file.write_all("\n".as_bytes()).ok();
        file.write_all(&util::encode_quality_scores(&read.forward.quality))
            .ok();
        file.write_all("\n".as_bytes()).ok();

        // write the reverse read if it's a PE read
        if read.reverse.is_some() {
            file.write_all(rev_header.as_bytes()).ok();
            file.write_all("\n".as_bytes()).ok();
            file.write_all(&read.reverse.as_ref().unwrap().sequence)
                .ok();
            file.write_all("\n".as_bytes()).ok();
            file.write_all("+".as_bytes()).ok();
            file.write_all("\n".as_bytes()).ok();
            file.write_all(&util::encode_quality_scores(
                &read.reverse.as_ref().unwrap().quality,
            ))
            .ok();
            file.write_all("\n".as_bytes()).ok();
        }
    }

    Ok(())
}
