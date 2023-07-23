/**
 * file: fastq.rs
 * desc: Write simulated reads to a FASTQ file.
 */
use std::fs;
use std::io::Write;

use crate::genome;
use crate::simulate::SimulatedRead;
use crate::util;

use shared::util::bytes_to_string;

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
            .replace(
                "{:sequence_id:}",
                &bytes_to_string(&read.forward_metadata.sequence_id),
            )
            .replace(
                "{:start_position:}",
                &read.forward_metadata.start_pos.to_string(),
            )
            .replace(
                "{:end_position:}",
                &read.forward_metadata.end_pos.to_string(),
            )
            .replace(
                "{:reverse_complement:}",
                if read.forward_metadata.is_reverse_complement {
                    "t"
                } else {
                    "f"
                },
            )
            .replace("{:pair:}", "1");

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
            let rev_header = header_format
                .replace("{:genome_id:}", &uuid.to_string())
                .replace("{:read_id:}", &read.id.to_string())
                .replace(
                    "{:sequence_id:}",
                    &bytes_to_string(&read.reverse_metadata.as_ref().unwrap().sequence_id),
                )
                .replace(
                    "{:start_position:}",
                    &read
                        .reverse_metadata
                        .as_ref()
                        .unwrap()
                        .start_pos
                        .to_string(),
                )
                .replace(
                    "{:end_position:}",
                    &read.reverse_metadata.as_ref().unwrap().end_pos.to_string(),
                )
                .replace(
                    "{:reverse_complement:}",
                    if read
                        .reverse_metadata
                        .as_ref()
                        .unwrap()
                        .is_reverse_complement
                    {
                        "t"
                    } else {
                        "f"
                    },
                )
                .replace("{:pair:}", "2");

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
