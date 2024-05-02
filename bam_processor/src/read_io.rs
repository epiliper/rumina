use crate::bottomhash::BottomHashMap;
use crate::processor;
use crate::Grouper;
use bam::BamReader;
use bam::BamWriter;
use bam::Record;
use bam::RecordWriter;
use std::collections::HashMap;
use std::fs::File;

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

pub struct ChunkProcessor<'a> {
    pub separator: &'a String,
    pub reads_to_output: &'a mut Vec<Record>,
    pub outfile: BamWriter<File>,
}

impl<'a> ChunkProcessor<'a> {
    // run grouping on pulled reads
    // add tags to Records
    // output them to list for writing to bam
    pub fn output_reads(&mut self) {
        println! {"Writing {} reads", self.reads_to_output.len()};
        self.reads_to_output
            .iter()
            .for_each(|x| self.outfile.write(x).unwrap());
        self.reads_to_output.clear();
    }

    pub fn group_reads(
        &mut self,
        bottomhash: &mut BottomHashMap,
        grouper: &mut Grouper,
    ) {
        let max = bottomhash.bottom_dict.keys().len();

        for (i, position) in bottomhash.bottom_dict.values_mut().enumerate() {
            print! {"\rProcessing bundle {} of {}", i, max};
            let bundle = position.shift_remove(&0).unwrap();
            let umis = bundle
                .keys()
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let processor = processor::Processor { umis: &umis };

            let mut counts: HashMap<&String, i32> = HashMap::new();

            for umi in &umis {
                counts.entry(umi).or_insert(bundle[umi].count);
            }

            let groupies = processor.main_grouper(counts);
            grouper.tag_records(groupies, bundle, self.reads_to_output);
        }
        self.output_reads();
        bottomhash.bottom_dict.clear();
    }

    pub fn pull_read(
        &mut self,
        read: &Record,
        bottomhash: &mut BottomHashMap,
        separator: &String,
    ) {
        if read.flag().is_mapped() && read.flag().is_reverse_strand() {
            bottomhash.update_dict(
                &(&read.calculate_end() + 1),
                0,
                &get_umi(&read, separator),
                &read,
            );

            // if bottomhash.bottom_dict.keys().last().unwrap() % chunksize == 0 {
            //     println! {"last pos {}", bottomhash.bottom_dict.keys().last().unwrap()};
            //     Self::group_reads(self, bottomhash, grouper);
            // }
        // otherwise, use it's first position to reference
        } else if read.flag().is_mapped() {
            bottomhash.update_dict(&(&read.start() + 1), 0, &get_umi(&read, separator), &read);
            // if bottomhash.bottom_dict.keys().last().unwrap() % chunksize == 0 {
            //     println! {"last pos {}", bottomhash.bottom_dict.keys().last().unwrap()};
            //     Self::group_reads(self, bottomhash, grouper);
            // }
        }
    }

    pub fn process_chunks(
        &mut self,
        input_file: BamReader<File>,
        mut bottomhash: BottomHashMap,
        mut grouper: Grouper,
    ) {
        let mut counter = 0;

        for r in input_file {
            let read = &r.unwrap();
            Self::pull_read(
                self,
                read,
                &mut bottomhash,
                self.separator,
            );

            counter += 1;
            if counter % 100_000 == 0 {
                print!{"\rRead in {counter} reads" }
            }
        }
        println!{};

        Self::group_reads(self, &mut bottomhash, &mut grouper);
    }
}
