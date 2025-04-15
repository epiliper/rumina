use bio::io::fastq;
use rust_htslib::bam;

pub trait Record {
    fn seq(&self) -> String;
    fn qual(&self) -> &[u8];
    fn get_umi(&self, separator: &String) -> String;
    fn qname(&self) -> &[u8];
}

pub type BamRecord = bam::Record;

impl Record for BamRecord {
    fn seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq().encoded.to_vec()) }
    }

    fn get_umi(&self, separator: &String) -> String {
        unsafe {
            std::str::from_utf8_unchecked(self.qname())
                .rsplit_once(separator)
                .expect("ERROR: failed to get UMI from read QNAME. Check --separator. Exiting.")
                .1
                .to_string()
        }
    }

    fn qual(&self) -> &[u8] {
        self.qual()
    }

    fn qname(&self) -> &[u8] {
        self.qname()
    }
}

pub type FastqRecord = fastq::Record;

impl Record for fastq::Record {
    fn seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.seq().to_vec()) }
    }

    fn get_umi(&self, separator: &String) -> String {
        self.id()
            .rsplit_once(separator)
            .expect("ERROR: failed to get UMI from read QNAME: Check --separator. Exiting.")
            .1
            .to_string()
    }

    fn qual(&self) -> &[u8] {
        self.qual()
    }

    fn qname(&self) -> &[u8] {
        self.id().as_bytes()
    }
}
