use colored::Colorize;
use std::fmt;

pub enum MergeResult {
    Discordant,
    NoMerge,
    Merge,
}

pub struct MergeReport {
    num_discordant: i32,
    num_unmerged: i32,
    num_merged: i32,
    pub num_inreads: i32,
    pub num_outreads: i32,
}

impl MergeReport {
    pub fn count(&mut self, merge_result: MergeResult) {
        match merge_result {
            MergeResult::Discordant => self.num_discordant += 1,
            MergeResult::NoMerge => self.num_unmerged += 1,
            MergeResult::Merge => self.num_merged += 1,
        }
    }

    pub fn new() -> Self {
        MergeReport {
            num_discordant: 0,
            num_unmerged: 0,
            num_merged: 0,
            num_inreads: 0,
            num_outreads: 0,
        }
    }
}

impl fmt::Display for MergeReport {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            {}: {}\n\
            ",
            "Discordant read pairs".yellow(),
            self.num_discordant,
            "Unmerged reads".yellow(),
            self.num_unmerged,
            "Merged read pairs".yellow(),
            self.num_merged,
            "Reads in".yellow(),
            self.num_inreads,
            "Reads out".yellow(),
            self.num_outreads
        )
    }
}
