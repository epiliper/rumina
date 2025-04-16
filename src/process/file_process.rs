use crate::args::Args;

pub trait FileProcess {
    fn init_from_args(args: &Args, bam_file_path: &String, bam_file_name: &String) -> Self;
    fn process(self);
}
