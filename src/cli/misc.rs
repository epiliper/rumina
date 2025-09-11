use crate::cli::DedupArgs;
use colored::Colorize;

const LOGO: &str = r#"

██████╗ ██╗   ██╗███╗   ███╗██╗███╗   ██╗ █████╗ 
██╔══██╗██║   ██║████╗ ████║██║████╗  ██║██╔══██╗
██████╔╝██║   ██║██╔████╔██║██║██╔██╗ ██║███████║
██╔══██╗██║   ██║██║╚██╔╝██║██║██║╚██╗██║██╔══██║
██║  ██║╚██████╔╝██║ ╚═╝ ██║██║██║ ╚████║██║  ██║
╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝

"#;

const DIVIDER: &str = "=========================";

pub fn print_logo() {
    println!("{}", LOGO.bright_blue())
}

pub fn print_init(args: &DedupArgs) {
    println!("{}", "INIT".purple());
    println!("{}", DIVIDER.purple());
    println!("{args}");
}

pub fn print_file_info(file_name: &String, cur_file_num: usize, num_files: usize) {
    println!(
        "{} {} {}: {}",
        format!("File {}", cur_file_num).cyan(),
        "of".cyan(),
        format!("{}", num_files).cyan(),
        file_name
    );
    println!("{}", DIVIDER.cyan());
}
