use std::path::PathBuf;
use std::process::{exit, Command};

// This file contains code to run test files packaged with the tool in its CARGO_MANIFEST_DIR.
// This file is unrelated to unit tests, which are found in the same file as the functionality they
// test.

fn get_test_dir() -> PathBuf {
    let mut tpath = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    tpath.push("test");
    tpath
}

fn get_test_file(filename: &str) -> PathBuf {
    let mut testpath = PathBuf::new();
    testpath.push(filename);
    testpath
}

fn run_tests(runner: &str) {
    let test_runner = get_test_file(runner);

    let status = Command::new("sh")
        .current_dir(get_test_dir())
        .arg(test_runner)
        .status()
        .expect("Failed to execute test file");

    if !status.success() {
        eprintln!("Test script failed with status: {status}");
        exit(status.code().unwrap_or(1));
    } else {
        println! {"Tests completed!"}
    }
}

pub fn run_dedup_tests() {
    run_tests("test_dedup.sh");
}

pub fn run_extract_tests() {
    run_tests("test_extract.sh");
}
