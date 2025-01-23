use std::env;
use std::path::Path;

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let version_path = Path::new(&manifest_dir).join("../VERSION");

    match std::fs::read_to_string(&version_path) {
        Ok(version) => {
            println!("cargo:rustc-env=PROJECT_VERSION={}", version.trim());
        }
        Err(_) => {
            println!("cargo:rustc-env=PROJECT_VERSION=0.0.0");
        }
    }
}
