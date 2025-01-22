use std::fs;
use std::path::Path;

fn main() {
    let version = fs::read_to_string("../VERSION").expect("Failed to read VERSION file");
    println!("cargo:rustc-env=CARGO_PKG_VERSION={}", version.trim());

    let cargo_toml_path = Path::new("Cargo.toml");
    let cargo_toml_content =
        fs::read_to_string(cargo_toml_path).expect("Failed to read Cargo.toml");

    let modified_cargo_toml = cargo_toml_content.replace(
        r#"version = "0.0.0""#,
        &format!(r#"version = "{}""#, version),
    );

    fs::write(cargo_toml_path, modified_cargo_toml).expect("Failed to write modified Cargo.toml");
}
