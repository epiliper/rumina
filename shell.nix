{ pkgs ? import <nixpkgs> {} }:
let 
  overrides = (builtins.fromTOML (builtins.readFile ./rust-toolchain.toml));
in

with pkgs; mkShell {
  buildInputs = [
    rustup
    llvmPackages.clang # importing so rustc can find the hts.h file
    libclang
    gnused
  ];

  RUSTC_VERSION = overrides.toolchain.channel;
}

