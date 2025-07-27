{
  stdenv, 
  fetchFromGitHub,
  htslib,
  rustup,
  clang,
  libclang,
  gnused,
  git
}:


stdenv.mkDerivation {
  pname = "rumina";
  version = "0.9.81";

  src = fetchFromGitHub {
    owner = "epiliper";
    repo = "rumina";
    rev = "0.9.81";
    sha256 = "sha256-XtFlLresxuTG7vPYAd9b2VAqCIkJHDlUiSkxS2sYx5A=";
  };
  
  RUSTC_VERSION = (builtins.fromTOML(builtins.readFile ./rust-toolchain.toml)).toolchain.channel;

  buildInputs = [ git htslib rustup clang libclang gnused ];

  buildPhase = ''
  export HOME=$(pwd)
  make
  '';

  installPhase = ''
  runHook preInstall
  mkdir -p $out/bin
  cp rumina $out/bin
  runHook postInstall
  '';
}
