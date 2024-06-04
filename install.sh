#!/usr/bin/env bash

cp * $HOME/.local/bin/rumina

cd $HOME/.local/bin/rumina

cd bam_processor
cargo build --release

cd ../multibam 
cargo build --release

cd ..
python3 -m venv python_env

source python_env/bin/activate 
python3.12 -m pip install -r requirements.txt

deactivate


