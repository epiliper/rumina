#!/usr/bin/env bash

# if [ ! -d $HOME/.local/bin/rumina ]; then
# 	mkdir -p $HOME/.local/bin/rumina
# fi

# cp * $HOME/.local/bin/rumina

# cd $HOME/.local/bin/rumina

cd bam_processor
cargo build --release

cd ../multibam 
cargo build --release

cd ..
python3 -m venv python_env

source python_env/bin/activate 
python3.12 -m pip install -r requirements.txt

deactivate

echo -e "$PWD/python_env/bin/activate\n" > $HOME/.local/bin/rumina
echo -e "python3.12 activate $PWD/main.py "$@"" >> $HOME/.local/bin/rumina
echo -e "deactivate" >> $HOME/.local/bin/rumina




