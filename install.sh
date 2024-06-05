#!/usr/bin/env bash

cd bam_processor
cargo build --release

cd ../multibam 
cargo build --release

cd ..
python3 -m venv python_env

source python_env/bin/activate 
python3.12 -m pip install -r requirements.txt

deactivate

home_dir=$PWD

echo "#!/usr/bin/env bash" > $HOME/.cargo/bin/rumina
echo "source $home_dir/python_env/bin/activate" >> $HOME/.cargo/bin/rumina
echo "python3.12 $home_dir/main.py \"\$@\"" >> $HOME/.cargo/bin/rumina
echo "deactivate" >> $HOME/.cargo/bin/rumina
chmod +x $HOME/.cargo/bin/rumina




