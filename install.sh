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

echo "#!/usr/bin/env bash" > $HOME/.local/bin/rumina
echo "source $home_dir/python_env/bin/activate\n" >> $HOME/.local/bin/rumina
echo "python3.12 $home_dir/main.py \"\$@\"" >> $HOME/.local/bin/rumina
echo "deactivate" >> $HOME/.local/bin/rumina
chmod +x $HOME/.local/bin/rumina




