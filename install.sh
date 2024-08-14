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
install_dir="$HOME/.cargo/bin"

if [ ! -d "$install_dir" ]; then
	mkdir -p "$install_dir"
fi

# echo "#!/usr/bin/env bash" > "$install_dir"/rumina
echo "#!/usr/bin/env bash" > "$install_dir"/rumina
echo "source $home_dir/python_env/bin/activate" >> "$install_dir"/rumina
echo "python3.12 $home_dir/main.py \"\$@\"" >> "$install_dir"/rumina
echo "deactivate" >> "$install_dir"/rumina
chmod +x "$install_dir"/rumina
echo "-------------------------"
echo "RUMINA installed successfully!\nTry running "rumina -h" to see the help screen.\nRUMINA installation dir: $home_dir"
