#!/usr/bin/env bash

mkdir -p bin

cd RUMINA
cargo build --release
mv target/release/RUMINA ../bin/
cd ..

python3 -m venv python_env
source python_env/bin/activate
python3 -m pip install -r requirements.txt
deactivate

home_dir=$PWD
install_dir="$HOME/.cargo/bin"

if [ ! -d "$install_dir" ]; then
	mkdir -p "$install_dir"
fi

echo "#!/usr/bin/env bash" > "$install_dir"/rumina
echo "source $home_dir/python_env/bin/activate" >> "$install_dir"/rumina
echo "python3 $home_dir/main.py \"\$@\"" >> "$install_dir"/rumina
echo "deactivate" >> "$install_dir"/rumina
chmod +x "$install_dir"/rumina
