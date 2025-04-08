build:
	cargo build --release && cp target/release/rumina ~/.cargo/bin
test:
	cargo test -- --nocapture
