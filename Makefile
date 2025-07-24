build:
	VERSION=$$(git tag | tail -n 1); \
	sed -i .bak "s/0.0.0/$$VERSION/" Cargo.toml; \
	cargo build --release && cp target/release/rumina ~/.cargo/bin

test:
	cargo test -- --show-output
