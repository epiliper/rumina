build:
	VERSION=$$(git tag | tail -n 1); \
	cargo build --release && mv target/release/rumina .

test:
	cargo test -- --show-output
