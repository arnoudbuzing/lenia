#!/bin/bash
# Build the Rust Lenia library and deploy to the paclet
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
RUST_DIR="$PROJECT_DIR/lenia-rs"
PACLET_LIB_DIR="$PROJECT_DIR/Lenia/LibraryResources/MacOSX-ARM64"

echo "=== Building Rust Lenia library ==="
cd "$RUST_DIR"
cargo build --release

echo "=== Deploying to paclet ==="
mkdir -p "$PACLET_LIB_DIR"
cp "$RUST_DIR/target/release/liblenia_rs.dylib" "$PACLET_LIB_DIR/"

echo "=== Done ==="
ls -lh "$PACLET_LIB_DIR/liblenia_rs.dylib"
