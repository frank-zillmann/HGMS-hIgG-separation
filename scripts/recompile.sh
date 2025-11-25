#!/bin/bash
# Recompile script for FS3 (no CMake reconfiguration)
# Usage: ./scripts/recompile.sh

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
INSTALL_DIR="$PROJECT_ROOT/install"

echo "=========================================="
echo "HGMS_hIgG_separation Recompile Script"
echo "=========================================="

if [ ! -d "$BUILD_DIR" ]; then
    echo "Error: Build directory does not exist!"
    echo "Please run ./scripts/build.sh first."
    exit 1
fi

cd "$BUILD_DIR"

# Build
echo "Building..."
cmake --build . --parallel $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Install
echo ""
echo "Installing to $INSTALL_DIR..."
cmake --install .

echo ""
echo "=========================================="
echo "Recompile completed successfully!"
echo "=========================================="
