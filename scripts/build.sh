#!/bin/bash
# Build script for HGMS_hIgG_separation
# Usage: ./scripts/build.sh [CMake options...]
# Example: ./scripts/build.sh -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/custom/path

set -e  # Exit on error

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Default directories
BUILD_DIR="$PROJECT_ROOT/build"
DEFAULT_INSTALL_DIR="$PROJECT_ROOT/install"

echo "=========================================="
echo "HGMS_hIgG_separation Build Script"
echo "=========================================="
echo "Project Root: $PROJECT_ROOT"
echo "Build Directory: $BUILD_DIR"
echo "CMake Args: $@"
echo "=========================================="

# Remove old build directory
if [ -d "$BUILD_DIR" ]; then
    echo "Removing old build directory..."
    rm -rf "$BUILD_DIR"
fi

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo ""
echo "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$DEFAULT_INSTALL_DIR" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    "$@" # Pass additional CMake args, overriding defaults if provided

# Build
echo ""
echo "Building..."
cmake --build . --parallel $(($(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) / 2))

# Install
echo ""
echo "Installing..."
cmake --install .

echo ""
echo "=========================================="
echo "Build completed successfully!"
echo "=========================================="
