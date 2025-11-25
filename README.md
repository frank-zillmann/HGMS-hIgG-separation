# HGMS-hIgG-separation

**Created by Frank Zillmann during his Bachelor's Thesis: High-Performance Process Simulations for Magnetic Separation of Biomolecules**

## Overview

Simulates a multi-step separation process for human Immunglobulin G purification using magnetic nanoparticles (MNPs) and Rotor-Stator High Gradient Magnetic Separation. The simulation includes:

- Convection and dispersion in many unit operations (pipes, chamber, volumes)
- Magnetic capture and slurry formation in rotor-stator process chamber
- Multi-stage loading, washing, and elution operations
- Chemical reactions (buffer systems, pH calculations)
- Langmuir adsorption kinetics of hIgG on MNPs

## Quick Start

```bash
# Clone with FS3 framework submodule
git clone --recursive https://github.com/frank-zillmann/HGMS-hIgG-separation.git
cd HGMS-hIgG-separation

# Build in Release mode (default)
./scripts/build.sh

# Build in Debug mode
./scripts/build.sh -DCMAKE_BUILD_TYPE=Debug

# Quick recompile after code changes
./scripts/recompile.sh

# Run simulation
./install/bin/HGMS_hIgG_separation_main_run

# Run hyperparameter search
./install/bin/HGMS_hIgG_separation_hyperparameter_grid_search
```

### VS Code Tasks

Press **`Ctrl + Shift + B`** (or **`Cmd + Shift + B`** on macOS) and select:

- **CMake: Release Build** - Optimized build (default)
- **CMake: Debug Build** - Debug build with AddressSanitizer
- **CMake: Recompile** - Quick rebuild without CMake reconfiguration
- **Run: Normal Run** - Execute main simulation
- **Run: Hyperparameter Search** - Execute hyperparameter grid search

## Built with FS³
This project uses [FS³ - Fast and Flexible Framework for Simple Simulations of Separation-Processes](https://github.com/frank-zillmann/FS3). See its README for setup and build instructions.