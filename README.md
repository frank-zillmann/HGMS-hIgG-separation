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

# Build
./scripts/build.sh

# Run simulation
./install/bin/HGMS_hIgG_separation_main_run

# Run hyperparameter search
./install/bin/HGMS_hIgG_separation_hyperparameter_grid_search
```

## Built with FS³
This project uses [FS³](https://github.com/frank-zillmann/FS3) - a Fast and Flexible Framework for Simple Simulations of Separation-Processes.

See [FS3 README](external/FS3/README.md) for detailed setup instructions.
