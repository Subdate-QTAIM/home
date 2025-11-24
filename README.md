# SUbstrate Bader Descriptors for biocATalytic reactivity Estimation (Subdate) computational workflow

A comprehensive collection of Python scripts and configuration files for 
automating computational chemistry workflows, including molecular structure 
generation, quantum chemistry calculations, and data processing.

## Table of Contents
- [Overview](#overview)
- [Directory Structure](#directory-structure)
- [Scripts and Files](#scripts-and-files)
- [Configuration Files](#configuration-files)
- [Quick Start](#quick-start)
- [System Requirements](#system-requirements)
- [Troubleshooting](#troubleshooting)
- [Citing](#citing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Overview

This repository provides automated tools for computational chemistry research:
- SMILES to 3D structure conversion with force field optimization
- Batch processing of molecular structures
- ORCA quantum chemistry calculations with parallel execution
- Molecular file manipulation and linking
- Wavefunction analysis with Multiwfn
- File format conversion utilities

## Directory Structure

Source code files are located in the `src/` directory, and usage examples can be found in the `examples/` directory.
```bash
├── src/ # Source code directory
│ ├── SMILES_to_XYZ.py # SMILES to 3D structure conversion
│ ├── Ligand_generator.py # Batch SMILES processing
│ ├── ORCA_launcher_multi.py # Parallel ORCA calculations
│ ├── ORCA_2mkl_launcher.py # GBW to Molden conversion
│ ├── XYZ_file_binder.py # Molecular linking tool
│ ├── opt_bader.inp # ORCA input template
│ └── wfn_commands.txt # Multiwfn analysis script
├── examples/ # Example files and usage cases
└── README.md # This file
```

## Scripts and Files

All scripts are located in the `src/` directory.

### 1. SMILES_to_XYZ.py

Converts SMILES strings to 3D XYZ coordinates using OpenBabel's UFF force field minimization.

**Features:**
- Validates SMILES strings using OpenBabel's parser
- Generates 3D coordinates with proper geometry
- Performs energy minimization with UFF force field
- Handles error checking and reporting

**Usage:**
```bash
python src/SMILES_to_XYZ.py "SMILES_STRING" output_filename.xyz 
```
### 2. Ligand_generator.py

Batch processes SMILES strings from a text file using SMILES_to_XYZ.py.

**Features:**
- Processes multiple SMILES strings from a file
- Creates organized output folders
- Generates safe filenames from SMILES strings
- Provides progress reporting

**Usage:**
```bash
python src/Ligand_generator.py input_smiles.txt folder_name
```
### 3. ORCA_launcher_multi.py

Organizes XYZ files and runs ORCA calculations with parallel execution.

**Features:**
- Automatically creates folders for each molecule
- Copies input templates to calculation directories
- Allows to run multiple calculations in parallel
- Monitors calculation progress and status
- Generates detailed reports

**Usage:**
```bash
python src/ORCA_launcher_multi.py path_to_input_template path_to_orca_binary [num_parallel_flows]
```
### 4. ORCA_2mkl_launcher.py

Processes ORCA .gbw files to generate Molden format files.

**Features:**
- Recursively finds .gbw files in directories
- Runs orca_2mkl conversion with Molden output
- Handles multiple files efficiently

**Usage:**
```bash
python src/ORCA_2mkl_launcher.py path_to_orca_2mkl_binary parent_folder
```
### 5. XYZ_file_binder.py

Links two XYZ files by removing specified atoms and creating new bonds.

**Features:**
- Removes user-defined atoms (e.g., linker atoms)
- Creates bonds between previously connected atoms
- Performs constrained UFF minimization
- Maintains structural integrity of specified fragments

**Usage:**
```bash
python src/XYZ_file_binder.py file1.xyz file2.xyz atom_to_remove constraint_coefficient
```
### CONFIGURATION FILES

### 6. opt_bader.inp

ORCA input template for QTAIM analysis with relativistic corrections.

**Key Features (please choose wisely for your system):**
- TPSS functional with ZORA relativistic correction
- def2-TZVP basis set with SARC for heavy elements
- Tight optimization criterion
- Numerical frequencies

### 7. wfn_commands.txt

Multiwfn script for wavefunction analysis.

**Analysis includes:**
- Electron density properties
- Bond critical points
- Atomic charges and volumes

## Quick Start

**Prerequisites:**
- Python 3.6 or greater (compatible with all required packages). <ins>The code was created and tested on Python 3.8.19 version</ins>
- OpenBabel 2.4.1: Install via conda: ```bash conda install -c conda-forge openbabel==2.4.1```
  * [OpenBabel Anaconda package](https://anaconda.org/channels/conda-forge/packages/openbabel/overview)
- OpenMPI 4.1.1: Download from [OpenMPI](https://www.open-mpi.org/software/ompi/v4.1/) *- strictly recommended for faster computing both for ORCA and DFTB+*
- Multiwfn: Download from [Multiwfn](http://sobereva.com/multiwfn/)
- Additional Python packages (install via conda):
  * ```bash conda install -c conda-forge dftbplus==22.1``` - [DFTB+ Anaconda package](https://anaconda.org/channels/conda-forge/packages/dftbplus/overview) <ins>**nompi_h0154332_100**  build is extremely recommended</ins>
  * ```bash conda install -c conda-forge plotly==XXX``` - [Plotly Anaconda package](https://anaconda.org/channels/conda-forge/packages/plotly/overview) Version 5.11.0 is recommended
  * ```bash conda install -c conda-forge tsne``` - [TSNE Anaconda package](https://anaconda.org/channels/conda-forge/packages/tsne/overview)

## Basic Workflow:

1. Generate 3D structures from SMILES

2. Run ORCA calculations

3. Convert output gbw files into molden input files

4. Analyze wavefunctions

## System Requirements

### Computational Requirements

This workflow is designed for high-performance computing environments:

* **Minimum**: 20-30 nodes with 64GB RAM

* **Recommended**: 120+ nodes with 100GB+ RAM

* **Operating System**: Ubuntu 22.04 (the whole pipeline was developed, run, and tested on Ubuntu 22.04 only)

[!TIP]
The *ab initio* metaMD runtime and number of simultaneous calculaions should be chosen wisely in agreement with user's computational resources. Insufficient resources may cause "Out of Memory" crashes and very long runs of individual systems

## Troubleshooting
### Common Issues:
+ OpenBabel import errors: Ensure OpenBabel Python bindings are properly installed
+ Conversion failures: Check SMILES string validity and file permissions
+ ORCA not found: Verify ORCA binary path and permissions
+ Parallel computational issues: ensure OpenMPI is installed correctly
+ Multiwfn errors: check **settings.ini** file, define corresponding amount of memory

## CITING
Please, cite

**Subdate workflow**,
*Solovev et al.,* "Chemical Space Exploration of Reactive Substrates for Biocatalytic Reactions" (under review)

## LICENSE

This project is licensed under the MIT License.

## ACKNOWLEDGMENTS

- OpenBabel team for molecular file conversion
- ORCA development team for quantum chemistry software
- Multiwfn developers for wavefunction analysis tools
- The computational chemistry community for methodologies and best practices

**NOTE:** This software is intended for research purposes. Users should verify 
results and methodologies for their specific applications.
