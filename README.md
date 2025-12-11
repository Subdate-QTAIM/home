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
│ ├── Bader_PCA_analysis.py # QTAIM data PCA analysis and visualization
│ ├── opt_bader.inp # ORCA input template
│ └── wfn_commands.txt # Multiwfn analysis script
├── examples/ # Example files and usage cases
└── README.md # This file
```

## Scripts and configuraion files

All scripts are located in the `src/` directory.
> **Note:**
All Python scripts should be run from the directory where they are located!

### 1. SMILES_to_XYZ.py

Converts SMILES strings to 3D XYZ coordinates using OpenBabel's UFF force field minimization.

**Features:**
- Validates SMILES strings using OpenBabel's parser
- Generates 3D coordinates with proper geometry
- Performs energy minimization with UFF force field
- Handles error checking and reporting

**Usage:**
```bash
python SMILES_to_XYZ.py n1ccccc1 pyridine.xyz 
```
### 2. Ligand_generator.py

Batch processes SMILES strings from a text file using SMILES_to_XYZ.py.
> **Note:**
Please, make sure the SMILES_to_XYZ.py file located in the same directory as well!

**Features:**
- Processes multiple SMILES strings from a file
- Creates organized output folders
- Generates safe filenames from SMILES strings
- Provides progress reporting

**Usage:**
```bash
python Ligand_generator.py examples/phosphate_molecules_examples generated_mols_dir
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
python ORCA_launcher_multi.py examples/opt_bader.inp path_to_orca_binary [num_parallel_flows]
```
### 4. ORCA_2mkl_launcher.py

Processes ORCA .gbw files to generate Molden format files.

**Features:**
- Recursively finds .gbw files in directories
- Runs orca_2mkl conversion with Molden output
- Handles multiple files efficiently

**Usage:**
```bash
python ORCA_2mkl_launcher.py path_to_orca_2mkl_binary path_to_parent_folder_with_subdirs_containing_gbw_files
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
python XYZ_file_binder.py examples/A17_para_theozyme.xyz examples/Au_nitro.xyz Au 999
```

### 6. Bader_PCA_analysis.py

Performs principal component analysis (PCA) on QTAIM-derived Bader descriptors and assigns compounds to computational rounds using a grid-based distribution strategy.

**Features:**
- Loads QTAIM data from CSV files with European decimal format (comma as decimal separator)
- Performs MinMax scaling followed by PCA transformation
- Projects compounds onto a user-defined n_rows × m_cols grid in the PCA space
- Distributes compounds across rounds by selecting one compound from each populated grid cell per round
- Generates interactive Plotly visualizations (HTML) and static images (PNG) showing PCA results colored by round assignment
- Outputs detailed text files with round assignments and grid cell contents

**Usage:**
```bash
python Bader_PCA_analysis.py examples/QTAIM_data_example.csv 5 5
```

### CONFIGURATION FILES

### 7. opt_bader.inp

ORCA input template for QTAIM analysis with relativistic corrections.

**Key Features (please choose wisely for your system):**
- TPSS functional with ZORA relativistic correction
- def2-TZVP basis set with SARC for heavy elements
- Tight optimization criterion
- Numerical frequencies

### 8. wfn_commands.txt

Multiwfn script for wavefunction analysis.

**Analysis includes:**
- Electron density properties
- Bond critical points
- Atomic charges and volumes

## Examples files
The ```examples/``` directory contains several theozyme structures in XYZ format used in our original study, including models for rH BChe, S198C rHBChe, and the A17 catalytic antibody (see ```A17_para_theozyme.xyz```, ```BuChe_Cys_meta_theozyme.xyz```, ```BuChe_Ser_ortho_theozyme.xyz``` etc). This directory also includes examples of functional groups with modified linker atoms.
- For completeness, typical input files (```dftb_in.hsd``` for DFTB+ and ```plumed.dat``` for metadynamics) are provided.
- Additional files:
  * ```QTAIM_data_example.csv```: example dataset from a QTAIM analysis, used to generate the chemical space population plots 
  * ```Example_pca_visualization.png/.html```: typical 2D visualization output from the principal component analysis (PCA) of the QTAIM data.
  * ```Example_rounds_assignment.txt```:  corresponding file showing the assigned *Rounds* distribution from the data analysis
  * ```PMe3_Pd_PhBr.zip```: a ready-to-run archive for testing your DFTB+ installation. It contains input files for a DFTB+/metadynamics simulation of a small phosphine ligand system (PMe₃-Pd-PhBr), along with the required XYZ structure.

## Quick Start

**Prerequisites:**
- Python 3.6 or greater (compatible with all required packages). <ins>The code was created and tested on Python 3.8.19 version</ins>
- OpenBabel 3.1.1: Install via conda: ```conda install -c conda-forge openbabel==3.1.1```
  * [OpenBabel Anaconda package](https://anaconda.org/channels/conda-forge/packages/openbabel/overview)
- OpenMPI 4.1.1: Download from [OpenMPI website](https://www.open-mpi.org/software/ompi/v4.1/) *- strictly recommended for faster computing both for ORCA and DFTB+*
- ORCA ver. 5.0.3: Download from [ORCA official site](https://orcaforum.kofo.mpg.de/app.php/portal)
- Multiwfn: Download from [Multiwfn official site](http://sobereva.com/multiwfn/)
- Additional Python packages (install via conda):
  * ```conda install -c conda-forge dftbplus==22.1=nompi_h0154332_100``` - [DFTB+ Anaconda package](https://anaconda.org/channels/conda-forge/packages/dftbplus/overview) <ins>nompi_h0154332_100  build is extremely recommended</ins>. You can also define the number of threads using ```bash export OMP_NUM_THREADS=n``` command
  * ```conda install -c conda-forge plotly==5.24.1``` - [Plotly Anaconda package](https://anaconda.org/channels/conda-forge/packages/plotly/overview) Version 5.24.1 is recommended
  * ```conda install -c conda-forge tsne``` - [TSNE Anaconda package](https://anaconda.org/channels/conda-forge/packages/tsne/overview)
  * ```conda install -c conda-forge pandas==1.3.5``` - [Pandas Anaconda package](https://anaconda.org/channels/conda-forge/packages/pandas/overview)
  * ```conda install -c conda-forge scikit-learn=1.0.2``` - [Sklearn Anaconda package](https://anaconda.org/channels/anaconda/packages/scikit-learn/overview)
  * ```pip install -U kaleido``` - [Kaleido library](https://pypi.org/project/kaleido/1.0.0rc0/)
Conda environment yaml file could be found in ```/src/subdate.yaml```

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

> **Tip:**
The *ab initio* metaMD runtime and number of simultaneous calculations should be chosen wisely in agreement with user's computational resources. Insufficient resources may cause "Out of Memory" crashes and very long runtimes for individual systems

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
