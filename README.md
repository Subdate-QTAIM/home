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
├── src/ # Source codSMILES_to_XYZ.py # SMILES to 3D structure conversion
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
> **Tip:**
[num_parallel_flows] - number of sumultaneous running queues. The number of threads for single calulation is defined on ORCA input file (opt_bader.inp)

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

### 7. Running DFTB+ Metadynamics Simulations
This section describes how to set up and run DFTB+ metadynamics simulations for theozyme systems, using the round-based distribution strategy generated by ```Bader_PCA_analysis.py```.
#### Prerequisites
Before running simulations, ensure you have:
- DFTB+ installed (recommended: version 22.1 nompi build via conda)
- Adequate computational resources (see [System Requirements](#system-requirements))
#### Workflow Overview
- Round Assignment: Systems are grouped into rounds based on chemical space distribution
- System Preparation: Each theozyme structure must be manually validated and prepared
- Configuration: DFTB+ and Plumed input files must be customized per system
- Execution: Simulations run sequentially by round, with parallel execution within rounds

#### Step-by-Step Instructions
1.  **Understanding Round Assignments**\
  The ```Bader_PCA_analysis.py``` script generates round assignments in a format similar to Example_rounds_assignment.txt:
  ```bash
  Round 1: Compound 9, Compound 8, Compound 5, Compound 45, ...
  Round 2: Compound 7, Compound 44, Compound 38, ...
  ...
  ```
  Key points:
  * Compounds are distributed across rounds to maximize chemical diversity
  * Each round should be completed before moving to the next
  * Within a round, systems can be run in parallel based on available resources
2.  **System Preparation and Validation**\
  For each system in a round:
  * Create a dedicated directory:
    ```bash
    mkdir Round1_Compound9
    cd Round1_Compound9
    ```
  * Copy the theozyme XYZ file (e.g., from ```examples/```) and verify its adequacy:
    - Check atom connectivity and bond orders
    - Validate geometry (no unrealistic bond lengths/angles)
    - Confirm the active site is properly oriented
  * Copy configuration files:
    ```bash
    cp ../examples/dftb_in.hsd 
    cp ../examples/plumed.dat 
    ```
3. **Customizing DFTB+ Input (```dftb_in.hsd```)**\
Edit ```dftb_in.hsd``` for each system:
```bash
Geometry = xyzFormat {
<<< "Put the name of your XYZ file here"
}

Driver = VelocityVerlet {

 Plumed = yes # if needed

 MovedAtoms = 1:m n:-1 #define your constrained/active atoms here (see Plumed/dftbplus documentation for details)
 TimeStep [fs] = 0.5
 Steps = 40000 #trajectory length
 KeepStationary = Yes
 MDRestartFrequency = 10
# ConvergentForcesOnly = Yes
 
 Thermostat = Berendsen { 
    Temperature [K] = 298.15 #define temperature here
    Timescale[fs] = 100.0 #choose wisely for your system!
 }

Hamiltonian = xTB {
  Method = "GFN1-xTB" 	#computational method
  
  Solvation = GeneralizedBorn { # GBSA(water) or define your XTB1-2 specific solvent parameters
  .....

Charge = -1 # total charge of the system
....
```
 * Critical checks:
   - System charge: Calculate and set the correct ```Charge``` value
   - XYZ filename: Update ```Geometry``` section with your file
   - Constraints: Identify atoms to fix (e.g., enzyme backbone)
   - Trajectory length: Adjust ```Steps``` based on system size (typically 1-5 million steps)
   - Solvation: Enable/configure if simulating in implicit solvent
4. **Configuring Plumed Input (```plumed.dat```)**\
Edit plumed.dat to define collective variables (CVs):
```bash
# Define CVs for metadynamics
d1: DISTANCE ATOMS=10,20
a1: ANGLE ATOMS=15,20,25
t1: TORSION ATOMS=5,10,15,20

# Metadynamics bias
metad: METAD ARG=d1,a1,t1 PACE=500 HEIGHT=1.2 SIGMA=0.1,0.2,0.3 FILE=HILLS

# Print CV values
PRINT ARG=d1,a1,t1,metad.bias STRIDE=100 FILE=COLVAR
```
* Validation required:
  - CV definitions: Ensure atom indices match your system
  - CV relevance: Choose CVs that describe the reaction coordinate
  - Metadynamics parameters: Adjust ```PACE```, ```HEIGHT```, ```SIGMA``` appropriately
  - File outputs: Verify output filenames and frequencies

5. **Setting Environment Variables**\
Before execution, set thread and memory parameters:
```bash
# Set number of threads (adjust X based on available cores)
export OMP_NUM_THREADS=16

# Set RAM per thread (adjust XX based on system size)
export OMP_STACKSIZE=4G
```
6. **Execution Command**\
Run DFTB+ with real-time execution monitoring:
```bash
dftb+ | tee target_system_output.log
```
7. **Monitoring and Management**
* Check progress: ```target_system_output -f output.log```
* Monitor COLVAR: ```tail -f COLVAR```
* Check ```HILLS``` file: Ensures Gaussian hills are being deposited
* Verify disk space: Trajectories can be large (10-100 GB)

#### Performance Considerations
Resource Requirements

| System Size	Recommended | Threads	| Estimated Runtime	| Memory/Thread |
| --- | --- | --- |--- |
| Small (50 atoms) |	8-12 |	3-7 days |	2-4 GB |
| Medium (200 atoms) |	12-16 |	1-2 weeks |	4-6 GB |
| Large (500+ atoms, theozymes) |	16-20 |	2-3 weeks |	6-8 GB |

* Parallel Execution Strategy
  - Within a round: Run multiple systems in parallel based on available cores
  - Across rounds: Complete all systems in Round N before starting Round N+1
  - Resource allocation: Balance between parallel systems and per-system threads
* Example for 64-core node:
  - Run 4 systems from Round 1 simultaneously
  - Allocate 16 threads to each system
  - Monitor memory usage to avoid swapping
#### Troubleshooting
* Common Issues:
  - "Out of Memory" errors:
    1. Reduce ```OMP_NUM_THREADS```
    2. Increase ```OMP_STACKSIZE```
    3. Check system size vs. available RAM
  - Plumed initialization failures:
    1. Verify Plumed is compiled with DFTB+
    2. Check ```plumed.dat``` syntax
    3. Ensure CV atom indices are valid
  - Slow convergence:
    1. Adjust metadynamics parameters (```PACE```, ```HEIGHT```)
    2. Verify force field parameters
    3. Reduce timestep if needed
#### Best Practices
 - Start small: Test with short trajectories (10,000 steps) before full runs
 - Checkpoint regularly: Use Plumed restart capabilities
 - Monitor convergence: Watch CV distributions in COLVAR file
 - Backup data: Regularly archive important files
 - Document changes: Keep notes of all parameter modifications

#### Notes
* Theozyme systems are computationally expensive: Plan for 2-3 weeks runtime on 15-20 CPU threads
* Manual validation is critical: Do not skip structure checking steps
* Round-based execution: Follow the assignment order for efficient chemical space exploration
* Resource management: Adjust the number of simultaneous runs based on your available computational resources\

For examples of complete setup, refer to the ```PMe3_Pd_PhBr.zip``` archive in the ```examples/``` directory.

### CONFIGURATION FILES

### 8. opt_bader.inp

ORCA input template for QTAIM analysis with relativistic corrections.

**Key Features (please choose wisely for your system):**
- TPSS functional with ZORA relativistic correction
- def2-TZVP basis set with SARC for heavy elements
- Tight optimization criterion
- Numerical frequencies

### 9. wfn_commands.txt

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
- ORCA ver. 5.0.3: Download from [ORCA official site](https://orcaforum.kofo.mpg.de/app.php/dlext/?cat=20)
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
