COMPUTATIONAL CHEMISTRY AUTOMATION SUITE
=========================================

A comprehensive collection of Python scripts and configuration files for 
automating computational chemistry workflows, including molecular structure 
generation, quantum chemistry calculations, and data processing.

OVERVIEW
--------
This repository provides automated tools for computational chemistry research:
- SMILES to 3D structure conversion with force field optimization
- Batch processing of molecular structures
- ORCA quantum chemistry calculations with parallel execution
- Molecular file manipulation and linking
- Wavefunction analysis with Multiwfn
- File format conversion utilities

SCRIPTS AND FILES
-----------------

1. SMILES_to_XYZ.py
-------------------
Converts SMILES strings to 3D XYZ coordinates using OpenBabel's UFF force field minimization.

Features:
- Validates SMILES strings using OpenBabel's parser
- Generates 3D coordinates with proper geometry
- Performs energy minimization with UFF force field
- Handles error checking and reporting

Usage:
python SMILES_to_XYZ.py "SMILES_STRING" output_filename.xyz

2. Ligand_generator.py
----------------------
Batch processes SMILES strings from a text file using SMILES_to_XYZ.py.

Features:
- Processes multiple SMILES strings from a file
- Creates organized output folders
- Generates safe filenames from SMILES strings
- Provides progress reporting

Usage:
python Ligand_generator.py input_smiles.txt folder_name

3. ORCA_launcher_multi.py
-------------------------
Organizes XYZ files and runs ORCA calculations with parallel execution.

Features:
- Automatically creates folders for each molecule
- Copies input templates to calculation directories
- Runs multiple calculations in parallel
- Monitors calculation progress and status
- Generates detailed reports

Usage:
python ORCA_launcher_multi.py path_to_input_template path_to_orca_binary [num_parallel_flows]

4. ORCA_2mkl_launcher.py
------------------------
Processes ORCA .gbw files to generate Molden format files.

Features:
- Recursively finds .gbw files in directories
- Runs orca_2mkl conversion with Molden output
- Handles multiple files efficiently

Usage:
python ORCA_2mkl_launcher.py path_to_orca_2mkl_binary parent_folder

5. XYZ_file_binder.py
---------------------
Links two XYZ files by removing specified atoms and creating new bonds.

Features:
- Removes user-defined atoms (e.g., linker atoms)
- Creates bonds between previously connected atoms
- Performs constrained UFF minimization
- Maintains structural integrity of specified fragments

Usage:
python XYZ_file_binder.py file1.xyz file2.xyz atom_to_remove constraint_coefficient

CONFIGURATION FILES
-------------------

6. opt_bader.inp
----------------
ORCA input template for QTAIM analysis with relativistic corrections.

Key Features:
- TPSS functional with ZORA relativistic correction
- def2-TZVP basis set with SARC for heavy elements
- Tight optimization and numerical frequencies
- Bader analysis compatible settings

7. wfn_commands.txt
-------------------
Multiwfn script for wavefunction analysis.

Analysis Includes:
- Electron density properties
- Bond critical points
- Atomic charges and volumes
- Molecular surface properties

QUICK START
-----------

Prerequisites:
- Python 3.6+
- OpenBabel 2.4.1 with Python bindings
- ORCA quantum chemistry package
- Multiwfn for wavefunction analysis

Basic Workflow:

1. Generate 3D structures from SMILES:
   python Ligand_generator.py molecules.txt output_structures

2. Run ORCA calculations:
   python ORCA_launcher_multi.py opt_bader.inp /path/to/orca 4

3. Convert output files:
   python ORCA_2mkl_launcher.py /path/to/orca_2mkl calculation_results

4. Analyze wavefunctions:
   ./wfn_commands.sh molecule.wfn

FILE STRUCTURE
--------------
.
├── SMILES_to_XYZ.py          # SMILES to 3D structure conversion
├── Ligand_generator.py       # Batch SMILES processing
├── ORCA_launcher_multi.py    # Parallel ORCA calculations
├── ORCA_2mkl_launcher.py     # GBW to Molden conversion
├── XYZ_file_binder.py        # Molecular linking tool
├── opt_bader.inp            # ORCA input template
├── wfn_commands.txt         # Multiwfn analysis script
└── README.md              # This file

CONFIGURATION
-------------

ORCA Settings:
Modify opt_bader.inp for different:
- Density functionals
- Basis sets
- Relativistic methods
- Convergence criteria

Multiwfn Settings:
Adjust wfn_commands.txt for specific:
- Analysis types
- Output formats
- Calculation parameters

OUTPUT FILES
------------
- XYZ files: 3D molecular structures
- ORCA output: Calculation results and properties
- Molden files: Visualization and analysis
- Multiwfn output: Electronic structure analysis

TROUBLESHOOTING
---------------

Common Issues:
1. OpenBabel import errors: Ensure OpenBabel Python bindings are properly installed
2. ORCA not found: Verify ORCA binary path and permissions
3. Memory issues: Adjust %pal nprocs in ORCA input files
4. Conversion failures: Check SMILES string validity and file permissions

Debug Mode:
Most scripts provide verbose output. For detailed debugging, run with:
python script.py --debug

LICENSE
-------
This project is licensed under the MIT License.

ACKNOWLEDGMENTS
---------------
- OpenBabel team for molecular file conversion
- ORCA development team for quantum chemistry software
- Multiwfn developers for wavefunction analysis tools
- The computational chemistry community for methodologies and best practices

NOTE: This software is intended for research purposes. Users should verify 
results and methodologies for their specific applications.
