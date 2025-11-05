#!/usr/bin/env python3
"""
SMILES to XYZ Converter with MMFF99 Minimization using OpenBabel C++ API

This script converts SMILES strings to 3D XYZ coordinates using OpenBabel's
UFF force field minimization without using pybel.

Requirements:
- openbabel 2.4.1 (python bindings)

Usage:
    python smiles_to_xyz_obapi.py "SMILES_STRING" output_filename.xyz
"""

import sys
import os
from openbabel import openbabel as ob

class SMILESConversionError(Exception):
    """Custom exception for SMILES conversion errors"""
    pass

class ForceFieldError(Exception):
    """Custom exception for force field minimization errors"""
    pass

class OBConversionError(Exception):
    """Custom exception for OpenBabel conversion errors"""
    pass

def validate_smiles(smiles_string):
    """
    Validate the SMILES string format using OpenBabel's built-in parser
    
    Args:
        smiles_string (str): Input SMILES string
        
    Returns:
        bool: True if valid, raises SMILESConversionError otherwise
    """
    if not smiles_string or not isinstance(smiles_string, str):
        raise SMILESConversionError("SMILES string cannot be empty or non-string")
    
    if len(smiles_string.strip()) == 0:
        raise SMILESConversionError("SMILES string cannot be empty")
    
    # Use OpenBabel's own parser to validate the SMILES
    try:
        conv = ob.OBConversion()
        mol = ob.OBMol()
        
        if not conv.SetInFormat("smi"):
            raise SMILESConversionError("OpenBabel SMILES parser not available")
        
        # Try to read the SMILES string - if it fails, it's invalid
        if not conv.ReadString(mol, smiles_string):
            raise SMILESConversionError(f"OpenBabel failed to parse SMILES: {smiles_string}")
        
        # Basic sanity check - molecule should have at least one atom
        if mol.NumAtoms() == 0:
            raise SMILESConversionError(f"SMILES string contains no atoms: {smiles_string}")
            
        return True
        
    except Exception as e:
        if isinstance(e, SMILESConversionError):
            raise
        raise SMILESConversionError(f"SMILES validation failed: {str(e)}")

def setup_force_field(mol, forcefield_type="UFF"):
    """
    Set up force field for molecular mechanics minimization
    
    Args:
        mol: OpenBabel OBMol object
        forcefield_type: Type of force field to use
        
    Returns:
        ob.OBForceField: Initialized force field object
        
    Raises:
        ForceFieldError: If force field setup fails
    """
    try:
        # Find the force field
        ff = ob.OBForceField.FindForceField(forcefield_type)
        if not ff:
            raise ForceFieldError(f"Force field '{forcefield_type}' not available")
        
        # Setup the force field with the molecule
        if not ff.Setup(mol):
            raise ForceFieldError(f"Failed to set up {forcefield_type} force field for the molecule")
        
        return ff
        
    except Exception as e:
        raise ForceFieldError(f"Force field setup failed: {str(e)}")

def generate_3d_coordinates(mol):
    """
    Generate 3D coordinates for a molecule using OpenBabel's builder
    
    Args:
        mol: OpenBabel OBMol object
        
    Returns:
        bool: True if successful
        
    Raises:
        OBConversionError: If 3D coordinate generation fails
    """
    try:
        # Create a builder for generating 3D coordinates
        builder = ob.OBBuilder()
        
        # Build 3D coordinates
        success = builder.Build(mol)
        if not success:
            raise OBConversionError("Failed to generate 3D coordinates")
        
        return True
        
    except Exception as e:
        raise OBConversionError(f"3D coordinate generation failed: {str(e)}")

def smiles_to_xyz(smiles_string, output_filename, max_iterations=1000):
    """
    Convert SMILES to XYZ with MMFF99 minimization using OpenBabel C++ API
    
    Args:
        smiles_string (str): Input SMILES string
        output_filename (str): Output XYZ filename
        max_iterations (int): Maximum iterations for minimization
        
    Returns:
        bool: True if successful, raises appropriate exceptions otherwise
        
    Raises:
        SMILESConversionError: If SMILES parsing fails
        ForceFieldError: If force field minimization fails
        OBConversionError: If OpenBabel conversion fails
        IOError: If file writing fails
    """
    try:
        # Validate input SMILES using OpenBabel's parser
        validate_smiles(smiles_string)
        
        # Create OpenBabel objects
        conv = ob.OBConversion()
        mol = ob.OBMol()
        
        # Set input format to SMILES
        if not conv.SetInFormat("smi"):
            raise OBConversionError("Failed to set input format to SMILES")
        
        # Read SMILES string
        if not conv.ReadString(mol, smiles_string):
            raise SMILESConversionError(f"Failed to parse SMILES string: {smiles_string}")
        
        # Add hydrogens
        mol.AddHydrogens()
        
        # Generate 3D coordinates
        generate_3d_coordinates(mol)
        
        # Set up force field for minimization
        forcefield = setup_force_field(mol, "UFF")
        
        # Perform energy minimization
        forcefield.ConjugateGradients(max_iterations)
        forcefield.GetCoordinates(mol)
        
        # Set output format to XYZ
        if not conv.SetOutFormat("xyz"):
            raise OBConversionError("Failed to set output format to XYZ")
        
        # Write to file
        if not conv.WriteFile(mol, output_filename):
            raise IOError(f"Failed to write to file: {output_filename}")
        
        # Print success information
        print(f"Successfully converted SMILES to XYZ: {output_filename}")
        print(f"Molecule: {smiles_string}")
        print(f"Atoms: {mol.NumAtoms()}")
        print(f"Bonds: {mol.NumBonds()}")
        
        # Optional: Print energy information
        energy = forcefield.Energy()
        print(f"Final energy: {energy:.4f} kcal/mol")
        
        return True
        
    except IOError as e:
        raise IOError(f"File operation error: {str(e)}")
    except Exception as e:
        # Re-raise already handled exceptions
        if isinstance(e, (SMILESConversionError, ForceFieldError, OBConversionError)):
            raise
        else:
            raise OBConversionError(f"Unexpected OpenBabel error: {str(e)}")

def check_openbabel_availability():
    """
    Check if OpenBabel is properly installed and available
    
    Returns:
        bool: True if OpenBabel is available
        
    Raises:
        ImportError: If OpenBabel cannot be imported
    """
    try:
        # Test basic OpenBabel functionality
        test_mol = ob.OBMol()
        test_conv = ob.OBConversion()
        
        # Check if UFF force field is available
        ff = ob.OBForceField.FindForceField("UFF")
        if not ff:
            raise ImportError("UFF force field not available in OpenBabel installation")
        
        return True
        
    except AttributeError as e:
        raise ImportError(f"OpenBabel installation appears incomplete: {str(e)}")
    except Exception as e:
        raise ImportError(f"OpenBabel check failed: {str(e)}")

def main():
    """Main function to handle command line arguments and execution"""
    
    # Check OpenBabel availability
    try:
        check_openbabel_availability()
    except ImportError as e:
        print(f"OpenBabel Error: {e}")
        print("Please ensure OpenBabel 2.4.1 is properly installed")
        sys.exit(1)
    
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python smiles_to_xyz_obapi.py \"SMILES_STRING\" output_filename.xyz")
        print("Example: python smiles_to_xyz_obapi.py \"CCO\" ethanol.xyz")
        sys.exit(1)
    
    smiles_string = sys.argv[1]
    output_filename = sys.argv[2]
    
    try:
        # Check if output directory exists
        output_dir = os.path.dirname(output_filename)
        if output_dir and not os.path.exists(output_dir):
            raise IOError(f"Output directory does not exist: {output_dir}")
        
        # Check file extension
        if not output_filename.lower().endswith('.xyz'):
            print("Warning: Output filename does not have .xyz extension")
        
        # Perform conversion
        success = smiles_to_xyz(smiles_string, output_filename)
        
        if success:
            print("Conversion completed successfully!")
            
    except SMILESConversionError as e:
        print(f"SMILES Conversion Error: {e}")
        sys.exit(2)
    except ForceFieldError as e:
        print(f"Force Field Error: {e}")
        sys.exit(3)
    except OBConversionError as e:
        print(f"OpenBabel Conversion Error: {e}")
        sys.exit(4)
    except IOError as e:
        print(f"I/O Error: {e}")
        sys.exit(5)
    except Exception as e:
        print(f"Unexpected Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(6)

if __name__ == "__main__":
    main()
