#!/usr/bin/env python3
"""
Script to link two XYZ files by removing user-defined atoms and creating bonds
between atoms previously connected to the removed atoms.

Features:
- Loads two XYZ files using OpenBabel C++ bindings only (no pybel)
- Removes user-specified atoms (e.g., Au)
- Creates bonds between atoms that were connected to removed atoms
- Performs UFF minimization with constraints on the first molecule
- Saves the united structure as XYZ file

Usage: python link_xyz_files.py file1.xyz file2.xyz atom_to_remove constraint_coefficient
"""

import sys
import os
from openbabel import openbabel as ob

def load_molecule_from_xyz(filename):
    """
    Load a molecule from XYZ file using OpenBabel.
    
    Args:
        filename (str): Path to XYZ file
    
    Returns:
        ob.OBMol: Loaded molecule object
    
    Raises:
        RuntimeError: If file cannot be read or is empty
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Input file not found: {filename}")
    
    mol = ob.OBMol()
    conv = ob.OBConversion()
    
    if not conv.SetInFormat("xyz"):
        raise RuntimeError("XYZ format not supported by OpenBabel")
    
    if not conv.ReadFile(mol, filename):
        raise RuntimeError(f"Failed to read molecule from {filename}")
    
    if mol.NumAtoms() == 0:
        raise RuntimeError(f"No atoms found in {filename}")
    
    return mol

def save_molecule_to_xyz(mol, filename):
    """
    Save a molecule to XYZ file using OpenBabel.
    
    Args:
        mol (ob.OBMol): Molecule to save
        filename (str): Output filename
    """
    conv = ob.OBConversion()
    if not conv.SetOutFormat("xyz"):
        raise RuntimeError("XYZ output format not supported")
    
    if not conv.WriteFile(mol, filename):
        raise RuntimeError(f"Failed to write molecule to {filename}")

def get_atom_neighbors(atom):
    """
    Get all neighbors of an atom.
    
    Args:
        atom (ob.OBAtom): Atom object
    
    Returns:
        list: List of neighbor atom indices
    """
    neighbors = []
    for neighbor in ob.OBAtomAtomIter(atom):
        neighbors.append(neighbor.GetIdx())
    return neighbors

def setup_uff_minimization(mol, constrained_atoms, constraint_coeff):
    """
    Set up UFF minimization with constraints.
    
    Args:
        mol (ob.OBMol): Molecule to minimize
        constrained_atoms (list): List of atom indices to constrain
        constraint_coeff (float): Constraint strength coefficient
    
    Returns:
        bool: True if minimization was successful
    """
    ff = ob.OBForceField.FindForceField("UFF")
    if not ff:
        raise RuntimeError("UFF force field not available")
    
    # Set up the force field without constraints first
    if not ff.Setup(mol):
        raise RuntimeError("Failed to set up force field")
    
    # Apply constraints to specific atoms
    for atom_idx in constrained_atoms:
        # In OpenBabel 2.4.1, we need to set constraints differently
        # We'll use the FixAtom method which is available in older versions
        atom = mol.GetAtom(atom_idx)
        if atom:
            # Set the atom as fixed (fully constrained)
            # The constraint_coeff is used to determine if we want full or partial constraints
            if constraint_coeff >= 1000:  # Strong constraint - fix atom completely
                ff.SetFixAtom(atom_idx)
            else:  # Weaker constraint - not directly supported in older versions, so we'll fix completely
                ff.SetFixAtom(atom_idx)
    
    return ff

def link_xyz_files(file1, file2, atom_to_remove, constraint_coeff):
    """
    Link two XYZ files by removing specified atoms and creating new bonds.
    
    Args:
        file1 (str): Path to first XYZ file
        file2 (str): Path to second XYZ file
        atom_to_remove (str): Element symbol of atoms to remove (e.g., 'Au')
        constraint_coeff (float): Constraint coefficient for UFF minimization
    
    Returns:
        ob.OBMol: United molecule object
    
    Raises:
        ValueError: If connection points are invalid
        RuntimeError: If force field setup fails
    """
    
    # Load molecules
    mol1 = load_molecule_from_xyz(file1)
    mol2 = load_molecule_from_xyz(file2)
    
    # Get atomic number of atom to remove
    atomic_num_to_remove = ob.GetAtomicNum(atom_to_remove)
    if atomic_num_to_remove == 0:
        raise ValueError(f"Invalid element symbol: {atom_to_remove}")
    
    # Store connection points: (molecule_id, atom_idx)
    connection_points = []
    
    # Process first molecule - find atoms to remove and their neighbors
    atoms_to_remove_mol1 = []
    for i in range(1, mol1.NumAtoms() + 1):
        atom = mol1.GetAtom(i)
        if atom.GetAtomicNum() == atomic_num_to_remove:
            atoms_to_remove_mol1.append(i)
            # Store neighbors of the atom to be removed
            neighbors = get_atom_neighbors(atom)
            for neighbor_idx in neighbors:
                connection_points.append((1, neighbor_idx))
    
    # Process second molecule - find atoms to remove and their neighbors
    atoms_to_remove_mol2 = []
    for i in range(1, mol2.NumAtoms() + 1):
        atom = mol2.GetAtom(i)
        if atom.GetAtomicNum() == atomic_num_to_remove:
            atoms_to_remove_mol2.append(i)
            # Store neighbors of the atom to be removed
            neighbors = get_atom_neighbors(atom)
            for neighbor_idx in neighbors:
                connection_points.append((2, neighbor_idx))
    
    # Validate connection points
    if len(connection_points) != 2:
        raise ValueError(f"Expected exactly 2 connection points, found {len(connection_points)}. "
                        f"Check if each molecule has exactly one {atom_to_remove} atom with one bond.")
    
    # Create united molecule
    united_mol = ob.OBMol()
    
    # Map old atom indices to new indices
    atom_map = {}  # Maps (molecule_id, old_idx) to new_idx
    
    # Add atoms from first molecule (excluding removed atoms)
    for i in range(1, mol1.NumAtoms() + 1):
        if i not in atoms_to_remove_mol1:
            old_atom = mol1.GetAtom(i)
            new_atom = united_mol.NewAtom()
            new_atom.SetAtomicNum(old_atom.GetAtomicNum())
            new_atom.SetVector(old_atom.GetVector())
            atom_map[(1, i)] = new_atom.GetIdx()
    
    # Add atoms from second molecule (excluding removed atoms)
    for i in range(1, mol2.NumAtoms() + 1):
        if i not in atoms_to_remove_mol2:
            old_atom = mol2.GetAtom(i)
            new_atom = united_mol.NewAtom()
            new_atom.SetAtomicNum(old_atom.GetAtomicNum())
            new_atom.SetVector(old_atom.GetVector())
            atom_map[(2, i)] = new_atom.GetIdx()
    
    # Add bonds from first molecule
    for i in range(mol1.NumBonds()):
        bond = mol1.GetBond(i)
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx not in atoms_to_remove_mol1 and end_idx not in atoms_to_remove_mol1:
            if (1, begin_idx) in atom_map and (1, end_idx) in atom_map:
                united_mol.AddBond(atom_map[(1, begin_idx)], 
                                  atom_map[(1, end_idx)], 
                                  bond.GetBondOrder())
    
    # Add bonds from second molecule
    for i in range(mol2.NumBonds()):
        bond = mol2.GetBond(i)
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx not in atoms_to_remove_mol2 and end_idx not in atoms_to_remove_mol2:
            if (2, begin_idx) in atom_map and (2, end_idx) in atom_map:
                united_mol.AddBond(atom_map[(2, begin_idx)], 
                                  atom_map[(2, end_idx)], 
                                  bond.GetBondOrder())
    
    # Create new bond between connection points
    mol1_id, atom1_idx = connection_points[0]
    mol2_id, atom2_idx = connection_points[1]
    
    if (mol1_id, atom1_idx) in atom_map and (mol2_id, atom2_idx) in atom_map:
        united_mol.AddBond(atom_map[(mol1_id, atom1_idx)], 
                          atom_map[(mol2_id, atom2_idx)], 
                          1)  # Single bond
    
    # Get constrained atoms (from first molecule)
    constrained_atoms = []
    for (mol_id, old_idx), new_idx in atom_map.items():
        if mol_id == 1:  # Only constrain atoms from first molecule
            constrained_atoms.append(new_idx)
    
    # Perform UFF minimization with constraints
    ff = setup_uff_minimization(united_mol, constrained_atoms, constraint_coeff)
    
    # Perform energy minimization
    ff.ConjugateGradients(1000, 1.0e-6)  # 1000 steps, convergence threshold 1.0e-6
    ff.GetCoordinates(united_mol)
    
    return united_mol

def main():
    """Main function to handle command line arguments and execute the script."""
    
    if len(sys.argv) != 5:
        print("Usage: python link_xyz_files.py file1.xyz file2.xyz atom_to_remove constraint_coefficient")
        print("Example: python link_xyz_files.py mol1.xyz mol2.xyz Au 100.0")
        sys.exit(1)
    
    try:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        atom_to_remove = sys.argv[3]
        constraint_coeff = float(sys.argv[4])
        
        # Validate constraint coefficient
        if constraint_coeff <= 0:
            raise ValueError("Constraint coefficient must be positive")
        
        # Process the molecules
        united_mol = link_xyz_files(file1, file2, atom_to_remove, constraint_coeff)
        
        # Generate output filename
        base1 = os.path.splitext(os.path.basename(file1))[0]
        base2 = os.path.splitext(os.path.basename(file2))[0]
        output_filename = f"{base1}_{base2}.xyz"
        
        # Write output file
        save_molecule_to_xyz(united_mol, output_filename)
        
        print(f"Successfully created united molecule: {output_filename}")
        print(f"Total atoms: {united_mol.NumAtoms()}")
        print(f"Total bonds: {united_mol.NumBonds()}")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()