#!/usr/bin/env python3
"""
Script to link two XYZ files by removing user-defined atoms and creating bonds
between atoms previously connected to the removed atoms.

Features:
- Loads two XYZ files using OpenBabel C++ bindings only (no pybel)
- Removes user-specified atoms (e.g., Au)
- Aligns coordinates of linker atoms and transposes second fragment
- Rotates second fragment to find least sterically hindered conformation
- Creates bonds between atoms that were connected to removed atoms
- Performs UFF minimization with constraints on the first molecule
- Saves the united structure as XYZ file

Usage: python link_xyz_files.py file1.xyz file2.xyz atom_to_remove constraint_coefficient
"""

import sys
import os
import math
import copy
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

def calculate_distance(coord1, coord2):
    """
    Calculate Euclidean distance between two coordinates.
    
    Args:
        coord1 (ob.vector3): First coordinate
        coord2 (ob.vector3): Second coordinate
    
    Returns:
        float: Distance between coordinates
    """
    dx = coord1.GetX() - coord2.GetX()
    dy = coord1.GetY() - coord2.GetY()
    dz = coord1.GetZ() - coord2.GetZ()
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def add_vectors(vec1, vec2):
    """
    Add two OpenBabel vector3 objects.
    
    Args:
        vec1 (ob.vector3): First vector
        vec2 (ob.vector3): Second vector
    
    Returns:
        ob.vector3: Sum of the two vectors
    """
    result = ob.vector3()
    result.SetX(vec1.GetX() + vec2.GetX())
    result.SetY(vec1.GetY() + vec2.GetY())
    result.SetZ(vec1.GetZ() + vec2.GetZ())
    return result

def subtract_vectors(vec1, vec2):
    """
    Subtract two OpenBabel vector3 objects (vec1 - vec2).
    
    Args:
        vec1 (ob.vector3): First vector
        vec2 (ob.vector3): Second vector
    
    Returns:
        ob.vector3: Difference of the two vectors
    """
    result = ob.vector3()
    result.SetX(vec1.GetX() - vec2.GetX())
    result.SetY(vec1.GetY() - vec2.GetY())
    result.SetZ(vec1.GetZ() - vec2.GetZ())
    return result

def scale_vector(vec, scale):
    """
    Scale an OpenBabel vector3 by a scalar.
    
    Args:
        vec (ob.vector3): Vector to scale
        scale (float): Scaling factor
    
    Returns:
        ob.vector3: Scaled vector
    """
    result = ob.vector3()
    result.SetX(vec.GetX() * scale)
    result.SetY(vec.GetY() * scale)
    result.SetZ(vec.GetZ() * scale)
    return result

def dot_product(vec1, vec2):
    """
    Calculate dot product of two vectors.
    
    Args:
        vec1 (ob.vector3): First vector
        vec2 (ob.vector3): Second vector
    
    Returns:
        float: Dot product
    """
    return vec1.GetX() * vec2.GetX() + vec1.GetY() * vec2.GetY() + vec1.GetZ() * vec2.GetZ()

def cross_product(vec1, vec2):
    """
    Calculate cross product of two vectors.
    
    Args:
        vec1 (ob.vector3): First vector
        vec2 (ob.vector3): Second vector
    
    Returns:
        ob.vector3: Cross product
    """
    result = ob.vector3()
    result.SetX(vec1.GetY() * vec2.GetZ() - vec1.GetZ() * vec2.GetY())
    result.SetY(vec1.GetZ() * vec2.GetX() - vec1.GetX() * vec2.GetZ())
    result.SetZ(vec1.GetX() * vec2.GetY() - vec1.GetY() * vec2.GetX())
    return result

def vector_length(vec):
    """
    Calculate length of a vector.
    
    Args:
        vec (ob.vector3): Vector
    
    Returns:
        float: Vector length
    """
    return math.sqrt(vec.GetX()**2 + vec.GetY()**2 + vec.GetZ()**2)

def normalize_vector(vec):
    """
    Normalize a vector to unit length.
    
    Args:
        vec (ob.vector3): Vector to normalize
    
    Returns:
        ob.vector3: Normalized vector
    """
    length = vector_length(vec)
    if length == 0:
        return vec
    return scale_vector(vec, 1.0 / length)

def rotate_around_axis(point, axis_start, axis_direction, angle):
    """
    Rotate a point around an axis by a given angle.
    
    Args:
        point (ob.vector3): Point to rotate
        axis_start (ob.vector3): Point on the rotation axis
        axis_direction (ob.vector3): Direction vector of the rotation axis
        angle (float): Rotation angle in radians
    
    Returns:
        ob.vector3: Rotated point
    """
    # Normalize axis direction
    u = normalize_vector(axis_direction)
    
    # Translate point to origin relative to axis
    p = subtract_vectors(point, axis_start)
    
    # Rotation matrix components
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    
    # Rodrigues' rotation formula
    dot = dot_product(u, p)
    cross = cross_product(u, p)
    
    # Calculate rotated point
    rotated = add_vectors(
        add_vectors(
            scale_vector(p, cos_angle),
            scale_vector(cross, sin_angle)
        ),
        scale_vector(u, dot * (1 - cos_angle))
    )
    
    # Translate back
    return add_vectors(rotated, axis_start)

def calculate_inter_fragment_distance(mol1, mol2, connection_point_mol1, connection_point_mol2):
    """
    Calculate minimum distance between atoms of two fragments (excluding connection points).
    
    Args:
        mol1 (ob.OBMol): First molecule
        mol2 (ob.OBMol): Second molecule
        connection_point_mol1 (int): Atom index in mol1 that connects to mol2
        connection_point_mol2 (int): Atom index in mol2 that connects to mol1
    
    Returns:
        float: Minimum distance between atoms of different fragments
    """
    min_distance = float('inf')
    
    for i in range(1, mol1.NumAtoms() + 1):
        if i == connection_point_mol1:
            continue  # Skip connection point
        
        atom1 = mol1.GetAtom(i)
        coord1 = atom1.GetVector()
        
        for j in range(1, mol2.NumAtoms() + 1):
            if j == connection_point_mol2:
                continue  # Skip connection point
            
            atom2 = mol2.GetAtom(j)
            coord2 = atom2.GetVector()
            
            distance = calculate_distance(coord1, coord2)
            if distance < min_distance:
                min_distance = distance
    
    return min_distance

def optimize_fragment_rotation(mol1, mol2, connection_point_mol1, connection_point_mol2, rotation_axis_atom=None):
    """
    Rotate the second fragment around the connection point to maximize inter-fragment distance.
    
    Args:
        mol1 (ob.OBMol): First molecule (fixed)
        mol2 (ob.OBMol): Second molecule (to be rotated)
        connection_point_mol1 (int): Atom index in mol1 that connects to mol2
        connection_point_mol2 (int): Atom index in mol2 that connects to mol1
        rotation_axis_atom (int): Optional atom index in mol2 to define rotation axis
    
    Returns:
        ob.OBMol: Optimally rotated second molecule
    """
    print("Optimizing fragment rotation to minimize steric clashes...")
    
    # Get connection point coordinates
    cp1 = mol1.GetAtom(connection_point_mol1).GetVector()
    cp2 = mol2.GetAtom(connection_point_mol2).GetVector()
    
    # Define rotation axis
    if rotation_axis_atom is not None:
        # Use bond between connection point and specified atom as rotation axis
        axis_atom = mol2.GetAtom(rotation_axis_atom)
        axis_direction = subtract_vectors(axis_atom.GetVector(), cp2)
    else:
        # Find a suitable rotation axis (perpendicular to the connection vector)
        # First, find any atom in mol2 that's not the connection point
        for i in range(1, mol2.NumAtoms() + 1):
            if i != connection_point_mol2:
                axis_atom = mol2.GetAtom(i)
                axis_direction = subtract_vectors(axis_atom.GetVector(), cp2)
                if vector_length(axis_direction) > 0.1:  # Ensure it's not too close
                    break
        else:
            # Fallback: use arbitrary axis
            axis_direction = ob.vector3(1.0, 0.0, 0.0)
    
    # Normalize axis direction
    axis_direction = normalize_vector(axis_direction)
    
    # Make a copy of mol2 for rotation trials
    mol2_working = ob.OBMol(mol2)
    
    # Try different rotation angles and find the one with maximum inter-fragment distance
    best_angle = 0.0
    best_distance = calculate_inter_fragment_distance(mol1, mol2_working, connection_point_mol1, connection_point_mol2)
    
    print(f"Initial minimum inter-fragment distance: {best_distance:.3f} Å")
    
    n_steps = 36  # 10-degree steps
    for step in range(n_steps):
        angle = 2.0 * math.pi * step / n_steps
        
        # Create a fresh copy for this rotation
        mol2_trial = ob.OBMol(mol2)
        
        # Rotate all atoms in mol2 around the connection point
        for i in range(1, mol2_trial.NumAtoms() + 1):
            atom = mol2_trial.GetAtom(i)
            if i != connection_point_mol2:  # Don't rotate the connection point itself
                current_coord = atom.GetVector()
                rotated_coord = rotate_around_axis(current_coord, cp2, axis_direction, angle)
                atom.SetVector(rotated_coord)
        
        # Calculate inter-fragment distance for this rotation
        distance = calculate_inter_fragment_distance(mol1, mol2_trial, connection_point_mol1, connection_point_mol2)
        
        if distance > best_distance:
            best_distance = distance
            best_angle = angle
            mol2_working = ob.OBMol(mol2_trial)  # Keep the best rotated version
    
    print(f"Optimal rotation angle: {math.degrees(best_angle):.1f}°")
    print(f"Best minimum inter-fragment distance: {best_distance:.3f} Å")
    
    return mol2_working

def align_and_transpose_fragments(mol1, mol2, atoms_to_remove_mol1, atoms_to_remove_mol2):
    """
    Align linker atoms and transpose second fragment closer to first fragment.
    
    Args:
        mol1 (ob.OBMol): First molecule
        mol2 (ob.OBMol): Second molecule  
        atoms_to_remove_mol1 (list): Atom indices to remove from mol1
        atoms_to_remove_mol2 (list): Atom indices to remove from mol2
    
    Returns:
        tuple: (modified_mol2, connection_points) where connection_points are
               the atoms that will form the new bond
    """
    
    # Find linker atoms (neighbors of atoms to be removed)
    linker_atoms_mol1 = []
    for atom_idx in atoms_to_remove_mol1:
        atom = mol1.GetAtom(atom_idx)
        neighbors = get_atom_neighbors(atom)
        linker_atoms_mol1.extend(neighbors)
    
    linker_atoms_mol2 = []
    for atom_idx in atoms_to_remove_mol2:
        atom = mol2.GetAtom(atom_idx)
        neighbors = get_atom_neighbors(atom)
        linker_atoms_mol2.extend(neighbors)
    
    if len(linker_atoms_mol1) != 1 or len(linker_atoms_mol2) != 1:
        raise ValueError(f"Expected exactly one linker atom per molecule. Found: {len(linker_atoms_mol1)} in mol1, {len(linker_atoms_mol2)} in mol2")
    
    # Get the linker atoms
    linker1_idx = linker_atoms_mol1[0]
    linker2_idx = linker_atoms_mol2[0]
    
    linker1 = mol1.GetAtom(linker1_idx)
    linker2 = mol2.GetAtom(linker2_idx)
    
    # Get coordinates
    linker1_coord = linker1.GetVector()
    linker2_coord = linker2.GetVector()
    
    # Calculate translation vector to align linker atoms
    translation_vector = ob.vector3()
    translation_vector.SetX(linker1_coord.GetX() - linker2_coord.GetX())
    translation_vector.SetY(linker1_coord.GetY() - linker2_coord.GetY())
    translation_vector.SetZ(linker1_coord.GetZ() - linker2_coord.GetZ())
    
    # Apply a small additional translation to avoid exact overlap
    # Move linker2 slightly away from linker1 along the bond direction
    bond_distance = 1.5  # Reasonable bond length in Angstroms
    current_distance = calculate_distance(linker1_coord, linker2_coord)
    
    if current_distance > 0:
        # Normalize and scale translation
        scale_factor = (current_distance - bond_distance) / current_distance
        translation_vector = scale_vector(translation_vector, scale_factor)
    
    # Apply translation to all atoms in mol2
    for i in range(1, mol2.NumAtoms() + 1):
        atom = mol2.GetAtom(i)
        current_coord = atom.GetVector()
        new_coord = add_vectors(current_coord, translation_vector)
        atom.SetVector(new_coord)
    
    # Now optimize the rotation of mol2 to minimize steric clashes
    mol2 = optimize_fragment_rotation(mol1, mol2, linker1_idx, linker2_idx)
    
    connection_points = [(1, linker1_idx), (2, linker2_idx)]
    return mol2, connection_points

def setup_uff_minimization_with_constraints(mol, constrained_atoms, constraint_coeff):
    """
    Set up UFF minimization with proper constraints.
    
    Args:
        mol (ob.OBMol): Molecule to minimize
        constrained_atoms (list): List of atom indices to constrain (1-based indexing)
        constraint_coeff (float): Constraint strength coefficient
    
    Returns:
        ob.OBForceField: Force field object
    
    Raises:
        RuntimeError: If force field setup fails
    """
    ff = ob.OBForceField.FindForceField("UFF")
    if not ff:
        raise RuntimeError("UFF force field not available")
    
    print(f"Setting up constraints for {len(constrained_atoms)} atoms")
    
    # First, set up the molecule in the force field
    if not ff.Setup(mol):
        raise RuntimeError("Failed to set up force field")
    
    # Apply constraints - OpenBabel uses 1-based indexing for constraints
    constraints = ob.OBFFConstraints()
    
    # In OpenBabel Python bindings, AddAtomConstraint only takes the atom index
    # The constraint strength is handled differently
    for atom_idx in constrained_atoms:
        # Verify the atom exists
        atom = mol.GetAtom(atom_idx)
        if atom:
            # Set the atom as fixed - OpenBabel Python API only supports fixed atoms
            constraints.AddAtomConstraint(atom_idx)
            print(f"  Constraining atom {atom_idx} (element {atom.GetAtomicNum()})")
        else:
            print(f"  Warning: Atom index {atom_idx} not found in molecule")
    
    # Apply constraints to force field
    if constrained_atoms:
        ff.SetConstraints(constraints)
        print(f"Applied fixed constraints to {len(constrained_atoms)} atoms")
        print("Note: OpenBabel Python API only supports fixed atom constraints (infinite strength)")
    
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
    
    print(f"Loaded molecules: {mol1.NumAtoms()} atoms in {file1}, {mol2.NumAtoms()} atoms in {file2}")
    
    # Get atomic number of atom to remove
    atomic_num_to_remove = ob.GetAtomicNum(atom_to_remove)
    if atomic_num_to_remove == 0:
        raise ValueError(f"Invalid element symbol: {atom_to_remove}")
    
    # Find atoms to remove in both molecules
    atoms_to_remove_mol1 = []
    for i in range(1, mol1.NumAtoms() + 1):
        atom = mol1.GetAtom(i)
        if atom.GetAtomicNum() == atomic_num_to_remove:
            atoms_to_remove_mol1.append(i)
    
    atoms_to_remove_mol2 = []
    for i in range(1, mol2.NumAtoms() + 1):
        atom = mol2.GetAtom(i)
        if atom.GetAtomicNum() == atomic_num_to_remove:
            atoms_to_remove_mol2.append(i)
    
    print(f"Atoms to remove: {len(atoms_to_remove_mol1)} from mol1, {len(atoms_to_remove_mol2)} from mol2")
    
    # Validate we have exactly one atom to remove in each molecule
    if len(atoms_to_remove_mol1) != 1 or len(atoms_to_remove_mol2) != 1:
        raise ValueError(f"Expected exactly one {atom_to_remove} atom in each molecule. "
                        f"Found {len(atoms_to_remove_mol1)} in first, {len(atoms_to_remove_mol2)} in second.")
    
    # Align coordinates and transpose second fragment
    mol2, connection_points = align_and_transpose_fragments(
        mol1, mol2, atoms_to_remove_mol1, atoms_to_remove_mol2
    )
    
    print("Fragments aligned, transposed, and optimally rotated")
    
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
    
    print(f"United molecule created with {united_mol.NumAtoms()} atoms")
    
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
        new_bond_begin = atom_map[(mol1_id, atom1_idx)]
        new_bond_end = atom_map[(mol2_id, atom2_idx)]
        united_mol.AddBond(new_bond_begin, new_bond_end, 1)  # Single bond
        print(f"Created new bond between atoms {new_bond_begin} and {new_bond_end}")
    
    print(f"Total bonds in united molecule: {united_mol.NumBonds()}")
    
    # Get constrained atoms (from first molecule) - using NEW indices in united_mol
    constrained_atoms = []
    for (mol_id, old_idx), new_idx in atom_map.items():
        if mol_id == 1:  # Only constrain atoms from first molecule
            constrained_atoms.append(new_idx)
    
    print(f"Constraining {len(constrained_atoms)} atoms from first molecule")
    
    # Perform UFF minimization with constraints
    ff = setup_uff_minimization_with_constraints(united_mol, constrained_atoms, constraint_coeff)
    
    # Store coordinates before minimization for comparison
    coords_before = []
    for i in range(1, united_mol.NumAtoms() + 1):
        atom = united_mol.GetAtom(i)
        coords_before.append((atom.GetVector().GetX(), atom.GetVector().GetY(), atom.GetVector().GetZ()))
    
    # Perform energy minimization
    print("Starting energy minimization...")
    ff.ConjugateGradients(10000, 1.0e-6)  # 100000 steps, convergence threshold 1.0e-6
    ff.GetCoordinates(united_mol)
    print("Energy minimization completed")
    
    # Check how much constrained atoms moved
    max_movement = 0.0
    constrained_movements = []
    for atom_idx in constrained_atoms:
        atom = united_mol.GetAtom(atom_idx)
        coords_after = (atom.GetVector().GetX(), atom.GetVector().GetY(), atom.GetVector().GetZ())
        movement = calculate_distance(
            ob.vector3(coords_before[atom_idx-1][0], coords_before[atom_idx-1][1], coords_before[atom_idx-1][2]),
            ob.vector3(coords_after[0], coords_after[1], coords_after[2])
        )
        constrained_movements.append((atom_idx, movement))
        max_movement = max(max_movement, movement)
    
    print(f"Maximum movement of constrained atoms: {max_movement:.6f} Å")
    
    # Report constraint effectiveness
    if max_movement > 0.01:
        print("WARNING: Fixed atoms moved significantly. Constraints may not be working properly.")
    else:
        print("Fixed atoms maintained their positions well.")
    
    return united_mol

def main():
    """Main function to handle command line arguments and execute the script."""
    
    if len(sys.argv) != 5:
        print("Usage: python link_xyz_files.py file1.xyz file2.xyz atom_to_remove constraint_coefficient")
        print("Example: python link_xyz_files.py mol1.xyz mol2.xyz Au 100.0")
        print("Note: OpenBabel Python API only supports fixed atom constraints (infinite strength)")
        print("The constraint_coefficient parameter is accepted but not used for constraint strength.")
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
        output_filename = f"{base1}_{base2}_linked.xyz"
        
        # Write output file
        save_molecule_to_xyz(united_mol, output_filename)
        
        print(f"Successfully created united molecule: {output_filename}")
        print(f"Total atoms: {united_mol.NumAtoms()}")
        print(f"Total bonds: {united_mol.NumBonds()}")
        print(f"We STRONGLY recommend you to check the obtained geometry manually before any other operations!")
        
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
