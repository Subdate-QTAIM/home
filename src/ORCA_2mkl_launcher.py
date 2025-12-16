#!/usr/bin/env python3
"""
Simple script to walk through directories and run orca_2mkl on .gbw files.

Usage: python run_orca_2mkl.py path_to_orca_2mkl_binary [parent_folder]
Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl
Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl /path/to/parent/folder
"""

import os
import sys
import subprocess
from pathlib import Path

def find_gbw_files(parent_folder):
    """
    Find all .gbw files in subdirectories of the parent folder.
    
    Args:
        parent_folder (str): Path to parent folder
    
    Returns:
        list: List of tuples (directory_path, gbw_file_path)
    """
    parent_path = Path(parent_folder)
    
    if not parent_path.exists():
        print(f"Error: Parent folder not found: {parent_folder}")
        return []
    
    if not parent_path.is_dir():
        print(f"Error: Path is not a directory: {parent_folder}")
        return []
    
    gbw_files = []
    
    for gbw_file in parent_path.rglob("*.gbw"):
        if gbw_file.is_file():
            gbw_files.append((gbw_file.parent, gbw_file))
    
    print(f"Found {len(gbw_files)} .gbw file(s) in {parent_folder}")
    return gbw_files

def run_orca_2mkl(orca_2mkl_path, directory, gbw_file):
    """
    Run orca_2mkl on a specific .gbw file.
    
    Args:
        orca_2mkl_path (str): Path to orca_2mkl binary
        directory (Path): Directory containing the .gbw file
        gbw_file (Path): Path to the .gbw file
    
    Returns:
        tuple: (success, message)
    """
    try:
        # Change to the directory containing the .gbw file
        original_cwd = os.getcwd()
        os.chdir(directory)
        
        # Get the base name without extension for the 'opt' argument
        base_name = gbw_file.stem
        
        # Run orca_2mkl opt -molden
        cmd = [orca_2mkl_path, base_name, "-molden"]
        
        print(f"Running: {orca_2mkl_path} {base_name} -molden in {directory}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        # Return to original directory
        os.chdir(original_cwd)
        
        if result.returncode == 0:
            return True, f"Success: {base_name}.gbw"
        else:
            return False, f"Failed (return code {result.returncode}): {result.stderr}"
            
    except subprocess.TimeoutExpired:
        os.chdir(original_cwd)
        return False, "Timeout after 5 minutes"
    except Exception as e:
        # Ensure we return to original directory even if error occurs
        if 'original_cwd' in locals():
            os.chdir(original_cwd)
        return False, f"Error: {str(e)}"

def validate_orca_2mkl_path(orca_2mkl_path):
    """
    Validate the orca_2mkl binary path.
    
    Args:
        orca_2mkl_path (str): Path to orca_2mkl binary
    
    Returns:
        bool: True if valid, False otherwise
    """
    orca_2mkl_binary = Path(orca_2mkl_path)
    
    if not orca_2mkl_binary.exists():
        print(f"Error: orca_2mkl binary not found: {orca_2mkl_path}")
        return False
    
    if not os.access(orca_2mkl_path, os.X_OK):
        print(f"Error: orca_2mkl binary is not executable: {orca_2mkl_path}")
        return False
    
    return True

def main():
    if len(sys.argv) not in [2, 3]:
        print("Usage: python run_orca_2mkl.py path_to_orca_2mkl_binary [parent_folder]")
        print("Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl")
        print("Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl /path/to/parent/folder")
        print("\nIf parent_folder is not provided, current directory will be used.")
        sys.exit(1)
    
    orca_2mkl_path = sys.argv[1]
    
    # Validate orca_2mkl binary
    if not validate_orca_2mkl_path(orca_2mkl_path):
        sys.exit(1)
    
    # Determine parent folder
    if len(sys.argv) == 3:
        parent_folder = Path(sys.argv[2])
    else:
        parent_folder = Path.cwd()
    
    # Validate parent folder
    if not parent_folder.exists():
        print(f"Error: Parent folder not found: {parent_folder}")
        sys.exit(1)
    
    if not parent_folder.is_dir():
        print(f"Error: Path is not a directory: {parent_folder}")
        sys.exit(1)
    
    print(f"Using parent folder: {parent_folder}")
    
    # Find all .gbw files
    gbw_files = find_gbw_files(parent_folder)
    
    if not gbw_files:
        print("No .gbw files found. Exiting.")
        sys.exit(0)
    
    # Process each .gbw file
    successful = 0
    failed = 0
    
    print(f"\nProcessing {len(gbw_files)} .gbw file(s)...")
    print("=" * 60)
    
    for directory, gbw_file in gbw_files:
        success, message = run_orca_2mkl(orca_2mkl_path, directory, gbw_file)
        
        if success:
            print(f"✓ {message}")
            successful += 1
        else:
            print(f"✗ {message}")
            failed += 1
    
    # Summary
    print("=" * 60)
    print(f"Summary: {successful} successful, {failed} failed")
    print(f"All files processed in: {parent_folder}")
    
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
