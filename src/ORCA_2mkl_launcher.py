#!/usr/bin/env python3
"""
Script to walk through directories, run orca_2mkl on .gbw files, 
and then run Multiwfn on the generated .molden.input files.

Usage: python run_orca_2mkl.py path_to_orca_2mkl_binary path_to_Multiwfn_binary [parent_folder]
Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl /home/user/Multiwfn/Multiwfn
Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl /home/user/Multiwfn/Multiwfn /path/to/parent/folder
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

def find_molden_files(parent_folder):
    """
    Find all .molden.input files in subdirectories of the parent folder.
    
    Args:
        parent_folder (str): Path to parent folder
    
    Returns:
        list: List of tuples (directory_path, molden_file_path)
    """
    parent_path = Path(parent_folder)
    
    if not parent_path.exists():
        return []
    
    molden_files = []
    
    # Search for both patterns: *.molden.input and molden.input
    for molden_file in parent_path.rglob("*.molden.input"):
        if molden_file.is_file():
            molden_files.append((molden_file.parent, molden_file))
    
    # Also search for just molden.input (if created with different options)
    for molden_file in parent_path.rglob("molden.input"):
        if molden_file.is_file():
            # Check if we haven't already added this file
            if not any(mf[1] == molden_file for mf in molden_files):
                molden_files.append((molden_file.parent, molden_file))
    
    print(f"Found {len(molden_files)} .molden.input file(s) in {parent_folder}")
    return molden_files

def create_wfn_script(multiwfn_path, output_dir):
    """
    Create the wfn_commands.sh script with the correct Multiwfn path.
    
    Args:
        multiwfn_path (str): Path to Multiwfn binary
        output_dir (Path): Directory where to create the script
    
    Returns:
        Path: Path to the created script
    """
    script_content = f"""#!/bin/bash

export KMP_STACKSIZE=2000M

"{multiwfn_path}" "$1" << EOF | tee "${{1%.wfn}}.scr"
2
2
3
4
5
7
0
-4
4
6
0
-10
17
1
1
3
-4
3
7
2
1
-10
q
EOF
"""
    
    # Use absolute path for output directory
    output_dir_abs = output_dir.resolve()
    script_path = output_dir_abs / "wfn_commands.sh"
    
    print(f"DEBUG: Creating script at: {script_path}")
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make the script executable
    script_path.chmod(0o755)
    
    return script_path

def run_orca_2mkl(orca_2mkl_path, directory, gbw_file):
    """
    Run orca_2mkl on a specific .gbw file.
    
    Args:
        orca_2mkl_path (str): Path to orca_2mkl binary
        directory (Path): Directory containing the .gbw file
        gbw_file (Path): Path to the .gbw file
    
    Returns:
        tuple: (success, message, molden_created, molden_path)
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
        
        # Check if .molden.input was created
        molden_file = directory / f"{base_name}.molden.input"
        molden_created = molden_file.exists()
        
        if result.returncode == 0:
            if molden_created:
                return True, f"Success: {base_name}.gbw -> {base_name}.molden.input created", True, molden_file
            else:
                # Try alternative naming
                alt_molden_file = directory / "molden.input"
                if alt_molden_file.exists():
                    return True, f"Success: {base_name}.gbw -> molden.input created", True, alt_molden_file
                else:
                    return True, f"Success but no .molden.input created for {base_name}.gbw", False, None
        else:
            return False, f"Failed (return code {result.returncode}): {result.stderr}", False, None
            
    except subprocess.TimeoutExpired:
        os.chdir(original_cwd)
        return False, "Timeout after 5 minutes", False, None
    except Exception as e:
        # Ensure we return to original directory even if error occurs
        if 'original_cwd' in locals():
            os.chdir(original_cwd)
        return False, f"Error: {str(e)}", False, None

def run_multiwfn(script_path, directory, molden_file):
    """
    Run the wfn_commands.sh script on a .molden.input file.
    
    Args:
        script_path (Path): Path to wfn_commands.sh script
        directory (Path): Directory containing the .molden.input file
        molden_file (Path): Path to the .molden.input file
    
    Returns:
        tuple: (success, message)
    """
    try:
        # Change to the directory containing the .molden.input file
        original_cwd = os.getcwd()
        os.chdir(directory)
        
        # Get absolute path to the script (from original location)
        script_abs_path = script_path.resolve()
        
        # Get just the filename for the molden file
        molden_filename = molden_file.name
        
        print(f"Running: bash {script_abs_path} {molden_filename}")
        print(f"Current directory: {Path.cwd()}")
        print(f"Script exists: {script_abs_path.exists()}")
        print(f"Molden file exists: {(directory / molden_filename).exists()}")
        print(f"Script path: {script_abs_path}")
        print(f"Molden file: {molden_filename}")
        
        cmd = ["bash", str(script_abs_path), molden_filename]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout for Multiwfn
        )
        
        # Return to original directory
        os.chdir(original_cwd)
        
        print(f"STDOUT: {result.stdout[:500]}...")  # Print first 500 chars of output
        if result.stderr:
            print(f"STDERR: {result.stderr[:500]}...")  # Print first 500 chars of error
        
        if result.returncode == 0:
            return True, f"Success: Multiwfn processed {molden_filename}"
        else:
            return False, f"Failed (return code {result.returncode}): {result.stderr}"
            
    except subprocess.TimeoutExpired:
        os.chdir(original_cwd)
        return False, "Timeout after 10 minutes"
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

def validate_multiwfn_path(multiwfn_path):
    """
    Validate the Multiwfn binary path.
    
    Args:
        multiwfn_path (str): Path to Multiwfn binary
    
    Returns:
        bool: True if valid, False otherwise
    """
    multiwfn_binary = Path(multiwfn_path)
    
    if not multiwfn_binary.exists():
        print(f"Error: Multiwfn binary not found: {multiwfn_path}")
        return False
    
    if not os.access(multiwfn_path, os.X_OK):
        print(f"Error: Multiwfn binary is not executable: {multiwfn_path}")
        return False
    
    return True

def main():
    if len(sys.argv) not in [3, 4]:
        print("Usage: python run_orca_2mkl.py path_to_orca_2mkl_binary path_to_Multiwfn_binary [parent_folder]")
        print("Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl /home/user/Multiwfn/Multiwfn")
        print("Example: python run_orca_2mkl.py /home/user/orca/orca_2mkl /home/user/Multiwfn/Multiwfn /path/to/parent/folder")
        print("\nIf parent_folder is not provided, current directory will be used.")
        sys.exit(1)
    
    orca_2mkl_path = sys.argv[1]
    multiwfn_path = sys.argv[2]
    
    # Validate binaries
    if not validate_orca_2mkl_path(orca_2mkl_path):
        sys.exit(1)
    
    if not validate_multiwfn_path(multiwfn_path):
        sys.exit(1)
    
    # Determine parent folder - ALWAYS convert to absolute path
    if len(sys.argv) == 4:
        parent_folder = Path(sys.argv[3]).resolve()  # Use resolve() to get absolute path
    else:
        parent_folder = Path.cwd().resolve()  # Use resolve() to get absolute path
    
    # Validate parent folder
    if not parent_folder.exists():
        print(f"Error: Parent folder not found: {parent_folder}")
        sys.exit(1)
    
    if not parent_folder.is_dir():
        print(f"Error: Path is not a directory: {parent_folder}")
        sys.exit(1)
    
    print(f"Using parent folder: {parent_folder}")
    
    # Create the wfn_commands.sh script in the parent folder
    print(f"Creating wfn_commands.sh script with Multiwfn path: {multiwfn_path}")
    wfn_script_path = create_wfn_script(multiwfn_path, parent_folder)
    print(f"Script created at: {wfn_script_path}")
    print(f"Script actually exists: {wfn_script_path.exists()}")
    
    # Find all .gbw files
    gbw_files = find_gbw_files(parent_folder)
    
    if not gbw_files:
        print("No .gbw files found. Exiting.")
        sys.exit(0)
    
    # Step 1: Process each .gbw file with orca_2mkl
    print(f"\nStep 1: Processing {len(gbw_files)} .gbw file(s) with orca_2mkl...")
    print("=" * 80)
    
    orca_successful = 0
    orca_failed = 0
    molden_created_count = 0
    molden_paths = []  # Track paths to created molden files
    
    for directory, gbw_file in gbw_files:
        success, message, molden_created, molden_path = run_orca_2mkl(orca_2mkl_path, directory, gbw_file)
        
        if success:
            print(f"✓ {message}")
            orca_successful += 1
            if molden_created and molden_path:
                molden_created_count += 1
                molden_paths.append(molden_path)
        else:
            print(f"✗ {message}")
            orca_failed += 1
    
    # Step 2: Find and process .molden.input files with Multiwfn
    print("\n" + "=" * 80)
    print(f"Step 2: Processing .molden.input files with Multiwfn...")
    
    # Now search for .molden.input files AFTER processing all .gbw files
    molden_files = find_molden_files(parent_folder)
    
    # Also use the paths we tracked during creation
    if molden_paths and not molden_files:
        print(f"Using tracked molden files from step 1")
        molden_files = [(path.parent, path) for path in molden_paths]
    
    if not molden_files:
        print("No .molden.input files found for Multiwfn processing.")
    else:
        multiwfn_successful = 0
        multiwfn_failed = 0
        
        for directory, molden_file in molden_files:
            success, message = run_multiwfn(wfn_script_path, directory, molden_file)
            
            if success:
                print(f"✓ {message}")
                multiwfn_successful += 1
            else:
                print(f"✗ {message}")
                multiwfn_failed += 1
        
        # Multiwfn summary
        print("\n" + "=" * 80)
        print(f"Multiwfn Summary: {multiwfn_successful} successful, {multiwfn_failed} failed")
    
    # Overall summary
    print("=" * 80)
    print("OVERALL SUMMARY")
    print("-" * 80)
    print(f"orca_2mkl: {orca_successful} successful, {orca_failed} failed")
    print(f".molden.input files created: {molden_created_count}")
    print(f".molden.input files found for Multiwfn: {len(molden_files)}")
    if molden_files:
        print(f"Multiwfn: {multiwfn_successful} successful, {multiwfn_failed} failed")
    print(f"All files processed in: {parent_folder}")
    print(f"wfn_commands.sh script at: {wfn_script_path}")
    
    if orca_failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
