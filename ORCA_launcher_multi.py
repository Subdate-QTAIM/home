#!/usr/bin/env python3
"""
Script to organize XYZ files and run ORCA calculations with parallel execution.

Features:
- Finds all XYZ files in current directory
- Moves each XYZ file to its own folder named after the file
- Renames XYZ files to "opt.xyz" in their respective folders
- Copies a user-provided input template to each folder
- Runs ORCA calculations with user-defined parallel flows
- Provides detailed reporting on calculation status

Usage: python run_orca_calculations.py path_to_input_template path_to_orca_binary [num_parallel_flows]
"""

import os
import sys
import shutil
import subprocess
import time
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import queue

def find_xyz_files():
    """
    Find all XYZ files in the current directory.
    
    Returns:
        list: List of Path objects for XYZ files
    """
    current_dir = Path.cwd()
    xyz_files = list(current_dir.glob("*.xyz"))
    
    if not xyz_files:
        print("No XYZ files found in the current directory.")
        return []
    
    print(f"Found {len(xyz_files)} XYZ file(s):")
    for xyz_file in xyz_files:
        print(f"  - {xyz_file.name}")
    
    return xyz_files

def create_folders_and_move_files(xyz_files):
    """
    Create folders for each XYZ file and move files into them.
    Renames each XYZ file to "opt.xyz" in its folder.
    
    Args:
        xyz_files (list): List of Path objects for XYZ files
    
    Returns:
        dict: Dictionary mapping folder names to new file paths
    """
    folder_mapping = {}
    
    for xyz_file in xyz_files:
        # Get folder name (without extension)
        folder_name = xyz_file.stem
        folder_path = Path(folder_name)
        
        # Create folder if it doesn't exist
        folder_path.mkdir(exist_ok=True)
        
        # Move XYZ file to folder and rename to opt.xyz
        destination = folder_path / "opt.xyz"
        shutil.move(str(xyz_file), str(destination))
        
        folder_mapping[folder_name] = destination
        print(f"Moved and renamed {xyz_file.name} to {folder_name}/opt.xyz")
    
    return folder_mapping

def copy_input_template(folder_mapping, input_template_path):
    """
    Copy input template to each folder.
    
    Args:
        folder_mapping (dict): Dictionary mapping folder names to file paths
        input_template_path (str): Path to input template file
    
    Returns:
        bool: True if successful, False otherwise
    """
    input_path = Path(input_template_path)
    
    if not input_path.exists():
        print(f"Error: Input template file not found: {input_template_path}")
        return False
    
    if not input_path.is_file():
        print(f"Error: Input template path is not a file: {input_template_path}")
        return False
    
    for folder_name in folder_mapping.keys():
        folder_path = Path(folder_name)
        destination = folder_path / input_path.name
        
        try:
            shutil.copy2(str(input_path), str(destination))
            print(f"Copied input template to {folder_name}/")
        except Exception as e:
            print(f"Error copying input template to {folder_name}: {e}")
            return False
    
    return True

def check_orca_termination(out_file_path):
    """
    Check if ORCA calculation terminated normally.
    
    Args:
        out_file_path (Path): Path to ORCA output file
    
    Returns:
        tuple: (bool, str) - (True if normal termination, status message)
    """
    if not out_file_path.exists():
        return False, "Output file not found"
    
    try:
        with open(out_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Check for normal termination signals
        if "ORCA TERMINATED NORMALLY" in content:
            return True, "Normal termination"
        elif "****ORCA TERMINATED NORMALLY****" in content:
            return True, "Normal termination"
        elif "TOTAL RUN TIME" in content:
            return True, "Completed (run time found)"
        
        # Check for error signals
        if "ORCA finished by error termination" in content:
            return False, "Error termination"
        elif "ABNORMAL TERMINATION" in content:
            return False, "Abnormal termination"
        elif "SEVERE ERROR" in content:
            return False, "Severe error"
        elif "not converged" in content.lower():
            return False, "Convergence failed"
        
        # If we can't determine, check if file ends with common ORCA patterns
        if content.strip().endswith(('***', '***\n', '****', '****\n')):
            return True, "Likely normal termination"
        
        return False, "Unknown status - cannot determine termination"
        
    except Exception as e:
        return False, f"Error reading output file: {e}"

def run_single_orca_calculation(folder_name, orca_binary_path):
    """
    Run a single ORCA calculation in a specific folder.
    
    Args:
        folder_name (str): Name of the folder containing the calculation
        orca_binary_path (str): Path to ORCA binary
    
    Returns:
        tuple: (folder_name, bool, str) - (folder name, True if successful, status message)
    """
    folder_path = Path(folder_name)
    
    if not folder_path.exists():
        return (folder_name, False, f"Folder not found: {folder_name}")
    
    # Find input file
    input_files = list(folder_path.glob("*.inp"))
    if not input_files:
        # Try to find any file that might be an input file
        all_files = list(folder_path.iterdir())
        input_files = [f for f in all_files if f.is_file() and f.suffix in ['.inp', '.txt', '.in']]
    
    if not input_files:
        return (folder_name, False, "No input file found")
    
    input_file = input_files[0]
    
    # Prepare output filename - unique name based on folder name
    output_file = folder_path / f"{folder_name}.out"
    
    try:
        thread_id = threading.current_thread().name
        print(f"[FLOW {thread_id}] Starting calculation: {folder_name}")
        
        # Run ORCA command from the folder path without changing global directory
        cmd = [orca_binary_path, str(input_file.name)]
        
        # Ensure output directory exists
        folder_path.mkdir(exist_ok=True)
        
        with open(output_file, 'w') as outfile:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                cwd=str(folder_path)  # Run the process in the specific folder
            )
            
            # Wait for process to complete without streaming output to screen
            stdout, stderr = process.communicate()
            
            # Write all output to file
            if stdout:
                outfile.write(stdout)
            if stderr:
                outfile.write(f"\nSTDERR:\n{stderr}\n")
        
        return_code = process.returncode
        
        # Check if output file was created and has content
        if not output_file.exists():
            return (folder_name, False, "Output file was not created")
        
        if output_file.stat().st_size == 0:
            return (folder_name, False, "Output file is empty")
        
        # Check termination status
        success, message = check_orca_termination(output_file)
        
        if success:
            print(f"[FLOW {thread_id}] ✓ Completed: {folder_name}")
            return (folder_name, True, f"Success: {message}")
        else:
            print(f"[FLOW {thread_id}] ✗ Failed: {folder_name} - {message}")
            return (folder_name, False, f"Failed: {message} (return code: {return_code})")
            
    except FileNotFoundError as e:
        error_msg = f"File not found error: {e}. Check if ORCA binary is accessible."
        print(f"[FLOW {thread_id}] ✗ Error: {folder_name} - {error_msg}")
        return (folder_name, False, error_msg)
    except Exception as e:
        error_msg = f"Exception during calculation: {e}"
        print(f"[FLOW {thread_id}] ✗ Error: {folder_name} - {error_msg}")
        return (folder_name, False, error_msg)

def run_parallel_orca_calculations(folder_names, orca_binary_path, num_parallel_flows):
    """
    Run ORCA calculations in parallel with specified number of flows.
    
    Args:
        folder_names (list): List of folder names to process
        orca_binary_path (str): Path to ORCA binary
        num_parallel_flows (int): Number of parallel calculations to run
    
    Returns:
        dict: Dictionary mapping folder names to (success, message) tuples
    """
    print(f"\nRunning {len(folder_names)} calculations with {num_parallel_flows} parallel flows...")
    
    results = {}
    
    # Use ThreadPoolExecutor for parallel execution
    with ThreadPoolExecutor(max_workers=num_parallel_flows, 
                          thread_name_prefix="Flow") as executor:
        # Submit all tasks and create a mapping of future to folder name
        future_to_folder = {}
        for folder_name in folder_names:
            future = executor.submit(run_single_orca_calculation, folder_name, orca_binary_path)
            future_to_folder[future] = folder_name
        
        # Monitor progress and collect results
        completed = 0
        total = len(folder_names)
        
        print(f"Submitted {total} calculations to {num_parallel_flows} parallel flows")
        
        # Process results as they complete
        for future in as_completed(future_to_folder):
            folder_name = future_to_folder[future]
            completed += 1
            
            try:
                # Get the result from the future
                folder_name, success, message = future.result()
                results[folder_name] = (success, message)
                
                status_icon = "✓" if success else "✗"
                print(f"[PROGRESS] {status_icon} {completed}/{total}: {folder_name} - {message}")
                
            except Exception as e:
                print(f"[ERROR] Exception in calculation for {folder_name}: {e}")
                results[folder_name] = (False, f"Execution error: {e}")
    
    return results

def validate_parameters(input_template_path, orca_binary_path, num_parallel_flows):
    """
    Validate all input parameters.
    
    Args:
        input_template_path (str): Path to input template file
        orca_binary_path (str): Path to ORCA binary
        num_parallel_flows (int): Number of parallel flows
    
    Returns:
        bool: True if all parameters are valid
    """
    # Validate input template
    input_path = Path(input_template_path)
    if not input_path.exists():
        print(f"Error: Input template file not found: {input_template_path}")
        return False
    
    # Validate ORCA binary
    orca_path = Path(orca_binary_path)
    if not orca_path.exists():
        print(f"Error: ORCA binary not found: {orca_binary_path}")
        return False
    
    if not os.access(orca_binary_path, os.X_OK):
        print(f"Error: ORCA binary is not executable: {orca_binary_path}")
        return False
    
    # Validate number of parallel flows
    if num_parallel_flows < 1:
        print(f"Error: Number of parallel flows must be at least 1")
        return False
    
    if num_parallel_flows > 16:  # Reasonable upper limit
        print(f"Warning: Using {num_parallel_flows} parallel flows may be excessive")
        response = input("Do you want to continue? (y/n): ")
        if response.lower() != 'y':
            return False
    
    return True

def main():
    """Main function to handle command line arguments and execute the script."""
    
    if len(sys.argv) not in [3, 4]:
        print("Usage: python run_orca_calculations.py path_to_input_template path_to_orca_binary [num_parallel_flows]")
        print("Example: python run_orca_calculations.py template.inp /home/user/orca/orca")
        print("Example: python run_orca_calculations.py template.inp /home/user/orca/orca 4")
        sys.exit(1)
    
    input_template_path = sys.argv[1]
    orca_binary_path = sys.argv[2]
    
    # Set default number of parallel flows to 1 if not specified
    if len(sys.argv) == 4:
        try:
            num_parallel_flows = int(sys.argv[3])
        except ValueError:
            print("Error: Number of parallel flows must be an integer")
            sys.exit(1)
    else:
        num_parallel_flows = 1
    
    # Validate parameters
    if not validate_parameters(input_template_path, orca_binary_path, num_parallel_flows):
        sys.exit(1)
    
    try:
        # Step 1: Find XYZ files
        print("Step 1: Finding XYZ files...")
        xyz_files = find_xyz_files()
        
        if not xyz_files:
            print("No XYZ files to process. Exiting.")
            sys.exit(0)
        
        # Step 2: Create folders and move files
        print("\nStep 2: Creating folders and moving files...")
        folder_mapping = create_folders_and_move_files(xyz_files)
        
        # Step 3: Copy input template
        print("\nStep 3: Copying input template...")
        if not copy_input_template(folder_mapping, input_template_path):
            print("Error copying input template. Exiting.")
            sys.exit(1)
        
        # Step 4: Run ORCA calculations in parallel
        print("\nStep 4: Running ORCA calculations...")
        print(f"Total calculations to run: {len(folder_mapping)}")
        print(f"Number of parallel flows: {num_parallel_flows}")
        
        folder_names = list(folder_mapping.keys())
        results = run_parallel_orca_calculations(folder_names, orca_binary_path, num_parallel_flows)
        
        # Step 5: Generate report
        print("\n" + "="*80)
        print("CALCULATION SUMMARY REPORT")
        print("="*80)
        
        successful = [name for name, (success, _) in results.items() if success]
        failed = [name for name, (success, _) in results.items() if not success]
        
        print(f"Successful calculations: {len(successful)}/{len(results)}")
        for name in successful:
            print(f"  ✓ {name}")
        
        print(f"\nFailed calculations: {len(failed)}/{len(results)}")
        for name in failed:
            status = results[name][1]
            print(f"  ✗ {name} - {status}")
        
        if failed:
            print(f"\nSome calculations failed. Check the individual output files for details.")
            sys.exit(1)
        else:
            print(f"\nAll calculations completed successfully!")
        
    except KeyboardInterrupt:
        print("\n\nCalculation interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()