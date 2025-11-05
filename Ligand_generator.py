#!/usr/bin/env python3
"""
API for batch processing SMILES strings from a text file
using the existing SMILES_to_XYZ.py script.

Usage:
    python API_script.py input_smiles.txt folder_name
"""

import os
import sys
import subprocess
import re
from pathlib import Path

def create_output_folder(folder_name):
    """Create output folder if it doesn't exist"""
    try:
        Path(folder_name).mkdir(parents=True, exist_ok=True)
        print(f"Created output folder: {folder_name}")
        return True
    except Exception as e:
        print(f"Error creating folder {folder_name}: {e}")
        return False

def read_smiles_from_file(input_file):
    """Read SMILES strings from input file"""
    smiles_list = []
    try:
        with open(input_file, 'r') as f:
            for line in f:
                # Skip empty lines and comments
                line = line.strip()
                if line and not line.startswith('#'):
                    smiles_list.append(line)
        return smiles_list
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        return []

def generate_safe_filename(smiles, index):
    """Generate a safe filename from SMILES string"""
    # Remove or replace characters that might be problematic in filenames
    safe_name = re.sub(r'[\\/*?:"<>|]', "_", smiles)
    # Limit length to avoid filesystem issues
    if len(safe_name) > 50:
        safe_name = safe_name[:50] + f"_{index}"
    return f"{safe_name}.xyz"

def process_smiles_batch(input_file, output_folder):
    """Process all SMILES strings in the input file"""
    # Read SMILES from file
    smiles_list = read_smiles_from_file(input_file)
    if not smiles_list:
        print("No valid SMILES strings found in the input file.")
        return False
    
    print(f"Found {len(smiles_list)} SMILES strings to process.")
    
    # Create output folder
    if not create_output_folder(output_folder):
        return False
    
    # Process each SMILES string
    success_count = 0
    for i, smiles in enumerate(smiles_list, 1):
        try:
            print(f"Processing {i}/{len(smiles_list)}: {smiles}")
            
            # Generate output filename
            output_filename = generate_safe_filename(smiles, i)
            output_path = os.path.join(output_folder, output_filename)
            
            # Run the conversion script
            cmd = [
                sys.executable,  # Use the same Python interpreter
                "SMILES_to_XYZ.py",
                f'{smiles}',  # Quote the SMILES string
                output_path
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"  Success: {output_filename}")
                success_count += 1
            else:
                print(f"  Failed: {result.stderr}")
                
        except Exception as e:
            print(f"  Error processing {smiles}: {e}")
    
    print(f"\nProcessing complete. Successfully converted {success_count}/{len(smiles_list)} molecules.")
    return success_count > 0

def main():
    """Main function"""
    if len(sys.argv) != 3:
        print("Usage: python API_script.py input_smiles.txt folder_name")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    
    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Input file not found: {input_file}")
        sys.exit(1)
    
    # Check if SMILES_to_XYZ.py exists
    if not os.path.isfile("SMILES_to_XYZ.py"):
        print("SMILES_to_XYZ.py not found in the current directory.")
        sys.exit(1)
    
    # Process the batch
    success = process_smiles_batch(input_file, output_folder)
    
    if success:
        print(f"\nXYZ files have been saved to: {output_folder}")
    else:
        print("\nBatch processing completed with errors.")
        sys.exit(1)

if __name__ == "__main__":
    main()
