import os
import subprocess
import argparse
from pathlib import Path
import sys

def check_create_index(bam_path):
    """
    Checks if BAM index exists. If not, creates it.
    Returns True if index exists or was created, False on failure.
    """
    bam_path = Path(bam_path)
    # Check for .bam.bai or .bai
    if bam_path.with_suffix('.bam.bai').exists() or bam_path.with_suffix('.bai').exists():
        return True
    
    print(f"Index missing for {bam_path.name}. Indexing...")
    try:
        # Run samtools index
        cmd = ["samtools", "index", str(bam_path)]
        subprocess.check_call(cmd)
        print(f"Index created for {bam_path.name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to create index for {bam_path}: {e}")
        return False
    except FileNotFoundError:
        print("Error: samtools not found. Please install samtools.")
        return False

def run_yleaf(bam_path, output_base_dir):
    """
    Runs Yleaf on the BAM file.
    """
    bam_path = Path(bam_path)
    # Create output directory name based on file stem
    output_dir = Path(output_base_dir) / f"output_{bam_path.stem}"
    
    # Ensure output dir exists (Yleaf creates it, but we need it for log file if we want to put it there immediately, 
    # but Yleaf might delete it if force is used. So let's write log AFTER or store in variable)
    # Better: Write to a separate log file in the base dir or inside the output dir after creation.
    # Actually, Yleaf creates 'run.log' internally, but that only captures python logging.
    # To capture terminal output (stdout/stderr), we need to redirect.
    
    # Let's create the dir first
    output_dir.mkdir(parents=True, exist_ok=True)
    
    log_file_path = output_dir / "full_terminal_output.log"
    
    print(f"--- Processing {bam_path.name} ---")
    print(f"Logs will be saved to: {log_file_path}")
    
    cmd = [
        "Yleaf",
        "-bam", str(bam_path),
        "-o", str(output_dir)
        # -rg is auto-detected
        # -force is default True
    ]
    
    try:
        # Open log file
        with open(log_file_path, "w") as log_file:
            # Start process, piping stdout and stderr
            process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, # Merge stderr into stdout
                text=True, # Decode bytes to string
                bufsize=1 # Line buffered
            )
            
            # Read line by line
            for line in process.stdout:
                # Print to console
                print(line, end='')
                # Write to file
                log_file.write(line)
                
            # Wait for finish
            return_code = process.wait()
            
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, cmd)
            
        print(f"Successfully processed {bam_path.name}\n")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {bam_path.name}: {e}")
        print(f"Check log for details: {log_file_path}\n")

def main():
    parser = argparse.ArgumentParser(description="Batch process BAM files with Yleaf.")
    parser.add_argument("file_list", help="Text file containing list of BAM file paths (one per line)")
    parser.add_argument("-o", "--output_dir", default=".", help="Base directory for output folders")
    
    args = parser.parse_args()
    
    file_list_path = Path(args.file_list)
    if not file_list_path.exists():
        print(f"File list not found: {file_list_path}")
        sys.exit(1)
        
    with open(file_list_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        
    print(f"Found {len(lines)} files to process.")
    
    for file_path_str in lines:
        file_path = Path(file_path_str)
        
        if not file_path.exists():
            print(f"File not found, skipping: {file_path}")
            continue
            
        # 1. Check/Create Index
        if not check_create_index(file_path):
            print(f"Skipping {file_path.name} due to indexing failure.")
            continue
            
        # 2. Run Yleaf
        run_yleaf(file_path, args.output_dir)

if __name__ == "__main__":
    main()
