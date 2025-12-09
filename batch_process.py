import os
import subprocess
import argparse
from pathlib import Path
import sys
import shutil

print("DEBUG: Script started", flush=True)

# Add 'yleaf' to python path to import internal modules
sys.path.append(str(Path(__file__).parent))
from yleaf import summary_logger

def check_index_fast(bam_path):
    """
    Quickly verifies if the index is usable using samtools idxstats.
    Returns True if index is good, False otherwise.
    """
    try:
        # 10s timeout for idxstats - it should be instant
        subprocess.check_call(
            ["samtools", "idxstats", str(bam_path)], 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL,
            timeout=10
        )
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False

def ensure_index(bam_path):
    """
    Ensures a valid index exists. 
    1. Checks if current index works (fast).
    2. If not, attempts to re-index.
    3. Verifies new index.
    """
    bam_path = Path(bam_path)
    
    # 1. Fast check existing index
    if check_index_fast(bam_path):
        return True
        
    print(f"  [Index] Index missing or invalid for {bam_path.name}. Re-indexing...")
    
    # Delete potentially bad indices
    for ext in ['.bam.bai', '.bai', '.cram.crai', '.crai']:
        idx = bam_path.with_suffix(ext)
        if idx.exists():
            try: idx.unlink()
            except OSError: pass
            
    # 2. Re-create index
    try:
        cmd = ["samtools", "index", str(bam_path)]
        # We allow time for indexing (it's heavy), but we don't wait for Yleaf to hang
        subprocess.check_call(cmd, timeout=600)
        print(f"  [Index] Re-indexing complete.")
    except subprocess.TimeoutExpired:
        print(f"  [Index] Failed: Indexing timed out (>10m).")
        return False
    except subprocess.CalledProcessError:
        print(f"  [Index] Failed: 'samtools index' returned error.")
        return False
        
    # 3. Verify new index
    if check_index_fast(bam_path):
        return True
    else:
        print(f"  [Index] Failed: Index created but 'idxstats' still fails.")
        return False

def run_yleaf(bam_path, output_base_dir):
    bam_path = Path(bam_path)
    output_dir = Path(output_base_dir) / f"output_{bam_path.stem}"
    # Don't create output_dir here, let Yleaf do it.
    # output_dir.mkdir(parents=True, exist_ok=True) 
    
    # Write log outside the target dir to prevent deletion by Yleaf
    log_file_path = output_dir.with_suffix('.log')
    print(f"--- Processing {bam_path.name} ---")
    
    # Use sys.executable and relative path to Yleaf.py
    cmd = [sys.executable, "Yleaf/yleaf/Yleaf.py", "-bam", str(bam_path), "-o", str(output_dir)]
    
    # Set PYTHONPATH so Yleaf can find its package
    env = os.environ.copy()
    project_root = str(Path(__file__).parent.absolute())
    env["PYTHONPATH"] = project_root + os.pathsep + env.get("PYTHONPATH", "")

    def run_cmd_realtime():
        print(f"  [1/2] Starting Yleaf subprocess...", end='', flush=True)
        start_time = time.time()
        # Safety net timeout still exists, but we rely on ensure_index to catch bad files
        MAX_DURATION = 900 
        
        with open(log_file_path, "w") as log:
            # Pass the modified environment
            process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
            
            while process.poll() is None:
                time.sleep(2) 
                if time.time() - start_time > MAX_DURATION:
                    print(f"  [Monitor] Safety Timeout ({MAX_DURATION}s). Killing...", end='')
                    process.kill()
                    print(" Killed.")
                    return -1
            
            print(" Done.")
            return process.returncode

    import time 

    # Run Yleaf
    exit_code = run_cmd_realtime()
    
    # Move log file into output directory if it was created
    final_log_path = output_dir / "full_terminal_output.log"
    if output_dir.exists() and output_dir.is_dir():
        try:
            shutil.move(str(log_file_path), str(final_log_path))
            log_file_path = final_log_path # Update ref for error message
        except Exception: pass
    
    if exit_code == 0:
        print(f"Successfully processed {bam_path.name}\n")
    elif exit_code == -1:
        print(f"Skipping {bam_path.name} due to timeout.")
        summary_logger.log_failure(bam_path, "Skipped: Safety Timeout")
    else:
        print(f"Yleaf failed (code {exit_code}). See {log_file_path}")
        summary_logger.log_failure(bam_path, f"Failed: Yleaf Error {exit_code}")

def check_bam_integrity(bam_path):
    """
    Uses samtools quickcheck to verify BAM integrity.
    """
    print(f"Checking integrity of {bam_path.name}...")
    try:
        subprocess.check_call(["samtools", "quickcheck", str(bam_path)])
        return True
    except subprocess.CalledProcessError:
        print(f"ERROR: File is corrupted: {bam_path}")
        return False
    except FileNotFoundError:
        # samtools not found, assume ok
        return True

def main():
    parser = argparse.ArgumentParser(description="Batch process BAM files with Yleaf.")
    parser.add_argument("input_source", help="Text file with list of paths OR directory to scan (if -d is used)")
    parser.add_argument("-d", "--directory-mode", action="store_true", help="Treat input_source as a directory and scan recursively for .bam/.cram")
    parser.add_argument("-o", "--output_dir", default=".", help="Base directory for output folders")
    
    args = parser.parse_args()
    
    files_to_process = []
    
    if args.directory_mode:
        search_dir = Path(args.input_source)
        if not search_dir.exists():
            print(f"Directory not found: {search_dir}")
            sys.exit(1)
            
        print(f"Scanning {search_dir} for BAM/CRAM files...")
        # Recursive glob
        extensions = ['*.bam', '*.cram']
        for ext in extensions:
            files_to_process.extend(list(search_dir.rglob(ext)))
            
    else:
        file_list_path = Path(args.input_source)
        if not file_list_path.exists():
            print(f"File list not found: {file_list_path}")
            sys.exit(1)
            
        with open(file_list_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]
            files_to_process = [Path(line) for line in lines]
        
    print(f"Found {len(files_to_process)} files to process.")
    
    for file_path in files_to_process:
        
        if not file_path.exists():
            print(f"File not found, skipping: {file_path}")
            continue
            
        # 0. Integrity Check
        if not check_bam_integrity(file_path):
            print(f"Skipping corrupted file: {file_path.name}\n")
            summary_logger.log_failure(file_path, "Skipped: File Corrupted")
            continue

        # 1. Check/Ensure Index
        if not ensure_index(file_path):
            print(f"Skipping {file_path.name} due to invalid/missing index.")
            summary_logger.log_failure(file_path, "Skipped: Index Invalid")
            continue
            
        # 2. Run Yleaf
        run_yleaf(file_path, args.output_dir)

if __name__ == "__main__":
    main()
