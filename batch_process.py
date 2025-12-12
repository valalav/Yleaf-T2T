import os
import subprocess
import argparse
from pathlib import Path
import sys
import shutil
import time
import multiprocessing
from functools import partial

# Add 'yleaf' to python path to import internal modules
sys.path.append(str(Path(__file__).parent))
from yleaf import summary_logger, yleaf_constants, external_tools

def detect_reference_path(bam_path):
    """
    Detects reference genome path (T2T, hg38, hg19) from BAM/CRAM header.
    Returns Path object or None.
    """
    try:
        # Use external_tools to get header lines safely
        header_lines = external_tools.samtools_view_header(bam_path)
        
        LEN_HG19 = 249250621
        LEN_HG38 = 248956422
        LEN_T2T = 248387328
        
        for line in header_lines:
            if "SN:chr1" in line or "SN:1" in line:
                for part in line.split():
                    if part.startswith("LN:"):
                        length = int(part.split(":")[1])
                        if length == LEN_T2T: return yleaf_constants.T2T_FULL_GENOME
                        if length == LEN_HG38: return yleaf_constants.HG38_FULL_GENOME
                        if length == LEN_HG19: return yleaf_constants.HG19_FULL_GENOME
    except Exception:
        pass
    return None

def check_index_fast(bam_path, reference_path=None):
    """
    Quickly verifies if the index is usable using samtools idxstats.
    Returns True if index is good, False otherwise.
    """
    env = os.environ.copy()
    
    # For CRAM files, help samtools find the reference to avoid hanging on downloads
    if reference_path and str(bam_path).endswith('.cram'):
        env['SAMTOOLS_CRAM_REF'] = str(reference_path)
        # Some versions might check REF_PATH
        env['REF_PATH'] = str(reference_path)

    try:
        # 10s timeout for idxstats - it should be instant
        subprocess.check_call(
            ["samtools", "idxstats", str(bam_path)], 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL,
            timeout=10,
            env=env
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
    
    # Detect reference (vital for CRAM)
    ref_path = detect_reference_path(bam_path)
    
    # Check if index exists and is fresh
    if bam_path.suffix == '.cram':
        idx_path = bam_path.with_suffix('.cram.crai')
        alt_idx = bam_path.with_suffix('.crai')
    else:
        idx_path = bam_path.with_suffix('.bam.bai')
        alt_idx = bam_path.with_suffix('.bai')
        
    existing_index = None
    if idx_path.exists() and idx_path.stat().st_mtime >= bam_path.stat().st_mtime:
        existing_index = idx_path
    elif alt_idx.exists() and alt_idx.stat().st_mtime >= bam_path.stat().st_mtime:
        existing_index = alt_idx
        
    if existing_index:
        # If it's CRAM, trust the timestamp. Remote idxstats is flaky.
        if bam_path.suffix == '.cram':
            # print(f"  [Index] CRAM index found: {existing_index.name} (Skipping deep check)")
            return True
        
        # For BAM, fast check is cheap and reliable
        if check_index_fast(bam_path):
            return True
        
    print(f"  [Index] Index missing, outdated, or invalid for {bam_path.name}. Re-indexing...")
    if ref_path:
        print(f"  [Index] Using reference: {ref_path.name}")
    
    # Delete potentially bad indices
    for ext in ['.bam.bai', '.bai', '.cram.crai', '.crai']:
        idx = bam_path.with_suffix(ext)
        if idx.exists():
            try:
                idx.unlink()
            except OSError as e:
                print(f"  [Index] Warning: Could not remove old index {idx}: {e}")
            
    # 2. Re-create index
    try:
        cmd = ["samtools", "index"]
        
        # Prepare environment for CRAM reference resolution
        env = os.environ.copy()
        if bam_path.suffix == '.cram':
            if ref_path and ref_path.exists():
                env['SAMTOOLS_CRAM_REF'] = str(ref_path)
                # Older samtools versions might also use REF_PATH
                env['REF_PATH'] = str(ref_path)
            else:
                # If no ref found for CRAM, indexing will likely fail/hang
                print(f"  [Index] Warning: CRAM file but no reference detected/found. Indexing may fail.")
        
        cmd.append(str(bam_path))
        
        # We allow time for indexing (it's heavy), but we don't wait for Yleaf to hang
        subprocess.check_call(cmd, timeout=yleaf_constants.INDEX_TIMEOUT, env=env)
        print(f"  [Index] Re-indexing complete.")
    except subprocess.TimeoutExpired:
        print(f"  [Index] Failed: Indexing timed out (>10m).")
        return False
    except subprocess.CalledProcessError:
        print(f"  [Index] Failed: 'samtools index' returned error.")
        return False
        
    # 3. Verify new index
    # For CRAM, we skip validation to avoid infinite loops if idxstats fails
    if bam_path.suffix == '.cram':
        return True
        
    if check_index_fast(bam_path, ref_path):
        return True
    else:
        print(f"  [Index] Failed: Index created but 'idxstats' still fails.")
        return False

def run_yleaf(bam_path, output_base_dir, read_thresh=None, quality_thresh=None):
    bam_path = Path(bam_path)
    # Ensure base output directory exists
    Path(output_base_dir).mkdir(parents=True, exist_ok=True)
    
    output_dir = Path(output_base_dir) / f"output_{{bam_path.stem}}"
    
    # Write log outside the target dir to prevent deletion by Yleaf
    log_file_path = output_dir.with_suffix('.console.log')
    
    print(f"--- Processing {bam_path.name} ---")
    
    # Determine path to Yleaf.py relative to this script
    script_dir = Path(__file__).resolve().parent
    yleaf_script = script_dir / "yleaf" / "Yleaf.py"
    
    cmd = [sys.executable, str(yleaf_script), "-bam", str(bam_path), "-o", str(output_dir)]
    
    # Add optional arguments
    if read_thresh is not None:
        cmd.extend(["-r", str(read_thresh)])
    if quality_thresh is not None:
        cmd.extend(["-q", str(quality_thresh)])
    
    # Set PYTHONPATH so Yleaf can find its package
    # batch_process.py is in Yleaf/, so we add Yleaf/ to PYTHONPATH
    env = os.environ.copy()
    env["PYTHONPATH"] = str(script_dir) + os.pathsep + env.get("PYTHONPATH", "")

    def run_cmd_realtime():
        print(f"  [1/2] Starting Yleaf subprocess...", end='', flush=True)
        start_time = time.time()
        # Safety net timeout still exists, but we rely on ensure_index to catch bad files
        MAX_DURATION = yleaf_constants.BATCH_PROCESS_TIMEOUT 
        
        with open(log_file_path, "w") as log:
            # Pass the modified environment
            process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
            
            while process.poll() is None:
                time.sleep(2) 
                if time.time() - start_time > MAX_DURATION:
                    print(f"  [Monitor] Safety Timeout ({{MAX_DURATION}}s). Killing...", end='')
                    process.kill()
                    print(" Killed.")
                    return -1
            
            print(" Done.")
            return process.returncode

    # Run Yleaf
    exit_code = run_cmd_realtime()
    
    # We do NOT move the log file. It stays as .console.log next to the folder.
    
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

def process_single_file(file_path, output_dir, read_thresh=None, quality_thresh=None):
    """
    Worker function to process a single BAM/CRAM file.
    """
    if not file_path.exists():
        print(f"File not found, skipping: {file_path}")
        return

    # 0. Integrity Check
    if not check_bam_integrity(file_path):
        print(f"Skipping corrupted file: {file_path.name}\n")
        summary_logger.log_failure(file_path, "Skipped: File Corrupted")
        return

    # 1. Check/Ensure Index
    if not ensure_index(file_path):
        print(f"Skipping {file_path.name} due to invalid/missing index.")
        summary_logger.log_failure(file_path, "Skipped: Index Invalid")
        return
        
    # 2. Run Yleaf
    run_yleaf(file_path, output_dir, read_thresh, quality_thresh)

def main():
    parser = argparse.ArgumentParser(description="Batch process BAM files with Yleaf.")
    parser.add_argument("input_source", help="Text file with list of paths OR directory to scan (if -d is used)")
    parser.add_argument("-d", "--directory-mode", action="store_true", help="Treat input_source as a directory and scan recursively for .bam/.cram")
    parser.add_argument("-o", "--output_dir", default=".", help="Base directory for output folders")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of files to process in parallel (default: 1)")
    parser.add_argument("-r", "--reads_threshold", type=int, help="Minimum read depth (default: Yleaf default)")
    parser.add_argument("-q", "--quality_threshold", type=int, help="Minimum read quality (default: Yleaf default)")
    
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
        
    print(f"Found {len(files_to_process)} files to process. Using {args.threads} threads.")
    
    # Prepare partial function with arguments
    process_func = partial(
        process_single_file, 
        output_dir=args.output_dir,
        read_thresh=args.reads_threshold,
        quality_thresh=args.quality_threshold
    )
    
    if args.threads > 1:
        # Create output directory ahead of time to avoid race conditions
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        
        with multiprocessing.Pool(processes=args.threads) as pool:
            pool.map(process_func, files_to_process)
    else:
        for file_path in files_to_process:
            process_func(file_path)

if __name__ == "__main__":
    main()