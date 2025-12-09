import os
import sys
import subprocess
from pathlib import Path
import time
import multiprocessing

def process_file(bam_path_str):
    """
    Checks and fixes index for a single BAM/CRAM file.
    Returns: (path, status_message)
    """
    bam_path = Path(bam_path_str)
    
    # Determine expected index path
    if bam_path.suffix == '.cram':
        idx_path = bam_path.with_suffix('.cram.crai')
        alt_idx = bam_path.with_suffix('.crai')
    else:
        idx_path = bam_path.with_suffix('.bam.bai')
        alt_idx = bam_path.with_suffix('.bai')
        
    # Check if any valid index exists
    existing_idx = None
    if idx_path.exists() and idx_path.stat().st_size > 0:
        existing_idx = idx_path
    elif alt_idx.exists() and alt_idx.stat().st_size > 0:
        existing_idx = alt_idx
        
    needs_indexing = False
    reason = ""
    
    if not existing_idx:
        needs_indexing = True
        reason = "Missing"
    else:
        # Check timestamp
        try:
            bam_mtime = bam_path.stat().st_mtime
            idx_mtime = existing_idx.stat().st_mtime
            if idx_mtime < bam_mtime:
                needs_indexing = True
                reason = "Outdated"
                # Remove old index to be safe
                try: existing_idx.unlink()
                except (OSError, PermissionError): pass
        except FileNotFoundError:
            needs_indexing = True
            reason = "Stat failed"

    if needs_indexing:
        # Create index
        try:
            cmd = ["samtools", "index", str(bam_path)]
            # If multithreaded samtools is available, use it? No, we parallelize files.
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return str(bam_path), f"FIXED ({reason})"
        except subprocess.CalledProcessError:
            return None, f"FAILED ({reason})"
    
    return str(bam_path), "OK"

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pre_check_indexes.py <directory> [threads]")
        sys.exit(1)
        
    root_dir = Path(sys.argv[1])
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    
    print(f"Scanning {root_dir} for BAM/CRAM files...")
    files = list(root_dir.rglob("*.bam")) + list(root_dir.rglob("*.cram"))
    print(f"Found {len(files)} files. Processing with {threads} threads...")
    
    valid_files = []
    
    with multiprocessing.Pool(processes=threads) as pool:
        # Map process_file over all files
        results = pool.map(process_file, [str(f) for f in files])
        
        for path, status in results:
            if path:
                valid_files.append(path)
                if status != "OK":
                    print(f"[+] {Path(path).name}: {status}")
            else:
                # Error case (path is None), status contains error msg
                # We don't know name easily here unless we pass it back, but message has it? 
                # Actually process_file returns None as path if failed.
                print(f"[-] Indexing failed: {status}")

    # Save clean list
    out_list = "clean_file_list.txt"
    with open(out_list, 'w') as f:
        f.write('\n'.join(valid_files))
        
    print(f"\nDone. {len(valid_files)} valid files saved to {out_list}")

if __name__ == "__main__":
    main()