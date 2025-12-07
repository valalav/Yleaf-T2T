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
    output_dir = Path(output_base_dir) / f"output_{bam_path.stem}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    log_file_path = output_dir / "full_terminal_output.log"
    
    print(f"--- Processing {bam_path.name} ---")
    print(f"Logs will be saved to: {log_file_path}")
    
    cmd = [
        "Yleaf",
        "-bam", str(bam_path),
        "-o", str(output_dir)
    ]
    
    def execute_process():
        with open(log_file_path, "w") as log_file:
            process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )
            output_lines = []
            for line in process.stdout:
                print(line, end='')
                log_file.write(line)
                output_lines.append(line)
            
            return_code = process.wait()
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, cmd, output="".join(output_lines))

    try:
        execute_process()
        print(f"Successfully processed {bam_path.name}\n")
        
    except subprocess.CalledProcessError as e:
        # Check if error is index related
        err_text = e.output.lower() if e.output else ""
        if "index" in err_text or "bgzf" in err_text or "older than" in err_text:
            print(f"\n[Auto-Heal] Detected potential index issue for {bam_path.name}. Re-indexing...")
            
            # Try to remove old index
            bai = bam_path.with_suffix('.bam.bai')
            bai2 = bam_path.with_suffix('.bai')
            crai = bam_path.with_suffix('.cram.crai')
            crai2 = bam_path.with_suffix('.crai')
            
            for idx in [bai, bai2, crai, crai2]:
                if idx.exists():
                    try:
                        idx.unlink()
                        print(f"Removed old index: {idx}")
                    except OSError:
                        print(f"Could not remove index: {idx} (permission?)")

            # Re-create index using our helper
            if check_create_index(bam_path):
                print(f"[Auto-Heal] Re-indexing successful. Retrying Yleaf...")
                try:
                    execute_process()
                    print(f"Successfully processed {bam_path.name} (after fix)\n")
                    return # Success
                except subprocess.CalledProcessError:
                    print(f"[Auto-Heal] Retry failed. File might be truly corrupted.")
            else:
                print(f"[Auto-Heal] Re-indexing failed.")

        print(f"Error processing {bam_path.name}: Process exited with code {e.returncode}")
        print(f"Check log for details: {log_file_path}\n")

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
            continue

        # 1. Check/Create Index
        if not check_create_index(file_path):
            print(f"Skipping {file_path.name} due to indexing failure.")
            continue
            
        # 2. Run Yleaf
        run_yleaf(file_path, args.output_dir)

if __name__ == "__main__":
    main()
