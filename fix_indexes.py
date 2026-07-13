#!/usr/bin/env python3
"""
fix_indexes.py - Batch index checker and fixer for BAM/CRAM files.

Recursively scans directories for BAM/CRAM files, validates their indexes,
and creates/recreates indexes as needed. Optimized for parallel execution.

Usage:
    python3 fix_indexes.py /path/to/data -t 16
    python3 fix_indexes.py /path/to/data -t 32 --storage ssd
    python3 fix_indexes.py /path/to/data --dry-run
"""

import os
import sys
import argparse
import subprocess
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional, List, Tuple
import time
import shutil

# Reference genome paths (adjust to your system)
REFERENCE_GENOMES = {
    "hg19": "/home/valalav/wgs/WGSExtractv4/reference/genomes/hs37d5.fa",
    "hg38": "/home/valalav/wgs/WGSExtractv4/reference/genomes/hg38.fa",
    "t2t": "/home/valalav/wgs/WGSExtractv4/reference/genomes/chm13v2.0.fa",
}

# Directory containing reference genomes (for REF_PATH)
REFERENCE_DIR = "/home/valalav/wgs/WGSExtractv4/reference/genomes"

# Chromosome 1 lengths for reference detection
CHR1_LENGTHS = {
    249250621: "hg19",   # GRCh37
    248956422: "hg38",   # GRCh38
    248387328: "t2t",    # T2T-CHM13
}


@dataclass
class IndexResult:
    """Result of index check/fix operation."""
    file_path: Path
    status: str  # 'ok', 'created', 'recreated', 'failed', 'skipped'
    message: str
    duration: float = 0.0


def get_samtools_path() -> str:
    """Find samtools executable."""
    samtools = shutil.which("samtools")
    if not samtools:
        raise RuntimeError("samtools not found in PATH")
    return samtools


def get_cram_env() -> dict:
    """Get environment with REF_PATH set for CRAM operations."""
    env = os.environ.copy()
    env['REF_PATH'] = REFERENCE_DIR
    return env


def detect_reference(file_path: Path) -> Optional[str]:
    """
    Detect reference genome from BAM/CRAM header.
    Returns path to reference FASTA or None.
    """
    try:
        # Use REF_PATH to prevent samtools from trying to download references
        env = get_cram_env() if file_path.suffix.lower() == '.cram' else None

        result = subprocess.run(
            ["samtools", "view", "-H", str(file_path)],
            capture_output=True, text=True, timeout=60,
            env=env
        )
        if result.returncode != 0:
            return None

        for line in result.stdout.split('\n'):
            if line.startswith("@SQ") and ("SN:chr1\t" in line or "SN:1\t" in line):
                for part in line.split('\t'):
                    if part.startswith("LN:"):
                        length = int(part.split(':')[1])
                        ref_name = CHR1_LENGTHS.get(length)
                        if ref_name and ref_name in REFERENCE_GENOMES:
                            ref_path = REFERENCE_GENOMES[ref_name]
                            if Path(ref_path).exists():
                                return ref_path
                        break
    except Exception:
        pass
    return None


def get_index_path(file_path: Path) -> Path:
    """Get expected index path for BAM/CRAM file."""
    if file_path.suffix.lower() == '.cram':
        return file_path.with_suffix('.cram.crai')
    else:
        # BAM can have .bam.bai or .bai
        bam_bai = file_path.with_suffix('.bam.bai')
        if bam_bai.exists():
            return bam_bai
        return Path(str(file_path) + '.bai')


def check_index_valid(file_path: Path, reference: Optional[str] = None) -> bool:
    """
    Check if index is valid using samtools idxstats.
    Returns True if index works, False otherwise.
    """
    try:
        # For CRAM, always use REF_PATH to prevent network access
        if file_path.suffix.lower() == '.cram':
            env = get_cram_env()
        else:
            env = os.environ.copy()

        result = subprocess.run(
            ["samtools", "idxstats", str(file_path)],
            capture_output=True, timeout=120, env=env
        )
        return result.returncode == 0 and len(result.stdout) > 0
    except Exception:
        return False


def create_index(file_path: Path, reference: Optional[str] = None) -> Tuple[bool, str]:
    """
    Create index for BAM/CRAM file.
    Returns (success, message).
    """
    try:
        # For CRAM, always use REF_PATH directory to prevent network access
        if file_path.suffix.lower() == '.cram':
            env = get_cram_env()
        else:
            env = os.environ.copy()

        result = subprocess.run(
            ["samtools", "index", str(file_path)],
            capture_output=True, text=True, timeout=1800,  # 30 min timeout
            env=env
        )

        if result.returncode == 0:
            return True, "Index created successfully"
        else:
            return False, f"samtools error: {result.stderr.strip()}"
    except subprocess.TimeoutExpired:
        return False, "Timeout (>30 min)"
    except Exception as e:
        return False, f"Error: {str(e)}"


def remove_index(file_path: Path) -> None:
    """Remove all possible index files for a BAM/CRAM."""
    extensions = ['.bai', '.bam.bai', '.crai', '.cram.crai']
    base = str(file_path)

    for ext in extensions:
        # Try both patterns: file.bam.bai and file.bai
        for idx_path in [Path(base + ext), file_path.with_suffix(ext)]:
            if idx_path.exists():
                try:
                    idx_path.unlink()
                except OSError:
                    pass


def process_file(file_path: Path, force_recreate: bool = False, dry_run: bool = False) -> IndexResult:
    """
    Process a single BAM/CRAM file - check and fix its index.
    """
    start_time = time.time()
    file_path = Path(file_path)

    if not file_path.exists():
        return IndexResult(file_path, 'skipped', 'File not found')

    is_cram = file_path.suffix.lower() == '.cram'
    reference = detect_reference(file_path) if is_cram else None

    if is_cram and not reference:
        return IndexResult(
            file_path, 'skipped',
            'CRAM file but could not detect reference genome',
            time.time() - start_time
        )

    index_path = get_index_path(file_path)
    index_exists = index_path.exists()

    # Check if index is valid
    if index_exists and not force_recreate:
        if check_index_valid(file_path, reference):
            return IndexResult(
                file_path, 'ok',
                'Index valid',
                time.time() - start_time
            )
        else:
            status = 'recreated'
            action = 'Recreating invalid index'
    else:
        status = 'created' if not index_exists else 'recreated'
        action = 'Creating index' if not index_exists else 'Recreating index (forced)'

    if dry_run:
        return IndexResult(
            file_path, 'dry_run',
            f"Would {action.lower()}",
            time.time() - start_time
        )

    # Remove old index if exists
    if index_exists:
        remove_index(file_path)

    # Create new index
    success, message = create_index(file_path, reference)

    if success:
        # Verify the new index works
        if check_index_valid(file_path, reference):
            return IndexResult(
                file_path, status,
                f"{action} - OK",
                time.time() - start_time
            )
        else:
            return IndexResult(
                file_path, 'failed',
                f"{action} - created but validation failed",
                time.time() - start_time
            )
    else:
        return IndexResult(
            file_path, 'failed',
            f"{action} - {message}",
            time.time() - start_time
        )


def find_files(root_path: Path, extensions: List[str] = None) -> List[Path]:
    """Recursively find all BAM/CRAM files."""
    if extensions is None:
        extensions = ['.bam', '.cram']

    files = []
    for ext in extensions:
        files.extend(root_path.rglob(f'*{ext}'))

    return sorted(files)


def worker_init():
    """Initialize worker process - ignore SIGINT for clean shutdown."""
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def main():
    parser = argparse.ArgumentParser(
        description='Batch check and fix BAM/CRAM indexes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan directory with 16 parallel workers
  python3 fix_indexes.py /mnt/truenas-data/Data/DNA/wgs -t 16

  # Optimized for SSD (more workers)
  python3 fix_indexes.py /data/wgs -t 32 --storage ssd

  # Dry run - just show what would be done
  python3 fix_indexes.py /path/to/data --dry-run

  # Force recreate all indexes
  python3 fix_indexes.py /path/to/data -t 16 --force

  # Process only CRAM files
  python3 fix_indexes.py /path/to/data -t 16 --cram-only
        """
    )

    parser.add_argument('path', type=Path, help='Directory to scan recursively')
    parser.add_argument('-t', '--threads', type=int, default=None,
                        help='Number of parallel workers (default: auto based on storage type)')
    parser.add_argument('--storage', choices=['nas', 'ssd', 'hdd'], default='nas',
                        help='Storage type for optimal thread count (default: nas)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without making changes')
    parser.add_argument('--force', action='store_true',
                        help='Force recreate all indexes even if valid')
    parser.add_argument('--bam-only', action='store_true',
                        help='Process only BAM files')
    parser.add_argument('--cram-only', action='store_true',
                        help='Process only CRAM files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show progress for each file')

    args = parser.parse_args()

    if not args.path.exists():
        print(f"Error: Path not found: {args.path}")
        sys.exit(1)

    # Determine optimal thread count based on storage type
    if args.threads is None:
        cpu_count = multiprocessing.cpu_count()
        if args.storage == 'ssd':
            # SSD can handle more parallel I/O
            args.threads = min(cpu_count, 32)
        elif args.storage == 'hdd':
            # HDD is I/O bound, limit parallelism
            args.threads = min(cpu_count // 2, 8)
        else:  # NAS
            # NAS depends on network, moderate parallelism
            args.threads = min(cpu_count, 16)

    # Determine extensions to process
    extensions = []
    if args.bam_only:
        extensions = ['.bam']
    elif args.cram_only:
        extensions = ['.cram']
    else:
        extensions = ['.bam', '.cram']

    print(f"Scanning {args.path} for {', '.join(extensions)} files...")
    files = find_files(args.path, extensions)

    if not files:
        print("No files found.")
        sys.exit(0)

    print(f"Found {len(files)} files to process")
    print(f"Using {args.threads} parallel workers (storage: {args.storage})")
    if args.dry_run:
        print("DRY RUN - no changes will be made\n")
    print()

    # Process files in parallel
    results = {
        'ok': 0,
        'created': 0,
        'recreated': 0,
        'failed': 0,
        'skipped': 0,
        'dry_run': 0,
    }
    failed_files = []

    start_time = time.time()

    try:
        with ProcessPoolExecutor(max_workers=args.threads, initializer=worker_init) as executor:
            futures = {
                executor.submit(process_file, f, args.force, args.dry_run): f
                for f in files
            }

            completed = 0
            for future in as_completed(futures):
                completed += 1
                result = future.result()
                results[result.status] += 1

                if result.status == 'failed':
                    failed_files.append((result.file_path, result.message))

                if args.verbose or result.status in ('created', 'recreated', 'failed'):
                    status_icon = {
                        'ok': '✓',
                        'created': '+',
                        'recreated': '↻',
                        'failed': '✗',
                        'skipped': '-',
                        'dry_run': '?',
                    }.get(result.status, '?')

                    print(f"[{completed}/{len(files)}] {status_icon} {result.file_path.name}: {result.message} ({result.duration:.1f}s)")
                elif completed % 10 == 0:
                    print(f"Progress: {completed}/{len(files)} files processed...", end='\r')

    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Partial results:")

    # Print summary
    duration = time.time() - start_time
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total files:      {len(files)}")
    print(f"Already valid:    {results['ok']}")
    print(f"Created:          {results['created']}")
    print(f"Recreated:        {results['recreated']}")
    print(f"Failed:           {results['failed']}")
    print(f"Skipped:          {results['skipped']}")
    if args.dry_run:
        print(f"Would process:    {results['dry_run']}")
    print(f"Duration:         {duration:.1f}s ({duration/60:.1f} min)")
    print(f"Speed:            {len(files)/duration:.1f} files/sec")

    if failed_files:
        print("\nFailed files:")
        for path, msg in failed_files:
            print(f"  {path}: {msg}")

    sys.exit(1 if results['failed'] > 0 else 0)


if __name__ == '__main__':
    main()
