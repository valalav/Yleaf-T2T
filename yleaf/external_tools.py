#!/usr/bin/env python
"""
Cross-platform abstraction layer for external bioinformatics tools.
Provides Windows-compatible wrappers for samtools, bcftools, and other CLI tools.

This module replaces direct shell=True subprocess calls with cross-platform alternatives.
"""

import os
import platform
import subprocess
import logging
import shutil
from pathlib import Path
from typing import List, Optional, Union, Tuple

LOG: logging.Logger = logging.getLogger("yleaf_logger")

# Platform detection
IS_WINDOWS = platform.system() == "Windows"

# Tool paths (will be set during initialization)
_TOOL_PATHS = {
    "samtools": None,
    "bcftools": None,
    "graphviz": None,
}


def _get_bundled_bin_dir() -> Path:
    """Get the directory containing bundled Windows binaries."""
    src_folder = Path(__file__).absolute().parent
    if IS_WINDOWS:
        return src_folder / "bin" / "windows"
    return src_folder / "bin" / "linux"


def find_tool(tool_name: str) -> Optional[Path]:
    """
    Find the path to an external tool.

    Search order:
    1. Already configured path in _TOOL_PATHS
    2. Bundled binaries in yleaf/bin/
    3. System PATH

    Args:
        tool_name: Name of the tool (e.g., 'samtools', 'bcftools')

    Returns:
        Path to the tool executable, or None if not found
    """
    # Check cached path
    if _TOOL_PATHS.get(tool_name):
        return _TOOL_PATHS[tool_name]

    # Check bundled binaries
    bundled_dir = _get_bundled_bin_dir()
    if IS_WINDOWS:
        bundled_path = bundled_dir / f"{tool_name}.exe"
    else:
        bundled_path = bundled_dir / tool_name

    if bundled_path.exists():
        _TOOL_PATHS[tool_name] = bundled_path
        return bundled_path

    # Check system PATH
    system_path = shutil.which(tool_name)
    if system_path:
        _TOOL_PATHS[tool_name] = Path(system_path)
        return Path(system_path)

    return None


def set_tool_path(tool_name: str, path: Union[str, Path]) -> None:
    """Manually set the path for a tool."""
    _TOOL_PATHS[tool_name] = Path(path)


def run_command(
    args: List[str],
    capture_output: bool = False,
    stdout_file: Optional[Path] = None,
    check: bool = True,
    timeout: Optional[int] = None
) -> subprocess.CompletedProcess:
    """
    Run a command without shell=True for cross-platform compatibility.

    Args:
        args: Command and arguments as a list
        capture_output: Whether to capture stdout/stderr
        stdout_file: Optional file to redirect stdout to
        check: Whether to raise exception on non-zero return code
        timeout: Optional timeout in seconds

    Returns:
        subprocess.CompletedProcess object

    Raises:
        SystemExit: If check=True and command fails
    """
    LOG.debug(f"Running command: {' '.join(str(a) for a in args)}")

    try:
        stdout_target = subprocess.PIPE if capture_output else None

        if stdout_file:
            with open(stdout_file, 'w') as f:
                result = subprocess.run(
                    args,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    timeout=timeout,
                    text=True
                )
        else:
            # Use manual PIPE assignment to avoid conflict with capture_output arg
            stdout_target = subprocess.PIPE if capture_output else None
            
            result = subprocess.run(
                args,
                stdout=stdout_target,
                stderr=subprocess.PIPE,
                timeout=timeout,
                text=True
            )

        if check and result.returncode != 0:
            LOG.error(f"Command failed: {' '.join(str(a) for a in args)}")
            LOG.error(f"Error: {result.stderr}")
            raise SystemExit("Failed command execution")

        LOG.debug(f"Command completed successfully")
        return result

    except subprocess.TimeoutExpired:
        LOG.error(f"Command timed out: {' '.join(str(a) for a in args)}")
        raise SystemExit("Command timeout")
    except FileNotFoundError as e:
        LOG.error(f"Command not found: {args[0]}")
        raise SystemExit(f"Tool not found: {args[0]}")


# =============================================================================
# SAMTOOLS WRAPPERS
# =============================================================================

def samtools_view_header(bam_file: Path) -> List[str]:
    """
    Run 'samtools view -H' and return header lines.
    Replaces: samtools view -H file | grep @SQ

    Args:
        bam_file: Path to BAM/CRAM file

    Returns:
        List of @SQ header lines
    """
    tool = find_tool("samtools")
    if not tool:
        raise SystemExit("samtools not found")

    result = run_command([str(tool), "view", "-H", str(bam_file)], capture_output=True)

    # Python replacement for grep @SQ
    header_lines = [line for line in result.stdout.splitlines() if "@SQ" in line]
    return header_lines


def samtools_index(bam_file: Path, threads: int = 1) -> None:
    """
    Run 'samtools index'.

    Args:
        bam_file: Path to BAM/CRAM file
        threads: Number of threads to use
    """
    tool = find_tool("samtools")
    if not tool:
        raise SystemExit("samtools not found")

    run_command([str(tool), "index", f"-@{threads}", str(bam_file)])


def samtools_idxstats(bam_file: Path, output_file: Path) -> None:
    """
    Run 'samtools idxstats' and save to file.

    Args:
        bam_file: Path to BAM/CRAM file
        output_file: Path to output file
    """
    tool = find_tool("samtools")
    if not tool:
        raise SystemExit("samtools not found")

    run_command([str(tool), "idxstats", str(bam_file)], stdout_file=output_file)


def samtools_mpileup(
    bam_file: Path,
    output_file: Path,
    quality_thresh: int,
    bed_file: Optional[Path] = None,
    reference: Optional[Path] = None,
    region: Optional[str] = None
) -> None:
    """
    Run 'samtools mpileup' with options.
    Replaces: samtools mpileup -r region -l bed -AQ20q1 bam > output

    Args:
        bam_file: Path to BAM/CRAM file
        output_file: Path to output pileup file
        quality_thresh: Minimum base quality
        bed_file: Optional BED file for positions
        reference: Optional reference FASTA
        region: Optional region (e.g., 'chrY')
    """
    tool = find_tool("samtools")
    if not tool:
        raise SystemExit("samtools not found")

    args = [str(tool), "mpileup"]

    if region:
        args.extend(["-r", region])

    if bed_file:
        args.extend(["-l", str(bed_file)])

    if reference:
        args.extend(["-f", str(reference)])

    args.extend([f"-AQ{quality_thresh}q1", str(bam_file)])

    run_command(args, stdout_file=output_file)


def samtools_view_sort(
    input_file: Path,
    output_bam: Path,
    threads: int = 1,
    input_is_sam: bool = False
) -> None:
    """
    Run 'samtools view | samtools sort' as two separate commands.
    Replaces: samtools view -@ N -bS input | samtools sort -@ N -m 2G -o output

    Args:
        input_file: Path to input SAM/BAM file
        output_bam: Path to output sorted BAM file
        threads: Number of threads
        input_is_sam: Whether input is SAM format
    """
    tool = find_tool("samtools")
    if not tool:
        raise SystemExit("samtools not found")

    # Create temp BAM file
    temp_bam = output_bam.parent / f"temp_{output_bam.name}"

    # Step 1: Convert SAM to BAM if needed
    if input_is_sam:
        view_args = [str(tool), "view", f"-@{threads}", "-bS", str(input_file), "-o", str(temp_bam)]
        run_command(view_args)
        sort_input = temp_bam
    else:
        sort_input = input_file

    # Step 2: Sort
    sort_args = [str(tool), "sort", f"-@{threads}", "-m", "2G", "-o", str(output_bam), str(sort_input)]
    run_command(sort_args)

    # Cleanup temp file
    if input_is_sam and temp_bam.exists():
        os.remove(temp_bam)


# =============================================================================
# BCFTOOLS WRAPPERS
# =============================================================================

def bcftools_query(
    vcf_file: Path,
    format_string: str,
    output_file: Optional[Path] = None
) -> Optional[str]:
    """
    Run 'bcftools query' with format string.

    Args:
        vcf_file: Path to VCF file
        format_string: Format string for bcftools query
        output_file: Optional output file (if None, returns stdout)

    Returns:
        Query output as string if output_file is None
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    args = [str(tool), "query", "-f", format_string, str(vcf_file)]

    if output_file:
        run_command(args, stdout_file=output_file)
        return None
    else:
        result = run_command(args, capture_output=True)
        return result.stdout


def bcftools_query_samples(vcf_file: Path) -> List[str]:
    """
    Get list of samples from VCF file.
    Replaces: bcftools query -l vcf | wc -l

    Args:
        vcf_file: Path to VCF file

    Returns:
        List of sample names
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    result = run_command([str(tool), "query", "-l", str(vcf_file)], capture_output=True)
    samples = [s for s in result.stdout.strip().splitlines() if s]
    return samples


def bcftools_query_chromosomes(vcf_file: Path) -> List[str]:
    """
    Get unique chromosomes from VCF file.
    Replaces: bcftools query -f '%CHROM\n' vcf | uniq

    Args:
        vcf_file: Path to VCF file

    Returns:
        List of unique chromosome names (preserving order)
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    result = run_command([str(tool), "query", "-f", "%CHROM\n", str(vcf_file)], capture_output=True)

    # Python replacement for uniq (preserves order)
    seen = set()
    unique_chroms = []
    for chrom in result.stdout.strip().splitlines():
        if chrom and chrom not in seen:
            seen.add(chrom)
            unique_chroms.append(chrom)
    return unique_chroms


def bcftools_sort(
    input_vcf: Path,
    output_vcf: Path,
    compressed: bool = True
) -> None:
    """
    Sort VCF file.
    Replaces: bcftools sort -O z -o output input

    Args:
        input_vcf: Path to input VCF
        output_vcf: Path to output VCF
        compressed: Whether to output compressed VCF
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    output_format = "z" if compressed else "v"
    args = [str(tool), "sort", "-O", output_format, "-o", str(output_vcf), str(input_vcf)]
    run_command(args)


def bcftools_index(vcf_file: Path, force: bool = True) -> None:
    """
    Index VCF file.
    Replaces: bcftools index -f vcf

    Args:
        vcf_file: Path to VCF file
        force: Whether to overwrite existing index
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    args = [str(tool), "index"]
    if force:
        args.append("-f")
    args.append(str(vcf_file))
    run_command(args)


def bcftools_view_regions(
    input_vcf: Path,
    output_vcf: Path,
    bed_file: Path,
    compressed: bool = True
) -> None:
    """
    Filter VCF by regions from BED file.
    Replaces: bcftools view -O z -R bed vcf > output

    Args:
        input_vcf: Path to input VCF
        output_vcf: Path to output VCF
        bed_file: Path to BED file with regions
        compressed: Whether to output compressed VCF
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    output_format = "z" if compressed else "v"
    args = [str(tool), "view", "-O", output_format, "-R", str(bed_file), str(input_vcf), "-o", str(output_vcf)]
    run_command(args)


def bcftools_split(
    vcf_file: Path,
    output_dir: Path,
    compressed: bool = True
) -> None:
    """
    Split multi-sample VCF into individual files.
    Replaces: bcftools +split vcf -Oz -o dir

    Args:
        vcf_file: Path to input VCF
        output_dir: Path to output directory
        compressed: Whether to output compressed VCFs
    """
    tool = find_tool("bcftools")
    if not tool:
        raise SystemExit("bcftools not found")

    output_format = "Oz" if compressed else "Ov"
    args = [str(tool), "+split", str(vcf_file), f"-{output_format}", "-o", str(output_dir)]
    run_command(args)


# =============================================================================
# UTILITY FUNCTIONS (Python replacements for shell utilities)
# =============================================================================

def count_lines(text: str) -> int:
    """
    Count lines in text.
    Python replacement for: wc -l
    """
    return len(text.strip().splitlines()) if text.strip() else 0


def unique_lines(text: str) -> List[str]:
    """
    Get unique lines preserving order.
    Python replacement for: uniq
    """
    seen = set()
    result = []
    for line in text.splitlines():
        if line not in seen:
            seen.add(line)
            result.append(line)
    return result


def filter_lines(text: str, pattern: str) -> List[str]:
    """
    Filter lines containing pattern.
    Python replacement for: grep pattern
    """
    return [line for line in text.splitlines() if pattern in line]


def safe_remove(file_path: Path) -> None:
    """
    Safely remove a file if it exists.
    Cross-platform replacement for: rm file
    """
    try:
        if file_path.exists():
            os.remove(file_path)
            LOG.debug(f"Removed file: {file_path}")
    except OSError as e:
        LOG.warning(f"Failed to remove file {file_path}: {e}")


# =============================================================================
# INITIALIZATION
# =============================================================================

def check_tools() -> dict:
    """
    Check availability of required tools.

    Returns:
        Dict with tool names as keys and (found: bool, path: Optional[Path]) as values
    """
    tools = ["samtools", "bcftools"]
    status = {}

    for tool in tools:
        path = find_tool(tool)
        status[tool] = {
            "found": path is not None,
            "path": path
        }
        if path:
            LOG.debug(f"Found {tool} at: {path}")
        else:
            LOG.warning(f"{tool} not found")

    return status


def get_tool_versions() -> dict:
    """Get versions of available tools."""
    versions = {}

    for tool_name in ["samtools", "bcftools"]:
        tool = find_tool(tool_name)
        if tool:
            try:
                result = subprocess.run(
                    [str(tool), "--version"],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                first_line = result.stdout.splitlines()[0] if result.stdout else "unknown"
                versions[tool_name] = first_line
            except Exception:
                versions[tool_name] = "unknown"
        else:
            versions[tool_name] = "not found"

    return versions


# =============================================================================
# MINIMAP2 WRAPPERS
# =============================================================================

def minimap2_align(
    reference: Path,
    fastq_files: List[Path],
    output_sam: Path,
    threads: int = 1,
    preset: str = "sr"
) -> None:
    """
    Run minimap2 alignment.
    Replaces: minimap2 -ax sr -k14 -w7 -t threads reference fastq1 [fastq2] > output.sam

    Args:
        reference: Path to reference FASTA
        fastq_files: List of FASTQ files (1 for single-end, 2 for paired-end)
        output_sam: Path to output SAM file
        threads: Number of threads
        preset: Minimap2 preset (default: 'sr' for short reads)
    """
    tool = find_tool("minimap2")
    if not tool:
        raise SystemExit("minimap2 not found")

    args = [
        str(tool),
        "-ax", preset,
        "-k14", "-w7",
        "-t", str(threads),
        str(reference)
    ]

    # Add FASTQ files
    for fq in fastq_files:
        args.append(str(fq))

    run_command(args, stdout_file=output_sam)
