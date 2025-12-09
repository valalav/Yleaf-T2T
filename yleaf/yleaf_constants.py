"""
Yleaf constants and configuration module.

This module provides cross-platform path handling and configuration
for the Yleaf Y-chromosome haplogroup prediction tool.
"""

import platform
import logging
from pathlib import Path
from typing import Optional

# =============================================================================
# PLATFORM DETECTION
# =============================================================================

IS_WINDOWS = platform.system() == "Windows"
IS_LINUX = platform.system() == "Linux"
IS_MACOS = platform.system() == "Darwin"

# =============================================================================
# PATH CONSTANTS
# =============================================================================

SRC_FOLDER: Path = Path(__file__).absolute().parent
DATA_FOLDER: Path = SRC_FOLDER / "data"
CONFIG_PATH: Path = SRC_FOLDER / "config.txt"
CHR_NAMING_CONVENTION_FILE: Path = SRC_FOLDER / "chr_naming_convention.txt"
HG_PREDICTION_FOLDER: Path = DATA_FOLDER / "hg_prediction_tables"

# Bundled binaries directory
BIN_FOLDER: Path = SRC_FOLDER / "bin"
if IS_WINDOWS:
    BIN_PLATFORM_FOLDER: Path = BIN_FOLDER / "windows"
else:
    BIN_PLATFORM_FOLDER: Path = BIN_FOLDER / "linux"

# =============================================================================
# REFERENCE GENOME CONSTANTS
# =============================================================================

HG19: str = "hg19"
HG38: str = "hg38"
T2T: str = "t2t"

__HG19_FOLDER: Path = DATA_FOLDER / HG19
__HG38_FOLDER: Path = DATA_FOLDER / HG38
__T2T_FOLDER: Path = DATA_FOLDER / T2T

FULL_REF_FILE: str = "full_reference.fa"
Y_REF_FILE: str = "chrY.fa"
SNP_DATA_FILE: str = "snp_data.csv"
NEW_POSITION_FILE: str = "new_positions.txt"
OLD_POSITION_FILE: str = "old_positions.txt"
NEW_POSITION_BED_FILE: str = "new_positions.bed"
OLD_POSITION_BED_FILE: str = "old_positions.bed"
NEW_POSITION_ANCIENT_FILE: str = "new_positions_ancient.txt"
OLD_POSITION_ANCIENT_FILE: str = "old_positions_ancient.txt"
NEW_POSITION_ANCIENT_BED_FILE: str = "new_positions_ancient.bed"
OLD_POSITION_ANCIENT_BED_FILE: str = "old_positions_ancient.bed"

TREE_FILE: str = "tree.json"

# Default reference paths (will be overridden by config if available)
HG19_FULL_GENOME: Path = __HG19_FOLDER / FULL_REF_FILE
HG19_Y_CHROMOSOME: Path = __HG19_FOLDER / Y_REF_FILE
HG38_FULL_GENOME: Path = __HG38_FOLDER / FULL_REF_FILE
HG38_Y_CHROMOSOME: Path = __HG38_FOLDER / Y_REF_FILE
T2T_FULL_GENOME: Path = __T2T_FOLDER / FULL_REF_FILE
T2T_Y_CHROMOSOME: Path = __T2T_FOLDER / Y_REF_FILE

FASTQ_BAM_FILE_FOLDER: str = "bam_files"

# =============================================================================
# TOOL PATH CONFIGURATION
# =============================================================================

# These will be set from config or found automatically
SAMTOOLS_PATH: Optional[Path] = None
BCFTOOLS_PATH: Optional[Path] = None
GRAPHVIZ_PATH: Optional[Path] = None


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_path(name_: str, value_: str) -> Path:
    """
    Validate and return a Path from a config value.

    Args:
        name_: Configuration key name (for error messages)
        value_: Path string from config

    Returns:
        Validated Path object

    Raises:
        ValueError: If path is invalid or not a FASTA file
    """
    path = Path(value_)
    if not path.exists():
        if not path.parent.exists():
            raise ValueError(
                f"Cannot find provided config path ({path}) for: '{name_}'. "
                "Try to define an absolute path!"
            )
        else:
            # create an empty file at the location
            open(path, "w").close()
    if path.suffix not in [".fa", ".fasta", ".fna"]:
        raise ValueError("Please provide a fasta file. File ending in .fa, .fasta or .fna")
    return path


def get_tool_path(name_: str, value_: str) -> Optional[Path]:
    """
    Validate and return a Path to an executable tool.

    Args:
        name_: Tool name (for error messages)
        value_: Path string from config

    Returns:
        Path object if valid, None if empty or not found
    """
    if not value_ or value_.strip() == "":
        return None

    path = Path(value_)
    if path.exists():
        return path

    # On Windows, try adding .exe
    if IS_WINDOWS and not path.suffix:
        exe_path = Path(str(path) + ".exe")
        if exe_path.exists():
            return exe_path

    return None


def _load_config() -> None:
    """
    Load configuration from config.txt file.

    This function reads the config file and sets global path variables.
    It handles both reference genome paths and tool paths.
    """
    global HG19_FULL_GENOME, HG38_FULL_GENOME, T2T_FULL_GENOME
    global HG19_Y_CHROMOSOME, HG38_Y_CHROMOSOME, T2T_Y_CHROMOSOME
    global SAMTOOLS_PATH, BCFTOOLS_PATH, GRAPHVIZ_PATH

    if not CONFIG_PATH.exists():
        # Config file doesn't exist, use defaults
        return

    try:
        with open(CONFIG_PATH) as f:
            for line in f:
                line = line.strip()
                if not line or '=' not in line:
                    continue

                name, value = line.split('=', 1)
                value = value.strip()
                name = name.strip()

                if value == "":
                    continue

                # Reference genome paths
                if name == "full hg19 genome fasta location":
                    HG19_FULL_GENOME = get_path(name, value)
                elif name == "full hg38 genome fasta location":
                    HG38_FULL_GENOME = get_path(name, value)
                elif name == "hg19 chromosome Y fasta location":
                    HG19_Y_CHROMOSOME = get_path(name, value)
                elif name == "hg38 chromosome Y fasta location":
                    HG38_Y_CHROMOSOME = get_path(name, value)
                elif name == "full t2t genome fasta location":
                    T2T_FULL_GENOME = get_path(name, value)
                elif name == "t2t chromosome Y fasta location":
                    T2T_Y_CHROMOSOME = get_path(name, value)

                # Tool paths (new for Windows support)
                elif name == "samtools_path":
                    SAMTOOLS_PATH = get_tool_path(name, value)
                elif name == "bcftools_path":
                    BCFTOOLS_PATH = get_tool_path(name, value)
                elif name == "graphviz_path":
                    GRAPHVIZ_PATH = get_tool_path(name, value)

    except Exception as e:
        # Log warning but don't fail - use defaults
        logging.warning(f"Error reading config file: {e}")


def create_default_config() -> None:
    """
    Create a default config.txt file if it doesn't exist.

    This is useful for first-time setup on Windows.
    """
    if CONFIG_PATH.exists():
        return

    default_content = """# Yleaf Configuration File
#
# Reference genome locations (leave empty to use bundled references or download automatically)
full hg19 genome fasta location =
hg19 chromosome Y fasta location =
full hg38 genome fasta location =
hg38 chromosome Y fasta location =
full t2t genome fasta location =
t2t chromosome Y fasta location =

# Tool paths (leave empty to auto-detect from PATH or bundled binaries)
samtools_path =
bcftools_path =
graphviz_path =
"""

    try:
        with open(CONFIG_PATH, 'w') as f:
            f.write(default_content)
    except Exception as e:
        logging.warning(f"Could not create default config: {e}")


# =============================================================================
# INITIALIZATION
# =============================================================================

# Load configuration on module import
_load_config()
