import os
import sys
from pathlib import Path

# Files we are looking for
REF_MAP = {
    "hg19": ["hs37d5.fa", "hs37d5.fasta", "hg19.fa", "hg19.fasta"],
    "hg38": ["hg38.fa", "hg38.fasta", "GRCh38.fa", "GRCh38.fasta"],
    "t2t": ["chm13v2.0.fa", "chm13v2.0.fasta", "T2T-CHM13v2.0.fa"]
}

# Paths to search (recursive)
SEARCH_ROOTS = [
    Path.cwd(),
    Path.home() / "wgs",
    Path("/mnt/truenas-data/Data/DNA"),
    Path.home()
]

CONFIG_FILE = Path("yleaf/config.txt")

def find_file(names, search_paths):
    print(f"Searching for {names}...")
    for root in search_paths:
        if not root.exists(): continue
        
        # Try direct lookups in specific known subdirs first for speed
        wgsextract_ref = root / "WGSExtractv4" / "reference" / "genomes"
        if wgsextract_ref.exists():
            for name in names:
                f = wgsextract_ref / name
                if f.exists():
                    return f.absolute()

        # Shallow search in root
        for name in names:
            f = root / name
            if f.exists():
                return f.absolute()
                
        # Deep search (limited depth to avoid hanging)
        # Using 'rg' or 'find' command is faster, but let's use python walk with limit
        # Actually, let's just rely on the WGSExtract path or direct paths for now to be safe.
        # Users can add more logic here.
    return None

def main():
    if not CONFIG_FILE.exists():
        print(f"Config file not found at {CONFIG_FILE}. Run from Yleaf root.")
        sys.exit(1)

    found_refs = {}

    print("--- Auto-configuring References ---")
    
    for build, filenames in REF_MAP.items():
        path = find_file(filenames, SEARCH_ROOTS)
        if path:
            print(f"Found {build}: {path}")
            found_refs[build] = path
        else:
            print(f"Could not find reference for {build}")

    # Read existing config
    with open(CONFIG_FILE, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        key = line.split('=')[0].strip()
        
        # Update values if found
        if key == "full hg19 genome fasta location" and "hg19" in found_refs:
            new_lines.append(f"{key} = {found_refs['hg19']}\n")
        elif key == "hg19 chromosome Y fasta location" and "hg19" in found_refs:
            new_lines.append(f"{key} = {found_refs['hg19']}\n")
            
        elif key == "full hg38 genome fasta location" and "hg38" in found_refs:
            new_lines.append(f"{key} = {found_refs['hg38']}\n")
        elif key == "hg38 chromosome Y fasta location" and "hg38" in found_refs:
            new_lines.append(f"{key} = {found_refs['hg38']}\n")
            
        elif key == "full t2t genome fasta location" and "t2t" in found_refs:
            new_lines.append(f"{key} = {found_refs['t2t']}\n")
        elif key == "t2t chromosome Y fasta location" and "t2t" in found_refs:
            new_lines.append(f"{key} = {found_refs['t2t']}\n")
        else:
            new_lines.append(line)

    # Write back
    with open(CONFIG_FILE, 'w') as f:
        f.writelines(new_lines)
        
    print("--- Config updated ---")

if __name__ == "__main__":
    main()
