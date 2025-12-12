import os
import glob
from pathlib import Path
import pandas as pd

OUTPUT_DIR = "output_wgs_full"
DATA_DIR = "/mnt/truenas-data/Data/DNA/wgs"
OUT_LIST = "rerun_list.txt"

files_to_rerun = set()

print("Scanning results...")
for hg_file in Path(OUTPUT_DIR).rglob("hg_prediction.hg"):
    try:
        df = pd.read_csv(hg_file, sep="\t")
        if df.empty: continue
        row = df.iloc[0]
        hg = str(row["Hg"])
        
        # Criteria: NA or Paragroup (*)
        if "NA" in hg or "*" in hg or "Low_Y_Signal" in hg:
            sample_name = str(row["Sample_name"])
            print(f"Candidate: {sample_name} ({hg})")
            
            # Search logic
            found = False
            # 1. Exact match BAM
            for p in Path(DATA_DIR).rglob(f"{sample_name}.bam"):
                files_to_rerun.add(str(p))
                found = True
                break
            
            if not found:
                # 2. Exact match CRAM
                for p in Path(DATA_DIR).rglob(f"{sample_name}.cram"):
                    files_to_rerun.add(str(p))
                    found = True
                    break
            
            if not found:
                # 3. Partial match (folder name usually contains sample name)
                # This is slower but necessary if sample_name in bam is different from filename
                # Actually, Yleaf uses filename as sample name usually.
                print(f"  -> WARNING: Could not find source file for {sample_name}")

    except Exception as e:
        print(f"Error reading {hg_file}: {e}")

with open(OUT_LIST, "w") as f:
    f.write("\n".join(list(files_to_rerun)))

print(f"Found {len(files_to_rerun)} files to rerun. Saved to {OUT_LIST}")
