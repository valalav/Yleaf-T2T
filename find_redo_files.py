import os
from pathlib import Path

TARGETS = [
    "TSBC5520_bwa_aligned_chrYM",
    "TSAC3435.YM",
    "68284_bwa-mem_t2tl_chrYM",
    "68129_t2t_YM",
    "63864_bwa-mem_hg38_sorted_final_chrYM",
    "55477.YM",
    "69983_final_chrYM"
]

DATA_DIR = "/mnt/truenas-data/Data/DNA/wgs"
OUT_LIST = "redo_list.txt"

found_paths = []

print("Scanning for targets...")
# Use set to avoid duplicates if multiple matches
found_set = set()

for root, dirs, files in os.walk(DATA_DIR):
    for f in files:
        if f.endswith(".bam") or f.endswith(".cram"):
            for t in TARGETS:
                # Match filename stem or full name
                if t == f or t == Path(f).stem:
                    full_path = os.path.join(root, f)
                    if full_path not in found_set:
                        print(f"Found: {full_path}")
                        found_set.add(full_path)

with open(OUT_LIST, "w") as f:
    f.write("\n".join(list(found_set)))

print(f"Saved {len(found_set)} paths to {OUT_LIST}")
