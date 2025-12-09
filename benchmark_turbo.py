import time
import pysam
import pandas as pd
import os
from pathlib import Path
import subprocess

# Config
BAM_FILE = "/home/valalav/wgs/output/72969_t2t_chrYM.bam"
POS_FILE = "Yleaf/yleaf/data/t2t/new_positions.txt"
LIMIT_SNPS = 5000 # Test on first 5000 SNPs first

def legacy_method(bam_path, bed_path):
    start = time.time()
    print("Running Legacy (Samtools)...")
    output_pu = "temp_bench.pu"
    
    # Only using -r chrY for fair comparison (we optimized this already)
    cmd = f"samtools mpileup -r chrY -l {bed_path} -AQ10q1 {bam_path} > {output_pu}"
    subprocess.check_call(cmd, shell=True)
    
    # Read back
    df = pd.read_csv(output_pu, sep="\t", header=None)
    end = time.time()
    
    os.remove(output_pu)
    print(f"Legacy Time: {end - start:.4f} sec. Rows: {len(df)}")
    return end - start

def turbo_method(bam_path, positions_df):
    start = time.time()
    print("Running Turbo (Pysam count_coverage)...")
    
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = []
    
    # Iterate over SNPs directly
    for idx, row in positions_df.iterrows():
        pos_0based = row['pos'] - 1
        # count_coverage returns tuple of arrays (A, C, G, T)
        # quality_threshold=10 mimics -Q10
        try:
            counts = bam.count_coverage("chrY", start=pos_0based, stop=pos_0based+1, quality_threshold=10, read_callback="all")
            # counts is usually tuple of arrays of length 1. e.g. ([10], [0], [0], [5])
            # Sum them to get reads
            total_reads = sum([c[0] for c in counts])
            results.append(total_reads)
        except (ValueError, KeyError) as e:
            # Skip positions that cannot be read (e.g., out of range, missing contig)
            print(f"  Warning: Skipping position {row['pos']}: {e}")
            
    end = time.time()
    print(f"Turbo Time: {end - start:.4f} sec. Processed: {len(results)}")
    return end - start

def create_temp_bed(pos_file, limit):
    df = pd.read_csv(pos_file, sep="\t", header=None, names=["chr", "name", "hg", "pos", "mut", "anc", "der"])
    if limit:
        df = df.head(limit)
    
    # Fix chromosome name for T2T BAM (chry -> chrY)
    # Ideally we should detect this from BAM header like Yleaf does, but for bench we hardcode "chrY"
    df['chr'] = "chrY"
    
    bed_path = "temp_bench.bed"
    # Create 1-col bed file format for samtools -l (chr pos) or just pos list?
    # Yleaf uses: chr <TAB> pos
    df[['chr', 'pos']].to_csv(bed_path, sep="\t", index=False, header=False)
    return bed_path, df

if __name__ == "__main__":
    print("Preparing benchmark...")
    # Ensure index exists
    if not os.path.exists(BAM_FILE + ".bai"):
        pysam.index(BAM_FILE)
        
    bed_file, pos_df = create_temp_bed(POS_FILE, LIMIT_SNPS)
    
    try:
        t_legacy = legacy_method(BAM_FILE, bed_file)
        t_turbo = turbo_method(BAM_FILE, pos_df)
        
        print(f"\nSpeedup: {t_legacy / t_turbo:.2f}x")
        
    finally:
        if os.path.exists(bed_file):
            os.remove(bed_file)
