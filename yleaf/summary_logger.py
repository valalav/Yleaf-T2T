import os
import csv
import pandas as pd
from pathlib import Path
from datetime import datetime

# Path to the summary CSV file in the root of the Yleaf project (or current working dir)
# We assume the script is run from the root or we find the package root.
# For simplicity, we'll put it in the current working directory where the user runs the command, 
# or specifically in the Yleaf/ directory if we can detect it.

SUMMARY_FILENAME = "summary_table.csv"

def get_summary_path():
    # Try to find the 'Yleaf' source directory to keep the log central
    # But often it is better to keep it where the user is working.
    # Let's check if 'Yleaf' dir exists in current path (cloned repo structure)
    if os.path.isdir("Yleaf"):
        return Path("Yleaf") / SUMMARY_FILENAME
    return Path(SUMMARY_FILENAME)

def log_run(output_folder, input_file_path, ref_genome):
    output_path = Path(output_folder)
    hg_file = output_path / "hg_prediction.hg"
    
    if not hg_file.exists():
        return

    # Parse Prediction File
    data_dict = {}
    try:
        df = pd.read_csv(hg_file, sep="\t")
        if not df.empty:
            data_dict = df.iloc[0].to_dict()
    except Exception as e:
        print(f"Error reading HG prediction for summary: {e}")
        return

    # Extract info
    sample_name = data_dict.get("Sample_name", "Unknown")
    hg_raw = data_dict.get("Hg", "")
    marker = data_dict.get("Hg_marker", "")
    qc_score = data_dict.get("QC-score", 0)
    total_reads = data_dict.get("Total_reads", 0)
    
    # Process Haplogroup for cleanliness
    hg_clean = str(hg_raw).split("(x")[0].replace("*", "").strip() # remove exclusions and stars
    
    # Generate Link
    # Use full haplogroup name (e.g. I-Z2336) as YFull supports it and it's safer than just SNP
    link_target = hg_clean 
    yfull_link = f"https://www.yfull.com/tree/{link_target}/" if link_target else ""

    # Prepare Row
    row_data = {
        "Date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Sample Name": sample_name,
        "Input File": str(Path(input_file_path).absolute()),
        "Reference": ref_genome,
        "Haplogroup": hg_clean,
        "Terminal Marker": marker,
        "QC Score": qc_score,
        "Total Reads": total_reads,
        "YFull Link": yfull_link
    }

    fieldnames = ["Date", "Sample Name", "Input File", "Reference", "Haplogroup", 
                  "Terminal Marker", "QC Score", "Total Reads", "YFull Link"]

    csv_path = get_summary_path()
    file_exists = csv_path.exists()

    try:
        with open(csv_path, mode='a', newline='', encoding='utf-8') as f:
            # Use QUOTE_NONNUMERIC to quote all non-numeric fields.
            # This handles semicolons or commas inside fields safely.
            writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
            
            # If creating new file, write header
            if not file_exists:
                writer.writeheader()
            
            writer.writerow(row_data)
        print(f"Summary log updated: {csv_path}")
    except Exception as e:
        print(f"Failed to update summary log: {e}")
