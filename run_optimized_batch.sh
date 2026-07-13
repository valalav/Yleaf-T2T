#!/bin/bash
# run_optimized_batch.sh — Extracts chrY from heavy VCFs into a local folder and runs deepen.sh

set -euo pipefail

INPUT_DIR="${1:-/mnt/16tb/Biotech/VCF}"
LOG_FILE="/tmp/optimize_vcf.log"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEEPEN_SH="$SCRIPT_DIR/deepen.sh"
OUTPUT_BASE="$SCRIPT_DIR/vcf_chrY"

echo "=== Starting VCF Optimization and Deepening ===" > "$LOG_FILE"
echo "Logging to $LOG_FILE"

mkdir -p "$OUTPUT_BASE"

count=0
echo "Extracting chrY from VCFs in $INPUT_DIR to $OUTPUT_BASE..."
for DIR in "$INPUT_DIR"/*/; do
    VCF=$(ls "$DIR"*DeepVariant.vcf.gz 2>/dev/null | grep -v 'g\.vcf\.gz$' | head -n 1 || true)

    if [ -n "$VCF" ]; then
        BASENAME=$(basename "$VCF" .vcf.gz)
        SAMPLE_DIR=$(basename "$DIR")

        # Recreate directory structure needed by deepen.sh --batch
        TARGET_DIR="$OUTPUT_BASE/$SAMPLE_DIR"
        mkdir -p "$TARGET_DIR"

        OUT="$TARGET_DIR/${BASENAME}.chrY.vcf.gz"

        if [ ! -f "$OUT" ]; then
            echo "-> Extracting chrY: $(basename "$VCF")" | tee -a "$LOG_FILE"
            bcftools view -r chrY "$VCF" -O z -o "$OUT" 2>>"$LOG_FILE"
            bcftools index -t "$OUT" 2>>"$LOG_FILE"
        fi
        count=$((count + 1))
    fi
done

echo "=== Extraction finished. Processed $count files. ===" | tee -a "$LOG_FILE"

# Step 2: Run deepening batch
echo "=== Running deepen.sh on optimized chrY VCFs ===" | tee -a "$LOG_FILE"
cd "$SCRIPT_DIR"
export PATH=$(echo "$PATH" | sed -e 's|/home/valalav/_dna/wgs/venv/bin:||g')
bash "$DEEPEN_SH" --batch "$OUTPUT_BASE" --ext chrY.vcf.gz >> "$LOG_FILE" 2>&1

echo "=== All done. Results in reports/deepening_results.csv ===" | tee -a "$LOG_FILE"
