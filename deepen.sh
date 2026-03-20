#!/bin/bash
# deepen.sh — End-to-end: VCF/BAM → Yleaf → FTDNA deepening
#
# Single sample:
#   bash deepen.sh /path/to/sample.vcf.gz
#   bash deepen.sh /path/to/sample.bam
#
# Batch mode (directory of VCF/BAM files):
#   bash deepen.sh --batch /path/to/vcf_dir/
#   bash deepen.sh --batch /path/to/bam_dir/ --ext bam

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
YLEAF_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
YLEAF_PY="$SCRIPT_DIR/yleaf/Yleaf.py"
DEEPEN_JS="$SCRIPT_DIR/deepen_with_ftdna.js"
RESULTS_DIR="$YLEAF_DIR/results"
REPORTS_DIR="$YLEAF_DIR/reports"

# Load ref from config
REF=$(node -e "try{const c=JSON.parse(require('fs').readFileSync('$SCRIPT_DIR/deepen_config.json','utf8'));console.log(c.reference_genome.hg38||'')}catch(e){}" 2>/dev/null)

# Create directories
mkdir -p "$RESULTS_DIR" "$REPORTS_DIR"

usage() {
    echo "Usage:"
    echo "  bash deepen.sh <vcf.gz|bam|cram>       — single sample"
    echo "  bash deepen.sh --batch <dir> [--ext bam] — batch mode"
    echo ""
    echo "Reports saved to: $REPORTS_DIR/"
    exit 1
}

[ $# -eq 0 ] && usage

# --- Batch mode ---
if [ "$1" = "--batch" ]; then
    BATCH_DIR="${2:-}"
    EXT="${4:-vcf.gz}"  # after --ext
    [ "$3" = "--ext" ] 2>/dev/null && EXT="$4"
    [ -z "$BATCH_DIR" ] && usage
    [ ! -d "$BATCH_DIR" ] && { echo "ERROR: $BATCH_DIR is not a directory"; exit 1; }

    # Clear previous batch CSV
    rm -f "$REPORTS_DIR/deepening_results.csv"
    
    COUNT=0
    echo "=== Batch deepening: $BATCH_DIR (*.$EXT) ==="
    echo ""
    for DIR in "$BATCH_DIR"/*/; do
        INPUT=$(ls "$DIR"*."$EXT" 2>/dev/null | head -1)
        [ -z "$INPUT" ] && continue
        COUNT=$((COUNT + 1))
        echo "[$COUNT] $(basename "$DIR")"
        bash "$0" "$INPUT" 2>&1 | grep -E "^(Yleaf|Clean|Original|Deepened|Depth|FTDNA path|YFull path|  Reports|  CSV):" || true
        echo ""
    done
    echo "=== Done: $COUNT samples processed ==="
    echo "CSV: $REPORTS_DIR/deepening_results.csv"
    exit 0
fi

# --- Single sample mode ---
INPUT="$1"
[ ! -f "$INPUT" ] && { echo "ERROR: File not found: $INPUT"; exit 1; }

SAMPLE_ID=$(basename "$INPUT" | grep -oP '^\d{12}' || true)
BASENAME=$(basename "$INPUT" | sed -E 's/\.(DeepVariant_v1\.6\.vcf\.gz|vcf\.gz|vcf|bam|cram)$//')

echo "=== Sample: ${SAMPLE_ID:-$BASENAME} ==="

# --- Step 1: Get Yleaf prediction ---
PRED_FILE=""
for dir in $(ls -d "$RESULTS_DIR"/output_${SAMPLE_ID:-NOSAMPLE}* 2>/dev/null); do
    [ -f "$dir/hg_prediction.hg" ] && PRED_FILE="$dir/hg_prediction.hg" && break
done
[ -z "$PRED_FILE" ] && [ -f "$RESULTS_DIR/output_${BASENAME}/hg_prediction.hg" ] && \
    PRED_FILE="$RESULTS_DIR/output_${BASENAME}/hg_prediction.hg"

if [ -z "$PRED_FILE" ]; then
    echo "Running Yleaf..."
    OUTPUT_DIR="$RESULTS_DIR/output_${BASENAME}"
    mkdir -p "$OUTPUT_DIR"
    export PYTHONPATH="$SCRIPT_DIR:${PYTHONPATH:-}"
    cd "$YLEAF_DIR"
    python3 "$YLEAF_PY" -bam "$INPUT" -o "$OUTPUT_DIR" 2>&1 | tail -3
    PRED_FILE="$OUTPUT_DIR/hg_prediction.hg"
    [ ! -f "$PRED_FILE" ] && { echo "ERROR: Yleaf failed"; exit 1; }
else
    echo "(Yleaf cached)"
fi

RAW_HG=$(awk 'NR==2 {print $2}' "$PRED_FILE")
HG=$(echo "$RAW_HG" | sed 's/\*(x.*//')
echo "Yleaf:    $RAW_HG"
echo "Clean:    $HG"
echo ""

if echo "$HG" | grep -qE "^(NA|Low_Y)"; then
    echo "⚠ Low Y signal — cannot deepen"
    exit 0
fi

# --- Step 2: Deepen with FTDNA ---
echo "=== Deepening with FTDNA ==="
EXT_LOW=$(echo "$INPUT" | tr '[:upper:]' '[:lower:]')
if echo "$EXT_LOW" | grep -qE '\.(bam|cram)$'; then
    if [ -z "$REF" ]; then
        echo "ERROR: BAM/CRAM requires a reference genome."
        echo "Set reference_genome.hg38 in deepen_config.json"
        exit 1
    fi
    node "$DEEPEN_JS" -i "$INPUT" --ref "$REF" --hg "$HG" -v
else
    node "$DEEPEN_JS" -i "$INPUT" --hg "$HG" -v
fi
