#!/bin/bash

# Script to intelligently index BAM/CRAM files.
# Usage: ./smart_index.sh /path/to/scan

LOG_FILE="smart_index.log"

main() {
    # Target directory
    SEARCH_DIR="${1:-.}"

    # --- CONFIGURATION: REFERENCE GENOME PATHS ---
    # Paths inside the NAS/Container
    REF_T2T="/data/ref/chm13v2.0.fa"
    REF_HG38="/data/ref/hg38.fa"
    REF_HG19="/data/ref/hs37d5.fa"
    # ---------------------------------------------

    echo "Starting Smart Indexer at $(date)"
    echo "Scanning $SEARCH_DIR for BAM/CRAM files..."

    find "$SEARCH_DIR" -type f \( -name "*.bam" -o -name "*.cram" \) | while read -r FILE; do
        echo "---------------------------------------------------"
        echo "Checking: $FILE"
        
        EXT="${FILE##*.}"
        if [ "$EXT" == "bam" ]; then
            INDEX_FILE="${FILE}.bai"
        else
            INDEX_FILE="${FILE}.crai"
        fi

        NEED_INDEX=0

        # 1. Check Existence
        if [ ! -f "$INDEX_FILE" ]; then
            echo "Status: Index MISSING."
            NEED_INDEX=1
        else
            # 2. Check Timestamp (Index must be newer than data)
            if [ "$INDEX_FILE" -ot "$FILE" ]; then
                echo "Status: Index OUTDATED."
                NEED_INDEX=1
            else
                # 3. Check Validity (idxstats)
                SKIP_IDXSTATS=0
                export SAMTOOLS_CRAM_REF=""
                
                if [ "$EXT" == "cram" ]; then
                    HEADER=$(samtools view -H "$FILE" 2>/dev/null | grep "@SQ" | head -n 1)
                    LEN=$(echo "$HEADER" | grep -o "LN:[0-9]*" | cut -d: -f2)
                    
                    if [ "$LEN" == "248387328" ]; then
                        export SAMTOOLS_CRAM_REF="$REF_T2T"
                    elif [ "$LEN" == "248956422" ]; then
                        export SAMTOOLS_CRAM_REF="$REF_HG38"
                    elif [ "$LEN" == "249250621" ]; then
                        export SAMTOOLS_CRAM_REF="$REF_HG19"
                    fi
                    
                    # CRITICAL FIX: If we cannot find the reference locally, DO NOT run idxstats.
                    # idxstats will fail/hang without ref, causing false positive "CORRUPTED" status.
                    if [ -z "$SAMTOOLS_CRAM_REF" ] || [ ! -f "$SAMTOOLS_CRAM_REF" ]; then
                        echo "Status: Index exists (Skipping validation - no local ref)."
                        SKIP_IDXSTATS=1
                    fi
                fi
                
                if [ "$SKIP_IDXSTATS" -eq 0 ]; then
                    if ! timeout 10s samtools idxstats "$FILE" > /dev/null 2>&1; then
                        echo "Status: Index CORRUPTED (idxstats failed)."
                        NEED_INDEX=1
                    else
                        echo "Status: Index OK."
                    fi
                fi
            fi
        fi

        # Re-indexing
        if [ "$NEED_INDEX" -eq 1 ]; then
            echo "Action: Re-indexing..."
            
            [ -f "$INDEX_FILE" ] && rm "$INDEX_FILE"
            
            if [ "$EXT" == "cram" ]; then
                 # Logic to set ref again for indexing command
                 HEADER=$(samtools view -H "$FILE" 2>/dev/null | grep "@SQ" | head -n 1)
                 LEN=$(echo "$HEADER" | grep -o "LN:[0-9]*" | cut -d: -f2)
                 
                 REF_PATH=""
                 if [ "$LEN" == "248387328" ]; then REF_PATH="$REF_T2T"; fi
                 if [ "$LEN" == "248956422" ]; then REF_PATH="$REF_HG38"; fi
                 if [ "$LEN" == "249250621" ]; then REF_PATH="$REF_HG19"; fi
                 
                 if [ -n "$REF_PATH" ] && [ -f "$REF_PATH" ]; then
                     export SAMTOOLS_CRAM_REF="$REF_PATH"
                     echo "  -> Using reference: $REF_PATH"
                 else
                     echo "  -> Warning: Local reference not found. Indexing will rely on network/cache."
                 fi
            fi

            if samtools index "$FILE"; then
                echo "Result: SUCCESS."
            else
                echo "Result: FAILED."
            fi
        fi
    done

    echo "Batch indexing completed at $(date)"
}

# Run main and redirect output to log (and console)
main "$@" 2>&1 | tee "$LOG_FILE"
