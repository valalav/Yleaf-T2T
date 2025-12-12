# Changelog

All notable changes to the Yleaf project (Enhanced Version) will be documented in this file.

## [v3.2.1-enhanced] - 2025-12-12

### Critical Fixes & Database Updates
- **Database (T2T):** Removed marker `G-PF3346` (pos 23141720) from `new_positions.txt`.
    - *Reason:* This marker consistently showed Ancestral state in valid G-L1264+ samples, causing the prediction algorithm to stop early at G-CTS796. Removal allows traversing deeper branches.
- **Algorithm Sensitivity:** Lowered `DEFAULT_MAJORITY_THRESHOLD` from 90% to **75%**.
    - *Reason:* WGS data with low coverage (e.g. 8x) often has single-read errors. A ratio of 7/8 (87.5%) was previously rejected. 75% allows for 1 error in 4 reads, recovering many valid markers.
- **Algorithm Sensitivity:** Lowered `DEFAULT_READ_THRESHOLD` from 10 to **2**.
    - *Reason:* Essential for WGS data where key markers might have low coverage. Recovered ~200k markers for some samples.

### Tools & Scripts
- **`batch_process.py`:**
    - Added multi-threading support (`-t` argument).
    - Added Fast Fail logic for indices (checks validity before running Yleaf).
    - Fixed f-string formatting bug in logging.
    - Added auto-detection of Reference Genome for CRAM files to set `SAMTOOLS_CRAM_REF` environment variable automatically.
    - **Fix:** Skips `samtools idxstats` validation for CRAM files if local reference is not found or index is fresh, preventing infinite re-indexing loops on NAS.
- **`smart_index.sh`:**
    - New Bash script for running indexing directly on NAS/Docker containers to avoid network I/O bottlenecks.
    - Supports auto-detection of references for CRAMs.

### Core Yleaf Improvements
- **Refactoring:** Split monolithic `main()` function in `Yleaf.py` into modular helper functions for better maintainability.
- **Reporting:** Added detection for `Low_Y_Signal` (female/failed samples) in `summary_table.csv` when haplogroup is NA and derived markers < 1000.
- **Stability:** Removed unused `pysam` import that caused crashes on some systems.
- **Cross-Platform:** Replaced `shell=True` subprocess calls with `external_tools.py` wrapper for better security and Windows compatibility.

---

## [v3.2.0] - 2025-12-09
### Added
- Initial implementation of `batch_process.py`.
- `summary_logger.py` for aggregating results.
- `external_tools.py` abstraction layer.
