"""
Analysis Worker Thread for Yleaf GUI

This module provides a QThread-based worker for running Yleaf analysis
in the background, allowing the GUI to remain responsive.
"""

import sys
import logging
import traceback
from pathlib import Path
from typing import Dict, Any, Optional
import argparse

from PyQt6.QtCore import QThread, pyqtSignal


class LogHandler(logging.Handler):
    """Custom logging handler that emits signals for GUI updates."""

    def __init__(self, signal):
        super().__init__()
        self.signal = signal

    def emit(self, record):
        msg = self.format(record)
        self.signal.emit(msg)


class AnalysisWorker(QThread):
    """
    Worker thread for running Yleaf analysis.

    Signals:
        progress_updated(int, str): Emitted when progress changes (percent, message)
        log_message(str): Emitted for log messages
        analysis_finished(dict): Emitted when analysis completes successfully
        error_occurred(str): Emitted when an error occurs
    """

    progress_updated = pyqtSignal(int, str)
    log_message = pyqtSignal(str)
    analysis_finished = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)

    def __init__(self, settings: Dict[str, Any], parent=None):
        """
        Initialize the worker.

        Args:
            settings: Dictionary containing analysis settings:
                - bam_file: Path to BAM/CRAM file
                - vcf_file: Path to VCF file
                - output_folder: Output directory path
                - reference_genome: Reference genome (hg19/hg38/t2t or None)
                - threads: Number of threads
                - reads_threshold: Minimum reads threshold
                - quality_threshold: Minimum quality threshold
                - base_majority: Base majority percentage
                - prediction_quality: Prediction quality threshold
                - draw_haplogroups: Whether to draw haplogroup tree
                - ancient_dna: Whether sample is ancient DNA
                - use_old: Whether to use old prediction method
            parent: Parent QObject
        """
        super().__init__(parent)
        self.settings = settings
        self._cancelled = False

    def cancel(self):
        """Request cancellation of the analysis."""
        self._cancelled = True

    def run(self):
        """Run the analysis in a separate thread."""
        try:
            self._run_analysis()
        except Exception as e:
            error_msg = f"{str(e)}\n\n{traceback.format_exc()}"
            self.error_occurred.emit(error_msg)

    def _run_analysis(self):
        """Perform the actual analysis."""
        # Setup logging to capture Yleaf output
        self._setup_logging()

        self.progress_updated.emit(5, "Initializing analysis...")
        self.log_message.emit("Starting Yleaf analysis...")

        # Import yleaf modules
        try:
            from yleaf import Yleaf
            from yleaf import yleaf_constants
            from yleaf import external_tools
        except ImportError as e:
            self.error_occurred.emit(f"Failed to import Yleaf modules: {e}")
            return

        # Check for required tools
        self.progress_updated.emit(10, "Checking dependencies...")
        self.log_message.emit("Checking for samtools and bcftools...")

        tool_status = external_tools.check_tools()
        for tool, info in tool_status.items():
            if not info["found"]:
                self.error_occurred.emit(
                    f"Required tool not found: {tool}\n\n"
                    f"Please install {tool} or configure its path in config.txt"
                )
                return
            self.log_message.emit(f"Found {tool} at: {info['path']}")

        if self._cancelled:
            self.log_message.emit("Analysis cancelled by user.")
            return

        # Prepare output directory
        output_folder = Path(self.settings["output_folder"])
        output_folder.mkdir(parents=True, exist_ok=True)

        self.progress_updated.emit(15, "Preparing analysis...")

        # Build arguments namespace (simulating command-line args)
        args = self._build_args_namespace()

        if self._cancelled:
            self.log_message.emit("Analysis cancelled by user.")
            return

        # Run the main Yleaf workflow
        try:
            self.progress_updated.emit(20, "Running analysis...")

            if args.vcffile:
                self.log_message.emit(f"Processing VCF file: {args.vcffile}")
                self._run_vcf_analysis(args, output_folder)
            elif args.bamfile:
                self.log_message.emit(f"Processing BAM file: {args.bamfile}")
                self._run_bam_analysis(args, output_folder)
            elif args.cramfile:
                self.log_message.emit(f"Processing CRAM file: {args.cramfile}")
                self._run_cram_analysis(args, output_folder)

            if self._cancelled:
                self.log_message.emit("Analysis cancelled by user.")
                return

            # Run haplogroup prediction
            self.progress_updated.emit(70, "Predicting haplogroups...")
            self._run_prediction(args, output_folder)

            if self._cancelled:
                self.log_message.emit("Analysis cancelled by user.")
                return

            # Draw haplogroups if requested
            if args.draw_haplogroups:
                self.progress_updated.emit(85, "Drawing haplogroup tree...")
                self._draw_haplogroups(args, output_folder)

            # Generate HTML report
            self.progress_updated.emit(95, "Generating report...")
            self._generate_report(output_folder)

            # Collect results
            results = self._collect_results(output_folder)
            results["output_folder"] = str(output_folder)

            self.progress_updated.emit(100, "Analysis complete!")
            self.log_message.emit("Analysis completed successfully!")
            self.analysis_finished.emit(results)

        except SystemExit as e:
            self.error_occurred.emit(f"Analysis failed: {str(e)}")
        except Exception as e:
            self.error_occurred.emit(f"Unexpected error: {str(e)}\n\n{traceback.format_exc()}")

    def _setup_logging(self):
        """Setup logging to capture Yleaf messages."""
        # Get the yleaf logger
        yleaf_logger = logging.getLogger("yleaf_logger")

        # Add our custom handler
        handler = LogHandler(self.log_message)
        handler.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
        handler.setLevel(logging.INFO)
        yleaf_logger.addHandler(handler)

    def _build_args_namespace(self) -> argparse.Namespace:
        """Build argparse.Namespace from settings."""
        bam_path = self.settings.get("bam_file")
        vcf_path = self.settings.get("vcf_file")

        # Determine input type
        bamfile = None
        cramfile = None
        vcffile = None

        if bam_path:
            path = Path(bam_path)
            if path.suffix.lower() == ".cram":
                cramfile = path
            else:
                bamfile = path
        elif vcf_path:
            vcffile = Path(vcf_path)

        args = argparse.Namespace(
            bamfile=bamfile,
            cramfile=cramfile,
            vcffile=vcffile,
            cram_reference=None,
            fastq=None,
            output=str(self.settings["output_folder"]),
            reference_genome=self.settings.get("reference_genome"),
            ref_fasta=None,
            threads=self.settings.get("threads", 1),
            reads_treshold=self.settings.get("reads_threshold", 10),
            quality_thresh=self.settings.get("quality_threshold", 20),
            base_majority=self.settings.get("base_majority", 90),
            prediction_quality=self.settings.get("prediction_quality", 0.95),
            draw_haplogroups=self.settings.get("draw_haplogroups", True),
            collapsed_draw_mode=False,
            ancient_DNA=self.settings.get("ancient_dna", False),
            use_old=self.settings.get("use_old", False),
            private_mutations=False,
            minor_allele_frequency=0.01,
            force=True,
            reanalyze=False,
        )

        return args

    def _run_vcf_analysis(self, args, output_folder: Path):
        """Run VCF-based analysis."""
        from yleaf.Yleaf import main_vcf
        main_vcf(args, output_folder)
        self.progress_updated.emit(60, "VCF processing complete...")

    def _run_bam_analysis(self, args, output_folder: Path):
        """Run BAM-based analysis."""
        from yleaf.Yleaf import main_bam_cram
        main_bam_cram(args, output_folder, is_bam=True)
        self.progress_updated.emit(60, "BAM processing complete...")

    def _run_cram_analysis(self, args, output_folder: Path):
        """Run CRAM-based analysis."""
        from yleaf.Yleaf import main_bam_cram, get_reference_path
        from yleaf import yleaf_constants

        # Try to auto-detect CRAM reference
        if args.cram_reference is None and args.reference_genome:
            ref_path = get_reference_path(args.reference_genome, is_full=True)
            if ref_path and ref_path.exists() and ref_path.stat().st_size > 100:
                args.cram_reference = ref_path
                self.log_message.emit(f"Auto-detected CRAM reference: {ref_path}")

        main_bam_cram(args, output_folder, is_bam=False)
        self.progress_updated.emit(60, "CRAM processing complete...")

    def _run_prediction(self, args, output_folder: Path):
        """Run haplogroup prediction."""
        from yleaf.Yleaf import predict_haplogroup, PREDICTION_OUT_FILE_NAME

        hg_out = output_folder / PREDICTION_OUT_FILE_NAME
        predict_haplogroup(
            output_folder,
            hg_out,
            args.use_old,
            args.prediction_quality,
            args.threads
        )
        self.log_message.emit("Haplogroup prediction complete.")

    def _draw_haplogroups(self, args, output_folder: Path):
        """Draw haplogroup tree."""
        from yleaf.Yleaf import draw_haplogroups, PREDICTION_OUT_FILE_NAME

        hg_out = output_folder / PREDICTION_OUT_FILE_NAME
        if hg_out.exists():
            draw_haplogroups(hg_out, args.collapsed_draw_mode)
            self.log_message.emit("Haplogroup tree generated.")
        else:
            self.log_message.emit("Warning: Prediction file not found, skipping tree drawing.")

    def _generate_report(self, output_folder: Path):
        """Generate HTML report."""
        try:
            from yleaf import html_report
            html_report.generate_html(output_folder)
            self.log_message.emit("HTML report generated.")
        except Exception as e:
            self.log_message.emit(f"Warning: Failed to generate HTML report: {e}")

    def _collect_results(self, output_folder: Path) -> Dict[str, Any]:
        """Collect results from output files."""
        from yleaf.Yleaf import PREDICTION_OUT_FILE_NAME

        results = {
            "samples": [],
            "output_folder": str(output_folder)
        }

        # Parse prediction file
        hg_file = output_folder / PREDICTION_OUT_FILE_NAME
        if hg_file.exists():
            try:
                with open(hg_file) as f:
                    lines = f.readlines()

                # Skip header if present
                start_idx = 0
                if lines and lines[0].startswith("sample"):
                    start_idx = 1

                for line in lines[start_idx:]:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        results["samples"].append({
                            "name": parts[0],
                            "haplogroup": parts[1],
                            "qc_score": parts[2] if len(parts) > 2 else "",
                            "markers": parts[3] if len(parts) > 3 else ""
                        })
            except Exception as e:
                self.log_message.emit(f"Warning: Could not parse results file: {e}")

        return results
