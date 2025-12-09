"""
Yleaf GUI Main Window

This module provides the main application window for the Yleaf GUI.
"""

import os
import sys
from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QSpinBox, QDoubleSpinBox,
    QProgressBar, QTextEdit, QFileDialog, QGroupBox, QTabWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QMessageBox, QStatusBar,
    QCheckBox, QSplitter
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont, QAction, QIcon

# Import workers
from yleaf_gui.workers.analysis_worker import AnalysisWorker


class YleafMainWindow(QMainWindow):
    """Main window for Yleaf GUI application."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Yleaf - Y-Chromosome Haplogroup Prediction")
        self.setMinimumSize(900, 700)

        # Worker thread
        self.worker: Optional[AnalysisWorker] = None

        # Create UI
        self._create_menu_bar()
        self._create_central_widget()
        self._create_status_bar()

        # Connect signals
        self._connect_signals()

    def _create_menu_bar(self):
        """Create the application menu bar."""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("&File")

        open_action = QAction("&Open BAM/VCF...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._open_file_dialog)
        file_menu.addAction(open_action)

        file_menu.addSeparator()

        exit_action = QAction("E&xit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Tools menu
        tools_menu = menubar.addMenu("&Tools")

        check_tools_action = QAction("&Check Dependencies...", self)
        check_tools_action.triggered.connect(self._check_dependencies)
        tools_menu.addAction(check_tools_action)

        # Help menu
        help_menu = menubar.addMenu("&Help")

        about_action = QAction("&About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    def _create_central_widget(self):
        """Create the main central widget with all controls."""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QVBoxLayout(central_widget)

        # Create splitter for resizable sections
        splitter = QSplitter(Qt.Orientation.Vertical)

        # Top section: Input and Settings
        top_widget = QWidget()
        top_layout = QVBoxLayout(top_widget)

        # File Selection Group
        file_group = self._create_file_selection_group()
        top_layout.addWidget(file_group)

        # Settings Group
        settings_group = self._create_settings_group()
        top_layout.addWidget(settings_group)

        # Control Buttons
        button_layout = self._create_control_buttons()
        top_layout.addLayout(button_layout)

        splitter.addWidget(top_widget)

        # Bottom section: Progress and Results
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout(bottom_widget)

        # Progress Group
        progress_group = self._create_progress_group()
        bottom_layout.addWidget(progress_group)

        # Results Tabs
        results_tabs = self._create_results_tabs()
        bottom_layout.addWidget(results_tabs)

        splitter.addWidget(bottom_widget)

        main_layout.addWidget(splitter)

    def _create_file_selection_group(self) -> QGroupBox:
        """Create the file selection group box."""
        group = QGroupBox("Input Files")
        layout = QGridLayout(group)

        # BAM/CRAM file
        layout.addWidget(QLabel("BAM/CRAM File:"), 0, 0)
        self.bam_input = QLineEdit()
        self.bam_input.setPlaceholderText("Select a BAM or CRAM file...")
        layout.addWidget(self.bam_input, 0, 1)
        self.bam_browse_btn = QPushButton("Browse...")
        self.bam_browse_btn.clicked.connect(lambda: self._browse_file("bam"))
        layout.addWidget(self.bam_browse_btn, 0, 2)

        # VCF file
        layout.addWidget(QLabel("VCF File:"), 1, 0)
        self.vcf_input = QLineEdit()
        self.vcf_input.setPlaceholderText("Select a VCF.GZ file...")
        layout.addWidget(self.vcf_input, 1, 1)
        self.vcf_browse_btn = QPushButton("Browse...")
        self.vcf_browse_btn.clicked.connect(lambda: self._browse_file("vcf"))
        layout.addWidget(self.vcf_browse_btn, 1, 2)

        # Output folder
        layout.addWidget(QLabel("Output Folder:"), 2, 0)
        self.output_input = QLineEdit()
        self.output_input.setPlaceholderText("Select output folder...")
        layout.addWidget(self.output_input, 2, 1)
        self.output_browse_btn = QPushButton("Browse...")
        self.output_browse_btn.clicked.connect(self._browse_output_folder)
        layout.addWidget(self.output_browse_btn, 2, 2)

        return group

    def _create_settings_group(self) -> QGroupBox:
        """Create the settings group box."""
        group = QGroupBox("Analysis Settings")
        layout = QGridLayout(group)

        # Reference genome
        layout.addWidget(QLabel("Reference Genome:"), 0, 0)
        self.reference_combo = QComboBox()
        self.reference_combo.addItems(["Auto-detect", "hg19", "hg38", "t2t"])
        layout.addWidget(self.reference_combo, 0, 1)

        # Threads
        layout.addWidget(QLabel("Threads:"), 0, 2)
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, os.cpu_count() or 4)
        self.threads_spin.setValue(min(4, os.cpu_count() or 1))
        layout.addWidget(self.threads_spin, 0, 3)

        # Read threshold
        layout.addWidget(QLabel("Min Reads:"), 1, 0)
        self.reads_spin = QSpinBox()
        self.reads_spin.setRange(1, 100)
        self.reads_spin.setValue(10)
        layout.addWidget(self.reads_spin, 1, 1)

        # Quality threshold
        layout.addWidget(QLabel("Min Quality:"), 1, 2)
        self.quality_spin = QSpinBox()
        self.quality_spin.setRange(10, 40)
        self.quality_spin.setValue(20)
        layout.addWidget(self.quality_spin, 1, 3)

        # Base majority
        layout.addWidget(QLabel("Base Majority %:"), 2, 0)
        self.majority_spin = QSpinBox()
        self.majority_spin.setRange(50, 99)
        self.majority_spin.setValue(90)
        layout.addWidget(self.majority_spin, 2, 1)

        # Prediction quality
        layout.addWidget(QLabel("Prediction Quality:"), 2, 2)
        self.pred_quality_spin = QDoubleSpinBox()
        self.pred_quality_spin.setRange(0.0, 1.0)
        self.pred_quality_spin.setSingleStep(0.05)
        self.pred_quality_spin.setValue(0.95)
        layout.addWidget(self.pred_quality_spin, 2, 3)

        # Options row
        layout.addWidget(QLabel("Options:"), 3, 0)
        options_layout = QHBoxLayout()

        self.draw_tree_check = QCheckBox("Draw Haplogroup Tree")
        self.draw_tree_check.setChecked(True)
        options_layout.addWidget(self.draw_tree_check)

        self.ancient_dna_check = QCheckBox("Ancient DNA")
        options_layout.addWidget(self.ancient_dna_check)

        self.use_old_check = QCheckBox("Use Old Method (ISOGG)")
        options_layout.addWidget(self.use_old_check)

        options_layout.addStretch()
        layout.addLayout(options_layout, 3, 1, 1, 3)

        return group

    def _create_control_buttons(self) -> QHBoxLayout:
        """Create the control buttons layout."""
        layout = QHBoxLayout()

        self.start_btn = QPushButton("Start Analysis")
        self.start_btn.setMinimumHeight(40)
        self.start_btn.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                font-weight: bold;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:disabled {
                background-color: #cccccc;
            }
        """)
        self.start_btn.clicked.connect(self._start_analysis)
        layout.addWidget(self.start_btn)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setMinimumHeight(40)
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.setStyleSheet("""
            QPushButton {
                background-color: #f44336;
                color: white;
                font-weight: bold;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #da190b;
            }
            QPushButton:disabled {
                background-color: #cccccc;
            }
        """)
        self.cancel_btn.clicked.connect(self._cancel_analysis)
        layout.addWidget(self.cancel_btn)

        return layout

    def _create_progress_group(self) -> QGroupBox:
        """Create the progress display group."""
        group = QGroupBox("Progress")
        layout = QVBoxLayout(group)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        # Log output
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setMaximumHeight(150)
        self.log_output.setFont(QFont("Consolas", 9))
        layout.addWidget(self.log_output)

        return group

    def _create_results_tabs(self) -> QTabWidget:
        """Create the results tab widget."""
        tabs = QTabWidget()

        # Results table tab
        results_widget = QWidget()
        results_layout = QVBoxLayout(results_widget)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(4)
        self.results_table.setHorizontalHeaderLabels([
            "Sample", "Haplogroup", "QC Score", "Markers"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        results_layout.addWidget(self.results_table)

        # Export buttons
        export_layout = QHBoxLayout()
        self.export_csv_btn = QPushButton("Export CSV")
        self.export_csv_btn.clicked.connect(self._export_csv)
        export_layout.addWidget(self.export_csv_btn)

        self.open_html_btn = QPushButton("Open HTML Report")
        self.open_html_btn.clicked.connect(self._open_html_report)
        export_layout.addWidget(self.open_html_btn)

        self.view_tree_btn = QPushButton("View Haplogroup Tree")
        self.view_tree_btn.clicked.connect(self._view_tree)
        export_layout.addWidget(self.view_tree_btn)

        export_layout.addStretch()
        results_layout.addLayout(export_layout)

        tabs.addTab(results_widget, "Results")

        # Details tab (raw output)
        details_widget = QWidget()
        details_layout = QVBoxLayout(details_widget)

        self.details_text = QTextEdit()
        self.details_text.setReadOnly(True)
        self.details_text.setFont(QFont("Consolas", 9))
        details_layout.addWidget(self.details_text)

        tabs.addTab(details_widget, "Details")

        return tabs

    def _create_status_bar(self):
        """Create the status bar."""
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

    def _connect_signals(self):
        """Connect additional signals."""
        pass

    # =========================================================================
    # Slot methods
    # =========================================================================

    def _browse_file(self, file_type: str):
        """Open file browser for input files."""
        if file_type == "bam":
            file_filter = "BAM/CRAM Files (*.bam *.cram);;All Files (*)"
            line_edit = self.bam_input
        else:
            file_filter = "VCF Files (*.vcf.gz *.vcf);;All Files (*)"
            line_edit = self.vcf_input

        file_path, _ = QFileDialog.getOpenFileName(
            self, f"Select {file_type.upper()} File", "", file_filter
        )

        if file_path:
            line_edit.setText(file_path)

    def _browse_output_folder(self):
        """Open folder browser for output."""
        folder = QFileDialog.getExistingDirectory(
            self, "Select Output Folder"
        )
        if folder:
            self.output_input.setText(folder)

    def _open_file_dialog(self):
        """Open file dialog from menu."""
        self._browse_file("bam")

    def _start_analysis(self):
        """Start the haplogroup analysis."""
        # Validate inputs
        bam_path = self.bam_input.text().strip()
        vcf_path = self.vcf_input.text().strip()
        output_path = self.output_input.text().strip()

        if not bam_path and not vcf_path:
            QMessageBox.warning(
                self, "Input Required",
                "Please select a BAM/CRAM or VCF file."
            )
            return

        if not output_path:
            QMessageBox.warning(
                self, "Output Required",
                "Please select an output folder."
            )
            return

        # Get settings
        reference = self.reference_combo.currentText()
        if reference == "Auto-detect":
            reference = None

        settings = {
            "bam_file": bam_path if bam_path else None,
            "vcf_file": vcf_path if vcf_path else None,
            "output_folder": output_path,
            "reference_genome": reference,
            "threads": self.threads_spin.value(),
            "reads_threshold": self.reads_spin.value(),
            "quality_threshold": self.quality_spin.value(),
            "base_majority": self.majority_spin.value(),
            "prediction_quality": self.pred_quality_spin.value(),
            "draw_haplogroups": self.draw_tree_check.isChecked(),
            "ancient_dna": self.ancient_dna_check.isChecked(),
            "use_old": self.use_old_check.isChecked(),
        }

        # Clear previous results
        self.log_output.clear()
        self.results_table.setRowCount(0)
        self.progress_bar.setValue(0)

        # Create and start worker
        self.worker = AnalysisWorker(settings)
        self.worker.progress_updated.connect(self._on_progress_updated)
        self.worker.log_message.connect(self._on_log_message)
        self.worker.analysis_finished.connect(self._on_analysis_finished)
        self.worker.error_occurred.connect(self._on_error_occurred)

        # Update UI state
        self.start_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.status_bar.showMessage("Analysis in progress...")

        # Start worker
        self.worker.start()

    def _cancel_analysis(self):
        """Cancel the running analysis."""
        if self.worker and self.worker.isRunning():
            self.worker.cancel()
            self.log_output.append("Cancelling analysis...")

    def _on_progress_updated(self, value: int, message: str):
        """Handle progress update signal."""
        self.progress_bar.setValue(value)
        self.status_bar.showMessage(message)

    def _on_log_message(self, message: str):
        """Handle log message signal."""
        self.log_output.append(message)
        # Auto-scroll to bottom
        scrollbar = self.log_output.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def _on_analysis_finished(self, results: dict):
        """Handle analysis completion."""
        self.start_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress_bar.setValue(100)
        self.status_bar.showMessage("Analysis complete!")

        # Populate results table
        if "samples" in results:
            for sample in results["samples"]:
                row = self.results_table.rowCount()
                self.results_table.insertRow(row)
                self.results_table.setItem(row, 0, QTableWidgetItem(sample.get("name", "")))
                self.results_table.setItem(row, 1, QTableWidgetItem(sample.get("haplogroup", "")))
                self.results_table.setItem(row, 2, QTableWidgetItem(str(sample.get("qc_score", ""))))
                self.results_table.setItem(row, 3, QTableWidgetItem(str(sample.get("markers", ""))))

        # Show completion message
        QMessageBox.information(
            self, "Analysis Complete",
            f"Haplogroup prediction completed successfully!\n\n"
            f"Results saved to: {results.get('output_folder', 'N/A')}"
        )

    def _on_error_occurred(self, error: str):
        """Handle error signal."""
        self.start_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.status_bar.showMessage("Error occurred")

        QMessageBox.critical(
            self, "Error",
            f"An error occurred during analysis:\n\n{error}"
        )

    def _check_dependencies(self):
        """Check and display dependency status."""
        from yleaf import external_tools

        status = external_tools.check_tools()
        versions = external_tools.get_tool_versions()

        message = "Dependency Status:\n\n"
        for tool, info in status.items():
            status_str = "OK" if info["found"] else "NOT FOUND"
            path_str = str(info["path"]) if info["path"] else "N/A"
            version_str = versions.get(tool, "unknown")
            message += f"{tool}:\n  Status: {status_str}\n  Path: {path_str}\n  Version: {version_str}\n\n"

        QMessageBox.information(self, "Dependency Check", message)

    def _show_about(self):
        """Show about dialog."""
        from yleaf_gui import __version__ as gui_version
        from yleaf import __version__ as yleaf_version

        QMessageBox.about(
            self, "About Yleaf",
            f"<h2>Yleaf</h2>"
            f"<p>Y-Chromosome Haplogroup Prediction Tool</p>"
            f"<p>GUI Version: {gui_version}</p>"
            f"<p>Core Version: {yleaf_version}</p>"
            f"<p></p>"
            f"<p>Developed at Erasmus Medical Center</p>"
            f"<p>Department of Genetic Identification</p>"
        )

    def _export_csv(self):
        """Export results to CSV file."""
        if self.results_table.rowCount() == 0:
            QMessageBox.warning(self, "No Results", "No results to export.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save CSV", "", "CSV Files (*.csv)"
        )

        if file_path:
            try:
                with open(file_path, 'w') as f:
                    # Header
                    headers = []
                    for col in range(self.results_table.columnCount()):
                        headers.append(self.results_table.horizontalHeaderItem(col).text())
                    f.write(",".join(headers) + "\n")

                    # Data
                    for row in range(self.results_table.rowCount()):
                        row_data = []
                        for col in range(self.results_table.columnCount()):
                            item = self.results_table.item(row, col)
                            row_data.append(item.text() if item else "")
                        f.write(",".join(row_data) + "\n")

                QMessageBox.information(
                    self, "Export Complete",
                    f"Results exported to:\n{file_path}"
                )
            except Exception as e:
                QMessageBox.critical(self, "Export Error", str(e))

    def _open_html_report(self):
        """Open the HTML report in default browser."""
        output_folder = self.output_input.text().strip()
        if not output_folder:
            QMessageBox.warning(self, "No Output", "No output folder specified.")
            return

        html_path = Path(output_folder) / "report.html"
        if html_path.exists():
            import webbrowser
            webbrowser.open(str(html_path))
        else:
            QMessageBox.warning(
                self, "Report Not Found",
                f"HTML report not found at:\n{html_path}"
            )

    def _view_tree(self):
        """View the haplogroup tree image."""
        output_folder = self.output_input.text().strip()
        if not output_folder:
            QMessageBox.warning(self, "No Output", "No output folder specified.")
            return

        # Look for tree image files
        for ext in [".svg", ".png", ".pdf"]:
            tree_path = Path(output_folder) / f"hg_tree_image{ext}"
            if tree_path.exists():
                import webbrowser
                webbrowser.open(str(tree_path))
                return

        QMessageBox.warning(
            self, "Tree Not Found",
            "Haplogroup tree image not found.\n"
            "Make sure 'Draw Haplogroup Tree' is enabled."
        )
