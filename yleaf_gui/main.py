#!/usr/bin/env python
"""
Yleaf GUI - Main Entry Point

This script launches the Yleaf graphical user interface.
"""

import sys
import os

# Ensure the yleaf package is importable
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)


def main():
    """Launch the Yleaf GUI application."""
    # Check for PyQt6
    try:
        from PyQt6.QtWidgets import QApplication
        from PyQt6.QtCore import Qt
    except ImportError:
        print("Error: PyQt6 is required for the GUI.")
        print("Install it with: pip install PyQt6")
        sys.exit(1)

    # Enable High DPI scaling
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"

    # Create application
    app = QApplication(sys.argv)
    app.setApplicationName("Yleaf")
    app.setOrganizationName("Erasmus MC")
    app.setOrganizationDomain("erasmusmc.nl")

    # Set application style
    app.setStyle("Fusion")

    # Import and create main window
    from yleaf_gui.main_window import YleafMainWindow

    window = YleafMainWindow()
    window.show()

    # Run event loop
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
