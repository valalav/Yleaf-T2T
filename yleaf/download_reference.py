#!/usr/bin/python
"""
Code for downloading the reference genome and extracting the specific y-chromomse data

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Extensively modified by: Bram van Wersch
"""

import os
from pathlib import Path
import urllib.request
import gzip
import shutil
import urllib.request
import logging

from yleaf import yleaf_constants

LOG: logging = logging.getLogger("yleaf_logger")


def main(
    choice: str
):
    if choice == yleaf_constants.HG19:
        reference_choice = [yleaf_constants.HG19]
    elif choice == yleaf_constants.HG38:
        reference_choice = [yleaf_constants.HG38]
    elif choice == yleaf_constants.T2T:
        reference_choice = [yleaf_constants.T2T]
    else:
        reference_choice = [yleaf_constants.HG19, yleaf_constants.HG38, yleaf_constants.T2T]

    # running downloader
    for dir_name in reference_choice:
        install_genome_files(dir_name)


def install_genome_files(
    reference_choice: str
):
    LOG.info(f"Starting with preparing {reference_choice}...")

    if reference_choice == yleaf_constants.HG19:
        ref_file = yleaf_constants.HG19_FULL_GENOME
        url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{reference_choice}/bigZips/{reference_choice}.fa.gz"
    elif reference_choice == yleaf_constants.HG38:
        ref_file = yleaf_constants.HG38_FULL_GENOME
        url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{reference_choice}/bigZips/{reference_choice}.fa.gz"
    else:
        # T2T (hs1)
        ref_file = yleaf_constants.T2T_FULL_GENOME
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"

    ref_gz_file = Path(str(ref_file) + ".gz")
    try:
        # Increase threshold to 1000 bytes to catch placeholders
        if os.path.getsize(ref_file) < 1000 and not ref_gz_file.exists():

            LOG.debug(f"Downloading the {reference_choice} genome...")
            print(f"Downloading {url} to {ref_gz_file}...")
            urllib.request.urlretrieve(url, ref_gz_file)
            
        if os.path.getsize(ref_file) < 1000:
            LOG.debug("Unpacking the downloaded archive...")
            print("Unpacking...")
            with gzip.open(ref_gz_file, 'rb') as f_in:
                with open(ref_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # Don't remove gz immediately if we want to keep cache? 
            # Original code removed it. Let's keep consistent.
            os.remove(ref_gz_file)

        if reference_choice == yleaf_constants.HG19:
            ychrom_file = yleaf_constants.HG19_Y_CHROMOSOME
        elif reference_choice == yleaf_constants.HG38:
            ychrom_file = yleaf_constants.HG38_Y_CHROMOSOME
        else:
            ychrom_file = yleaf_constants.T2T_Y_CHROMOSOME

        if os.path.getsize(ychrom_file) < 1000:
            LOG.debug("Writing Ychromosomal data")
            print("Extracting chrY...")
            get_ychrom_data(ref_file, ychrom_file, reference_choice)
            
    except KeyboardInterrupt:
        try:
            if ref_gz_file.exists(): os.remove(ref_gz_file)
            if ref_file.exists() and os.path.getsize(ref_file) < 1000: os.remove(ref_file)
        finally:
            raise


def get_ychrom_data(
    full_data_path: Path,
    yhcrom_file: Path,
    ref_type: str = ""
):
    # T2T uses 'chrY', others usually 'chrY' too.
    target_header = ">chrY\n"
    
    with open(yhcrom_file, "w") as fo:
        with open(full_data_path) as fi:
            record = False
            for line in fi:
                # Simple check for chrY header
                if line.startswith(">chrY") or line.startswith(">Y"):
                    record = True
                    fo.write(line)
                elif record:
                    if line.startswith(">"):
                        break
                    fo.write(line)
