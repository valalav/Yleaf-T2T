#!/usr/bin/env python

"""
Yleaf detection of Y-Haplogroups in Human DNA v3.1

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Extensively modified by: Bram van Wersch
"""

import argparse
import os
import sys
import logging
import shutil
import subprocess
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
from argparse import ArgumentParser
from pathlib import Path
from typing import Union, List, TextIO, Tuple, Dict, Set
from collections import defaultdict
import time
import datetime

from yleaf import __version__
from yleaf.tree import Tree
from yleaf import yleaf_constants, download_reference, html_report, summary_logger
from yleaf import external_tools

pd.options.mode.chained_assignment = None  # default='warn'

CACHED_POSITION_DATA: Union[Set[str], None] = None
CACHED_SNP_DATABASE: Union[Dict[str, List[Dict[str, str]]], None] = None
CACHED_REFERENCE_FILE: Union[List[str], None] = None
NUM_SET: Set[str] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"}
ACCEPTED_REF_BASES: Set[str] = {"A", "C", "G", "T"}

# path constants
PREDICTION_OUT_FILE_NAME: str = "hg_prediction.hg"
HAPLOGROUP_IMAGE_FILE_NAME: str = "hg_tree_image"

LOG: logging.Logger = logging.getLogger("yleaf_logger")


# =============================================================================
# Helper functions for reducing code duplication
# =============================================================================

def _load_marker_file(path_markerfile: Path) -> pd.DataFrame:
    """
    Load marker reference file and remove duplicate positions.

    Args:
        path_markerfile: Path to tab-separated marker file.

    Returns:
        DataFrame with columns: chr, marker_name, haplogroup, pos, mutation, anc, der
    """
    markerfile = pd.read_csv(
        path_markerfile, header=None, sep="\t",
        names=["chr", "marker_name", "haplogroup", "pos", "mutation", "anc", "der"],
        dtype={
            "chr": str, "marker_name": str, "haplogroup": str,
            "pos": int, "mutation": str, "anc": str, "der": str
        }
    )
    return markerfile.drop_duplicates(subset='pos', keep='first', inplace=False)


def _apply_quality_filters(
        df: pd.DataFrame,
        reads_thresh: int,
        base_majority: float,
        reads_col: str = 'reads'
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    """
    Apply standard Yleaf quality filters.
    
    Returns:
        df_passed: Markers passing all filters
        df_filtered: Markers failing filters (concatenated)
        stats: Dictionary of counts
    """
    stats = {}
    
    # 1. Zero reads
    # VCF logic sometimes has 'called_reads', pileup uses 'reads'
    zero_mask = df[reads_col] == 0
    df_zero = df[zero_mask].copy()
    df_zero["Description"] = "Position with zero reads"
    # Reset fields for zero reads
    cols_to_na = ["called_perc", "called_base", "state"]
    for c in cols_to_na:
        if c in df_zero.columns: df_zero[c] = "NA"
    
    stats["zero"] = len(df_zero)
    df = df[~zero_mask]

    # 2. Discordant (if bool_state exists)
    if "bool_state" in df.columns:
        discordant_mask = df["bool_state"] == False
        df_discordant = df[discordant_mask].copy()
        if "bool_state" in df_discordant.columns:
            df_discordant = df_discordant.drop(["bool_state"], axis=1)
        df_discordant["state"] = "NA"
        df_discordant["Description"] = "Discordant genotype"
        stats["discordant"] = len(df_discordant)
        df = df[~discordant_mask]
    else:
        df_discordant = pd.DataFrame()
        stats["discordant"] = 0

    # 3. Read Threshold
    low_reads_mask = df[reads_col] < reads_thresh
    df_low_reads = df[low_reads_mask].copy()
    df_low_reads["Description"] = "Below read threshold"
    stats["low_reads"] = len(df_low_reads)
    df = df[~low_reads_mask]

    # 4. Base Majority
    low_majority_mask = df["called_perc"] < base_majority
    df_low_majority = df[low_majority_mask].copy()
    df_low_majority["Description"] = "Below base majority"
    stats["low_majority"] = len(df_low_majority)
    df = df[~low_majority_mask]

    stats["passed"] = len(df)
    
    # Combine failed markers
    df_filtered = pd.concat([df_zero, df_low_reads, df_low_majority, df_discordant], axis=0, sort=True)
    # Ensure correct columns for fmf
    fmf_cols = ['chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der', reads_col,
                'called_perc', 'called_base', 'state', 'Description']
    # Select available columns
    avail_cols = [c for c in fmf_cols if c in df_filtered.columns]
    df_filtered = df_filtered[avail_cols]
    if reads_col != 'reads':
        df_filtered = df_filtered.rename(columns={reads_col: 'reads'})
        
    stats["total_filtered"] = len(df_filtered)
    
    # Clean up passed df
    if "bool_state" in df.columns:
        df = df.drop(["bool_state"], axis=1)
        
    return df, df_filtered, stats


def _generate_info_statistics(
        markerfile_len: int,
        stats: Dict[str, int],
        reads_thresh: int,
        base_majority: float
) -> List[str]:
    """Generate stats strings for .info file."""
    info = []
    info.append(f"Markers with zero reads: {stats['zero']}")
    info.append(f"Markers below the read threshold {{{reads_thresh}}}: {stats['low_reads']}")
    info.append(f"Markers below the base majority threshold {{{base_majority}}}: {stats['low_majority']}")
    info.append(f"Markers with discordant genotype: {stats['discordant']}")
    info.append(f"Markers without haplogroup information: {stats['total_filtered']}")
    info.append(f"Markers with haplogroup information: {stats['passed']}")
    return info


def _apply_quality_filters(
        df: pd.DataFrame,
        reads_thresh: float,
        base_majority: float,
        reads_column: str = "reads",
        called_reads_column: str = None
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    """
    Apply quality filters to marker dataframe.

    Filters applied in order:
    1. Zero reads - positions with no coverage
    2. Discordant genotypes - called base doesn't match anc/der
    3. Read threshold - coverage below minimum
    4. Base majority - allele frequency below minimum

    Args:
        df: DataFrame with marker data including called_base, called_perc, state columns.
        reads_thresh: Minimum number of reads required.
        base_majority: Minimum percentage for called allele.
        reads_column: Name of column containing read counts.
        called_reads_column: Name of column with called reads (if different from reads_column).

    Returns:
        Tuple of (df_passed, df_filtered, filter_stats) where filter_stats is a dict
        with counts for each filter category.
    """
    if called_reads_column is None:
        called_reads_column = reads_column

    filter_stats = {}

    # Zero reads filter
    index_belowzero = df[df[called_reads_column] == 0].index
    df_belowzero = df[df.index.isin(index_belowzero)].copy()
    df_belowzero["called_perc"] = "NA"
    df_belowzero["called_base"] = "NA"
    df_belowzero["state"] = "NA"
    df_belowzero["Description"] = "Position with zero reads"
    filter_stats["zero_reads"] = len(df_belowzero)

    df = df[~df.index.isin(index_belowzero)]

    # Discordant genotypes filter
    if "bool_state" in df.columns:
        bool_list_state = df[df["bool_state"] == False].index
        df_discordant = df[df.index.isin(bool_list_state)].copy()
        df_discordant = df_discordant.drop(["bool_state"], axis=1, errors='ignore')
        df_discordant["state"] = "NA"
        df_discordant["Description"] = "Discordant genotype"
        df = df[~df.index.isin(bool_list_state)]
    else:
        df_discordant = pd.DataFrame()
    filter_stats["discordant"] = len(df_discordant)

    # Read threshold filter
    df_readsthreshold = df[df[called_reads_column] < reads_thresh].copy()
    df_readsthreshold["Description"] = "Below read threshold"
    df = df[df[called_reads_column] >= reads_thresh]
    filter_stats["below_reads"] = len(df_readsthreshold)

    # Base majority filter
    df_basemajority = df[df["called_perc"] < base_majority].copy()
    df_basemajority["Description"] = "Below base majority"
    df = df[df["called_perc"] >= base_majority]
    filter_stats["below_majority"] = len(df_basemajority)

    # Combine filtered markers
    df_fmf = pd.concat(
        [df_belowzero, df_readsthreshold, df_basemajority, df_discordant],
        axis=0, sort=True
    )

    filter_stats["total_filtered"] = len(df_fmf)
    filter_stats["passed"] = len(df)

    return df, df_fmf, filter_stats


def _generate_info_statistics(
        markerfile_len: int,
        filter_stats: Dict[str, int],
        reads_thresh: float,
        base_majority: float,
        is_vcf: bool = False
) -> List[str]:
    """
    Generate info statistics list from filter results.

    Args:
        markerfile_len: Total number of valid markers.
        filter_stats: Dictionary with filter counts from _apply_quality_filters.
        reads_thresh: Read threshold used.
        base_majority: Base majority threshold used.
        is_vcf: Whether processing VCF (affects initial lines).

    Returns:
        List of info strings for writing to .info file.
    """
    info_list = []
    if is_vcf:
        info_list.append("Total of mapped reads: VCF")
        info_list.append("Total of unmapped reads: VCF")
    info_list.append(f"Valid markers: {markerfile_len}")
    info_list.append(f"Markers with zero reads: {filter_stats['zero_reads']}")
    info_list.append(
        f"Markers below the read threshold {{{reads_thresh}}}: {filter_stats['below_reads']}"
    )
    info_list.append(
        f"Markers below the base majority threshold {{{base_majority}}}: {filter_stats['below_majority']}"
    )
    info_list.append(f"Markers with discordant genotype: {filter_stats['discordant']}")
    info_list.append(f"Markers without haplogroup information: {filter_stats['total_filtered']}")
    info_list.append(f"Markers with haplogroup information: {filter_stats['passed']}")
    return info_list


def _write_tree_sorted_output(
        df_out: pd.DataFrame,
        outputfile: Path,
        columns: List[str]
) -> None:
    """
    Write output file sorted by phylogenetic tree order.

    Args:
        df_out: DataFrame with marker results.
        outputfile: Path to output file.
        columns: Column names to write.
    """
    lst_df = df_out.values.tolist()
    mappable_df: Dict[str, List] = {}
    for lst in lst_df:
        hg = lst[3]  # haplogroup column index
        if hg not in mappable_df:
            mappable_df[hg] = []
        mappable_df[hg].append(lst)

    tree = Tree(
        yleaf_constants.DATA_FOLDER / yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.TREE_FILE
    )
    with open(outputfile, "w") as f:
        f.write('\t'.join(columns + ["depth"]) + "\n")
        for node_key in tree.node_mapping:
            if node_key not in mappable_df:
                continue
            depth = tree.get(node_key).depth
            for lst in mappable_df[node_key]:
                f.write('\t'.join(map(str, lst)) + f"\t{depth}\n")


# =============================================================================
# Main processing classes and functions
# =============================================================================


class MyFormatter(logging.Formatter):
    """
    Copied from MultiGeneBlast (my code)
    """

    def __init__(self, fmt, starttime=time.time()):
        logging.Formatter.__init__(self, fmt)
        self._start_time = starttime

    def format(self, record):
        """
        Overwrite of the format function that prints the passed time and adds
        current time to the existing format
        :See: logging.Formatter.format()
        """
        record.passedTime = "{:.3f}".format(time.time() - self._start_time)
        record.currentTime = datetime.datetime.now().time()
        return super(MyFormatter, self).format(record)


def run_vcf(
        path_markerfile: Path,
        base_out_folder: Path,
        args: argparse.Namespace,
        sample_vcf_file: Path,
):
    """
    Process a single VCF file to extract Y-chromosome haplogroup markers.

    This function reads marker positions from a reference file, filters the VCF
    for relevant Y-chromosome positions, and extracts genotype information for
    haplogroup prediction.

    Args:
        path_markerfile: Path to the marker reference file containing known
            Y-chromosome SNP positions and their haplogroup associations.
        base_out_folder: Base output directory for results.
        args: Namespace containing command-line arguments including:
            - reads_treshold: Minimum read count threshold
            - base_majority: Minimum base percentage for acceptance
            - reference_genome: Reference genome version (hg19/hg38/t2t)
        sample_vcf_file: Path to the input VCF file to process.

    Returns:
        None. Results are written to output files in sample-specific folders.

    Side Effects:
        - Creates output folder for the sample
        - Writes .out file with valid markers
        - Writes .fmf file with filtered markers
        - Writes .info file with processing statistics
    """
    LOG.debug("Starting with extracting haplogroups...")
    markerfile = _load_marker_file(path_markerfile)

    sample_vcf_folder = base_out_folder / (sample_vcf_file.name.replace(".vcf.gz", ""))
    safe_create_dir(sample_vcf_folder, args.force)

    sample_vcf_file_txt = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".txt"))
    external_tools.bcftools_query(
        sample_vcf_file,
        '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n',
        output_file=sample_vcf_file_txt
    )

    pileupfile = pd.read_csv(sample_vcf_file_txt,
                             dtype=str, header=None, sep="\t")

    # remove sample_vcf_file_txt
    external_tools.safe_remove(sample_vcf_file_txt)

    pileupfile.columns = ['chr', 'pos', 'refbase', 'altbase', 'reads']
    pileupfile['pos'] = pileupfile['pos'].astype(int)

    pileupfile['altbase'] = pileupfile['altbase'].str.split(',')
    pileupfile['reads'] = pileupfile['reads'].str.split(',')

    pileupfile['ref_reads'] = pileupfile['reads'].apply(lambda x: x[0])
    pileupfile['alt_reads'] = pileupfile['reads'].apply(lambda x: x[1:])

    pileupfile['alt_reads_dict'] = pileupfile.apply(lambda row: dict(zip(row['altbase'], row['alt_reads'])), axis=1)
    pileupfile['alt_reads_dict'] = pileupfile['alt_reads_dict'].apply(lambda x: {k: int(v) for k, v in x.items()})
    pileupfile['highest_alt_reads'] = pileupfile['alt_reads_dict'].apply(lambda x: max(x.values()) if len(x) > 0 else 0)
    pileupfile['highest_alt_reads_base'] = pileupfile['alt_reads_dict'].apply(
        lambda x: max(x, key=x.get) if len(x) > 0 else 'NA')
    pileupfile['total_reads'] = pileupfile.apply(lambda row: int(row['ref_reads']) + row['highest_alt_reads'], axis=1)
    pileupfile['called_ref_perc'] = pileupfile.apply(
        lambda row: round((int(row['ref_reads']) / row['total_reads']) * 100, 1) if row['total_reads'] > 0 else 0,
        axis=1)
    pileupfile['called_alt_perc'] = pileupfile.apply(
        lambda row: round((row['highest_alt_reads'] / row['total_reads']) * 100, 1) if row['total_reads'] > 0 else 0,
        axis=1)

    pileupfile['called_base'] = pileupfile.apply(
        lambda row: row['refbase'] if row['called_ref_perc'] >= row['called_alt_perc'] else row[
            'highest_alt_reads_base'], axis=1)
    pileupfile['called_perc'] = pileupfile.apply(
        lambda row: row['called_ref_perc'] if row['called_ref_perc'] >= row['called_alt_perc'] else row[
            'called_alt_perc'], axis=1).astype(float)
    pileupfile['called_reads'] = pileupfile.apply(
        lambda row: row['ref_reads'] if row['called_ref_perc'] >= row['called_alt_perc'] else row['highest_alt_reads'],
        axis=1).astype(int)

    intersect_pos = np.intersect1d(pileupfile['pos'], markerfile['pos'])
    markerfile = markerfile.loc[markerfile['pos'].isin(intersect_pos)]
    markerfile = markerfile.sort_values(by=['pos'])
    pileupfile = pileupfile.loc[pileupfile['pos'].isin(intersect_pos)]

    pileupfile = pileupfile.drop(['chr'], axis=1)
    df = pd.merge(markerfile, pileupfile, on='pos')

    df['state'] = df.apply(
        lambda row: 'A' if row['called_base'] == row['anc'] else 'D' if row['called_base'] == row['der'] else 'NA',
        axis=1)
    df['bool_state'] = df.apply(
        lambda row: True if (row['called_base'] == row['anc'] or row['called_base'] == row['der']) else False, axis=1)

    markerfile_len = len(markerfile)

    # valid markers from positionsfile.txt
    general_info_list = ["Total of mapped reads: VCF", "Total of unmapped reads: VCF"]
    general_info_list += ["Valid markers: " + str(markerfile_len)]

    reads_thresh = int(args.reads_treshold)
    base_majority = float(args.base_majority)

    # Use helper for filtering
    df_out, df_fmf, stats = _apply_quality_filters(
        df, 
        reads_thresh=reads_thresh, 
        base_majority=base_majority, 
        reads_col='called_reads'
    )
    
    # Rename 'called_reads' to 'reads' for output consistency if needed
    if 'called_reads' in df_out.columns:
        df_out['reads'] = df_out['called_reads']
    
    # Select final columns
    out_columns = ['chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der', 'reads',
                   'called_perc', 'called_base', 'state']
    df_out = df_out[out_columns]

    # Generate stats
    general_info_list.extend(_generate_info_statistics(markerfile_len, stats, reads_thresh, base_majority))

    write_info_file(sample_vcf_folder, general_info_list)

    use_old = args.use_old
    outputfile = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".out"))
    fmf_output = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".fmf"))

    if use_old:
        df_out = df_out.sort_values(by=['haplogroup'], ascending=True)
        df_fmf.to_csv(fmf_output, sep="\t", index=False)
        df_out.to_csv(outputfile, sep="\t", index=False)
        return

    df_fmf.to_csv(fmf_output, sep="\t", index=False)

    # Write output sorted by phylogenetic tree order
    _write_tree_sorted_output(df_out, outputfile, out_columns)

    LOG.info(f"Finished extracting genotypes for {sample_vcf_file.name}")


def main_vcf_split(
        position_bed_file: Path,
        base_out_folder: Path,
        args: argparse.Namespace,
        vcf_file: Path
):
    """
    Split a multi-sample VCF file and process each sample individually.

    This function handles VCF files containing multiple samples by:
    1. Sorting the VCF file
    2. Filtering by Y-chromosome regions from the BED file
    3. Splitting into individual sample VCF files
    4. Processing each sample through the haplogroup extraction pipeline

    Args:
        position_bed_file: Path to BED file containing Y-chromosome marker positions.
        base_out_folder: Base output directory for results.
        args: Namespace containing command-line arguments.
        vcf_file: Path to the multi-sample VCF file.

    Returns:
        None. Results are written to sample-specific output folders.

    Side Effects:
        - Creates filtered_vcf_files subdirectory
        - Writes sorted and filtered VCF files
        - Processes each sample and writes haplogroup results
    """
    # first sort the vcf file
    sorted_vcf_file = base_out_folder / (vcf_file.name.replace(".vcf.gz", ".sorted.vcf.gz"))
    try:
        external_tools.bcftools_sort(vcf_file, sorted_vcf_file, compressed=True)
    except SystemExit:
        LOG.error(f"Failed to sort the vcf file {vcf_file.name}. Skipping...")
        return None

    # next index the vcf file
    external_tools.bcftools_index(sorted_vcf_file, force=True)

    # get chromosome annotation using bcftools query (Python replacement for | uniq)
    chromosomes = external_tools.bcftools_query_chromosomes(sorted_vcf_file)
    chry = [x for x in chromosomes if "y" in x.lower()]
    if len(chry) == 0:
        LOG.error("Unable to find Y-chromosome in the vcf file. Exiting...")
        raise SystemExit("Unable to find Y-chromosome in the vcf file.")
    elif len(chry) > 1:
        LOG.error("Multiple Y-chromosome annotations found in the vcf file. Exiting...")
        raise SystemExit("Multiple Y-chromosome annotations found in the vcf file.")
    else:
        # make new position_bed_file with correct chrY annotation
        new_position_bed_file = base_out_folder / (vcf_file.name.replace(".vcf.gz", "temp_position_bed.bed"))
        with open(position_bed_file, "r") as f:
            with open(new_position_bed_file, "w") as f2:
                for line in f:
                    line = line.replace("chrY", chry[0])
                    f2.write(line)

    # filter the vcf file using the reference bed file
    filtered_vcf_file = base_out_folder / "filtered_vcf_files" / (
        sorted_vcf_file.name.replace(".sorted.vcf.gz", ".filtered.vcf.gz"))
    external_tools.bcftools_view_regions(sorted_vcf_file, filtered_vcf_file, new_position_bed_file, compressed=True)

    # remove temp_position_bed.bed
    external_tools.safe_remove(new_position_bed_file)

    # remove sorted.vcf.gz and sorted.vcf.gz.csi
    external_tools.safe_remove(sorted_vcf_file)
    external_tools.safe_remove(Path(str(sorted_vcf_file) + ".csi"))

    # check number of samples in the vcf file (Python replacement for | wc -l)
    samples = external_tools.bcftools_query_samples(filtered_vcf_file)
    num_samples = len(samples)

    if num_samples > 1:
        # split the vcf file into separate files for each sample
        split_vcf_folder = base_out_folder / (vcf_file.name.replace(".vcf.gz", "_split"))
        safe_create_dir(split_vcf_folder, args.force)
        external_tools.bcftools_split(filtered_vcf_file, split_vcf_folder, compressed=True)
        sample_vcf_files = get_files_with_extension(split_vcf_folder, '.vcf.gz')
    elif num_samples == 1:
        sample_vcf_files = [filtered_vcf_file]
    else:
        LOG.error("No samples found in the vcf file. Exiting...")
        raise SystemExit("No samples found in the vcf file.")

    return sample_vcf_files


def main_vcf(
        args: argparse.Namespace,
        base_out_folder: Path
):
    """
    Main entry point for VCF file processing workflow.

    Orchestrates the processing of VCF files for Y-chromosome haplogroup
    prediction. Handles both single-sample and multi-sample VCF files,
    routing them to appropriate processing functions.

    Args:
        args: Namespace containing command-line arguments including:
            - vcffile: Path to input VCF file(s)
            - reference_genome: Reference genome version
            - use_old: Whether to use old prediction method
            - ancient_DNA: Whether processing ancient DNA samples
            - threads: Number of parallel threads
        base_out_folder: Base output directory for all results.

    Returns:
        None. Processing results are written to output files.

    Side Effects:
        - Creates output directory structure
        - Writes filtered VCF files
        - Writes haplogroup prediction results
    """
    position_bed_file = get_position_bed_file(args.reference_genome, args.use_old, args.ancient_DNA)
    path_markerfile = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA)

    safe_create_dir(base_out_folder / "filtered_vcf_files", args.force)

    files = get_files_with_extension(args.vcffile, '.vcf.gz')

    if not args.reanalyze:
        with multiprocessing.Pool(processes=args.threads) as p:
            sample_vcf_files = p.map(partial(main_vcf_split, position_bed_file, base_out_folder, args), files)

        sample_vcf_files = [x for x in sample_vcf_files if x is not None]
        sample_vcf_files = [item for sublist in sample_vcf_files for item in sublist]

        with multiprocessing.Pool(processes=args.threads) as p:
            p.map(partial(run_vcf, path_markerfile, base_out_folder, args), sample_vcf_files)

    else:
        with multiprocessing.Pool(processes=args.threads) as p:
            p.map(partial(run_vcf, path_markerfile, base_out_folder, args), files)


def _resolve_reference_genome(args: argparse.Namespace) -> str:
    """Determine the reference genome version to use."""
    if args.reference_genome:
        return args.reference_genome

    LOG.info("No reference genome specified. Attempting to detect from input file header...")
    
    input_file = None
    if args.bamfile: input_file = args.bamfile
    elif args.cramfile: input_file = args.cramfile
    
    if input_file:
        detected = detect_reference_genome(input_file)
        LOG.info(f"Detected reference genome: {detected}")
        return detected
        
    LOG.error("Reference genome must be specified for VCF/FastQ inputs or if detection fails.")
    sys.exit(1)


def _setup_custom_reference(args: argparse.Namespace):
    """Apply manual reference FASTA override if specified."""
    if args.ref_fasta:
        LOG.info(f"Overriding reference FASTA with: {args.ref_fasta}")
        # Update constants dynamically
        if args.reference_genome == yleaf_constants.HG19:
            yleaf_constants.HG19_FULL_GENOME = args.ref_fasta
            yleaf_constants.HG19_Y_CHROMOSOME = args.ref_fasta
        elif args.reference_genome == yleaf_constants.HG38:
            yleaf_constants.HG38_FULL_GENOME = args.ref_fasta
            yleaf_constants.HG38_Y_CHROMOSOME = args.ref_fasta
        elif args.reference_genome == yleaf_constants.T2T:
            yleaf_constants.T2T_FULL_GENOME = args.ref_fasta
            yleaf_constants.T2T_Y_CHROMOSOME = args.ref_fasta


def _process_input_files(args: argparse.Namespace, out_folder: Path):
    """Dispatch processing based on input file type."""
    if args.fastq:
        main_fastq(args, out_folder)
    elif args.bamfile:
        main_bam_cram(args, out_folder, True)
    elif args.cramfile:
        # Auto-fill CRAM reference from config if not provided
        if args.cram_reference is None:
            ref_path = get_reference_path(args.reference_genome, is_full=True)
            if ref_path and ref_path.exists() and ref_path.stat().st_size > 100:
                LOG.info(f"Auto-detected CRAM reference: {ref_path}")
                args.cram_reference = ref_path
            else:
                raise ValueError("Please specify a reference genome for the CRAM file (-cr) or configure it in config.txt.")
        main_bam_cram(args, out_folder, False)
    elif args.vcffile:
        main_vcf(args, out_folder)
    else:
        LOG.error("Please specify either a bam, a cram, a fastq, or a vcf file")
        raise ValueError("Please specify either a bam, a cram, a fastq, or a vcf file")


def _run_post_processing(args: argparse.Namespace, out_folder: Path, start_time: float):
    """Run prediction, visualization, and reporting."""
    hg_out = out_folder / PREDICTION_OUT_FILE_NAME
    predict_haplogroup(out_folder, hg_out, args.use_old, args.prediction_quality, args.threads)
    
    if args.draw_haplogroups:
        draw_haplogroups(hg_out, args.collapsed_draw_mode)

    # Generate HTML Report
    try:
        html_report.generate_html(out_folder)
    except Exception as e:
        LOG.error(f"Failed to generate HTML report: {e}")

    # Update Summary Log
    try:
        # Determine input file path for logging
        log_input = "Unknown"
        if args.bamfile: log_input = args.bamfile
        elif args.cramfile: log_input = args.cramfile
        elif args.fastq: log_input = args.fastq 
        elif args.vcffile: log_input = args.vcffile
        
        duration = time.time() - start_time
        summary_logger.log_run(out_folder, log_input, args.reference_genome, duration=duration)
    except Exception as e:
        LOG.error(f"Failed to update summary log: {e}")


def main():
    start_time = time.time()
    print("Erasmus MC Department of Genetic Identification\nYleaf: software tool for human Y-chromosomal "
          f"phylogenetic analysis and haplogroup inference v{__version__}")
    logo()

    args = get_arguments()
    out_folder = Path(args.output)
    safe_create_dir(out_folder, args.force)
    setup_logger(out_folder)

    LOG.info(f"Running Yleaf with command: {' '.join(sys.argv)}")

    args.reference_genome = _resolve_reference_genome(args)
    _setup_custom_reference(args)
    
    # make sure the reference genome is present before doing something else, if not present it is downloaded
    check_reference(args.reference_genome)

    _process_input_files(args, out_folder)
    _run_post_processing(args, out_folder, start_time)

    LOG.info("Done!")


def logo():
    print(r"""

                   |
                  /|\          
                 /\|/\    
                \\\|///   
                 \\|//  
                  |||   
                  |||    
                  |||    

        """)


def get_arguments() -> argparse.Namespace:
    parser = ArgumentParser()

    parser.add_argument("-fastq", "--fastq", required=False,
                        help="Use raw FastQ files", metavar="PATH", type=check_file)
    parser.add_argument("-bam", "--bamfile", required=False,
                        help="input BAM file", metavar="PATH", type=check_file)
    parser.add_argument("-cram", "--cramfile", required=False,
                        help="input CRAM file", metavar="PATH", type=check_file)
    parser.add_argument("-cr", "--cram_reference", required=False,
                        help="Reference genome for the CRAM file. Required when using CRAM files.",
                        metavar="PATH", type=check_file)
    parser.add_argument("-vcf", "--vcffile", required=False,
                        help="input VCF file (.vcf.gz)", metavar="PATH", type=check_file)
    parser.add_argument("-ra", "--reanalyze", required=False,
                        help="reanalyze (skip filtering and splitting) the vcf file", action="store_true")
    parser.add_argument("-force", "--force", action="store_true", default=True,
                        help="Delete files without asking (Default: True)")
    parser.add_argument("-no-force", "--interactive", dest="force", action="store_false",
                        help="Ask before deleting files")
    parser.add_argument("-rg", "--reference_genome",
                        help="The reference genome build to be used. If no reference is available "
                             "they will be downloaded. If you added references in your config.txt file these"
                             " will be used instead as reference or the location will be used to download the "
                             "reference if those files are missing or empty.",
                        choices=[yleaf_constants.HG19, yleaf_constants.HG38, yleaf_constants.T2T], required=False)
    parser.add_argument("-rf", "--ref-fasta", required=False, metavar="PATH", type=check_file,
                        help="Manually specify path to the Reference Genome FASTA file (overrides config)")
    parser.add_argument("-o", "--output", required=True,
                        help="Folder name containing outputs", metavar="STRING")
    parser.add_argument("-r", "--reads_treshold",
                        help=f"The minimum number of reads for each base. (default={yleaf_constants.DEFAULT_READ_THRESHOLD})",
                        type=int, required=False,
                        default=yleaf_constants.DEFAULT_READ_THRESHOLD, metavar="INT")
    parser.add_argument("-q", "--quality_thresh",
                        help=f"Minimum quality for each read, integer between 10 and 40. [10-40] (default={yleaf_constants.DEFAULT_QUALITY_THRESHOLD})",
                        type=int, required=False, metavar="INT", default=yleaf_constants.DEFAULT_QUALITY_THRESHOLD)
    parser.add_argument("-b", "--base_majority",
                        help=f"The minimum percentage of a base result for acceptance, integer between 50 and 99."
                             f" [50-99] (default={yleaf_constants.DEFAULT_MAJORITY_THRESHOLD})",
                        type=int, required=False, metavar="INT", default=yleaf_constants.DEFAULT_MAJORITY_THRESHOLD)
    parser.add_argument("-t", "--threads", dest="threads",
                        help="The number of processes to use when running Yleaf.",
                        type=int, default=1, metavar="INT")
    parser.add_argument("-pq", "--prediction_quality", type=float, required=False,
                        default=yleaf_constants.DEFAULT_MIN_SCORE, metavar="FLOAT",
                        help=f"The minimum quality of the prediction (QC-scores) for it to be accepted. [0-1] (default={yleaf_constants.DEFAULT_MIN_SCORE})")

    # arguments for prediction
    parser.add_argument("-old", "--use_old", dest="use_old",
                        help="Add this value if you want to use the old prediction method of Yleaf (version 2.3). This"
                             " version only uses the ISOGG tree and slightly different prediction criteria.",
                        action="store_true")

    # arguments for drawing haplo group trees
    parser.add_argument("-dh", "--draw_haplogroups", help="Draw the predicted haplogroups in the haplogroup tree.",
                        action="store_true")
    parser.add_argument("-hc", "--collapsed_draw_mode",
                        help="Add this flag to compress the haplogroup tree image and remove all uninformative "
                             "haplogroups from it.", action="store_true")

    # arguments for ancient DNA samples
    parser.add_argument("-aDNA", "--ancient_DNA",
                        help="Add this flag if the sample is ancient DNA. This will ignore all G > A and C > T mutations.",
                        action="store_true")

    # arguments for private mutations
    parser.add_argument("-p", "--private_mutations",
                        help="Add this flag to search for private mutations. These are variations that are not"
                             " considered in the phylogenetic tree and thus not used for haplogroup prediction, "
                             "however can be informative and differentiate individuals within the same haplogroup "
                             "prediction.",
                        action="store_true")

    parser.add_argument("-maf", "--minor_allele_frequency", help="Maximum rate of minor allele for it to be considered"
                                                                 " as a private mutation. (default=0.01)",
                        default=0.01, type=float, metavar="FLOAT")

    args = parser.parse_args()
    return args


def check_file(
        path: str
) -> Path:
    """Check for the presence of a file and return a Path object"""
    object_path = Path(path)
    if not object_path.exists():
        raise argparse.ArgumentTypeError("Path to provided file/dir does not exist")
    return object_path


def setup_logger(
        out_folder: Path
):
    """Setup logging"""
    LOG.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    start_time = time.time()
    formatter = MyFormatter('%(levelname)s %(currentTime)s (%(passedTime)s s) - %(message)s',
                            starttime=start_time)
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    LOG.addHandler(handler)

    file_handler = logging.FileHandler(filename=out_folder / "run.log")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    LOG.addHandler(file_handler)

    LOG.debug("Logger created")


def safe_create_dir(
        folder: Path,
        force: bool
):
    """Create the given folder. If the folder is already present delete if the user agrees."""
    if folder.is_dir():
        while True and not force:
            LOG.warning("Folder " + str(folder) + " already exists, would you like to remove it?")
            choice = input("y/n: ")
            if str(choice).upper() == "Y":
                break
            elif str(choice).upper() == "N":
                sys.exit(0)
            else:
                print("Please type y/Y or n/N")
        shutil.rmtree(folder)
        os.mkdir(folder)
    else:
        try:
            os.mkdir(folder)
        except OSError:
            print("Failed to create directory. Exiting...")
            raise


def main_fastq(
        args: argparse.Namespace,
        out_folder: Path
):
    """
    Process FASTQ files by aligning with minimap2 and converting to BAM.

    Note: This function requires minimap2 and uses shell pipes.
    FASTQ processing is NOT supported on Windows in the current version.
    For Windows users, please pre-process FASTQ files to BAM format using
    a Linux/WSL environment before running Yleaf.
    """
    # Check Windows compatibility
    if external_tools.IS_WINDOWS:
        LOG.error("FASTQ processing is not supported on Windows.")
        LOG.error("Please convert your FASTQ files to BAM format using WSL or a Linux system.")
        raise SystemExit("FASTQ processing not available on Windows")

    files = get_files_with_extension(args.fastq, '.fastq')
    files += get_files_with_extension(args.fastq, '.fastq.gz')
    reference = get_reference_path(args.reference_genome, True)
    bam_folder = out_folder / yleaf_constants.FASTQ_BAM_FILE_FOLDER
    try:
        os.mkdir(bam_folder)
    except FileExistsError:
        LOG.debug(f"Directory {bam_folder} already exists, using existing directory")
    except OSError as e:
        LOG.warning(f"Could not create directory {bam_folder}: {e}")
    LOG.info("Creating bam files from fastq files...")

    # For paired end fastq files (with _R1 and _R2 in the file name) and gzipped fastq files with .gz extension align the pairs together
    for fastq_file in files:
        if "_R1" in str(fastq_file) and Path(str(fastq_file).replace("_R1", "_R2")).exists():
            fastq_file2 = Path(str(fastq_file).replace("_R1", "_R2"))
            LOG.info(f"Starting with running for {fastq_file} and {fastq_file2}")
            sam_file = bam_folder / "temp_fastq_sam.sam"
            # Use external_tools for cross-platform compatibility (no shell=True)
            external_tools.minimap2_align(
                reference=Path(reference),
                fastq_files=[fastq_file, fastq_file2],
                output_sam=sam_file,
                threads=args.threads
            )
            bam_file = bam_folder / (fastq_file.name.rsplit("_R1", 1)[0] + ".bam")
            external_tools.samtools_view_sort(
                input_file=sam_file,
                output_bam=bam_file,
                threads=args.threads,
                input_is_sam=True
            )
            external_tools.samtools_index(bam_file, threads=args.threads)
            os.remove(sam_file)
        elif "_R1" not in str(fastq_file) and "_R2" not in str(fastq_file) and ".gz" in str(fastq_file):
            LOG.info(f"Starting with running for {fastq_file}")
            sam_file = bam_folder / "temp_fastq_sam.sam"
            # Use external_tools for cross-platform compatibility (no shell=True)
            external_tools.minimap2_align(
                reference=Path(reference),
                fastq_files=[fastq_file],
                output_sam=sam_file,
                threads=args.threads
            )
            bam_file = bam_folder / (fastq_file.name.rsplit(".", 1)[0] + ".bam")
            external_tools.samtools_view_sort(
                input_file=sam_file,
                output_bam=bam_file,
                threads=args.threads,
                input_is_sam=True
            )
            external_tools.samtools_index(bam_file, threads=args.threads)
            os.remove(sam_file)
    args.bamfile = bam_folder
    main_bam_cram(args, out_folder, True)


def main_bam_cram(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool
):
    if args.bamfile is not None:
        files = get_files_with_extension(args.bamfile, '.bam')
    elif args.cramfile is not None:
        if args.cram_reference is None:
            raise ValueError("Please specify a reference genome for the CRAM file.")
        files = get_files_with_extension(args.cramfile, '.cram')
    else:
        print("Please specify either (a) bam or a cram file(s)")
        return

    with multiprocessing.Pool(processes=args.threads) as p:
        p.map(partial(run_bam_cram, args, base_out_folder, is_bam), files)


def run_bam_cram(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool,
        input_file: Path
):
    LOG.info(f"Starting with running for {input_file}")
    output_dir = base_out_folder / input_file.name.rsplit(".", 1)[0]
    safe_create_dir(output_dir, args.force)
    general_info_list = samtools(output_dir, input_file, is_bam, args)
    write_info_file(output_dir, general_info_list)
    if args.private_mutations:
        find_private_mutations(output_dir, input_file, args, is_bam)
    LOG.debug(f"Finished running for {input_file.name}")
    print()


def call_command(
        command_str: str,
        stdout_location: TextIO = None
):
    """
    Call a command on the command line and make sure to exit if it fails.

    .. deprecated::
        This function uses shell=True which is a security risk and not cross-platform.
        Use external_tools module functions instead:
        - external_tools.minimap2_align() for minimap2
        - external_tools.samtools_view_sort() for samtools view | sort
        - external_tools.samtools_index() for samtools index
        - external_tools.run_command() for other commands (without shell)

    Args:
        command_str: Shell command string to execute
        stdout_location: Optional file handle for stdout redirection

    Raises:
        SystemExit: If command fails with non-zero return code

    Warning:
        On Windows, shell=True requires commands to be Windows-compatible.
    """
    LOG.debug(f"Started running the following command: {command_str}")

    # Check if running on Windows
    if external_tools.IS_WINDOWS:
        LOG.warning("call_command() with shell=True may not work on Windows. "
                   "Consider using external_tools module instead.")

    if stdout_location is None:
        process = subprocess.Popen(command_str, stderr=subprocess.PIPE, shell=True)
    else:
        process = subprocess.Popen(command_str, stderr=subprocess.PIPE, stdout=stdout_location, shell=True)
    # blocking call
    stdout, stderr = process.communicate()
    # will only fail if returncode is not 0
    if process.returncode != 0:
        LOG.error(f"Call: '{command_str}' failed. Reason given: '{stderr.decode('utf-8')}'")
        raise SystemExit("Failed command execution")
    LOG.debug(f"Finished running the command {command_str}")


def get_files_with_extension(
        path: Union[str, Path],
        ext: str
) -> List[Path]:
    """Get all files with a certain extension from a path. The path can be a file or a dir."""
    filtered_files = []
    path = Path(path)  # to be sure
    if path.is_dir():
        for file in path.iterdir():
            if str(file)[-len(ext):] == ext:
                filtered_files.append(file)
        return filtered_files
    else:
        return [path]


def detect_reference_genome(input_path: Path) -> str:
    """
    Detects reference genome (hg19, hg38, t2t) by reading BAM/CRAM header.
    Checks chromosome lengths.
    """
    try:
        # Cross-platform replacement for: samtools view -H | grep @SQ
        header_lines = external_tools.samtools_view_header(input_path)
    except SystemExit:
        LOG.warning("Failed to read header for auto-detection")
        raise ValueError("Could not read input file header.")
    
    # Known lengths of chr1
    LEN_HG19 = 249250621
    LEN_HG38 = 248956422
    LEN_T2T = 248387328
    
    for line in header_lines:
        if "SN:chr1" in line or "SN:1" in line: # Check chr1
            parts = line.split()
            for part in parts:
                if part.startswith("LN:"):
                    length = int(part.split(":")[1])
                    if length == LEN_T2T:
                        return yleaf_constants.T2T
                    elif length == LEN_HG38:
                        return yleaf_constants.HG38
                    elif length == LEN_HG19:
                        return yleaf_constants.HG19
    
    # Fallback: check chrY if chr1 is missing (e.g. targeted bam)
    # T2T chrY: ? (Wait, user file only has chrY?)
    # If user file is ONLY chrY, we might not see chr1 in header IF it was subsetted without keeping full header.
    # But usually 'samtools view -H' keeps all @SQ lines even if reads are subsetted. 
    # Let's check user's previous output for 'samtools view -H'. 
    # It showed: @SQ SN:chr1 LN:248387328. So chr1 IS present in header!
    
    raise ValueError("Could not identify reference genome from header. Please specify -rg manually.")


def check_reference(
        requested_version: str,
):
    reference_file = get_reference_path(requested_version, True)
    if os.path.getsize(reference_file) < 100:
        LOG.info(f"No reference genome version was found. Downloading the {requested_version} reference genome. This "
                 f"should be a one time thing.")
        download_reference.main(requested_version)
        LOG.info("Finished downloading the reference genome.")


def get_reference_path(
        requested_version: str,
        is_full: bool
) -> Union[Path, None]:
    if is_full:
        if requested_version == yleaf_constants.HG19:
            reference_file = yleaf_constants.HG19_FULL_GENOME
        elif requested_version == yleaf_constants.HG38:
            reference_file = yleaf_constants.HG38_FULL_GENOME
        else:
            reference_file = yleaf_constants.T2T_FULL_GENOME
    else:
        if requested_version == yleaf_constants.HG19:
            reference_file = yleaf_constants.HG19_Y_CHROMOSOME
        elif requested_version == yleaf_constants.HG38:
            reference_file = yleaf_constants.HG38_Y_CHROMOSOME
        else:
            reference_file = yleaf_constants.T2T_Y_CHROMOSOME
    return reference_file


def samtools(
        output_folder: Path,
        path_file: Path,
        is_bam_pathfile: bool,
        args: argparse.Namespace,
) -> List[str]:
    outputfile = output_folder / (output_folder.name + ".out")
    fmf_output = output_folder / (output_folder.name + ".fmf")
    pileupfile = output_folder / "temp_haplogroup_pileup.pu"
    reference = args.cram_reference if not is_bam_pathfile else None

    if is_bam_pathfile:
        if not any([Path(str(path_file) + ".bai").exists(), Path(str(path_file).rstrip(".bam") + '.bai').exists()]):
            external_tools.samtools_index(path_file, threads=args.threads)
    else:
        if not any([Path(str(path_file) + ".crai").exists(), Path(str(path_file).rstrip(".cram") + '.crai').exists()]):
            external_tools.samtools_index(path_file, threads=args.threads)
    header, mapped, unmapped = chromosome_table(path_file, output_folder, output_folder.name)

    position_file = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA)

    bed = output_folder / "temp_position_bed.bed"
    write_bed_file(bed, position_file, header)

    execute_mpileup(bed, path_file, pileupfile, args.quality_thresh, reference, region=header)

    general_info_list = ["Total of mapped reads: " + str(mapped), "Total of unmapped reads: " + str(unmapped)]

    extract_haplogroups(position_file, args.reads_treshold, args.base_majority,
                        pileupfile, fmf_output, outputfile, is_bam_pathfile, args.use_old, general_info_list)

    os.remove(pileupfile)
    os.remove(bed)

    LOG.debug("Finished extracting haplogroups")
    return general_info_list


def chromosome_table(
        path_file: Path,
        path_folder: Path,
        file_name: str
) -> Tuple[str, int, int]:
    output = path_folder / (file_name + '.chr')
    tmp_output = path_folder / "tmp.txt"
    external_tools.samtools_idxstats(path_file, tmp_output)
    df_chromosome = pd.read_csv(tmp_output, sep="\t", header=None)

    total_reads = sum(df_chromosome[2])

    unmapped = df_chromosome[df_chromosome[0].str.contains("Y")][3].values[0]

    df_chromosome["perc"] = (df_chromosome[2] / total_reads) * 100
    df_chromosome = df_chromosome.round(decimals=2)
    df_chromosome['perc'] = df_chromosome['perc'].astype(str) + '%'
    df_chromosome = df_chromosome.drop(columns=[1, 3])
    df_chromosome.columns = ['chr', 'reads', 'perc']
    df_chromosome.to_csv(output, index=None, sep="\t")

    os.remove(tmp_output)

    if 'Y' in df_chromosome["chr"].values:
        return "Y", total_reads, unmapped
    elif 'chrY' in df_chromosome["chr"].values:
        return "chrY", total_reads, unmapped
    else:
        LOG.error("Unable to find Y-chromosomal data in the provided files. Exiting...")
        raise SystemExit("No Y-chromosomal data")


def get_position_file(
        reference_name: str,
        use_old: bool,
        ancient_DNA: bool,
) -> Path:
    if use_old:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_ANCIENT_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_FILE
    else:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_ANCIENT_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_FILE
    return position_file


def get_position_bed_file(
        reference_name: str,
        use_old: bool,
        ancient_DNA: bool,
) -> Path:
    if use_old:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_ANCIENT_BED_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_BED_FILE
    else:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_ANCIENT_BED_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_BED_FILE
    return position_file


def write_bed_file(
        bed: Path,
        markerfile: Path,
        header: str
):
    mf = pd.read_csv(markerfile, sep="\t", header=None)
    mf = mf[[0, 3]]
    mf[0] = header
    mf.to_csv(str(bed), sep="\t", index=False, header=False)


def execute_mpileup(
        bed: Union[Path, None],
        bam_file: Path,
        pileupfile: Path,
        quality_thresh: float,
        reference: Union[Path, None],
        region: str = None
):
    # Cross-platform mpileup using external_tools wrapper
    external_tools.samtools_mpileup(
        bam_file=bam_file,
        output_file=pileupfile,
        quality_thresh=int(quality_thresh),
        bed_file=bed,
        reference=reference,
        region=region
    )


def extract_haplogroups(
        path_markerfile: Path,
        reads_thresh: float,
        base_majority: int,
        path_pileupfile: Path,
        fmf_output: Path,
        outputfile: Path,
        is_bam_file: bool,
        use_old: bool,
        general_info_list: List[str]
):
    """
    Extract Y-chromosome haplogroup markers from pileup data.

    This is the core function for processing sequencing data and extracting
    Y-chromosome markers. It applies quality filters to determine which markers
    pass thresholds and which are filtered out.

    Args:
        path_markerfile: Path to reference file containing marker positions
            and their haplogroup associations.
        reads_thresh: Minimum number of reads required for a position to pass.
        base_majority: Minimum percentage of reads that must support the
            called allele (0-100).
        path_pileupfile: Path to samtools mpileup output file.
        fmf_output: Path to write filtered/failed markers (.fmf file).
        outputfile: Path to write valid markers (.out file).
        is_bam_file: True if processing BAM/CRAM, False if processing VCF.
        use_old: Whether to use old ISOGG-based marker positions.
        general_info_list: List to append processing statistics for logging.

    Returns:
        None. Results are written to output files.

    Output Files:
        - outputfile (.out): Tab-separated file with columns:
            haplogroup, pos, marker_name, mutation, anc, der, state
        - fmf_output (.fmf): Same format with additional 'Description' column
            explaining why marker was filtered.

    Filter Reasons:
        - "Position with zero reads": No coverage at position
        - "Reads below threshold": Coverage below reads_thresh
        - "Base majority below threshold": Allele frequency below base_majority
        - "Discordant genotype": Called allele doesn't match ancestral/derived
    """
    LOG.debug("Starting with extracting haplogroups...")
    markerfile = _load_marker_file(path_markerfile)

    # Load pileup file
    pileupfile = pd.read_csv(path_pileupfile, header=None, sep="\t",
                             dtype={0: str, 1: int, 2: str, 3: int, 4: str, 5: str},
                             on_bad_lines='skip')

    pileupfile.columns = ['chr', 'pos', 'refbase', 'reads', 'align', 'quality']

    if not is_bam_file:
        ref_base = pileupfile["refbase"].values
        read_results = pileupfile["align"].values
        new_read_results = list(map(replace_with_bases, ref_base, read_results))
        pileupfile["align"] = new_read_results

    intersect_pos = np.intersect1d(pileupfile['pos'], markerfile['pos'])
    markerfile = markerfile.loc[markerfile['pos'].isin(intersect_pos)]
    markerfile = markerfile.sort_values(by=['pos'])
    pileupfile = pileupfile.loc[pileupfile['pos'].isin(intersect_pos)]

    pileupfile = pileupfile.drop(['chr'], axis=1)
    df = pd.merge(markerfile, pileupfile, on='pos')

    markerfile_len = len(markerfile)

    # valid markers from positionsfile.txt
    general_info_list.append("Valid markers: " + str(markerfile_len))

    # Calculate statistics for non-zero reads
    non_zero_mask = df["reads"] > 0
    df_nonzero = df[non_zero_mask].copy()
    
    # Initialize columns with defaults
    df["called_perc"] = 0.0
    df["called_base"] = "NA"
    df["state"] = "NA"
    df["bool_state"] = False

    if not df_nonzero.empty:
        freq_dict = get_frequency_table(df_nonzero.values)
        df_freq_table = pd.DataFrame.from_dict(freq_dict, orient='index')
        df_freq_table.columns = ["A", "T", "G", "C", "+", "-"]
        df_freq_table = df_freq_table.drop(['+', '-'], axis=1)
        
        list_col_indices = np.argmax(df_freq_table.values, axis=1)
        called_base = df_freq_table.columns[list_col_indices]
        total_count_bases = np.sum(df_freq_table.values, axis=1)
        max_count_bases = np.max(df_freq_table, axis=1)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            called_perc = (max_count_bases / total_count_bases) * 100
            called_perc = called_perc.fillna(0).replace([np.inf, -np.inf], 0).round(1)
            
        bool_anc = np.equal(np.array(called_base), df_nonzero["anc"].values)
        bool_der = np.equal(np.array(called_base), df_nonzero["der"].values)
        
        bool_list_anc = np.where(bool_anc, 'A', 'D')
        bool_list_der = np.where(bool_der, 'D', 'A')
        bool_list_state = np.equal(bool_list_anc, bool_list_der)
        
        df_nonzero["called_perc"] = np.array(called_perc, dtype=float)
        df_nonzero["called_base"] = called_base
        df_nonzero["state"] = bool_list_anc
        df_nonzero["bool_state"] = bool_list_state
        
        df.update(df_nonzero)

    # Apply filters using helper
    reads_thresh = int(reads_thresh)
    base_majority = float(base_majority)
    
    df_out, df_fmf, stats = _apply_quality_filters(
        df,
        reads_thresh=reads_thresh,
        base_majority=base_majority,
        reads_col='reads'
    )
    
    out_columns = ['chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der', 'reads',
                   'called_perc', 'called_base', 'state']
    df_out = df_out[out_columns]

    # Generate stats
    general_info_list.extend(_generate_info_statistics(markerfile_len, stats, reads_thresh, base_majority))

    if use_old:
        df_out = df_out.sort_values(by=['haplogroup'], ascending=True)
        df_out = df_out[
            ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc", "called_base",
             "state"]]
        df_fmf.to_csv(fmf_output, sep="\t", index=False)
        df_out.to_csv(outputfile, sep="\t", index=False)
        return

    out_columns = ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads",
                   "called_perc", "called_base", "state"]
    df_out = df_out[out_columns]
    df_fmf.to_csv(fmf_output, sep="\t", index=False)

    # Write output sorted by phylogenetic tree order
    _write_tree_sorted_output(df_out, outputfile, out_columns)


def replace_with_bases(
        base: str,
        read_result: str
) -> str:
    return read_result.replace(",", base[0]).replace(".", base[0])


def get_frequency_table(
        mpileup: List[str]
) -> Dict[str, List[int]]:
    frequency_table = {}
    for i in mpileup:
        fastadict = get_frequencies(i[9])
        frequency_table[i[3]] = list(fastadict.values())
    return frequency_table


def get_frequencies(
        sequence: str
) -> Dict[str, int]:
    """
    Parse samtools pileup sequence string and count base frequencies.

    This function parses the complex pileup format from samtools mpileup
    and extracts counts for each nucleotide and indel type.

    Args:
        sequence: Pileup sequence string from samtools mpileup output.
            Uses standard pileup encoding:
            - A,T,G,C: nucleotides
            - ^: start of read (followed by mapping quality char)
            - $: end of read
            - +N: insertion of N bases (followed by inserted sequence)
            - -N: deletion of N bases (followed by deleted sequence)
            - *: placeholder for deleted base

    Returns:
        Dict mapping nucleotide/event to count:
        {"A": int, "T": int, "G": int, "C": int, "-": int, "+": int}
        Note: "-" includes both deletions and * placeholders.

    Example:
        >>> get_frequencies("AAATTT^~G+2AC-1G")
        {"A": 3, "T": 3, "G": 1, "C": 0, "-": 1, "+": 1}
    """
    fastadict = {"A": 0, "T": 0, "G": 0, "C": 0, "-": 0, "+": 0, "*": 0}
    sequence = sequence.upper()
    index = 0
    while index < len(sequence):
        char = sequence[index]
        if char in fastadict:
            fastadict[char] += 1
            index += 1
        elif char == "^":
            index += 2
        elif char in {"-", "+"}:
            index += 1
            digit, index = find_digit(sequence, index)
            index += digit
        else:
            index += 1
    fastadict["-"] += fastadict["*"]
    del fastadict["*"]
    return fastadict


def find_digit(
        sequence: str,
        index: int
) -> Tuple[int, int]:
    # first is always a digit
    nr = [sequence[index]]
    index += 1
    while True:
        char = sequence[index]
        # this seems to be faster than isdigit()
        if char in NUM_SET:
            nr.append(char)
            index += 1
            continue
        return int(''.join(nr)), index


def write_info_file(
        folder: Path,
        general_info_list: List[str]
):
    try:
        with open(folder / (folder.name + ".info"), "a") as f:
            for marker in general_info_list:
                f.write(marker)
                f.write("\n")
    except IOError:
        LOG.warning("Failed to write .info file")


def find_private_mutations(
        output_folder: Path,
        path_file: Path,
        args: argparse.Namespace,
        is_bam: bool
):
    # identify mutations not part of haplogroups that are annotated in dbsnp or differ from the reference genome
    LOG.debug("Starting with extracting private mutations...")
    snp_reference_file = get_reference_path(args.reference_genome, False)
    snp_database_file = yleaf_constants.DATA_FOLDER / args.reference_genome / yleaf_constants.SNP_DATA_FILE

    # run mpileup
    pileup_file = output_folder / "temp_private_mutation_pileup.pu"
    execute_mpileup(None, path_file, pileup_file, args.quality_thresh, snp_reference_file if not is_bam else None)

    LOG.debug("Loading reference files")

    position_file = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA)
    filter_positions = load_filter_data(position_file)
    known_snps = load_snp_database_file(snp_database_file, args.minor_allele_frequency)
    ychrom_reference = load_reference_file(snp_reference_file)

    LOG.debug("Finding private mutations...")
    private_mutations = []
    confirmed_private_mutations = []
    with open(pileup_file) as f:
        for index, line in enumerate(f):
            try:
                chrom, position, ref_base, count, aligned, quality = line.strip().split()
            except ValueError:
                LOG.warning(f"failed to read line {index} of pileupfile")
                continue
            if chrom != "chrY":
                continue
            # not enough reads
            if int(count) < args.reads_treshold:
                continue
            if not is_bam:
                aligned = replace_with_bases(ref_base, aligned)
            if position in filter_positions:
                continue

            # not annotated in dbsnp
            if position not in known_snps and snp_reference_file is not None:
                freq_dict = get_frequencies(aligned)
                actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
                called_percentage = round(allele_count / sum(freq_dict.values()) * 100, 2)
                # not a high enough base majority measured, cannot be sure of real allele
                if called_percentage < args.base_majority:
                    continue
                ref_base = ychrom_reference[int(position) - 1]

                # do not match against insertions or repeat regions (lower case)
                if ref_base not in ACCEPTED_REF_BASES or actual_allele not in ACCEPTED_REF_BASES:
                    continue
                if ref_base == actual_allele:
                    continue
                private_mutations.append(f"{chrom}\t{position}\t-\t{ref_base}->{actual_allele}\t{ref_base}\t"
                                         f"{actual_allele}\t{allele_count}\t{called_percentage}\tNA\n")
            elif position in known_snps and snp_database_file is not None:
                freq_dict = get_frequencies(aligned)
                actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
                called_percentage = round(allele_count / sum(freq_dict.values()) * 100, 2)
                # not a high enough base majority measured, cannot be sure of real allele
                if called_percentage < args.base_majority:
                    continue
                possible_minor_alleles = {dct["minor_allele"] for dct in known_snps[position]}
                # measured allele is not the dbsnp allele
                if actual_allele not in possible_minor_alleles:
                    continue

                matched_pos_dct = [dct for dct in known_snps[position]
                                   if dct["minor_allele"] == actual_allele][0]
                rs, major_allele, minor_allele, frequency = matched_pos_dct.values()
                confirmed_private_mutations.append(f"{chrom}\t{position}\t{rs}\t{major_allele}->{actual_allele}\t"
                                                   f"{major_allele}\t{actual_allele}\t{allele_count}\t"
                                                   f"{called_percentage}\t{frequency}\n")
    os.remove(pileup_file)
    with open(output_folder / f"{output_folder.name}.pmu", "w") as f:
        f.write(f"chrom\tposition\trn_no\tmutation\treference\tdetected\treads\tcalled_percentage\t"
                f"minor allele frequency\n")
        f.write(''.join(confirmed_private_mutations))
        f.write(''.join(private_mutations))
    LOG.debug("Finished extracting private mutations")


def load_filter_data(
        path: Path
) -> Set[str]:
    global CACHED_POSITION_DATA
    if CACHED_POSITION_DATA is None:
        CACHED_POSITION_DATA = set()
        with open(path) as f:
            for line in f:
                CACHED_POSITION_DATA.add(line.strip().split("\t")[3])
    return CACHED_POSITION_DATA


def load_snp_database_file(
        path: Path,
        minor_allele_frequency: float
) -> Union[Dict[str, List[Dict[str, str]]], None]:
    global CACHED_SNP_DATABASE

    if path is not None:
        if CACHED_SNP_DATABASE is None:
            CACHED_SNP_DATABASE = defaultdict(list)
            with open(path) as f:
                f.readline()
                for line in f:
                    rs, position, major_allele, minor_allele, frequency = line.strip().split(",")
                    if float(frequency) > minor_allele_frequency:
                        continue
                    CACHED_SNP_DATABASE[position].append({"rs": rs, "major_allele": major_allele,
                                                          "minor_allele": minor_allele, "frequency": frequency})
            CACHED_SNP_DATABASE = dict(CACHED_SNP_DATABASE)
        return CACHED_SNP_DATABASE
    else:
        return None


def load_reference_file(
        path: Path
) -> Union[List[str], None]:
    global CACHED_REFERENCE_FILE
    if path is not None:
        if CACHED_REFERENCE_FILE is None:
            CACHED_REFERENCE_FILE = []
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        continue
                    CACHED_REFERENCE_FILE.extend(line.strip())
        return CACHED_REFERENCE_FILE
    else:
        return None


def predict_haplogroup(
        path_file: Path,
        output: Path,
        use_old: bool,
        prediction_quality: float,
        threads: int,
):
    if use_old:
        script = yleaf_constants.SRC_FOLDER / "old_predict_haplogroup.py"
        # Cross-platform Python script execution
        result = external_tools.run_command(
            [sys.executable, str(script), "-i", str(path_file), "-o", str(output)],
            check=True
        )
    else:
        from yleaf import predict_haplogroup
        namespace = argparse.Namespace(input=path_file, outfile=output,
                                       minimum_score=prediction_quality, threads=threads)
        predict_haplogroup.main(namespace)


def draw_haplogroups(
        haplogroup_file: Path,
        collapsed_draw_mode: bool
):
    # make sure that it is only imported if requested by user
    from yleaf import haplogroup_tree_image
    namespace = argparse.Namespace(input=haplogroup_file, collapse_mode=collapsed_draw_mode,
                                   outfile=haplogroup_file.parent / HAPLOGROUP_IMAGE_FILE_NAME)
    haplogroup_tree_image.main(namespace)


if __name__ == "__main__":
    main()
