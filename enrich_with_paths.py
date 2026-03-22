#!/usr/bin/env python3
"""
enrich_with_paths.py — Post-processing enrichment for Yleaf haplogroup predictions.

Queries the haplo server (ftdna_haplo, port 9003) to add dual-tree paths
(YFull + FTDNA) to Yleaf prediction results.

Usage:
    python3 enrich_with_paths.py -i hg_prediction.hg -o hg_prediction_enriched.hg
    python3 enrich_with_paths.py -i hg_prediction.hg  # overwrites in-place with .enriched backup

Dependencies:
    - haplo server running on localhost:9003 (pm2 start ecosystem.config.js)
    - requests library (pip install requests)
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

try:
    import requests
except ImportError:
    print("ERROR: 'requests' library required. Install: pip install requests", file=sys.stderr)
    sys.exit(1)

LOG = logging.getLogger("enrich_paths")
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

DEFAULT_HAPLO_SERVER = "http://localhost:9003"
API_TIMEOUT = 5  # seconds per request


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Enrich Yleaf haplogroup predictions with dual-tree paths (YFull + FTDNA)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input .hg prediction file from Yleaf")
    parser.add_argument("-o", "--output", help="Output enriched file (default: <input>.enriched.hg)")
    parser.add_argument("--server", default=DEFAULT_HAPLO_SERVER, help=f"Haplo server URL (default: {DEFAULT_HAPLO_SERVER})")
    parser.add_argument("--short-path", action="store_true", default=True,
                        help="Show only last N nodes of path (default: True)")
    parser.add_argument("--path-depth", type=int, default=5,
                        help="Number of tail nodes to show in short path (default: 5)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser.parse_args()


def check_server(server_url: str) -> bool:
    """Check if haplo server is available."""
    try:
        resp = requests.get(f"{server_url}/api/health", timeout=3)
        data = resp.json()
        if data.get("status") in ("ok", "degraded"):
            LOG.info(f"Haplo server OK at {server_url}")
            return True
        LOG.warning(f"Haplo server degraded: {data}")
        return True  # degraded is still usable
    except Exception as e:
        LOG.error(f"Haplo server unreachable at {server_url}: {e}")
        return False


def shorten_path(path_string: str, depth: int = 5) -> str:
    """Shorten a full path string to last N nodes."""
    if not path_string:
        return ""
    nodes = path_string.split(" > ")
    if len(nodes) <= depth:
        return path_string
    return "... > " + " > ".join(nodes[-depth:])


def query_haplogroup(server_url: str, haplogroup: str) -> Dict:
    """Query haplo server for dual-tree path information."""
    # Clean haplogroup name: remove asterisks and ancestral annotations
    clean_hg = haplogroup.split("*")[0].strip()
    if not clean_hg or clean_hg == "NA":
        return {}

    try:
        resp = requests.get(f"{server_url}/api/search/{clean_hg}", timeout=API_TIMEOUT)
        if resp.status_code == 404:
            LOG.debug(f"  {clean_hg}: not found in haplo server")
            return {}
        resp.raise_for_status()
        return resp.json()
    except requests.Timeout:
        LOG.warning(f"  {clean_hg}: timeout ({API_TIMEOUT}s)")
        return {}
    except Exception as e:
        LOG.warning(f"  {clean_hg}: error querying server: {e}")
        return {}


def extract_paths(data: Dict, path_depth: int = 5) -> Tuple[str, str, str, str]:
    """Extract FTDNA and YFull names and paths from haplo server response.
    
    Returns: (hg_ftdna, hg_yfull, path_ftdna, path_yfull)
    """
    hg_ftdna = ""
    hg_yfull = ""
    path_ftdna = ""
    path_yfull = ""

    ftdna = data.get("ftdnaDetails")
    if ftdna and ftdna.get("path"):
        path_nodes = ftdna["path"].get("nodes", [])
        if path_nodes:
            hg_ftdna = path_nodes[-1].get("name", "")
        path_ftdna = shorten_path(ftdna["path"].get("string", ""), path_depth)

    yfull = data.get("yfullDetails")
    if yfull and yfull.get("path"):
        # Use canonicalId if available (more accurate)
        hg_yfull = yfull.get("canonicalId", "")
        if not hg_yfull:
            path_nodes = yfull["path"].get("nodes", [])
            if path_nodes:
                hg_yfull = path_nodes[-1].get("name", "")
        path_yfull = shorten_path(yfull["path"].get("string", ""), path_depth)

    return hg_ftdna, hg_yfull, path_ftdna, path_yfull


def enrich_file(input_path: str, output_path: str, server_url: str, path_depth: int = 5) -> None:
    """Read Yleaf prediction file, enrich with haplo server data, write output."""
    in_file = Path(input_path)
    if not in_file.exists():
        LOG.error(f"Input file not found: {input_path}")
        sys.exit(1)

    # Read input
    with open(in_file) as f:
        lines = f.readlines()

    if not lines:
        LOG.error("Empty input file")
        sys.exit(1)

    # Parse header
    header = lines[0].rstrip("\n")
    header_fields = header.split("\t")

    # New columns to add
    new_columns = ["Hg_FTDNA", "Hg_YFull", "Path_FTDNA", "Path_YFull"]
    enriched_header = header + "\t" + "\t".join(new_columns)

    # Process data lines
    enriched_lines = [enriched_header + "\n"]
    total = 0
    enriched_count = 0
    skipped = 0

    for line in lines[1:]:
        line = line.rstrip("\n")
        if not line.strip():
            continue

        fields = line.split("\t")
        total += 1

        # Hg is the second column (index 1)
        hg = fields[1] if len(fields) > 1 else ""

        # Skip NA and Low_Y_Signal
        if not hg or hg == "NA" or hg.startswith("Low_Y_Signal"):
            enriched_lines.append(line + "\t\t\t\t\n")
            skipped += 1
            continue

        # Query haplo server
        data = query_haplogroup(server_url, hg)
        if data:
            hg_ftdna, hg_yfull, path_ftdna, path_yfull = extract_paths(data, path_depth)
            enriched_lines.append(f"{line}\t{hg_ftdna}\t{hg_yfull}\t{path_ftdna}\t{path_yfull}\n")
            enriched_count += 1
            LOG.info(f"  {fields[0]}: {hg} → FTDNA:{hg_ftdna or 'N/A'} | YFull:{hg_yfull or 'N/A'}")
        else:
            enriched_lines.append(line + "\t\t\t\t\n")
            LOG.info(f"  {fields[0]}: {hg} → no cross-tree data")

    # Write output
    out_file = Path(output_path)
    with open(out_file, "w") as f:
        f.writelines(enriched_lines)

    LOG.info(f"\n{'='*60}")
    LOG.info(f"Results: {enriched_count}/{total} enriched, {skipped} skipped (NA/Low_Y_Signal)")
    LOG.info(f"Output: {out_file}")


def main():
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        in_path = Path(args.input)
        output_path = str(in_path.parent / (in_path.stem + ".enriched" + in_path.suffix))

    LOG.info(f"Input:  {args.input}")
    LOG.info(f"Output: {output_path}")
    LOG.info(f"Server: {args.server}")
    LOG.info(f"Path depth: {args.path_depth} nodes")

    # Check server
    if not check_server(args.server):
        LOG.error("Cannot proceed without haplo server. Start it with:")
        LOG.error("  cd ~/ystr-matcher/ftdna_haplo && pm2 start ecosystem.config.js")
        sys.exit(1)

    # Enrich
    enrich_file(args.input, output_path, args.server, args.path_depth)


if __name__ == "__main__":
    main()
