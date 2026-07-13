#!/usr/bin/env python3
"""
Enrich Yleaf SNP positions with FTDNA tree data.

Creates SEPARATE enriched files (does NOT modify originals):
  - new_positions_enriched.txt  (hg38 only, FTDNA uses hg38 coordinates)
  - tree_enriched.json          (merged YFull + FTDNA hierarchy)

Usage:
  python3 enrich_from_ftdna.py [--ftdna PATH] [--dry-run]

Default FTDNA path: /home/valalav/_dna/ystr-matcher/ftdna_haplo/data/get.json
"""

import json
import argparse
import logging
import sys
from pathlib import Path
from collections import defaultdict

LOG = logging.getLogger("enrich_ftdna")
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Paths
YLEAF_DATA = Path(__file__).parent / "data"
HG38_POSITIONS = YLEAF_DATA / "hg38" / "new_positions.txt"
HG38_ENRICHED = YLEAF_DATA / "hg38" / "new_positions_enriched.txt"
TREE_FILE = YLEAF_DATA / "hg_prediction_tables" / "tree.json"
TREE_ENRICHED = YLEAF_DATA / "hg_prediction_tables" / "tree_enriched.json"

DEFAULT_FTDNA = Path("/home/valalav/_dna/ystr-matcher/ftdna_haplo/data/get.json")


def load_existing_positions(path: Path) -> dict:
    """Load existing new_positions.txt into {SNP_NAME_UPPER: full_line}."""
    existing = {}
    if not path.exists():
        LOG.warning(f"File not found: {path}")
        return existing

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                snp_name = parts[1].upper()
                existing[snp_name] = line.strip()
    LOG.info(f"Loaded {len(existing)} existing SNPs from {path.name}")
    return existing


def load_ftdna_tree(path: Path) -> dict:
    """Load FTDNA get.json and return allNodes dict."""
    LOG.info(f"Loading FTDNA tree from {path} ({path.stat().st_size / 1024 / 1024:.1f} MB)...")
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    nodes = data.get("allNodes", {})
    LOG.info(f"Loaded {len(nodes)} FTDNA nodes")
    return nodes


def extract_ftdna_snps(nodes: dict) -> dict:
    """
    Extract SNPs with valid hg38 positions from FTDNA nodes.
    Returns {SNP_NAME_UPPER: (position, ancestral, derived, haplogroup_name)}
    """
    snps = {}
    neg_count = 0
    no_pos = 0

    for node_id, node in nodes.items():
        hg_name = node.get("name", "")
        for v in node.get("variants", []):
            variant_name = v.get("variant", "")
            pos = v.get("position", 0)
            anc = v.get("ancestral", "")
            der = v.get("derived", "")

            if not variant_name:
                continue
            if not pos:
                no_pos += 1
                continue
            if pos < 0:
                neg_count += 1
                continue
            if not anc or not der:
                continue

            snps[variant_name.upper()] = (str(pos), anc, der, hg_name)

    LOG.info(f"Extracted {len(snps)} SNPs with valid hg38 positions")
    LOG.info(f"  Skipped: {neg_count} negative positions, {no_pos} missing positions")
    return snps


def enrich_positions(
    existing: dict, ftdna_snps: dict, output_path: Path, dry_run: bool = False
) -> int:
    """
    Merge existing positions with FTDNA SNPs → write enriched file.
    Returns count of new SNPs added.
    """
    new_count = 0
    new_lines = []

    for snp_name, (pos, anc, der, hg_name) in ftdna_snps.items():
        if snp_name in existing:
            continue  # Already in Yleaf

        mutation = f"{anc}->{der}"
        # Use the original case from FTDNA for the SNP name
        # Find original case name
        orig_name = snp_name  # fallback
        for v_name in [snp_name]:  # already uppercased
            orig_name = v_name
            break

        line = f"chry\t{orig_name}\t{hg_name}\t{pos}\t{mutation}\t{anc}\t{der}"
        new_lines.append(line)
        new_count += 1

    if dry_run:
        LOG.info(f"[DRY RUN] Would add {new_count} new SNPs to {output_path.name}")
        return new_count

    # Write: first all existing lines, then new ones
    with open(output_path, "w", encoding="utf-8") as f:
        for line in existing.values():
            f.write(line + "\n")
        for line in new_lines:
            f.write(line + "\n")

    total = len(existing) + new_count
    LOG.info(f"Written {total} total SNPs to {output_path.name}")
    LOG.info(f"  Existing: {len(existing)}, New from FTDNA: {new_count}")
    return new_count


def load_existing_tree(path: Path) -> dict:
    """Load existing Yleaf tree.json."""
    if not path.exists():
        LOG.warning(f"Tree file not found: {path}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_ftdna_tree_dict(nodes: dict) -> dict:
    """
    Build an Yleaf-compatible tree dict from FTDNA flat graph.
    Format: {"parent_name": ["child1_name", "child2_name", "parent_name*"]}
    """
    tree = {}

    # Build name→node and id→name maps
    id_to_name = {}
    for node_id, node in nodes.items():
        name = node.get("name", "")
        if name:
            id_to_name[node.get("haplogroupId", int(node_id))] = name

    for node_id, node in nodes.items():
        name = node.get("name", "")
        if not name:
            continue

        children_ids = node.get("children", [])
        child_names = []
        for cid in children_ids:
            cname = id_to_name.get(cid)
            if cname:
                child_names.append(cname)

        if child_names:
            # Add star paragroup
            star_name = f"{name}*"
            child_names.append(star_name)
            tree[star_name] = []
            tree[name] = child_names
        else:
            tree[name] = []

    LOG.info(f"Built FTDNA tree with {len(tree)} nodes")
    return tree


def merge_trees(
    existing_tree: dict, ftdna_tree: dict, output_path: Path, dry_run: bool = False
) -> int:
    """
    Merge FTDNA tree into existing Yleaf tree.
    Only adds nodes/branches that don't exist yet.
    Returns count of new nodes added.
    """
    merged = dict(existing_tree)
    new_count = 0

    for node_name, children in ftdna_tree.items():
        if node_name not in merged:
            merged[node_name] = children
            new_count += 1
        else:
            # Node exists — add any missing children
            existing_children = set(merged[node_name])
            for child in children:
                if child not in existing_children:
                    merged[node_name].append(child)
                    # Make sure child node also exists
                    if child not in merged:
                        merged[child] = ftdna_tree.get(child, [])
                        new_count += 1

    if dry_run:
        LOG.info(f"[DRY RUN] Would add {new_count} new tree nodes to {output_path.name}")
        return new_count

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(merged, f, indent=2, ensure_ascii=False)

    LOG.info(f"Written merged tree to {output_path.name}")
    LOG.info(f"  Existing nodes: {len(existing_tree)}, New from FTDNA: {new_count}")
    LOG.info(f"  Total nodes: {len(merged)}")
    return new_count


def main():
    parser = argparse.ArgumentParser(description="Enrich Yleaf with FTDNA tree data")
    parser.add_argument(
        "--ftdna", type=Path, default=DEFAULT_FTDNA,
        help=f"Path to FTDNA get.json (default: {DEFAULT_FTDNA})"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Only show what would be done, don't write files"
    )
    parser.add_argument(
        "--positions-only", action="store_true",
        help="Only enrich positions, skip tree merge"
    )
    parser.add_argument(
        "--tree-only", action="store_true",
        help="Only merge tree, skip positions"
    )
    args = parser.parse_args()

    if not args.ftdna.exists():
        LOG.error(f"FTDNA file not found: {args.ftdna}")
        sys.exit(1)

    # Load FTDNA
    ftdna_nodes = load_ftdna_tree(args.ftdna)

    # Enrich positions (hg38)
    if not args.tree_only:
        existing = load_existing_positions(HG38_POSITIONS)
        ftdna_snps = extract_ftdna_snps(ftdna_nodes)
        new_snps = enrich_positions(existing, ftdna_snps, HG38_ENRICHED, args.dry_run)
        LOG.info(f"Positions: +{new_snps} new SNPs")

    # Merge tree
    if not args.positions_only:
        existing_tree = load_existing_tree(TREE_FILE)
        ftdna_tree = build_ftdna_tree_dict(ftdna_nodes)
        new_nodes = merge_trees(existing_tree, ftdna_tree, TREE_ENRICHED, args.dry_run)
        LOG.info(f"Tree: +{new_nodes} new nodes")

    LOG.info("Done! Enriched files created (originals untouched):")
    if not args.tree_only:
        LOG.info(f"  Positions: {HG38_ENRICHED}")
    if not args.positions_only:
        LOG.info(f"  Tree:      {TREE_ENRICHED}")
    LOG.info("")
    LOG.info("To USE enriched data with Yleaf, swap files:")
    LOG.info(f"  cp {HG38_POSITIONS} {HG38_POSITIONS}.bak")
    LOG.info(f"  cp {HG38_ENRICHED} {HG38_POSITIONS}")
    LOG.info(f"  cp {TREE_FILE} {TREE_FILE}.bak")
    LOG.info(f"  cp {TREE_ENRICHED} {TREE_FILE}")


if __name__ == "__main__":
    main()
