"""
YFull SNP Integrator for Yleaf (v2).

Maps SNPs from YFull CSV to branches via YFull tree JSON.
Tree structure is NOT modified — only new_positions.txt gets new markers.

Usage:
    python3 -m yleaf.yfull_integrator \\
        /path/to/yfull_snp.csv \\
        --yfull-tree /path/to/current_tree.json \\
        --positions yleaf/data/hg38/new_positions.txt
"""

import csv
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def build_snp_to_branch_map(yfull_tree_path):
    """Build SNP name → branch mapping from YFull tree JSON.

    YFull tree is recursive: each node has 'id' (branch name),
    'snps' (comma-separated SNP names with aliases), 'children'.
    Example snps: "YP2787, FGC25750/YP2704/V2222, L1284"

    For each node, we split snps by comma, then by slash (aliases),
    and map every individual SNP name to the node's id (branch).

    Returns:
        Dict mapping snp_name → branch_id (e.g. "BY122562" → "G-Y219854").
    """
    with open(yfull_tree_path, "r") as f:
        data = json.load(f)

    snp_map = {}

    def walk(node):
        node_id = node.get("id", "")
        snps_str = node.get("snps", "")
        if node_id and snps_str:
            # Parse: "YP2787, FGC25750/YP2704/V2222, L1284"
            for group in snps_str.split(","):
                group = group.strip()
                if not group:
                    continue
                for alias in group.split("/"):
                    alias = alias.strip()
                    if alias:
                        snp_map[alias] = node_id
        for child in node.get("children", []):
            walk(child)

    walk(data)
    return snp_map


def load_existing_snps(positions_path):
    """Load existing SNP names from Yleaf new_positions.txt.

    Returns:
        Set of snp_names.
    """
    snp_names = set()
    with open(positions_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                snp_names.add(parts[1])
    return snp_names


def filter_and_map(csv_path, existing_snps, snp_branch_map):
    """Filter YFull CSV and map each SNP to its branch via YFull tree.

    Criteria:
        1. confirmed == "Yes"
        2. build38 is not empty
        3. anc and der are not empty
        4. snp_id not already in Yleaf
        5. snp_id found in YFull tree (has a branch)

    Args:
        csv_path: Path to yfull_snp_p1-p658.csv.
        existing_snps: Set of SNP names already in Yleaf.
        snp_branch_map: Dict from build_snp_to_branch_map().

    Returns:
        Tuple of (filtered_list, skipped_dict).
    """
    filtered = []
    skipped = {
        "no_confirmed": 0,
        "no_build38": 0,
        "no_alleles": 0,
        "duplicate": 0,
        "not_in_tree": 0,
    }

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            snp_id = row.get("snp_id", "").strip()
            confirmed = row.get("confirmed", "").strip()
            build38 = row.get("build38", "").strip()
            anc = row.get("anc", "").strip()
            der = row.get("der", "").strip()

            if confirmed != "Yes":
                skipped["no_confirmed"] += 1
                continue
            if not build38:
                skipped["no_build38"] += 1
                continue
            if not anc or not der:
                skipped["no_alleles"] += 1
                continue
            if snp_id in existing_snps:
                skipped["duplicate"] += 1
                continue

            # Map SNP to branch via YFull tree
            branch = snp_branch_map.get(snp_id)
            if not branch:
                skipped["not_in_tree"] += 1
                continue

            filtered.append({
                "snp_id": snp_id,
                "branch": branch,
                "build38": build38,
                "anc": anc,
                "der": der,
            })

    logger.info(
        f"Filtered {len(filtered)} SNPs mapped to tree branches. "
        f"Skipped: {skipped}"
    )
    return filtered, skipped


def convert_to_yleaf_format(filtered_snps):
    """Convert filtered SNPs to Yleaf new_positions.txt format.

    Output format (TSV): chry  snp_name  haplogroup  pos  mutation  ref  alt
    """
    lines = []
    for snp in filtered_snps:
        line = (
            f"chry\t{snp['snp_id']}\t{snp['branch']}\t{snp['build38']}\t"
            f"{snp['anc']}>{snp['der']}\t{snp['anc']}\t{snp['der']}"
        )
        lines.append(line)
    return lines


def merge_positions(existing_path, new_lines, output_path=None):
    """Merge new SNP lines into existing new_positions.txt.

    Deduplicates by snp_name (column 2).
    """
    if output_path is None:
        output_path = existing_path

    existing_lines = []
    existing_names = set()
    with open(existing_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped:
                parts = stripped.split("\t")
                if len(parts) >= 2:
                    existing_names.add(parts[1])
                existing_lines.append(stripped)

    added = 0
    for new_line in new_lines:
        parts = new_line.split("\t")
        if len(parts) >= 2 and parts[1] not in existing_names:
            existing_lines.append(new_line)
            existing_names.add(parts[1])
            added += 1

    with open(output_path, "w") as f:
        for line in existing_lines:
            f.write(line + "\n")

    logger.info(f"Merged: {len(existing_lines)} total, {added} new added")
    return len(existing_lines), added


def integrate_yfull(csv_path, positions_path, yfull_tree_path,
                     output_positions=None, dry_run=False):
    """Main integration: map CSV SNPs to YFull tree branches, merge positions.

    Tree.json is NOT modified — branches come from YFull JSON only.
    """
    logging.basicConfig(level=logging.INFO)

    logger.info("=== YFull SNP Integration v2 ===")
    logger.info(f"CSV: {csv_path}")
    logger.info(f"YFull tree: {yfull_tree_path}")
    logger.info(f"Positions: {positions_path}")

    # Step 1: Build SNP → branch map from YFull tree
    logger.info("Step 1: Building SNP→branch map from YFull tree...")
    snp_branch_map = build_snp_to_branch_map(yfull_tree_path)
    logger.info(f"  Mapped {len(snp_branch_map)} SNPs to branches")

    # Step 2: Load existing Yleaf SNPs
    logger.info("Step 2: Loading existing Yleaf positions...")
    existing_snps = load_existing_snps(positions_path)
    logger.info(f"  Existing: {len(existing_snps)} SNPs")

    # Step 3: Filter CSV and map to branches
    logger.info("Step 3: Filtering and mapping CSV SNPs...")
    filtered, skipped = filter_and_map(csv_path, existing_snps, snp_branch_map)
    logger.info(f"  Result: {len(filtered)} new SNPs with branches")

    # Count unique branches
    branches = set(s["branch"] for s in filtered)

    stats = {
        "yfull_tree_snps": len(snp_branch_map),
        "existing_yleaf_snps": len(existing_snps),
        "new_snps_mapped": len(filtered),
        "unique_branches": len(branches),
        "skipped": skipped,
    }

    if dry_run:
        logger.info("DRY RUN: No files written.")
        logger.info(f"Stats: {stats}")
        return stats

    # Step 4: Convert and merge
    logger.info("Step 4: Converting and merging positions...")
    new_lines = convert_to_yleaf_format(filtered)
    total, added = merge_positions(
        positions_path, new_lines, output_positions
    )
    stats["total_positions"] = total
    stats["positions_added"] = added

    logger.info("=== Integration Complete ===")
    logger.info(f"Stats: {stats}")
    return stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Integrate YFull SNPs into Yleaf (v2: tree-mapped)"
    )
    parser.add_argument("csv_path", help="Path to yfull_snp CSV")
    parser.add_argument(
        "--yfull-tree", required=True,
        help="Path to YFull current_tree.json (REQUIRED)"
    )
    parser.add_argument(
        "--positions",
        default="yleaf/data/hg38/new_positions.txt",
        help="Path to Yleaf new_positions.txt",
    )
    parser.add_argument(
        "--output-positions", default=None,
        help="Output path for merged positions (default: overwrite)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show stats without writing files",
    )
    args = parser.parse_args()

    stats = integrate_yfull(
        csv_path=args.csv_path,
        positions_path=args.positions,
        yfull_tree_path=args.yfull_tree,
        output_positions=args.output_positions,
        dry_run=args.dry_run,
    )
    print(f"\nResult: {json.dumps(stats, indent=2)}")
