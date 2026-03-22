"""
YFull Deepener — post-processing tool for Yleaf predictions.

Takes an Yleaf prediction, finds deeper branches in YFull tree,
looks up SNP addresses from YFull CSV, checks them in BAM via samtools.

Usage:
    python3 yfull_deepener.py \
        --prediction "G-Y219854" \
        --bam /path/to/sample.bam \
        --yfull-tree /path/to/current_tree.json \
        --csv /path/to/yfull_snp.csv \
        --samtools /path/to/samtools
"""

import csv
import json
import subprocess
import sys
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_yfull_tree(tree_path):
    """Load YFull tree and build id→node lookup + parent map."""
    with open(tree_path) as f:
        root = json.load(f)

    nodes = {}       # id → node
    parent_map = {}  # id → parent_id
    snp_map = {}     # id → list of SNP names

    def walk(node, parent_id=None):
        nid = node.get("id", "")
        if nid:
            nodes[nid] = node
            parent_map[nid] = parent_id
            # Parse SNPs
            snps_str = node.get("snps", "")
            snp_list = []
            if snps_str:
                for group in snps_str.split(","):
                    for alias in group.strip().split("/"):
                        alias = alias.strip()
                        if alias:
                            snp_list.append(alias)
            snp_map[nid] = snp_list
        for ch in node.get("children", []):
            walk(ch, nid)

    walk(root)
    return nodes, parent_map, snp_map


def get_children_recursive(nodes, branch_id, max_depth=None):
    """Get all descendant branches. No depth limit by default."""
    result = []

    def walk(node, depth):
        if max_depth is not None and depth > max_depth:
            return
        for ch in node.get("children", []):
            ch_id = ch.get("id", "")
            if ch_id and not ch_id.endswith("*"):
                result.append((ch_id, depth))
            walk(ch, depth + 1)

    if branch_id in nodes:
        walk(nodes[branch_id], 1)
    return result


def load_csv_addresses(csv_path):
    """Load SNP → (build38, anc, der) from YFull CSV.
    Only confirmed SNPs with build38 addresses.
    """
    snp_addr = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            snp_id = row.get("snp_id", "").strip()
            confirmed = row.get("confirmed", "").strip()
            build38 = row.get("build38", "").strip()
            anc = row.get("anc", "").strip()
            der = row.get("der", "").strip()
            if confirmed == "Yes" and build38 and anc and der and snp_id:
                snp_addr[snp_id] = {
                    "pos": int(build38),
                    "anc": anc,
                    "der": der,
                }
    return snp_addr


def load_yleaf_positions(positions_path):
    """Load SNP → (pos, anc, der) from Yleaf new_positions.txt.
    Format: chry\tSNP_NAME\tHAPLOGROUP\tPOS\tMUTATION\tANC\tDER
    Used as fallback when SNP is not in CSV.
    """
    snp_addr = {}
    with open(positions_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                snp_name = parts[1]
                pos = parts[3]
                anc = parts[5]
                der = parts[6]
                if pos and anc and der:
                    try:
                        snp_addr[snp_name] = {
                            "pos": int(pos),
                            "anc": anc,
                            "der": der,
                        }
                    except ValueError:
                        pass
    return snp_addr


def check_positions_in_bam(bam_path, positions, samtools="samtools",
                           min_depth=3):
    """Check specific chrY positions in BAM via samtools mpileup.

    Args:
        bam_path: Path to BAM file.
        positions: List of (pos, anc, der, snp_name, branch_id).
        samtools: Path to samtools binary.
        min_depth: Minimum read depth to trust.

    Returns:
        List of dicts with check results.
    """
    if not positions:
        return []

    # Create temp BED for mpileup
    import tempfile
    bed_lines = []
    for pos, anc, der, snp_name, branch_id in positions:
        # BED is 0-based, mpileup position is 1-based
        bed_lines.append(f"chrY\t{pos-1}\t{pos}")

    # Deduplicate and sort
    bed_lines = sorted(set(bed_lines))

    tmp_bed = tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    )
    tmp_bed.write("\n".join(bed_lines) + "\n")
    tmp_bed.close()

    # Run samtools mpileup
    cmd = [
        samtools, "mpileup",
        "-r", "chrY",
        "-l", tmp_bed.name,
        "-AQ20q1",
        bam_path,
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
        if proc.returncode != 0:
            logger.error(f"samtools error: {proc.stderr[:200]}")
            return []
    finally:
        Path(tmp_bed.name).unlink(missing_ok=True)

    # Parse mpileup output → pos → (ref_base, depth, bases)
    pileup = {}
    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 5:
            pos = int(parts[1])
            ref = parts[2]
            depth = int(parts[3])
            bases = parts[4]
            pileup[pos] = (ref, depth, bases)

    # Check each position
    results = []
    for pos, anc, der, snp_name, branch_id in positions:
        if pos in pileup:
            ref, depth, bases = pileup[pos]
            # Count derived allele
            der_count = bases.upper().count(der.upper())
            anc_count = bases.upper().count(anc.upper())
            # Also count reference matches (. and ,)
            if ref.upper() == anc.upper():
                anc_count += bases.count(".") + bases.count(",")
            elif ref.upper() == der.upper():
                der_count += bases.count(".") + bases.count(",")

            if depth >= min_depth:
                if der_count > anc_count:
                    state = "DERIVED"
                elif anc_count > der_count:
                    state = "ANCESTRAL"
                else:
                    state = "UNCERTAIN"
            else:
                state = "LOW_DEPTH"

            results.append({
                "branch": branch_id,
                "snp": snp_name,
                "pos": pos,
                "mutation": f"{anc}>{der}",
                "depth": depth,
                "der_count": der_count,
                "anc_count": anc_count,
                "state": state,
            })
        else:
            results.append({
                "branch": branch_id,
                "snp": snp_name,
                "pos": pos,
                "mutation": f"{anc}>{der}",
                "depth": 0,
                "der_count": 0,
                "anc_count": 0,
                "state": "NO_DATA",
            })

    return results


def deepen(prediction, bam_path, yfull_tree_path, csv_path,
           samtools="samtools", max_depth=None, min_depth=3,
           positions_path=None):
    """Main deepening function.

    1. Find children of prediction in YFull tree
    2. For each child, find SNPs with addresses in CSV
    3. Check those positions in BAM
    4. Report which deeper branches are derived
    """
    logger.info(f"=== YFull Deepener ===")
    logger.info(f"Prediction: {prediction}")
    logger.info(f"BAM: {bam_path}")

    # Step 1: Load YFull tree
    logger.info("Loading YFull tree...")
    nodes, parent_map, snp_names = load_yfull_tree(yfull_tree_path)
    logger.info(f"  {len(nodes)} branches")

    # Find prediction in tree
    if prediction not in nodes:
        # Try stripping *(x...) suffix
        base = prediction.split("*")[0]
        if base in nodes:
            prediction = base
        else:
            logger.error(f"  Prediction '{prediction}' not found in YFull tree!")
            return None

    # Step 2: Get children
    children = get_children_recursive(nodes, prediction, max_depth)
    logger.info(f"  {len(children)} descendant branches (depth 1-{max_depth})")

    if not children:
        logger.info("  No deeper branches found — already at terminal!")
        return {"prediction": prediction, "deepened": None, "checks": []}

    # Step 3: Load CSV addresses + Yleaf positions fallback
    logger.info("Loading CSV addresses...")
    csv_addrs = load_csv_addresses(csv_path)
    logger.info(f"  {len(csv_addrs)} SNPs with addresses from CSV")

    yleaf_addrs = {}
    if positions_path:
        logger.info("Loading Yleaf positions (fallback)...")
        yleaf_addrs = load_yleaf_positions(positions_path)
        logger.info(f"  {len(yleaf_addrs)} SNPs from Yleaf positions")

    # Step 4: Match children SNPs to addresses (CSV first, then Yleaf)
    positions_to_check = []
    branches_with_snps = 0
    csv_hits = 0
    yleaf_hits = 0
    for branch_id, depth in children:
        branch_snps = snp_names.get(branch_id, [])
        found = False
        for snp in branch_snps:
            addr = None
            if snp in csv_addrs:
                addr = csv_addrs[snp]
                csv_hits += 1
            elif snp in yleaf_addrs:
                addr = yleaf_addrs[snp]
                yleaf_hits += 1
            if addr:
                positions_to_check.append((
                    addr["pos"], addr["anc"], addr["der"],
                    snp, branch_id
                ))
                found = True
        if found:
            branches_with_snps += 1

    logger.info(
        f"  {len(positions_to_check)} SNPs to check "
        f"across {branches_with_snps} branches "
        f"(CSV: {csv_hits}, Yleaf: {yleaf_hits})"
    )

    if not positions_to_check:
        logger.info("  No SNPs with addresses found for deeper branches")
        return {"prediction": prediction, "deepened": None, "checks": []}

    # Step 5: Check in BAM
    logger.info(f"Checking {len(positions_to_check)} positions in BAM...")
    results = check_positions_in_bam(
        bam_path, positions_to_check, samtools, min_depth
    )

    # Step 6: Aggregate by branch
    branch_results = {}
    for r in results:
        bid = r["branch"]
        if bid not in branch_results:
            branch_results[bid] = {"derived": 0, "ancestral": 0,
                                    "total": 0, "snps": []}
        branch_results[bid]["total"] += 1
        if r["state"] == "DERIVED":
            branch_results[bid]["derived"] += 1
        elif r["state"] == "ANCESTRAL":
            branch_results[bid]["ancestral"] += 1
        branch_results[bid]["snps"].append(r)

    # Step 7: Transitive deepening (bottom-up)
    # If a descendant is derived, all ancestors up to prediction are too.
    # Build child→depth lookup
    child_depth = {bid: d for bid, d in children}

    # Find all directly derived branches
    derived_set = set()
    for branch_id, br in branch_results.items():
        if br["derived"] > 0 and br["derived"] >= br["ancestral"]:
            derived_set.add(branch_id)

    # For each derived branch, mark all ancestors up to prediction
    inferred_derived = set()
    for dbranch in derived_set:
        cur = dbranch
        while cur and cur != prediction:
            inferred_derived.add(cur)
            cur = parent_map.get(cur)

    # Find deepest derived branch (by depth, prefer directly tested)
    deepened = None
    derived_path = []
    all_derived = derived_set | inferred_derived
    for branch_id, depth in sorted(children, key=lambda x: x[1]):
        if branch_id in all_derived:
            derived_path.append(branch_id)
            deepened = branch_id

    # Print results
    print(f"\n{'='*70}")
    print(f"Prediction: {prediction}")
    print(f"Deeper branches checked: {len(branch_results)}")
    if inferred_derived - derived_set:
        print(f"Inferred derived (transitive): "
              f"{inferred_derived - derived_set}")
    print(f"{'='*70}")
    print(f"{'Branch':<30} {'D':>3} {'A':>3} {'Total':>5}  Status")
    print("-" * 70)

    for branch_id, depth in sorted(children, key=lambda x: x[1]):
        br = branch_results.get(branch_id)
        if br:
            if branch_id in derived_set:
                status = "✅ DERIVED"
            elif branch_id in inferred_derived:
                status = "✅ DERIVED (inferred)"
            else:
                status = "❌ ancestral"
            print(
                f"{'  ' * (depth-1)}{branch_id:<{30 - 2*(depth-1)}} "
                f"{br['derived']:>3} {br['ancestral']:>3} "
                f"{br['total']:>5}  {status}"
            )
        elif branch_id in inferred_derived:
            # Branch had no SNPs in CSV but is inferred derived
            print(
                f"{'  ' * (depth-1)}{branch_id:<{30 - 2*(depth-1)}} "
                f"  -   -     -  ✅ DERIVED (inferred from child)"
            )

    if deepened:
        print(f"\n🔥 DEEPENED: {prediction} → {deepened}")
        if len(derived_path) > 1:
            print(f"   Path: {' → '.join(derived_path)}")
    else:
        print(f"\n— No deeper branches confirmed")

    return {
        "prediction": prediction,
        "deepened": deepened,
        "derived_path": derived_path,
        "branch_results": branch_results,
    }


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="YFull post-Yleaf deepener")
    parser.add_argument("--prediction", required=True,
                        help="Yleaf haplogroup prediction (e.g. G-Y219854)")
    parser.add_argument("--bam", required=True, help="BAM file path")
    parser.add_argument("--yfull-tree", required=True,
                        help="YFull current_tree.json")
    parser.add_argument("--csv", required=True,
                        help="YFull yfull_snp CSV with addresses")
    parser.add_argument("--positions", default=None,
                        help="Yleaf new_positions.txt (fallback for SNPs "
                             "not in CSV)")
    parser.add_argument("--samtools",
                        default="/home/valalav/.local/bin/samtools")
    parser.add_argument("--max-depth", type=int, default=None,
                        help="Max descent depth (default: unlimited)")
    parser.add_argument("--min-depth", type=int, default=3)
    args = parser.parse_args()

    result = deepen(
        prediction=args.prediction,
        bam_path=args.bam,
        yfull_tree_path=args.yfull_tree,
        csv_path=args.csv,
        positions_path=args.positions,
        samtools=args.samtools,
        max_depth=args.max_depth,
        min_depth=args.min_depth,
    )

    # Iterative deepening: keep going from deepened result
    if result and result.get("deepened"):
        iteration = 1
        while True:
            new_pred = result["deepened"]
            logger.info(f"\n=== Iteration {iteration + 1}: "
                        f"re-deepening from {new_pred} ===")
            next_result = deepen(
                prediction=new_pred,
                bam_path=args.bam,
                yfull_tree_path=args.yfull_tree,
                csv_path=args.csv,
                positions_path=args.positions,
                samtools=args.samtools,
                max_depth=args.max_depth,
                min_depth=args.min_depth,
            )
            if next_result and next_result.get("deepened"):
                result = next_result
                iteration += 1
            else:
                break
        if iteration > 1:
            print(f"\n🏁 FINAL (after {iteration} iterations): "
                  f"{args.prediction} → {result['deepened']}")
