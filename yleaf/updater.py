import json
import csv
import os
import sys
import urllib.request
import gzip
from pathlib import Path

print("Starting updater...", flush=True)

# Constants
YFULL_TREE_URL = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json"
YBROWSE_T2T_VCF_URL = "https://ybrowse.org/gbrowse2/gff/CP086569.1/snps_CP086569.1_concat_sorted.vcf.gz"
YBROWSE_HG38_URL = "https://ybrowse.org/gbrowse2/gff/snps_hg38.csv"

def download_file(url, output_path):
    print(f"Downloading {url} to {output_path}...", flush=True)
    try:
        req = urllib.request.Request(
            url, 
            data=None, 
            headers={
                'User-Agent': 'Mozilla/5.0'
            }
        )
        with urllib.request.urlopen(req) as response, open(output_path, 'wb') as out_file:
            data = response.read()
            out_file.write(data)
        print("Download complete.", flush=True)
        return True
    except Exception as e:
        print(f"Error downloading {url}: {e}", flush=True)
        return False

def load_yfull_tree(json_path):
    print("Loading YFull tree...", flush=True)
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

def parse_yfull_node_map(node, parent_hg, result_map):
    hg_name = node.get('id')
    snps = node.get('snps', '').replace(' ', '').split(',')
    if hg_name:
        for snp in snps:
            if snp:
                result_map[snp] = hg_name
    for child in node.get('children', []):
        parse_yfull_node_map(child, hg_name, result_map)

def build_yleaf_tree_dict(node, tree_dict):
    """
    Recursively builds Yleaf adjacency list from YFull tree.
    Format: {"Parent": ["Child1", "Child2", "Parent*"]}
    """
    hg_name = node.get('id')
    children = node.get('children', [])
    
    # Check explicitly for None, because empty string "" is a valid name for Root in YFull
    if hg_name is not None:
        child_names = [child.get('id') for child in children if child.get('id') is not None]
        
        # If there are children, Yleaf usually adds a "star" version of the parent 
        # to represent the paragroup (belonging to parent but not any known child).
        # e.g. "I-Z2336": ["I-Y3866", "I-Z2336*"]
        if child_names:
            if hg_name: # Only add star group if name is not empty
                star_name = f"{hg_name}*"
                child_names.append(star_name)
                # IMPORTANT: The star node must also exist in the dict as a leaf
                tree_dict[star_name] = []
                
            tree_dict[hg_name] = child_names
        else:
            # No children -> Leaf node
            tree_dict[hg_name] = []
            
        for child in children:
            build_yleaf_tree_dict(child, tree_dict)
    # else block removed because logic is now inside if hg_name is not None

def convert_ybrowse_vcf_gz(vcf_path, snp_to_hg_map, output_txt_path):
    """
    Parses VCF.gz for T2T.
    VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO...
    INFO usually has attributes. But ID is usually the SNP name.
    """
    print(f"Processing VCF.gz {vcf_path} -> {output_txt_path}", flush=True)
    count_written = 0
    count_skipped = 0
    
    with gzip.open(vcf_path, 'rt', encoding='utf-8') as f_in, open(output_txt_path, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            if line.startswith("#"): continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5: continue
            
            # CHROM POS ID REF ALT
            # 0     1   2  3   4
            
            # chrom = parts[0] # CP086569.1
            pos = parts[1]
            snp_name = parts[2]
            ref = parts[3]
            alt = parts[4]
            
            if snp_name == ".": continue
            
            # Clean name if multiple IDs
            # Sometimes VCF ID is "rs123;Z2336"
            
            found_hg = None
            used_snp_name = snp_name
            
            for sub_name in snp_name.split(';'):
                if sub_name in snp_to_hg_map:
                    found_hg = snp_to_hg_map[sub_name]
                    used_snp_name = sub_name
                    break
            
            if not found_hg:
                count_skipped += 1
                continue
            
            mutation = f"{ref}->{alt}"
            chrom = "chry"
            
            f_out.write(f"{chrom}\t{used_snp_name}\t{found_hg}\t{pos}\t{mutation}\t{ref}\t{alt}\n")
            count_written += 1
            
    print(f"Written (T2T): {count_written}, Skipped: {count_skipped}", flush=True)

def normalize_ybrowse_haplogroup(raw_hg):
    """
    Normalize YBrowse YCC_haplogroup to Yleaf-compatible format.
    Examples:
      'R-P312' -> 'R-P312'
      'R1b'    -> 'R1b'
      'R1b-P312 (not listed)' -> 'R-P312' (strip parenthetical)
      'unknown' -> None
      'not listed' -> None
    """
    if not raw_hg:
        return None
    raw_hg = raw_hg.strip().strip('"')
    
    # Skip useless values
    if raw_hg.lower() in ('unknown', 'not listed', '.', ''):
        return None
    
    # Remove parenthetical annotations: "R1b-DF21 (not listed)" -> "R1b-DF21"
    if '(' in raw_hg:
        raw_hg = raw_hg[:raw_hg.index('(')].strip()
    
    # After stripping, check again
    if not raw_hg or raw_hg.lower() in ('unknown', 'not listed'):
        return None
    
    return raw_hg


def convert_ybrowse_csv(csv_path, snp_to_hg_map, output_txt_path):
    """
    Process YBrowse hg38 CSV with dual mapping:
    1. Primary: YFull tree SNP->haplogroup map (precise)
    2. Fallback: YBrowse's own YCC_haplogroup column (broader coverage)
    
    Returns set of new haplogroup names found via fallback (for tree enrichment).
    """
    print(f"Processing CSV {csv_path} -> {output_txt_path}", flush=True)
    count_yfull = 0
    count_ybrowse = 0
    count_skipped = 0
    new_haplogroups = set()  # haplogroups from YBrowse not in YFull tree
    seen_positions = set()   # dedup by position
    
    with open(csv_path, 'r', encoding='utf-8') as f_in, open(output_txt_path, 'w', encoding='utf-8') as f_out:
        reader = csv.DictReader(f_in)
        reader.fieldnames = [x.lower() for x in reader.fieldnames]
        
        for row in reader:
            # Columns: Name, ID, allele_anc, allele_der, start
            snp_name = row.get('name') or row.get('id')
            if not snp_name: continue
            snp_name = snp_name.strip().strip('"')
            
            pos = row.get('start') or row.get('position')
            if not pos: continue
            pos = pos.strip().strip('"')
            
            # Dedup by SNP name + position
            key = f"{snp_name}_{pos}"
            if key in seen_positions:
                continue
            seen_positions.add(key)
            
            # allele_anc / allele_der for hg38 YBrowse
            ref = row.get('allele_anc') or row.get('ref')
            alt = row.get('allele_der') or row.get('alt')
            if not ref or not alt: continue
            ref = ref.strip().strip('"')
            alt = alt.strip().strip('"')
            
            # Skip non-standard alleles (insertions, deletions, complex)
            if len(ref) > 1 or len(alt) > 1:
                continue
            if ref not in 'ACGT' or alt not in 'ACGT':
                continue
            
            # === Dual mapping ===
            # Primary: YFull tree map (precise haplogroup names)
            found_hg = snp_to_hg_map.get(snp_name)
            source = 'yfull'
            
            # Fallback: YBrowse's own YCC_haplogroup
            if not found_hg:
                ycc_hg = row.get('ycc_haplogroup', '')
                found_hg = normalize_ybrowse_haplogroup(ycc_hg)
                if found_hg:
                    source = 'ybrowse'
                    new_haplogroups.add(found_hg)
            
            if not found_hg:
                count_skipped += 1
                continue
            
            mutation = f"{ref}->{alt}"
            chrom = "chry"
            f_out.write(f"{chrom}\t{snp_name}\t{found_hg}\t{pos}\t{mutation}\t{ref}\t{alt}\n")
            
            if source == 'yfull':
                count_yfull += 1
            else:
                count_ybrowse += 1
            
    total = count_yfull + count_ybrowse
    print(f"Written (hg38): {total} total ({count_yfull} YFull + {count_ybrowse} YBrowse fallback), Skipped: {count_skipped}", flush=True)
    print(f"New haplogroups from YBrowse: {len(new_haplogroups)}", flush=True)
    return new_haplogroups


def main():
    script_dir = Path(__file__).parent  # yleaf/
    base_dir = script_dir / "data"
    temp_dir = script_dir.parent / "temp_update"
    temp_dir.mkdir(exist_ok=True)
    
    # 1. YFull Tree
    yfull_json_path = temp_dir / "current_tree.json"
    # Always redownload if needed or check size
    if not yfull_json_path.exists():
        download_file(YFULL_TREE_URL, yfull_json_path)
    
    yfull_data = load_yfull_tree(yfull_json_path)
    snp_map = {}
    parse_yfull_node_map(yfull_data, "Root", snp_map)
    print(f"Loaded {len(snp_map)} SNPs from YFull Tree.", flush=True)
    
    # 2. Update Tree JSON (Yleaf format)
    # Yleaf expects specific root name sometimes? Or just valid graph.
    # The original tree starts with "ROOT (Y-Chromosome "Adam")". 
    # Let's map YFull's root (which might be "Root") to that, or just use YFull's structure.
    # Yleaf code loads tree.json and uses keys.
    
    yleaf_tree_dict = {}
    build_yleaf_tree_dict(yfull_data, yleaf_tree_dict)
    
    # Fix Root Name for Yleaf compatibility
    # YFull uses empty string "" or "Root"
    if "" in yleaf_tree_dict:
        print("Renaming empty root node...", flush=True)
        yleaf_tree_dict['ROOT (Y-Chromosome "Adam")'] = yleaf_tree_dict.pop("")
    elif "Root" in yleaf_tree_dict:
        print("Renaming 'Root' node...", flush=True)
        yleaf_tree_dict['ROOT (Y-Chromosome "Adam")'] = yleaf_tree_dict.pop("Root")
    
    # Save tree.json
    tree_out_path = base_dir / "hg_prediction_tables" / "tree.json"
    print(f"Updating tree structure: {tree_out_path}", flush=True)
    with open(tree_out_path, 'w', encoding='utf-8') as f:
        json.dump(yleaf_tree_dict, f, indent=4)
    
    # 3. T2T Update (VCF.gz)
    t2t_vcf_path = temp_dir / "t2t_snps.vcf.gz"
    if not t2t_vcf_path.exists():
        download_file(YBROWSE_T2T_VCF_URL, t2t_vcf_path)
    
    if t2t_vcf_path.exists():
        t2t_out_dir = base_dir / "t2t"
        t2t_out_dir.mkdir(parents=True, exist_ok=True)
        convert_ybrowse_vcf_gz(t2t_vcf_path, snp_map, t2t_out_dir / "new_positions.txt")
    
    # 4. hg38 Update (CSV) — with YBrowse fallback for deeper predictions
    hg38_csv_path = temp_dir / "hg38_snps.csv"
    # Always re-download for fresh data (updated weekly)
    if not hg38_csv_path.exists():
        download_file(YBROWSE_HG38_URL, hg38_csv_path)
    
    if hg38_csv_path.exists():
        hg38_out_dir = base_dir / "hg38"
        hg38_out_dir.mkdir(parents=True, exist_ok=True)
        new_hgs = convert_ybrowse_csv(hg38_csv_path, snp_map, hg38_out_dir / "new_positions.txt")
        
        # 5. Enrich tree.json with new haplogroups from YBrowse
        if new_hgs:
            enriched = 0
            for hg_name in new_hgs:
                if hg_name in yleaf_tree_dict:
                    continue  # already in tree
                
                # Determine parent: base letter (G, J, R, etc.)
                # For "G-Z31275" -> parent is "G"
                # For "R1b" -> parent is "R1"
                # For "R-P312" -> parent is "R"
                base = hg_name[0] if hg_name[0].isalpha() else None
                if not base:
                    continue
                
                # Find best parent in tree
                parent_found = False
                if '-' in hg_name:
                    # Deep haplogroup like G-Z31275
                    base_group = hg_name.split('-')[0]  # "G" from "G-Z31275"
                    # Try exact base group first, then just the letter
                    for candidate in [base_group, base]:
                        if candidate in yleaf_tree_dict:
                            if hg_name not in yleaf_tree_dict[candidate]:
                                yleaf_tree_dict[candidate].append(hg_name)
                            yleaf_tree_dict[hg_name] = []  # leaf node
                            parent_found = True
                            enriched += 1
                            break
                else:
                    # Shallow like "R1b" — try parent "R1", then "R"
                    for i in range(len(hg_name) - 1, 0, -1):
                        candidate = hg_name[:i]
                        if candidate in yleaf_tree_dict:
                            if hg_name not in yleaf_tree_dict[candidate]:
                                yleaf_tree_dict[candidate].append(hg_name)
                            yleaf_tree_dict[hg_name] = []
                            parent_found = True
                            enriched += 1
                            break
                
                if not parent_found:
                    # Add under ROOT as last resort
                    root_key = 'ROOT (Y-Chromosome "Adam")'
                    if root_key in yleaf_tree_dict:
                        if hg_name not in yleaf_tree_dict[root_key]:
                            yleaf_tree_dict[root_key].append(hg_name)
                        yleaf_tree_dict[hg_name] = []
                        enriched += 1
            
            print(f"Enriched tree with {enriched} new haplogroups from YBrowse.", flush=True)
            
            # Re-save enriched tree
            with open(tree_out_path, 'w', encoding='utf-8') as f:
                json.dump(yleaf_tree_dict, f, indent=4)
            print(f"Tree nodes total: {len(yleaf_tree_dict)}", flush=True)
        
    print("Update complete.", flush=True)

if __name__ == "__main__":
    main()

