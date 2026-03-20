# Yleaf FTDNA Deepening Module

Extends [Yleaf](https://github.com/Erasmus-MC/Yleaf) Y-chromosome haplogroup predictions by walking the [FTDNA Y-DNA haplotree](https://www.familytreedna.com/public/y-dna-haplotree) to find deeper sub-branches using WGS variant data.

## What it does

Yleaf predicts haplogroups using ~65K known Y-SNP markers. This module takes that prediction and tries to go deeper by checking thousands of additional SNP positions from the FTDNA tree (768K+ positions). Results include:

- **FTDNA haplogroup** — deepest branch found in the FTDNA tree
- **YFull synonym** — equivalent haplogroup name in the YFull tree
- **Full paths** — complete haplogroup lineage from root to terminal branch (both FTDNA and YFull nomenclature)
- **Reports** — per-sample HTML (with clickable links), TXT, and batch CSV

Typical deepening: **+5 to +15 levels** below Yleaf's prediction.

## Requirements

- **Node.js** ≥ 18
- **bcftools** ≥ 1.15 (for VCF/BAM/CRAM processing)
- **samtools** ≥ 1.15 (for BAM/CRAM index access)
- **Yleaf** (for initial haplogroup prediction) — [Installation guide](https://github.com/Erasmus-MC/Yleaf)

### Reference Genome (for BAM/CRAM only)

BAM and CRAM files require a reference genome. Download hg38:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

Set the path in `deepen_config.json`:
```json
{
  "reference_genome": {
    "hg38": "/path/to/hg38.fa"
  }
}
```

### FTDNA Tree Data

The module requires two pre-built data files:

| File | Source | Description |
|------|--------|-------------|
| `yleaf/data/ftdna_tree.json` | Built from FTDNA API | Tree structure with SNP positions |
| `yleaf/data/yfull_snp_index.json` | Built from YFull tree | SNP→YFull haplogroup mapping |

Generate these with:
```bash
node prepare_ftdna_data.js
```

This fetches from FTDNA public API (`https://www.familytreedna.com/public/y-dna-haplotree/get`) and YFull tree data.

## Usage

### End-to-end (recommended)

Give a VCF or BAM file — the script automatically runs Yleaf, reads the prediction, and deepens:

```bash
# VCF (fastest, ~1 second)
bash deepen.sh /path/to/sample.vcf.gz

# BAM (needs reference, ~16 seconds)
bash deepen.sh /path/to/sample.bam

# CRAM
bash deepen.sh /path/to/sample.cram
```

### Batch mode

Process all samples in a directory:

```bash
# All VCF samples in a directory
bash deepen.sh --batch /path/to/vcf_samples/

# BAM samples
bash deepen.sh --batch /path/to/bam_samples/ --ext bam
```

### Direct deepening (skip Yleaf)

If you already know the haplogroup:

```bash
node deepen_with_ftdna.js -i sample.vcf.gz --hg R-L584 -v
```

### From prediction file

Deepen all predictions in a Yleaf `.hg` output file:

```bash
node deepen_with_ftdna.js -i sample.vcf.gz -p path/to/hg_prediction.hg
```

## Output

### Console

```
Original:      R-L584
Deepened:      R-FT409028 (FTDNA) = R-FT409028 (YFull)
Depth:         +12 (R-L584 > R-FTA62508 > ... > R-FT409028)
FTDNA path:    A-PR2921 > A-L1090 > ... > R-FT409028
YFull path:    A-PR2921 > A0-T > ... > R-FT409028
```

### Report files (saved to `reports/`)

| File | Format | Contents |
|------|--------|----------|
| `<sample>.html` | HTML | Styled report with clickable FTDNA/YFull links |
| `<sample>.txt` | Text | Plain text report |
| `deepening_results.csv` | CSV | Batch summary with all samples |

### CSV columns

```
Sample, Yleaf_Original, Deepened_FTDNA, Deepened_YFull, Depth_Gain, Chain, FTDNA_Path, YFull_Path
```

## Configuration

Edit `deepen_config.json`:

```json
{
  "reference_genome": {
    "hg38": "/path/to/hg38.fa",
    "hg19": null,
    "t2t": "/path/to/chm13v2.0.fa",
    "download_urls": {
      "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
      "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
      "t2t": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    }
  },
  "ftdna_data": {
    "ftdna_tree": "yleaf/data/ftdna_tree.json",
    "yfull_snp_index": "yleaf/data/yfull_snp_index.json"
  },
  "output_dir": "../reports"
}
```

## How it works

1. **Load FTDNA tree** (~768K positioned SNPs across ~135K branches)
2. **Extract variants** from input file:
   - VCF: single `bcftools query` call → in-memory Map (~0.1s)
   - BAM/CRAM: targeted `bcftools mpileup` using BED file of FTDNA positions (~16s)
3. **Walk the tree** from the Yleaf prediction downward:
   - At each node, check all children for derived SNPs
   - Move to the child with the highest derived ratio (≥50%)
   - Repeat until no deeper branch is confirmed
4. **Map to YFull** — cross-reference each FTDNA node with YFull nomenclature
5. **Generate reports** — HTML with clickable links, TXT, CSV

## Performance

| Input | Variants loaded | Time | Deepening |
|-------|-----------------|------|-----------|
| VCF (DeepVariant) | ~30K | **0.1s** | In-memory |
| BAM (30x WGS) | ~750K | **16s** | In-memory |

## Files

| File | Description |
|------|-------------|
| `deepen_with_ftdna.js` | Core deepening engine |
| `deepen.sh` | End-to-end wrapper (Yleaf → deepen → report) |
| `deepen_config.json` | Configuration (reference paths, output settings) |
| `prepare_ftdna_data.js` | Data preparation from FTDNA/YFull APIs |

## Credits

- [Yleaf](https://github.com/Erasmus-MC/Yleaf) — Erasmus MC, Y-chromosome haplogroup prediction
- [FTDNA Y-DNA Haplotree](https://www.familytreedna.com/public/y-dna-haplotree) — FamilyTreeDNA
- [YFull YTree](https://www.yfull.com/tree/) — YFull

## License

Same as the parent Yleaf project (see `LICENSE.txt`).
