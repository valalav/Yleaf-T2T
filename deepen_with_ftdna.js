#!/usr/bin/env node
/**
 * deepen_with_ftdna.js — Standalone FTDNA deepening for Yleaf predictions.
 * 
 * Optimized: extracts ALL chrY variants in ONE call, then walks tree in-memory.
 * Supports VCF(.gz), BAM, and CRAM input files.
 * 
 * Usage:
 *   node deepen_with_ftdna.js --input sample.vcf.gz --hg J-Z39973 [-v]
 *   node deepen_with_ftdna.js --input sample.bam --hg R-L584 [-v] [--ref ref.fa]
 *   node deepen_with_ftdna.js --input sample.vcf.gz -p hg_prediction.hg [-o out.hg]
 * 
 * Data preparation (one-time): node prepare_ftdna_data.js
 * Requires: bcftools (for VCF), samtools+bcftools (for BAM/CRAM)
 */

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// --- Configuration ---
const DATA_DIR = path.resolve(__dirname, 'yleaf/data');
const FTDNA_TREE_FILE = path.join(DATA_DIR, 'ftdna_tree.json');
const YFULL_INDEX_FILE = path.join(DATA_DIR, 'yfull_snp_index.json');
const MIN_DERIVED_RATIO = 0.5;

// --- Parse arguments ---
function parseArgs() {
    const args = { input: null, hg: null, prediction: null, output: null, ref: null, verbose: false };
    const argv = process.argv.slice(2);
    for (let i = 0; i < argv.length; i++) {
        switch (argv[i]) {
            case '--input': case '--vcf': case '-i': args.input = argv[++i]; break;
            case '--hg': args.hg = argv[++i]; break;
            case '--prediction': case '-p': args.prediction = argv[++i]; break;
            case '-o': case '--output': args.output = argv[++i]; break;
            case '--ref': case '-r': args.ref = argv[++i]; break;
            case '-v': case '--verbose': args.verbose = true; break;
            case '--help': case '-h':
                console.log(`Yleaf FTDNA Deepening — extend haplogroup predictions\n`);
                console.log(`Usage:`);
                console.log(`  node deepen_with_ftdna.js -i <vcf.gz|bam|cram> --hg <haplogroup> [-v]`);
                console.log(`  node deepen_with_ftdna.js -i <vcf.gz|bam|cram> -p <prediction.hg> [-o <output.hg>]`);
                console.log(`  --ref <ref.fa>  Reference genome (required for BAM/CRAM)\n`);
                console.log(`Data prep: node prepare_ftdna_data.js`);
                process.exit(0);
        }
    }
    if (!args.input) { console.error('ERROR: --input required (VCF, BAM, or CRAM)'); process.exit(1); }
    if (!args.hg && !args.prediction) { console.error('ERROR: --hg or -p required'); process.exit(1); }
    return args;
}

// --- Detect input file type ---
function detectInputType(filepath) {
    const ext = filepath.toLowerCase();
    if (ext.endsWith('.vcf') || ext.endsWith('.vcf.gz') || ext.endsWith('.bcf')) return 'vcf';
    if (ext.endsWith('.bam')) return 'bam';
    if (ext.endsWith('.cram')) return 'cram';
    return 'vcf'; // default
}

// --- Load data ---
function loadData() {
    if (!fs.existsSync(FTDNA_TREE_FILE)) {
        console.error(`ERROR: ${FTDNA_TREE_FILE} not found. Run: node prepare_ftdna_data.js`);
        process.exit(1);
    }
    const t0 = Date.now();
    console.error('Loading FTDNA tree...');
    const tree = JSON.parse(fs.readFileSync(FTDNA_TREE_FILE, 'utf8'));

    const nameIndex = {};
    const variantIndex = {};
    const parentIndex = {};
    for (const [nid, node] of Object.entries(tree)) {
        nameIndex[node.n] = nid;
        for (const v of node.v) {
            variantIndex[v.s.toUpperCase()] = nid;
        }
        for (const childId of (node.c || [])) {
            parentIndex[String(childId)] = nid;
        }
    }

    let yfullIndex = {};
    if (fs.existsSync(YFULL_INDEX_FILE)) {
        console.error('Loading YFull SNP index...');
        yfullIndex = JSON.parse(fs.readFileSync(YFULL_INDEX_FILE, 'utf8'));
    }

    console.error(`  Trees ready in ${((Date.now() - t0) / 1000).toFixed(1)}s`);
    return { tree, nameIndex, variantIndex, yfullIndex, parentIndex };
}

// --- Extract chrY variants from input file into memory ---
function extractChrYVariants(inputPath, inputType, refPath, data, verbose) {
    const t0 = Date.now();
    const variants = new Map();  // position → { ref, alt, gt }

    if (inputType === 'vcf') {
        // VCF: single bcftools query — all chrY variants
        if (verbose) console.error(`  Extracting chrY variants from VCF...`);
        try {
            const cmd = `bcftools query -r chrY -f '%POS\\t%REF\\t%ALT\\t[%GT]\\n' ${inputPath} 2>/dev/null`;
            const output = execSync(cmd, { encoding: 'utf8', maxBuffer: 100 * 1024 * 1024 });
            for (const line of output.split('\n')) {
                if (!line) continue;
                const parts = line.split('\t');
                const pos = parseInt(parts[0]);
                if (!isNaN(pos)) variants.set(pos, { ref: parts[1] || '', alt: parts[2] || '.', gt: parts[3] || '.' });
            }
        } catch (e) {
            console.error(`  WARNING: VCF extraction failed: ${e.message.split('\n')[0]}`);
        }
    } else {
        // BAM/CRAM: targeted mpileup — only FTDNA positions (fast!)
        if (!refPath) {
            console.error('ERROR: --ref <reference.fa> required for BAM/CRAM input');
            process.exit(1);
        }
        if (verbose) console.error(`  Building targeted positions from FTDNA tree...`);

        // Collect ALL unique positions from FTDNA tree
        const positions = new Set();
        for (const node of Object.values(data.tree)) {
            for (const v of node.v) {
                if (v.p > 0) positions.add(v.p);
            }
        }

        // Write positions to temp BED file
        const tmpBed = `/tmp/ftdna_positions_${process.pid}.bed`;
        const bedLines = [...positions].sort((a, b) => a - b).map(p => `chrY\t${p - 1}\t${p}`);
        fs.writeFileSync(tmpBed, bedLines.join('\n') + '\n');
        if (verbose) console.error(`  ${positions.size} target positions → ${tmpBed}`);

        try {
            // File-based pipeline (pipes fail on large datasets in execSync)
            if (verbose) console.error(`  Running targeted mpileup on BAM (~30s)...`);
            const tmpVcf = `/tmp/ftdna_mpileup_${process.pid}.vcf`;
            const tmpCalled = `/tmp/ftdna_called_${process.pid}.vcf`;

            // Step 1: mpileup → raw VCF
            execSync(`bcftools mpileup -r chrY -T ${tmpBed} -f ${refPath} --no-BAQ -d 500 ${inputPath} -o ${tmpVcf} 2>/dev/null`, { timeout: 300000 });
            
            // Step 2: call variants
            execSync(`bcftools call -m ${tmpVcf} -o ${tmpCalled} 2>/dev/null`, { timeout: 120000 });
            
            // Step 3: query genotypes
            const output = execSync(`bcftools query -f '%POS\\t%REF\\t%ALT\\t[%GT]\\n' ${tmpCalled} 2>/dev/null`, { encoding: 'utf8', maxBuffer: 100 * 1024 * 1024 });
            for (const line of output.split('\n')) {
                if (!line) continue;
                const parts = line.split('\t');
                const pos = parseInt(parts[0]);
                if (!isNaN(pos)) variants.set(pos, { ref: parts[1] || '', alt: parts[2] || '.', gt: parts[3] || '.' });
            }
            
            // Cleanup temp files
            try { fs.unlinkSync(tmpVcf); } catch (e) {}
            try { fs.unlinkSync(tmpCalled); } catch (e) {}
        } catch (e) {
            console.error(`  WARNING: BAM extraction failed: ${e.message.split('\n')[0]}`);
        }
        // Cleanup
        try { fs.unlinkSync(tmpBed); } catch (e) {}
    }

    const elapsed = ((Date.now() - t0) / 1000).toFixed(1);
    console.error(`  ${variants.size} chrY variants loaded in ${elapsed}s`);
    return variants;
}

// --- Check node derived status using in-memory variants ---
function checkNodeDerived(node, variants, verbose) {
    const snps = node.v.filter(v => v.p > 0);
    if (snps.length === 0) return { derived: 0, total: 0, ratio: 0, checked: false };

    let derived = 0, ancestral = 0, missing = 0;
    for (const v of snps) {
        const entry = variants.get(v.p);
        if (!entry) { missing++; continue; }

        const isDerived = entry.alt === v.d && entry.gt !== '0/0' && entry.gt !== '0|0';
        const isAncestral = (entry.gt === '0/0' || entry.gt === '0|0' || entry.alt === '.');

        if (isDerived) derived++;
        else if (isAncestral) ancestral++;

        if (verbose) {
            console.error(`    ${v.s} (${v.p}): ${isDerived ? 'DERIVED ✓' : isAncestral ? 'ancestral' : `? alt=${entry.alt} gt=${entry.gt}`}`);
        }
    }

    const checked = derived + ancestral;
    return { derived, ancestral, missing, total: snps.length, checked: checked > 0, ratio: checked > 0 ? derived / checked : 0 };
}

// --- Find node ---
function findNode(data, hg) {
    let nid = data.nameIndex[hg];
    if (nid) return nid;
    const snp = hg.includes('-') ? hg.split('-').slice(1).join('-') : hg;
    return data.variantIndex[snp.toUpperCase()] || null;
}

// --- Branch prefix ---
function getBranchPrefix(name) {
    return name ? name.split('-')[0].replace(/[0-9]/g, '') : '';
}

// --- Map to YFull with branch validation ---
function toYFull(data, ftdnaName) {
    const prefix = getBranchPrefix(ftdnaName);
    const snp = ftdnaName.includes('-') ? ftdnaName.split('-').slice(1).join('-') : ftdnaName;
    const yid = data.yfullIndex[snp.toUpperCase()];
    if (yid && getBranchPrefix(yid) === prefix) return yid;

    const nid = findNode(data, ftdnaName);
    if (nid) {
        for (const v of data.tree[nid].v) {
            const y = data.yfullIndex[v.s.toUpperCase()];
            if (y && getBranchPrefix(y) === prefix) return y;
        }
    }
    return null;
}

// --- Build full FTDNA path ---
function buildPath(data, nodeId) {
    const p = [];
    let cur = nodeId;
    while (cur) {
        p.unshift(data.tree[cur].n);
        cur = data.parentIndex[cur] || null;
    }
    return p;
}

// --- Deepen prediction (in-memory, zero I/O) ---
function deepen(data, startHg, variants, verbose) {
    const nid = findNode(data, startHg);
    if (!nid) {
        if (verbose) console.error(`  ${startHg}: not found in FTDNA tree`);
        return { original: startHg, deepened: null, yfull: null, deepPath: [], fullPath: [], depthGain: 0 };
    }

    let cur = nid;
    const deepPath = [data.tree[cur].n];
    if (verbose) console.error(`\nDeepening ${startHg} (FTDNA: ${data.tree[cur].n})`);

    let maxDepth = 30;
    while (maxDepth-- > 0) {
        const children = data.tree[cur].c || [];
        if (children.length === 0) break;
        if (verbose) console.error(`  ${children.length} children of ${data.tree[cur].n}...`);

        let bestChild = null, bestScore = { derived: 0 };
        for (const cid of children) {
            const child = data.tree[String(cid)];
            if (!child) continue;
            const score = checkNodeDerived(child, variants, verbose);
            if (verbose) console.error(`  ${child.n}: ${score.derived}/${score.total} (${(score.ratio*100).toFixed(0)}%)`);
            if (score.checked && score.ratio >= MIN_DERIVED_RATIO && score.derived > bestScore.derived) {
                bestChild = String(cid);
                bestScore = score;
            }
        }

        if (!bestChild) break;
        cur = bestChild;
        deepPath.push(data.tree[cur].n);
        if (verbose) console.error(`  → ${data.tree[cur].n} (${bestScore.derived}/${bestScore.total} ✓)`);
    }

    const deepenedName = data.tree[cur].n;

    // Build YFull path: map each FTDNA node in full path to YFull synonym
    const fullPath = buildPath(data, cur);
    const yfullPath = fullPath.map(ftdnaName => toYFull(data, ftdnaName) || ftdnaName);

    return {
        original: startHg,
        deepened: deepenedName,
        yfull: toYFull(data, deepenedName),
        deepPath,
        fullPath,
        yfullPath,
        depthGain: deepPath.length - 1
    };
}

// --- Load config ---
function loadConfig() {
    const configPath = path.join(__dirname, 'deepen_config.json');
    if (fs.existsSync(configPath)) return JSON.parse(fs.readFileSync(configPath, 'utf8'));
    return {};
}

// --- Report: TXT ---
function generateTXT(result, sampleName) {
    const ts = new Date().toISOString().replace('T', ' ').split('.')[0];
    let t = `FTDNA Deepening Report\nGenerated: ${ts}\nSample: ${sampleName}\n${'='.repeat(60)}\n\n`;
    t += `Original (Yleaf):  ${result.original}\n`;
    if (result.depthGain > 0) {
        t += `Deepened (FTDNA):   ${result.deepened}\n`;
        t += `Deepened (YFull):   ${result.yfull || 'N/A'}\n`;
        t += `Depth gain:         +${result.depthGain}\n`;
        t += `Chain:              ${result.deepPath.join(' > ')}\n\n`;
    } else {
        t += `Deepened:           — (at deepest known branch)\n\n`;
    }
    t += `FTDNA Path\n${'-'.repeat(40)}\n${result.fullPath.join(' > ')}\n\n`;
    t += `YFull Path\n${'-'.repeat(40)}\n${result.yfullPath.join(' > ')}\n`;
    return t;
}

// --- Report: HTML ---
function generateHTML(result, sampleName) {
    const ts = new Date().toISOString().replace('T', ' ').split('.')[0];
    const deep = result.depthGain > 0;
    const ftEnd = result.deepened || result.original;
    const yfEnd = result.yfull || result.original;
    const ftLink = `https://www.familytreedna.com/public/y-dna-haplotree/${ftEnd}`;
    const yfLink = `https://www.yfull.com/tree/${yfEnd}/`;
    const pathHTML = (arr, tree) => {
        const base = tree === 'ftdna' ? 'https://www.familytreedna.com/public/y-dna-haplotree/' : 'https://www.yfull.com/tree/';
        const suf = tree === 'yfull' ? '/' : '';
        return arr.map(n => `<a href="${base}${n}${suf}" target="_blank">${n}</a>`).join(' &gt; ');
    };
    return `<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>Deepening — ${sampleName}</title>
<style>
body{font-family:'Segoe UI',system-ui,sans-serif;max-width:900px;margin:2em auto;padding:0 1em;background:#f8f9fa;color:#222}
h1{color:#2c5282;border-bottom:3px solid #2c5282;padding-bottom:.3em}
.meta{color:#666;font-size:.9em;margin-bottom:1.5em}
.card{background:#fff;border-radius:8px;padding:1.2em 1.5em;margin:1em 0;box-shadow:0 1px 4px rgba(0,0,0,.1)}
.card h2{margin-top:0;color:#2d3748;font-size:1.1em}
.label{font-weight:600;color:#4a5568;min-width:140px;display:inline-block}
.gain{color:#38a169;font-weight:bold;font-size:1.1em}
.no-gain{color:#a0aec0}
.path{word-break:break-all;line-height:1.8;font-size:.85em}
.path a{color:#3182ce;text-decoration:none}
.path a:hover{text-decoration:underline}
table{border-collapse:collapse;width:100%}
td{padding:.4em .6em;vertical-align:top}
td:first-child{width:150px}
</style></head><body>
<h1>🧬 FTDNA Deepening Report</h1>
<div class="meta">Sample: <b>${sampleName}</b> · ${ts}</div>
<div class="card"><h2>Prediction</h2><table>
<tr><td class="label">Yleaf (original)</td><td>${result.original}</td></tr>
<tr><td class="label">FTDNA (deepened)</td><td>${deep ? `<a href="${ftLink}" target="_blank">${result.deepened}</a>` : '<span class="no-gain">—</span>'}</td></tr>
<tr><td class="label">YFull (synonym)</td><td>${deep ? `<a href="${yfLink}" target="_blank">${result.yfull || 'N/A'}</a>` : '<span class="no-gain">—</span>'}</td></tr>
<tr><td class="label">Depth gain</td><td>${deep ? `<span class="gain">+${result.depthGain}</span>` : '<span class="no-gain">0</span>'}</td></tr>
${deep ? `<tr><td class="label">Chain</td><td>${result.deepPath.join(' → ')}</td></tr>` : ''}
</table></div>
<div class="card"><h2>FTDNA Path</h2><div class="path">${pathHTML(result.fullPath, 'ftdna')}</div></div>
<div class="card"><h2>YFull Path</h2><div class="path">${pathHTML(result.yfullPath, 'yfull')}</div></div>
<div class="meta" style="margin-top:2em;text-align:center">
<a href="https://github.com/Erasmus-MC/Yleaf">Yleaf</a> ·
<a href="https://www.familytreedna.com/public/y-dna-haplotree">FTDNA</a> ·
<a href="https://www.yfull.com/tree/">YFull</a>
</div></body></html>`;
}

// --- CSV helpers ---
const CSV_HEADER = 'Sample,Yleaf_Original,Deepened_FTDNA,Deepened_YFull,Depth_Gain,Chain,FTDNA_Path,YFull_Path';
function toCSV(r, name) {
    const e = s => `"${(s||'').replace(/"/g, '""')}"`;
    return [e(name), e(r.original), e(r.deepened||''), e(r.yfull||''), r.depthGain||0,
        e(r.deepPath.join(' > ')), e(r.fullPath.join(' > ')), e(r.yfullPath.join(' > '))].join(',');
}

// --- Save reports ---
function saveReports(r, name, dir) {
    fs.mkdirSync(dir, { recursive: true });
    const base = path.join(dir, name);
    fs.writeFileSync(`${base}.txt`, generateTXT(r, name));
    fs.writeFileSync(`${base}.html`, generateHTML(r, name));
    console.error(`  Reports: ${base}.{txt,html}`);
    
    const csvFile = path.join(dir, 'deepening_results.csv');
    const needH = !fs.existsSync(csvFile);
    fs.appendFileSync(csvFile, (needH ? CSV_HEADER + '\n' : '') + toCSV(r, name) + '\n');
    console.error(`  CSV: ${csvFile}`);
}

// --- Main ---
function main() {
    const args = parseArgs();
    const config = loadConfig();
    const data = loadData();
    const inputType = detectInputType(args.input);
    
    if (!args.ref && (inputType === 'bam' || inputType === 'cram') && config.reference_genome) {
        args.ref = config.reference_genome.hg38;
        if (args.ref) console.error(`  Ref from config: ${args.ref}`);
    }

    const variants = extractChrYVariants(args.input, inputType, args.ref, data, args.verbose);
    const outDir = args.output || path.resolve(__dirname, config.output_dir || '../reports');
    const sample = path.basename(args.input).replace(/\.(DeepVariant_v1\.6\.vcf\.gz|vcf\.gz|vcf|bam|cram)$/i, '');

    if (args.hg) {
        const r = deepen(data, args.hg, variants, args.verbose);
        if (r.depthGain > 0) {
            console.log(`Original:      ${r.original}`);
            console.log(`Deepened:      ${r.deepened} (FTDNA) = ${r.yfull || 'N/A'} (YFull)`);
            console.log(`Depth:         +${r.depthGain} (${r.deepPath.join(' > ')})`);
            console.log(`FTDNA path:    ${r.fullPath.join(' > ')}`);
            console.log(`YFull path:    ${r.yfullPath.join(' > ')}`);
        } else {
            console.log(`${r.original}: no deepening found`);
            if (r.fullPath.length) {
                console.log(`FTDNA path:    ${r.fullPath.join(' > ')}`);
                console.log(`YFull path:    ${r.yfullPath.join(' > ')}`);
            }
        }
        saveReports(r, sample, outDir);
        return;
    }

    if (args.prediction) {
        const lines = fs.readFileSync(args.prediction, 'utf8').trim().split('\n');
        const csvFile = path.join(outDir, 'deepening_results.csv');
        fs.mkdirSync(outDir, { recursive: true });
        if (fs.existsSync(csvFile)) fs.unlinkSync(csvFile);

        for (let i = 1; i < lines.length; i++) {
            const fields = lines[i].split('\t');
            const hg = (fields[1] || '').split('*')[0].trim();
            if (!hg || hg === 'NA' || hg.startsWith('Low_Y_Signal')) continue;
            const r = deepen(data, hg, variants, args.verbose);
            const sn = fields[0] || `sample_${i}`;
            fs.appendFileSync(csvFile, (i === 1 ? CSV_HEADER + '\n' : '') + toCSV(r, sn) + '\n');
            if (r.depthGain > 0) {
                saveReports(r, sn, outDir);
                console.error(`  ${sn}: ${hg} → ${r.deepened} [+${r.depthGain}]`);
            } else {
                console.error(`  ${sn}: ${hg} → no deepening`);
            }
        }
        console.error(`\nBatch CSV: ${csvFile}`);
    }
}

main();
