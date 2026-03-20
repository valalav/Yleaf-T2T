#!/usr/bin/env node
/**
 * prepare_ftdna_data.js — Pre-process FTDNA get.json into lightweight format for Yleaf.
 * 
 * Strips unnecessary fields (countries, statistics, kitsCount, etc.)
 * Keeps only: name, children, variants[{variant, position, ancestral, derived}]
 * Also builds YFull SNP→canonical ID index from ytree.json
 * 
 * Usage:
 *   node prepare_ftdna_data.js --ftdna /path/to/get.json --yfull /path/to/ytree.json
 *   node prepare_ftdna_data.js   (reads paths from deepen_config.json)
 * 
 * Source data:
 *   FTDNA: curl -o get.json "https://www.familytreedna.com/public/y-dna-haplotree/get"
 *   YFull: from FTDNA Haplo Server data/ or YFull website
 * 
 * Output:
 *   yleaf/data/ftdna_tree.json   (~40MB)
 *   yleaf/data/yfull_snp_index.json  (~10MB)
 */

const fs = require('fs');
const path = require('path');

// --- Parse arguments ---
function getSourcePaths() {
    const args = process.argv.slice(2);
    let ftdna = null, yfull = null;
    
    for (let i = 0; i < args.length; i++) {
        if (args[i] === '--ftdna' && args[i+1]) ftdna = args[++i];
        if (args[i] === '--yfull' && args[i+1]) yfull = args[++i];
    }
    
    // Fallback to config
    if (!ftdna || !yfull) {
        const configPath = path.join(__dirname, 'deepen_config.json');
        if (fs.existsSync(configPath)) {
            const config = JSON.parse(fs.readFileSync(configPath, 'utf8'));
            if (!ftdna && config.ftdna_data && config.ftdna_data.get_json) ftdna = config.ftdna_data.get_json;
            if (!yfull && config.ftdna_data && config.ftdna_data.ytree_json) yfull = config.ftdna_data.ytree_json;
        }
    }
    
    if (!ftdna || !yfull) {
        console.error(`Usage: node prepare_ftdna_data.js --ftdna /path/to/get.json --yfull /path/to/ytree.json`);
        console.error(`\nOr set paths in deepen_config.json under ftdna_data.get_json and ftdna_data.ytree_json`);
        console.error(`\nTo download source data:`);
        console.error(`  FTDNA: curl -o get.json "https://www.familytreedna.com/public/y-dna-haplotree/get"`);
        console.error(`  YFull: available from FTDNA Haplo Server or YFull website`);
        process.exit(1);
    }
    
    return { ftdna, yfull };
}

const OUTPUT_DIR = path.resolve(__dirname, 'yleaf/data');
const { ftdna: FTDNA_SRC, yfull: YFULL_SRC } = getSourcePaths();

// --- FTDNA tree ---
function processFTDNA() {
    console.log(`Loading FTDNA from ${FTDNA_SRC}...`);
    if (!fs.existsSync(FTDNA_SRC)) {
        console.error(`  ERROR: File not found: ${FTDNA_SRC}`);
        console.error(`  Download: curl -o get.json "https://www.familytreedna.com/public/y-dna-haplotree/get"`);
        process.exit(1);
    }
    const data = JSON.parse(fs.readFileSync(FTDNA_SRC, 'utf8'));
    const nodes = data.allNodes;
    
    const slim = {};
    let totalVariants = 0;
    
    for (const [nid, node] of Object.entries(nodes)) {
        const variants = (node.variants || [])
            .filter(v => v.position && v.position !== 0)
            .map(v => ({
                s: v.variant,                    // snp name
                p: Math.abs(v.position),          // position (absolute)
                a: v.ancestral,                   // ancestral allele
                d: v.derived                      // derived allele
            }));
        
        totalVariants += variants.length;
        
        slim[nid] = {
            n: node.name,                         // name
            c: node.children || [],               // children IDs
            v: variants                           // variants (minimal)
        };
    }
    
    const outPath = path.join(OUTPUT_DIR, 'ftdna_tree.json');
    const json = JSON.stringify(slim);
    fs.writeFileSync(outPath, json);
    
    const sizeMB = (Buffer.byteLength(json) / 1024 / 1024).toFixed(1);
    console.log(`  ✅ ftdna_tree.json: ${Object.keys(slim).length} nodes, ${totalVariants} variants, ${sizeMB}MB`);
    return slim;
}

// --- YFull SNP index ---
function processYFull() {
    console.log(`Loading YFull from ${YFULL_SRC}...`);
    if (!fs.existsSync(YFULL_SRC)) {
        console.error(`  ERROR: File not found: ${YFULL_SRC}`);
        process.exit(1);
    }
    const data = JSON.parse(fs.readFileSync(YFULL_SRC, 'utf8'));
    
    // Build SNP name → canonical YFull ID index
    const snpIndex = {};  // { "Z45300": "J-Y94477", ... }
    let nodeCount = 0;
    
    function walkTree(node) {
        if (!node) return;
        nodeCount++;
        
        const canonicalId = node.id || node.name || '';
        
        // Index all SNP names for this node
        const snps = [];
        if (node.snps) {
            // YFull format: snps can be string or array
            const snpList = typeof node.snps === 'string' ? node.snps.split(',') : (Array.isArray(node.snps) ? node.snps : []);
            for (const snp of snpList) {
                const clean = snp.trim().replace(/\*/g, '').replace(/!/g, '');
                if (clean) {
                    snpIndex[clean.toUpperCase()] = canonicalId;
                    snps.push(clean);
                }
            }
        }
        
        // Also index by node ID itself
        if (canonicalId) {
            const snp = canonicalId.includes('-') ? canonicalId.split('-').slice(1).join('-') : '';
            if (snp) {
                snpIndex[snp.toUpperCase()] = canonicalId;
            }
        }
        
        // Recurse into children
        if (node.children && Array.isArray(node.children)) {
            for (const child of node.children) {
                walkTree(child);
            }
        }
    }
    
    // YFull tree can be a single root or array
    if (Array.isArray(data)) {
        for (const root of data) walkTree(root);
    } else if (data.children) {
        walkTree(data);
    } else {
        // Try known structures
        for (const key of Object.keys(data)) {
            if (typeof data[key] === 'object' && data[key].children) {
                walkTree(data[key]);
            }
        }
    }
    
    const outPath = path.join(OUTPUT_DIR, 'yfull_snp_index.json');
    const json = JSON.stringify(snpIndex);
    fs.writeFileSync(outPath, json);
    
    const sizeMB = (Buffer.byteLength(json) / 1024 / 1024).toFixed(1);
    console.log(`  ✅ yfull_snp_index.json: ${nodeCount} nodes, ${Object.keys(snpIndex).length} SNP mappings, ${sizeMB}MB`);
}

// --- Main ---
console.log(`Output directory: ${OUTPUT_DIR}\n`);

if (!fs.existsSync(OUTPUT_DIR)) {
    fs.mkdirSync(OUTPUT_DIR, { recursive: true });
}

processFTDNA();
processYFull();
console.log('\nDone! Data files ready for standalone deepening.');
