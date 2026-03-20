#!/usr/bin/env node
/**
 * prepare_ftdna_data.js — Pre-process FTDNA get.json into lightweight format for Yleaf.
 * 
 * Strips unnecessary fields (countries, statistics, kitsCount, etc.)
 * Keeps only: name, children, variants[{variant, position, ancestral, derived}]
 * Also builds YFull SNP→canonical ID index from ytree.json
 * 
 * Usage:
 *   node prepare_ftdna_data.js
 * 
 * Output:
 *   yleaf/data/ftdna_tree.json   (~15-20MB instead of 115MB)
 *   yleaf/data/yfull_snp_index.json  (~2-5MB)
 */

const fs = require('fs');
const path = require('path');
const os = require('os');

const FTDNA_SRC = path.resolve(os.homedir(), '_dna/ystr-matcher/ftdna_haplo/data/get.json');
const YFULL_SRC = path.resolve(os.homedir(), '_dna/ystr-matcher/ftdna_haplo/data/ytree.json');
const OUTPUT_DIR = path.resolve(__dirname, 'yleaf/data');

// --- FTDNA tree ---
function processFTDNA() {
    console.log(`Loading FTDNA from ${FTDNA_SRC}...`);
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
