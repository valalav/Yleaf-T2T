import os
import pandas as pd
from pathlib import Path
import argparse
import sys

def generate_html(output_folder):
    output_path = Path(output_folder)
    sample_name = None
    
    # 1. Find Prediction File
    hg_file = output_path / "hg_prediction.hg"
    prediction_info = {}
    if hg_file.exists():
        with open(hg_file, 'r') as f:
            header = f.readline().strip().split('\t')
            data = f.readline().strip().split('\t')
            if len(data) >= 2:
                prediction_info = dict(zip(header, data))
                sample_name = prediction_info.get("Sample_name", "Sample")

    if not sample_name:
        print("Could not find prediction file or sample name.")
        return

    # 2. Find Data Files
    sample_dir = output_path / sample_name
    out_file = sample_dir / f"{sample_name}.out"
    fmf_file = sample_dir / f"{sample_name}.fmf"
    
    if not out_file.exists():
        print(f"Data file not found: {out_file}")
        return

    # 3. Load Data
    try:
        df_out = pd.read_csv(out_file, sep="\t")
        df_fmf = pd.read_csv(fmf_file, sep="\t") if fmf_file.exists() else pd.DataFrame()
    except Exception as e:
        print(f"Error reading data files: {e}")
        return

    # 4. Process Data
    # Derived (Positive)
    positive_snps = df_out[df_out['state'] == 'D'].copy()
    # Ancestral (Negative)
    negative_snps = df_out[df_out['state'] == 'A'].copy()
    
    # Missing / Low Quality
    missing_snps = df_fmf.copy()
    
    # Add YFull Links
    def make_link(hg):
        if pd.isna(hg) or hg == "-" or str(hg).strip() == "": return hg
        clean_hg = str(hg).strip()
        # Extract SNP/Haplogroup for link (take last part after dash)
        link_target = clean_hg.split("-")[-1].replace("*", "")
        return f'<a href="https://www.yfull.com/tree/{link_target}/" target="_blank" class="text-decoration-none">{clean_hg}</a>'

    def format_full_haplogroup(hg_str):
        if not hg_str or pd.isna(hg_str):
            return "Unknown"
        
        # Check for exclusions pattern like "Main(xExcl1,Excl2)"
        if "(x" in hg_str and hg_str.endswith(")"):
            try:
                main_part, exclusions_part = hg_str.split("(x", 1)
                exclusions_part = exclusions_part.rstrip(")")
                excluded_branches = exclusions_part.split(",")
            except ValueError:
                # Fallback if split fails
                return make_link(hg_str)
        else:
            main_part = hg_str
            excluded_branches = []

        # Format Main Branch
        main_clean = main_part.strip()
        main_target = main_clean.split("-")[-1].replace("*", "")
        main_html = f'<a href="https://www.yfull.com/tree/{main_target}/" target="_blank" class="fs-4 fw-bold text-decoration-none">{main_clean}</a>'

        # Format Excluded Branches
        excl_html = ""
        if excluded_branches:
            excl_html = """
            <div class="mt-2">
                <small class="text-muted text-uppercase fw-bold" style="font-size: 0.75rem;">Excluded Downstream Branches (Not Found):</small>
                <div class="d-flex flex-wrap gap-1 mt-1">
            """
            for branch in excluded_branches:
                b_clean = branch.strip()
                b_target = b_clean.split("-")[-1].replace("*", "")
                excl_html += f'<a href="https://www.yfull.com/tree/{b_target}/" target="_blank" class="badge bg-light text-secondary border text-decoration-none fw-normal" title="Open {b_clean} on YFull">{b_clean}</a>'
            excl_html += "</div></div>"
        
        return main_html + excl_html

    # 5. HTML Template Construction
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Yleaf Report: {sample_name}</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css">
        <style>
            body {{ padding: 20px; background-color: #f8f9fa; }}
            .card {{ margin-bottom: 20px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border: none; }}
            .card-header {{ background: linear-gradient(45deg, #0d6efd, #0a58ca); color: white; }}
            .hg-title {{ color: #0d6efd; font-weight: bold; }}
            .table-container {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
            .badge-pos {{ background-color: #198754; }}
            .badge-neg {{ background-color: #dc3545; }}
            .badge-miss {{ background-color: #6c757d; }}
            .nav-tabs .nav-link.active {{ font-weight: bold; border-bottom: 3px solid #0d6efd; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="card">
                <div class="card-header">
                    <h2 class="mb-0 fs-3">🧬 Yleaf Analysis Report</h2>
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-md-7">
                            <h5 class="text-muted mb-1">Sample</h5>
                            <p class="fs-5 fw-bold mb-3">{sample_name}</p>
                            
                            <h5 class="text-muted mb-1">Predicted Haplogroup</h5>
                            <div class="mb-3">
                                {format_full_haplogroup(prediction_info.get('Hg', 'Unknown'))}
                            </div>
                            
                            <h5 class="text-muted mb-1">Terminal Marker</h5>
                            <p class="fs-5">{prediction_info.get('Hg_marker', 'N/A')}</p>
                        </div>
                        <div class="col-md-5">
                            <div class="card bg-light border-0">
                                <div class="card-body">
                                    <h6 class="card-title text-center mb-3">Quality Metrics</h6>
                                    <ul class="list-group list-group-flush bg-transparent">
                                        <li class="list-group-item d-flex justify-content-between align-items-center bg-transparent">
                                            QC Score
                                            <span class="badge bg-primary rounded-pill" style="font-size: 1rem;">{prediction_info.get('QC-score', 'N/A')}</span>
                                        </li>
                                        <li class="list-group-item d-flex justify-content-between align-items-center bg-transparent">
                                            Total Reads
                                            <span class="badge bg-secondary rounded-pill">{prediction_info.get('Total_reads', '0')}</span>
                                        </li>
                                        <li class="list-group-item d-flex justify-content-between align-items-center bg-transparent">
                                            Valid Markers
                                            <span class="badge bg-secondary rounded-pill">{prediction_info.get('Valid_markers', '0')}</span>
                                        </li>
                                    </ul>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <ul class="nav nav-tabs mb-3" id="myTab" role="tablist">
                <li class="nav-item" role="presentation">
                    <button class="nav-link active" id="pos-tab" data-bs-toggle="tab" data-bs-target="#pos" type="button" role="tab">
                        ✅ Positive SNPs <span class="badge badge-pos ms-2">{len(positive_snps)}</span>
                    </button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="neg-tab" data-bs-toggle="tab" data-bs-target="#neg" type="button" role="tab">
                        ❌ Negative SNPs <span class="badge badge-neg ms-2">{len(negative_snps)}</span>
                    </button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="miss-tab" data-bs-toggle="tab" data-bs-target="#miss" type="button" role="tab">
                        ⚠️ Missing/Filtered <span class="badge badge-miss ms-2">{len(missing_snps)}</span>
                    </button>
                </li>
            </ul>

            <div class="tab-content table-container mt-3" id="myTabContent">
                <!-- Positive SNPs -->
                <div class="tab-pane fade show active" id="pos" role="tabpanel">
                    {generate_table_html(positive_snps, "posTable", make_link)}
                </div>
                
                <!-- Negative SNPs -->
                <div class="tab-pane fade" id="neg" role="tabpanel">
                    {generate_table_html(negative_snps, "negTable", make_link)}
                </div>
                
                <!-- Missing SNPs -->
                <div class="tab-pane fade" id="miss" role="tabpanel">
                    {generate_table_html(missing_snps, "missTable", make_link, is_missing=True)}
                </div>
            </div>
            
            <footer class="mt-5 text-center text-muted">
                <small>Generated by Yleaf Report Module</small>
            </footer>
        </div>

        <script src="https://code.jquery.com/jquery-3.5.1.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.4/js/dataTables.bootstrap5.min.js"></script>
        <script>
            $(document).ready(function () {{
                $('#posTable').DataTable({{ "order": [[ 3, "desc" ]] }}); // Order by reads (col index 3)
                $('#negTable').DataTable({{ "order": [[ 0, "asc" ]] }});
                $('#missTable').DataTable({{ "order": [[ 0, "asc" ]] }});
            }});
        </script>
    </body>
    </html>
    """

    report_path = output_path / "report.html"
    with open(report_path, "w") as f:
        f.write(html_content)
    
    print(f"Report generated successfully: {report_path}")

def generate_table_html(df, table_id, link_func, is_missing=False):
    if df.empty:
        return "<p>No data available.</p>"
    
    # Select columns
    cols = ['marker_name', 'haplogroup', 'mutation', 'reads', 'called_perc', 'pos']
    if is_missing:
        cols.append('Description')
    
    # Apply links
    df_display = df.copy()
    if 'haplogroup' in df_display.columns:
        df_display['haplogroup'] = df_display['haplogroup'].apply(link_func)
    
    # Generate HTML
    table_html = f'<table id="{table_id}" class="table table-striped table-hover" style="width:100%">'
    table_html += "<thead><tr>"
    for col in cols:
        table_html += f"<th>{col}</th>"
    table_html += "</tr></thead><tbody>"
    
    for _, row in df_display.iterrows():
        table_html += "<tr>"
        for col in cols:
            val = row.get(col, '')
            table_html += f"<td>{val}</td>"
        table_html += "</tr>"
    
    table_html += "</tbody></table>"
    return table_html

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate HTML report for Yleaf.')
    parser.add_argument('output_folder', type=str, help='Path to Yleaf output folder')
    args = parser.parse_args()
    
    generate_html(args.output_folder)
