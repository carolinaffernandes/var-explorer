#!/usr/bin/env python3
import os
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from bs4 import BeautifulSoup
from datetime import datetime

# --- For the substitution graph:

VALID_SUBSTITUTIONS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

COSMIC_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

COSMIC_COLORS = {
    "C>A": "#03BCEE",
    "C>G": "#010101",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
}


def load_full_config(path="config.yaml"):
    if not os.path.exists(path):
        print(f"[CRITICAL ERROR] File {path} not found.")
        sys.exit(1)
    with open(path, "r") as f:
        return yaml.safe_load(f)


# --- Classification VEP-based ---

IMPACT_ORDER = ['MODIFIER', 'LOW', 'MODERATE', 'HIGH']
IMPACT_COLORS = {
    'HIGH': '#c0392b',
    'MODERATE': '#e67e22',
    'LOW': '#f1c40f',
    'MODIFIER': '#95a5a6'
}

def classify_impact(row):
    func = str(row.get('Func.refGene', '')).lower()
    exonic_func = str(row.get('ExonicFunc.refGene', '')).lower()

    if any(x in exonic_func for x in ['stopgain', 'stoploss', 'frameshift']) or 'splicing' in func:
        return 'HIGH'
    if 'nonsynonymous' in exonic_func:
        return 'MODERATE'
    if 'synonymous' in exonic_func:
        return 'LOW'
    return 'MODIFIER'


def classify_substitution(row):
    ref = str(row.get("Ref", "")).upper()
    alt = str(row.get("Alt", "")).upper()

    if len(ref) != 1 or len(alt) != 1:
        return None

    pair = f"{ref}>{alt}"

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    # Normalize to pyrimidine context
    if ref in ["A", "G"]:
        pair = f"{complement[ref]}>{complement[alt]}"

    return pair if pair in VALID_SUBSTITUTIONS else None


def get_natural_chr_order():
    return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


# Since I'm using files annotated by Annovar, this function does the parse

def parse_annovar(file_path):
    try:
        sep = ',' if file_path.endswith('.csv') else '\t'
        df = pd.read_csv(file_path, sep=sep, low_memory=False)

        df['Impact'] = df.apply(classify_impact, axis=1)
        df['Substitution'] = df.apply(classify_substitution, axis=1)

        relevant_vars = df[df['Impact'].isin(['HIGH', 'MODERATE'])].copy()

        res = {
            "total": len(df),
            "counts_func": df['Func.refGene'].value_counts().to_dict() if 'Func.refGene' in df.columns else {},
            "counts_exonic": df[df['ExonicFunc.refGene'] != '.'].get('ExonicFunc.refGene', pd.Series()).value_counts().to_dict() if 'ExonicFunc.refGene' in df.columns else {},
            "counts_impact": df['Impact'].value_counts().to_dict(),
            "counts_chr": {c: df['Chr'].value_counts().get(c, 0) for c in get_natural_chr_order()} if 'Chr' in df.columns else {},
            "counts_subst": df['Substitution'].value_counts().to_dict(),
            "high_mod_df": relevant_vars[['Gene.refGene', 'Impact', 'ExonicFunc.refGene', 'AAChange.refGene']]
        }

        res["gene_counts_for_plot"] = relevant_vars['Gene.refGene'].value_counts().head(10).to_dict()
        return res

    except Exception as e:
        print(f"[ERROR] Failed to read {file_path}: {e}")
        return None


# --- Function for plotting

def save_bar_plot(mean_dict, std_dict, title, outpath, color_map=None,
                  default_color='skyblue', horizontal=True,
                  sort_by_impact=False, cosmic_subst=False):

    if not mean_dict:
        return

    plt.figure(figsize=(7, 5))

    # --- Label order ---
    if cosmic_subst:
        labels = [l for l in COSMIC_ORDER if l in mean_dict]

    elif sort_by_impact:
        labels = [l for l in IMPACT_ORDER if l in mean_dict]

    elif horizontal:
        items = sorted(mean_dict.items(), key=lambda x: x[1], reverse=False)
        labels = [it[0] for it in items]

    else:
        labels = list(mean_dict.keys())

    values = [mean_dict.get(l, 0) for l in labels]

    # --- colors ---
    if cosmic_subst:
        colors = [COSMIC_COLORS[l] for l in labels]
    elif color_map:
        colors = [color_map.get(l, default_color) for l in labels]
    else:
        colors = default_color

    # --- error for groups ---
    yerr_val = None
    if std_dict:
        std_list = [std_dict.get(l, 0) for l in labels]
        lower_err = [min(v, s) for v, s in zip(values, std_list)]
        yerr_val = [lower_err, std_list]

    # --- plot ---
    if horizontal:
        plt.barh(labels, values, xerr=yerr_val, color=colors,
                 edgecolor='black', alpha=0.85, capsize=3)
    else:
        plt.bar(labels, values, yerr=yerr_val, color=colors,
                edgecolor='black', alpha=0.9, capsize=3)
        plt.xticks(rotation=45, fontsize=9)

    plt.title(title, fontweight='bold', fontsize=11)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# --- tsv export for data collection ---

def save_tsv(data_dict, std_dict, outpath, col_name="Category"):
    try:
        if std_dict:
            df = pd.DataFrame({
                col_name: list(data_dict.keys()),
                "Mean": list(data_dict.values()),
                "StdDev": [std_dict.get(k, 0) for k in data_dict.keys()]
            })
        else:
            df = pd.DataFrame({col_name: list(data_dict.keys()), "Count": list(data_dict.values())})

        df.to_csv(outpath, sep='\t', index=False)

    except Exception as e:
        print(f"[ERROR] Failed to save TSV {outpath}: {e}")


# --- HTML with BeautifulSoup ---

def generate_html_report(output_dir, results_list, report_title, filename, author, date_str):
    html_path = os.path.join(output_dir, filename)

    style = """
    <style>
        body { font-family: 'Segoe UI', sans-serif; padding: 25px; background-color: #f4f7f6; color: #333; }
        h1 { color: #2c3e50; border-bottom: 4px solid #3498db; padding-bottom: 10px; }
        .metadata { background: white; padding: 15px; border-radius: 8px; border-left: 5px solid #3498db; margin-bottom: 30px; }
        .sample-section { background: #fff; padding: 25px; border-radius: 12px; margin-bottom: 50px; }
        .plot-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0; }
        .plot-box { border: 1px solid #eee; padding: 10px; text-align: center; border-radius: 10px; background: #fafafa; }
        .plot-img { width: 100%; height: auto; }
        .stat-badge { display: inline-block; padding: 8px 20px; background: #3498db; color: white; border-radius: 25px; font-weight: bold; margin-bottom: 15px; }
        table { width: 100%; border-collapse: collapse; font-size: 12px; margin-top: 20px; }
        th { background: #34495e; color: white; padding: 10px; text-align: left; }
        td { padding: 8px; border-bottom: 1px solid #eee; }
        .impact-high { color: #c0392b; font-weight: bold; }
        .impact-moderate { color: #d35400; font-weight: bold; }
    </style>
    """

    soup = BeautifulSoup(f"<html><head><meta charset='utf-8'><title>{report_title}</title>{style}</head><body></body></html>", "html.parser")

    soup.body.append(BeautifulSoup(f"<h1>{report_title}</h1>", "html.parser"))
    soup.body.append(BeautifulSoup(f"<div class='metadata'><strong>Analyst:</strong> {author} | <strong>Date:</strong> {date_str}</div>", "html.parser"))

    for res in results_list:
        section = soup.new_tag("div", attrs={"class": "sample-section"})
        section.append(BeautifulSoup(f"<h2>Analysis: {res['sample']}</h2>", "html.parser"))
        section.append(BeautifulSoup(f"<div class='stat-badge'>Mean Variant Count: {res['total']}</div>", "html.parser"))

        grid = soup.new_tag("div", attrs={"class": "plot-grid"})

        for p_type, label in [
            ("Func", "Genomic Regions"),
            ("Exonic", "Exonic Function"),
            ("Impact", "Variant Impact"),
            ("Genes", "Top 10 Genes"),
            ("Chr", "Chromosomes"),
            ("Subst", "Substitution Types")
        ]:
            img_name = res.get(f"img_{p_type.lower()}")
            if img_name:
                box = soup.new_tag("div", attrs={"class": "plot-box"})
                box.append(BeautifulSoup(f"<span style='font-size:11px; font-weight:bold;'>{label}</span>", "html.parser"))
                box.append(soup.new_tag("img", src=f"plots/{img_name}", attrs={"class": "plot-img"}))
                grid.append(box)

        section.append(grid)

        if res.get("table_details"):
            section.append(BeautifulSoup("<h3>Recurrent High/Moderate Variants</h3>", "html.parser"))
            table = soup.new_tag("table")
            table.append(BeautifulSoup("<tr><th>Gene</th><th>Impact</th><th>Function</th><th>AA Change</th></tr>", "html.parser"))

            for item in res['table_details']:
                cl = "impact-high" if item['Impact'] == 'HIGH' else "impact-moderate"
                row = f"<tr><td>{item['Gene.refGene']}</td><td><span class='{cl}'>{item['Impact']}</span></td><td>{item['ExonicFunc.refGene']}</td><td>{item['AAChange.refGene']}</td></tr>"
                table.append(BeautifulSoup(row, "html.parser"))

            section.append(table)

        soup.body.append(section)

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(str(soup))


# --- PIPELINE ---

def run_pipeline(name, p_cfg, settings):
    print(f"\n>>> Running Pipeline: {name.upper()}")

    input_dir, out_dir = p_cfg['input_dir'], p_cfg['output_dir']

    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "tables"), exist_ok=True)

    files = sorted([f for f in os.listdir(input_dir) if f.endswith(('.txt', '.csv', '.tsv'))])

    raw_data = []

    for f in files:
        sample = f.split('.')[0]
        data = parse_annovar(os.path.join(input_dir, f))

        if data:
            data['sample'] = sample
            data['group'] = sample[:-1] if any(char.isdigit() for char in sample) else sample
            raw_data.append(data)

    final_results = []

    # --- GROUP MODE ---
    if 'groups' in name:
        groups = sorted(list(set([d['group'] for d in raw_data])))

        for g in groups:
            g_samples = [d for d in raw_data if d['group'] == g]
            if not g_samples:
                continue

            df_impact = pd.DataFrame([s['counts_impact'] for s in g_samples]).fillna(0)
            df_chr = pd.DataFrame([s['counts_chr'] for s in g_samples]).fillna(0)
            df_func = pd.DataFrame([s['counts_func'] for s in g_samples]).fillna(0)
            df_exonic = pd.DataFrame([s['counts_exonic'] for s in g_samples]).fillna(0)
            df_genes_plot = pd.DataFrame([s['gene_counts_for_plot'] for s in g_samples]).fillna(0)
            df_subst = pd.DataFrame([s['counts_subst'] for s in g_samples]).fillna(0)

            all_vars = pd.concat([s['high_mod_df'] for s in g_samples])

            consenso = (
                all_vars
                .groupby(['Gene.refGene', 'Impact', 'ExonicFunc.refGene', 'AAChange.refGene'])
                .size()
                .reset_index(name='Freq')
            )

            consenso = consenso[consenso['Freq'] >= 3].sort_values(by='Freq', ascending=False).head(20)

            avg_res = {
                "sample": f"Group {g}",
                "total": int(np.mean([s['total'] for s in g_samples])),
                "counts_func": df_func.mean().to_dict(), "std_func": df_func.std().to_dict(),
                "counts_exonic": df_exonic.mean().to_dict(), "std_exonic": df_exonic.std().to_dict(),
                "counts_impact": df_impact.mean().to_dict(), "std_impact": df_impact.std().to_dict(),
                "counts_genes": df_genes_plot.mean().to_dict(), "std_genes": df_genes_plot.std().to_dict(),
                "counts_chr": df_chr.mean().to_dict(), "std_chr": df_chr.std().to_dict(),
                "counts_subst": df_subst.mean().to_dict(), "std_subst": df_subst.std().to_dict(),
                "table_details": consenso.to_dict('records')
            }

            final_results.append(avg_res)

    # --- INDIVIDUAL MODE ---
    else:
        for r in raw_data:
            r['counts_genes'], r['counts_chr'] = r['gene_counts_for_plot'], r['counts_chr']
            r['std_func'] = r['std_exonic'] = r['std_impact'] = r['std_genes'] = r['std_chr'] = r['std_subst'] = None
            r['table_details'] = r['high_mod_df'].head(25).to_dict('records')
            final_results.append(r)

    # --- export and plotting ---
    for res in final_results:
        s_id = res['sample'].replace(" ", "_")

        configs = [
            ("Function", res['counts_func'], res['std_func'], 'mediumseagreen', True, "Region", False),
            ("Exonic", res['counts_exonic'], res['std_exonic'], 'cornflowerblue', True, "Function", False),
            ("Impact", res['counts_impact'], res['std_impact'], None, True, "Impact", True),
            ("Genes", res.get('counts_genes'), res.get('std_genes'), '#8e44ad', True, "Gene", False),
            ("Chr", res.get('counts_chr'), res.get('std_chr'), '#f39c12', False, "Chromosome", False),
            ("Substitution", res.get('counts_subst'), res.get('std_subst'), None, False, "Substitution", False)
        ]

        for p_t, mean_d, std_d, c, h, col_n, sort_imp in configs:
            img_file = f"{name}_{s_id}_{p_t}.png"
            tsv_file = f"{name}_{s_id}_{p_t}.tsv"

            plot_data = dict(sorted(mean_d.items(), key=lambda x: x[1], reverse=True)[:10]) if p_t == "Genes" else mean_d

            save_bar_plot(
                plot_data,
                std_d,
                f"{res['sample']} - {p_t}",
                os.path.join(out_dir, "plots", img_file),
                color_map=IMPACT_COLORS if sort_imp else None,
                default_color=c if c else 'skyblue',
                horizontal=h,
                sort_by_impact=sort_imp,
                cosmic_subst=(p_t == "Subst")
            )

            save_tsv(mean_d, std_d, os.path.join(out_dir, "tables", tsv_file), col_n)
            res[f"img_{p_t.lower()}"] = img_file

        if res.get('table_details'):
            pd.DataFrame(res['table_details']).to_csv(
                os.path.join(out_dir, "tables", f"{name}_{s_id}_Variants.tsv"),
                sep='\t', index=False
            )

    if final_results:
        generate_html_report(
            out_dir,
            final_results,
            p_cfg['report_title'],
            p_cfg['html_filename'],
            settings['author'],
            datetime.now().strftime("%Y-%m-%d")
        )


def main():
    cfg = load_full_config()

    for n, p in cfg.get('pipelines', {}).items():
        run_pipeline(n, p, cfg.get('settings', {}))


if __name__ == "__main__":
    main()

