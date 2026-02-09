#!/usr/bin/env python3
import os
import sys
import yaml
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from bs4 import BeautifulSoup
from datetime import datetime

#load yaml

def load_full_config(path="config.yaml"):
    if not os.path.exists(path):
        print(f"[ERRO CRÍTICO] Arquivo {path} não encontrado.")
        sys.exit(1)
    with open(path, "r") as f:
        return yaml.safe_load(f)

# initializing and recovering the vcf infos
def parse_vcf(vcf_path, dp_field="DP", af_field="AF"):
    snps, indels = 0, 0
    transitions, transversions = 0, 0
    dp_values, af_values = [], []
    ti_set = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
# possible vcf ou vcf.gz (good for genome)
    open_func = gzip.open if vcf_path.endswith(".gz") else open
    mode = "rt" if vcf_path.endswith(".gz") else "r"

    try:
        with open_func(vcf_path, mode) as vcf:
            for line in vcf:
                if line.startswith("#"): continue
                cols = line.strip().split("\t")
                if len(cols) < 10: continue #jumping the head

                ref, alt = cols[3].upper(), cols[4].upper()
                if len(ref) == 1 and len(alt) == 1:
                    snps += 1
                    if (ref, alt) in ti_set: transitions += 1
                    else: transversions += 1
                else: indels += 1

                fmt = dict(zip(cols[8].split(":"), cols[9].split(":")))
                if dp_field in fmt and fmt[dp_field].isdigit():
                    dp_values.append(int(fmt[dp_field]))
                if af_field in fmt:
                    try: af_values.append(float(fmt[af_field].split(",")[0]))
                    except: pass
    except Exception as e:
        print(f"[ERROR] Reading failure {vcf_path}: {e}")
        return None

    titv = transitions / transversions if transversions > 0 else 0.0
    return {"snps": snps, "indels": indels, "titv_ratio": round(titv, 2),
            "dp_values": dp_values, "af_values": af_values}

# histogram configs
def get_hist_max_y(values, xmin, xmax, bins=50):
    if not values: return 0
    counts, _ = np.histogram(values, bins=bins, range=(xmin, xmax))
    return counts.max()

def save_plot(values, title, outpath, xmin=0, xmax=None, ylim=None, color='skyblue'):
    plt.figure(figsize=(6, 4)) 
    histo_range = (xmin, xmax) if xmax is not None else None
    plt.hist(values, bins=50, range=histo_range, color=color, edgecolor='black', alpha=0.8, linewidth=0.5)

    if ylim and ylim > 0: plt.ylim(0, ylim)
    if xmax is not None: plt.xlim(xmin, xmax)

    plt.title(title, fontweight='bold', fontsize=10)
    plt.xlabel("Valor", fontsize=9); plt.ylabel("Contagem", fontsize=9)
    plt.xticks(fontsize=8); plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=100)
    plt.close()

# internal samplesheet based on input_dir
def auto_generate_samplesheet(input_dir, output_path):
    if not input_dir or not os.path.exists(input_dir): return False
    files = sorted([f for f in os.listdir(input_dir) if f.endswith(".vcf") or f.endswith(".vcf.gz")])
    if not files: return False
    
    data = []
    for f in files:
        sample_name = f.replace(".vcf.gz", "").replace(".vcf", "").split(".")[0]
        data.append({"sample": sample_name, "file": os.path.join(input_dir, f)})
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    pd.DataFrame(data).to_csv(output_path, sep='\t', index=False)
    return True

def generate_html_report(output_dir, df, report_title, filename, author, date_str):
    html_path = os.path.join(output_dir, filename)
    style = """
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 20px; background-color: #f4f7f6; color: #333; }
        h1 { color: #2c3e50; border-bottom: 4px solid #3498db; padding-bottom: 10px; }
        .metadata { background: white; padding: 15px; border-radius: 8px; border-left: 5px solid #3498db; margin-bottom: 25px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
        
        /* Grid de 4 colunas para os gráficos */
        .plot-grid { 
            display: grid; 
            grid-template-columns: repeat(4, 1fr); 
            gap: 15px; 
            margin-top: 20px;
        }
        
        .plot-box { 
            background: white; 
            border: 1px solid #ddd; 
            padding: 10px; 
            border-radius: 10px; 
            text-align: center; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        
        .plot-img { width: 100%; height: auto; border-radius: 4px; }
        .plot-label { font-size: 11px; font-weight: bold; color: #555; display: block; margin-bottom: 8px; }
        
        table { border-collapse: collapse; width: 100%; background: white; margin-bottom: 30px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; }
        th { background-color: #3498db; color: white; padding: 12px; text-align: left; }
        td { padding: 10px; border-bottom: 1px solid #eee; }
        tr:hover { background-color: #f1f1f1; }
    </style>
    """
    soup = BeautifulSoup(f"<html><head><meta charset='utf-8'><title>{report_title}</title>{style}</head><body></body></html>", "html.parser")
    soup.body.append(BeautifulSoup(f"<h1>{report_title}</h1>", "html.parser"))
    soup.body.append(BeautifulSoup(f"<div class='metadata'><strong>Analista:</strong> {author}<br><strong>Data de Emissão:</strong> {date_str}</div>", "html.parser"))
    
    # Tabela
    soup.body.append(BeautifulSoup("<h2>Métricas por Amostra</h2>", "html.parser"))
    cols_show = ['sample', 'snps', 'indels', 'titv_ratio', 'mean_DP', 'mean_AF']
    df_html = df[cols_show].to_html(index=False, border=0, classes='metrics-table')
    soup.body.append(BeautifulSoup(df_html, "html.parser"))

    # Gráficos em Grid
    soup.body.append(BeautifulSoup("<h2>Distribuições (Eixos Padronizados)</h2>", "html.parser"))
    grid_div = soup.new_tag("div", attrs={"class": "plot-grid"})
    
    for _, row in df.iterrows():
        # Adiciona DP e depois AF para cada amostra no grid
        for col_img, tipo in [('img_dp', 'DP'), ('img_af', 'AF')]:
            if col_img in df.columns and pd.notna(row[col_img]):
                box = soup.new_tag("div", attrs={"class": "plot-box"})
                label = soup.new_tag("span", attrs={"class": "plot-label"})
                label.string = f"{row['sample']} | {tipo}"
                
                img_tag = soup.new_tag("img", src=f"plots/{row[col_img]}", attrs={"class": "plot-img"})
                
                box.append(label)
                box.append(img_tag)
                grid_div.append(box)

    soup.body.append(grid_div)
    with open(html_path, "w", encoding="utf-8") as f: f.write(str(soup))

# executando
def run_pipeline(name, p_cfg, global_settings):
    print(f"\n>>> Iniciando Pipeline: {name.upper()}")
    out_dir = p_cfg['output_dir']
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    
    samples_path = os.path.join(out_dir, "samplesheet_internal.tsv")
    if not auto_generate_samplesheet(p_cfg['input_dir'], samples_path):
        print(f"[AVISO] Nenhum dado encontrado para {name}. Pulando..."); return

    df_samples = pd.read_csv(samples_path, sep='\t')
    
# caso não tenha limites definidos no config, o valor segundo valor é o default - impede de não rodar
    L_DP_MIN = global_settings.get('limit_dp_x_min', 0)
    L_DP_MAX = global_settings.get('limit_dp_x_max', 200)
    L_AF_MIN = global_settings.get('limit_af_x_min', 0.05)
    L_AF_MAX = global_settings.get('limit_af_x_max', 1.0)
    USER_Y_DP = global_settings.get('limit_dp_y_max', 0)
    USER_Y_AF = global_settings.get('limit_af_y_max', 0)

    # para o eixo Y, calcula-se automaticamente
    temp_metrics = []
    auto_max_y_dp, auto_max_y_af = 0, 0
    
    for _, row in df_samples.iterrows():
        m = parse_vcf(row['file'])
        if m:
            m['sample'] = row['sample']
            auto_max_y_dp = max(auto_max_y_dp, get_hist_max_y(m['dp_values'], L_DP_MIN, L_DP_MAX))
            auto_max_y_af = max(auto_max_y_af, get_hist_max_y(m['af_values'], L_AF_MIN, L_AF_MAX))
            temp_metrics.append(m)

    final_y_dp = USER_Y_DP if USER_Y_DP > 0 else (auto_max_y_dp * 1.1)
    final_y_af = USER_Y_AF if USER_Y_AF > 0 else (auto_max_y_af * 1.1)

    # plotagem -individual- de graficos para evitar sobrescrição no HTML 
    results = []
    for m in temp_metrics:
        sample = m['sample']
        img_dp_name = f"{name}_{sample}_DP.png"
        img_af_name = f"{name}_{sample}_AF.png"
        
        if m['dp_values']:
            save_plot(m['dp_values'], f"DP: {sample}", os.path.join(out_dir, "plots", img_dp_name), 
                      xmin=L_DP_MIN, xmax=L_DP_MAX, ylim=final_y_dp, color='royalblue')
            m['mean_DP'] = np.mean(m['dp_values'])
            m['img_dp'] = img_dp_name
            
        if m['af_values']:
            save_plot(m['af_values'], f"AF: {sample}", os.path.join(out_dir, "plots", img_af_name), 
                      xmin=L_AF_MIN, xmax=L_AF_MAX, ylim=final_y_af, color='salmon')
            m['mean_AF'] = np.mean(m['af_values'])
            m['img_af'] = img_af_name
            
        results.append(m)

    if results:
        df_final = pd.DataFrame(results).sort_values('sample')
        generate_html_report(out_dir, df_final, p_cfg['report_title'], p_cfg['html_filename'], 
                             global_settings['author'], datetime.now().strftime("%d/%m/%Y"))
        print(f"[SUCESSO] Relatório {name} gerado em: {out_dir}")

def main():
    full_cfg = load_full_config()
    pipelines = full_cfg.get('summary_pipelines', {})
    settings = full_cfg.get('summary_settings', {})
    
    for name, p_cfg in pipelines.items():
        run_pipeline(name, p_cfg, settings)

if __name__ == "__main__":
    main()
