import os, glob, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde, ks_2samp, wasserstein_distance

#config 
BASE_DIR = "./codigo/Diario-de-bordo"
GFF_PATH = os.path.join(BASE_DIR, "dados/Hsalinarum.gff")
TSS_GTF  = os.path.join(BASE_DIR, "dados/Hsalinarum_TSS.gtf")
GFA_DIR  = os.path.join(BASE_DIR, "ferramenta_analises_brutas/non-B_gfa")
NEW_TOOLS_DIR = os.path.join(BASE_DIR, "output_randomization_py/resultados")
OUT_DIR  = os.path.join(BASE_DIR, "resultados_densidade_nonb/plots_finais")
WINDOWS = [150, 250, 500]
FASTA_WINDOW = 250
os.makedirs(OUT_DIR, exist_ok=True)
stats_records = []

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 150})

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"]
CLASS_ORDER_PRIORITY = ["G-Quadruplex", "Z-DNA", "Cruciform / Inverted Repeat", "R-Loop", "Slipped DNA (STR)"]
GFA_MTYPE_LABELS = {
    "GQ": "G-Quadruplex", "IR": "Cruciform/Inverted Repeat", "SD": "Slipped DNA",
    "MR": "Mirror Repeat", "DR": "Direct Repeat", "Z": "Z-DNA", "APR": "A-Phased Repeat"}

TOOL_CLASS_KEYWORDS = [
    ("g4boost", "G-Quadruplex"), ("zseeker", "Z-DNA"), 
    ("qmrlfs", "R-Loop"), ("einverted", "Cruciform / Inverted Repeat")]

def classify_tool_folder(tool_name):
    tl = tool_name.lower()
    for kw, label in TOOL_CLASS_KEYWORDS:
        if kw in tl: return label
    return tool_name
def discover_extra_tools():
    if not os.path.isdir(NEW_TOOLS_DIR): return []
    entries = []
    for name in sorted(os.listdir(NEW_TOOLS_DIR)):
        full = os.path.join(NEW_TOOLS_DIR, name)
        if not os.path.isdir(full) or name.lower() == "gfa": continue
        m = re.match(r"^(real|vertical|horizontal)_(.+)$", name, re.IGNORECASE)
        clean_tool = m.group(2) if m else name
        entries.append({"folder": name, "tool": clean_tool, "root": full})
    return entries
def build_class_sources():
    sources = {}
    for f in glob.glob(os.path.join(GFA_DIR, "real_motifs_*.tsv")):
        m = re.search(r"real_motifs_(.+)\.tsv$", os.path.basename(f))
        if m:
            mtype = m.group(1)
            label = GFA_MTYPE_LABELS.get(mtype.upper(), mtype)
            sources.setdefault(label, []).append({"kind": "gfa", "mtype": mtype, "id": f"GFA({mtype})"})
    
    for entry in discover_extra_tools():
        label = classify_tool_folder(entry["tool"])
        sources.setdefault(label, []).append({"kind": "tool", "tool_root": entry["root"], "id": entry["tool"]})
    return sources
def ordered_class_names(sources):
    return [c for c in CLASS_ORDER_PRIORITY if c in sources] + sorted(c for c in sources if c not in CLASS_ORDER_PRIORITY)

# PARSERS 
def parse_gff(path):
    genes = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene": continue
            m = re.search(r'locus_tag=(VNG_RS[^;]+)', cols[8])
            if m: genes[m.group(1)] = {"start": int(cols[3]), "stop": int(cols[4]), "strand": cols[6]}
    return genes
def parse_tss(path):
    records = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or len((cols := line.rstrip("\n").split("\t"))) < 9: continue
            if 'tss_class "P"' in cols[8]:
                records.append({"pos": int(cols[3]), "strand": cols[6]})
    return records
def get_positions(provider, prefix):
    dists, gpos = [], []
    if provider["kind"] == "gfa":
        path = os.path.join(GFA_DIR, f"{prefix}_motifs_{provider['mtype']}.tsv")
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, sep="\t")
                if not df.empty:
                    mids = (df["Start"] + df["Stop"]) / 2
                    if prefix == "real": gpos = mids.tolist()
                    else: dists = (mids - FASTA_WINDOW).tolist()
            except: pass
    else:
        for root, _, files in os.walk(provider["tool_root"]):
            for fname in files:
                if prefix in os.path.join(root, fname).lower():
                    try:
                        df = pd.read_csv(os.path.join(root, fname), sep=None, engine="python", comment="#")
                        s_col, e_col = [c for c in df.columns if 'start' in str(c).lower()], [c for c in df.columns if 'end' in str(c).lower() or 'stop' in str(c).lower()]
                        if s_col and e_col:
                            mids = (pd.to_numeric(df[s_col[0]], errors="coerce") + pd.to_numeric(df[e_col[0]], errors="coerce")) / 2
                            dists.extend((mids.dropna() - FASTA_WINDOW).tolist())
                    except: pass
    return np.array(dists), np.array(gpos)
def calc_distances_and_bins(gpos, tss_list, genes):
    dists, pcts = [], []
    if len(gpos) > 0:
        tss_pos, tss_strand = np.array([t["pos"] for t in tss_list]), np.array([t["strand"] for t in tss_list])
        for m in gpos:
            i = np.argmin(np.abs(tss_pos - m))
            dists.append(m - tss_pos[i] if tss_strand[i] == "+" else tss_pos[i] - m)
            for tag, g in genes.items():
                if g["start"] <= m <= g["stop"]:
                    size = g["stop"] - g["start"]
                    pcts.append(((m - g["start"]) / size * 100) if g["strand"] == "+" else ((g["stop"] - m) / size * 100))
                    break
    return np.array(dists), np.array(pcts)

# calculo e plot
def calc_kde(positions, x_min, x_max, scale=1e4):
    grid = np.linspace(x_min, x_max, 200)
    if len(positions) < 5: return grid, np.zeros(200)
    try: return grid, gaussian_kde(positions, bw_method="scott")(grid) * scale
    except: return grid, np.zeros(200)
def plot_panel(ax, providers, window, tss_list, genes, color, is_bin=False):
    all_vert, all_hori = [], []
    real_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    
    for i, p in enumerate(providers):
        _, gpos_real = get_positions(p, "real")
        d_vert, _ = get_positions(p, "vertical")
        d_hori, _ = get_positions(p, "horizontal")
        
        d_real, p_real = calc_distances_and_bins(gpos_real, tss_list, genes)
        data = p_real if is_bin else d_real
        w_data = data[(data >= 0) & (data <= 100)] if is_bin else data[(data >= -window) & (data <= 0)]
        
        all_vert.extend(d_vert)
        all_hori.extend(d_hori)
        
        if len(w_data) > 0:
            c = real_colors[i % len(real_colors)]
            x, y = calc_kde(w_data, 0 if is_bin else -window, 100 if is_bin else 0, scale=1 if is_bin else 1e4)
            ax.plot(x, y, color=c, ls="-", lw=2, label=f"Real: {p['id']} (n={len(w_data)})")

    # Plota Ruído Agrupado
    if not is_bin:
        w_vert = np.array(all_vert)[(np.array(all_vert) >= -window) & (np.array(all_vert) <= 0)]
        w_hori = np.array(all_hori)[(np.array(all_hori) >= -window) & (np.array(all_hori) <= 0)]
        
        if len(w_vert) > 0:
            x, y = calc_kde(w_vert, -window, 0)
            ax.plot(x, y, color="gray", ls="--", lw=1.5, alpha=0.8, label=f"Vertical Pooled (n={len(w_vert)})")
        if len(w_hori) > 0:
            x, y = calc_kde(w_hori, -window, 0)
            ax.plot(x, y, color="gray", ls=":", lw=1.5, alpha=0.8, label=f"Horizontal Pooled (n={len(w_hori)})")
            
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Legenda movida para baixo para não cobrir o gráfico
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=1, frameon=False, fontsize=7)

def make_plots():
    tss_list = parse_tss(TSS_GTF)
    genes = parse_gff(GFF_PATH)
    sources = build_class_sources()
    class_names = ordered_class_names(sources)
    colors = {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(class_names)}
    
    fig = plt.figure(figsize=(18, 4.5 * len(class_names)))
    gs = gridspec.GridSpec(len(class_names), 4, wspace=0.3, hspace=0.7) # 4 colunas agora (inclui bins)
    
    for row, class_name in enumerate(class_names):
        provs = sources[class_name]
        # Colunas 0 a 2: Janelas Promotoras
        for col, window in enumerate(WINDOWS):
            ax = fig.add_subplot(gs[row, col])
            plot_panel(ax, provs, window, tss_list, genes, colors[class_name], is_bin=False)
            ax.set_title(f"{class_name} (-{window} bp)")
            ax.set_xlim(-window, 0)
            if col == 0: ax.set_ylabel("Densidade prob. (×10⁻⁴)")
        # Coluna 3: Bins do Corpo do Gene
        ax_bin = fig.add_subplot(gs[row, 3])
        plot_panel(ax_bin, provs, None, tss_list, genes, colors[class_name], is_bin=True)
        ax_bin.set_title(f"{class_name} (Corpo do Gene)")
        ax_bin.set_xlim(0, 100)
        ax_bin.set_xlabel("Progresso no Gene (%)")

    out_uni = os.path.join(OUT_DIR, "todas_classes_ferramentas_separadas.pdf")
    fig.suptitle("Distribuição posicional (Ferramentas Separadas + Bins do Gene)", fontsize=14, y=0.98)
    fig.savefig(out_uni, bbox_inches="tight")
    plt.close(fig)
    print(f"\n✓ Gráfico final gerado com sucesso: {out_uni}")
if __name__ == "__main__":
    make_plots()
