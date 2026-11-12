import os, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde

# ───────────────────────── config ─────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))   # coloque este script em Diario-de-bordo/
RAND_DIR = os.path.join(BASE_DIR, "output_randomization_py")
REF_DIR  = os.path.join(RAND_DIR, "00_reference")
VERT_DIR = os.path.join(RAND_DIR, "01_vertical_permutations")
HORI_DIR = os.path.join(RAND_DIR, "02_horizontal_permutations")
OUT_DIR  = os.path.join(BASE_DIR, "resultados_densidade_nonb/plots_finais")
PNG_DIR  = os.path.join(OUT_DIR, "paineis_png")

WINDOWS = [150, 250, 500]       # janelas upstream (-w a 0)
FLANK_WINDOW = 500              # janela bilateral (-500 a +500)
PROMOTER_TSS_POS = 500          # posição do TSS dentro dos promotores de 1000 pb (ASSUNÇÃO: centro)
MAX_NOISE_POINTS = 200000       # subamostra do ruído p/ KDE rápido (não muda a densidade)

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PNG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 150, "pdf.fonttype": 42})

CLASS_ORDER_PRIORITY = ["G-Quadruplex", "Z-DNA", "Cruciform / Inverted Repeat", "R-Loop", "Slipped DNA (STR)"]
# nome da pasta (sem o prefixo numérico) -> classe
FOLDER_CLASS = {
    "g4boost_g4": "G-Quadruplex", "gfa2_g4": "G-Quadruplex",
    "zseeker_zdna": "Z-DNA", "gfa2_zdna": "Z-DNA",
    "qmrlfs_rloop": "R-Loop",
    "einverted_cruciform": "Cruciform / Inverted Repeat", "gfa2_cruciform": "Cruciform / Inverted Repeat",
    "gfa2_inverted_repeats": "Cruciform / Inverted Repeat",
    "gfa2_triplex": "Triplex", "gfa2_slipped": "Slipped DNA", "gfa2_str": "Slipped DNA (STR)",
    "gfa2_apr": "A-Phased Repeat", "gfa2_mirror_repeats": "Mirror Repeat", "gfa2_direct_repeats": "Direct Repeat"}
REAL_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

# ───────────────────────── parsers ─────────────────────────
def read_fasta_lengths(path):
    lens, name = {}, None
    if not os.path.exists(path): return lens
    with open(path, errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]; lens[name] = 0
            elif name:
                lens[name] += len(line)
    return lens

def read_gene_ids(path):
    if not os.path.exists(path): return []
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]

def parse_einverted(path, gene_ids):
    """EMBOSS einverted (.out): ponto médio de cada bloco 'Score' (entre os dois braços)."""
    recs, name, coords = [], None, []
    def flush():
        if name and coords: recs.append((name, (min(coords) + max(coords)) / 2))
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = re.match(r"=+\s*(\S+)\s*=+\s*$", line.strip())
            if m:
                flush(); coords = []
                tok = m.group(1); ms = re.match(r"seq_(\d+)$", tok)
                if ms:
                    i = int(ms.group(1)) - 1
                    name = gene_ids[i] if 0 <= i < len(gene_ids) else tok
                else:
                    name = tok
                continue
            if line.startswith(": Score"):
                flush(); coords = []; continue
            m = re.match(r"^\s*(\d+)\s+[A-Za-z\-]+\s+(\d+)\s*$", line)
            if m: coords += [int(m.group(1)), int(m.group(2))]
    flush()
    return recs

def parse_positional(path, gene_ids):
    """Retorna [(nome_da_sequencia, ponto_medio)] para .bed/.gff/.csv/.tsv/.out."""
    recs, ext = [], os.path.splitext(path)[1].lower()
    try:
        if ext == ".out":
            return parse_einverted(path, gene_ids)
        if ext in (".bed", ".gff"):
            with open(path, errors="ignore") as fh:
                for line in fh:
                    if not line.strip() or line.startswith(("#", "track", "browser")): continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 3: continue
                    try: s, e = float(p[1]), float(p[2])
                    except ValueError:
                        if len(p) < 5: continue
                        try: s, e = float(p[3]), float(p[4])      # GFF padrão (9 colunas)
                        except ValueError: continue
                    recs.append((p[0], (s + e) / 2))
        elif ext in (".csv", ".tsv"):
            if os.path.getsize(path) == 0: return recs
            df = pd.read_csv(path, sep=None, engine="python", comment="#")
            if df.empty: return recs
            cols = {str(c).lower().strip(): c for c in df.columns}
            s_c = next((cols[c] for c in cols if "start" in c), None)
            e_c = next((cols[c] for c in cols if "end" in c or "stop" in c), None)
            if s_c is None or e_c is None: return recs
            n_c = next((cols[c] for c in cols if c in
                        ("sequence_name", "seqname", "seq_name", "name", "id", "locus_tag", "seqid", "chrom")),
                       df.columns[0])
            s = pd.to_numeric(df[s_c], errors="coerce"); e = pd.to_numeric(df[e_c], errors="coerce")
            for nm, a, b in zip(df[n_c].astype(str), s, e):
                if pd.notna(a) and pd.notna(b): recs.append((nm, (a + b) / 2))
    except Exception as ex:
        print(f"  [aviso] falha ao ler {path}: {ex}")
    return recs

def collect_sets(tool_dir, gene_ids, token, expect_perm):
    """{id_replica: [(nome, mid)]}. Entre variantes duplicadas (ex.: G4Boost vs tool-g4boost)
    de uma mesma réplica, fica a que tem mais posições."""
    sets = {}
    if not os.path.isdir(tool_dir): return sets
    for root, _, files in os.walk(tool_dir):
        for f in sorted(files):
            fl = f.lower()
            if token not in fl or "all_total" in fl or "g4scores" in fl: continue
            if fl.endswith((".fasta", ".fa", ".npy")): continue
            m = re.search(r"perm_?(\d+)", fl)
            if expect_perm and not m: continue
            if not expect_perm and m: continue
            key = m.group(1) if m else "real"
            recs = parse_positional(os.path.join(root, f), gene_ids)
            if len(recs) > len(sets.get(key, [])): sets[key] = recs
    return sets

# ───────────────────────── carga dos dados ─────────────────────────
def load_all():
    gene_ids = read_gene_ids(os.path.join(REF_DIR, "gene_ids.txt"))
    cds_len = read_fasta_lengths(os.path.join(REF_DIR, "real_cds_sequences.fasta"))
    folders = set()
    for d in (REF_DIR, VERT_DIR, HORI_DIR):
        root = os.path.join(d, "non_B_structures")
        if os.path.isdir(root): folders |= set(os.listdir(root))
    classes = {}
    print("\nResumo de motivos por ferramenta (promotor):")
    for folder in sorted(folders):
        tid = re.sub(r"^\d+_", "", folder)
        ref = os.path.join(REF_DIR, "non_B_structures", folder)
        real_sets = collect_sets(ref, gene_ids, "promoter", expect_perm=False)
        real = [mid - PROMOTER_TSS_POS for _, mid in real_sets.get("real", [])]
        vert, hori = [], []
        for s in collect_sets(os.path.join(VERT_DIR, "non_B_structures", folder), gene_ids, "promoter", True).values():
            vert += [mid - PROMOTER_TSS_POS for _, mid in s]
        for s in collect_sets(os.path.join(HORI_DIR, "non_B_structures", folder), gene_ids, "promoter", True).values():
            hori += [mid - PROMOTER_TSS_POS for _, mid in s]
        cds_sets = collect_sets(ref, gene_ids, "cds", expect_perm=False)
        pct = [mid / cds_len[nm] * 100 for nm, mid in cds_sets.get("real", [])
               if nm in cds_len and cds_len[nm] > 0 and 0 <= mid <= cds_len[nm]]
        print(f"  {folder:28s} real={len(real):6d}  vert={len(vert):8d}  hori={len(hori):8d}  corpo={len(pct):6d}")
        if not real and not pct: continue
        cname = FOLDER_CLASS.get(tid, tid)
        classes.setdefault(cname, []).append({
            "id": tid, "real_tss": np.array(real), "real_bin": np.array(pct),
            "vert": np.array(vert), "hori": np.array(hori)})
    return classes

def ordered_class_names(classes):
    return [c for c in CLASS_ORDER_PRIORITY if c in classes] + sorted(c for c in classes if c not in CLASS_ORDER_PRIORITY)

# ───────────────────────── cálculo e plot ─────────────────────────
def calc_kde(positions, x_min, x_max, scale=1e4):
    grid = np.linspace(x_min, x_max, 200)
    if len(positions) < 5: return grid, np.zeros(200)
    try: return grid, gaussian_kde(positions, bw_method="scott")(grid) * scale
    except Exception: return grid, np.zeros(200)

def pooled(provs, key, lo, hi):
    arr = np.concatenate([p[key] for p in provs]) if provs else np.array([])
    return arr[(arr >= lo) & (arr <= hi)]

def render_panel(ax, cname, provs, kind, window, show_ylabel=True):
    """kind: 'up' (-w..0), 'flank' (-w..+w), 'bin' (corpo do gene, 0..100%)."""
    if kind == "bin": lo, hi = 0, 100
    elif kind == "up": lo, hi = -window, 0
    else: lo, hi = -window, window
    for i, p in enumerate(provs):
        d = p["real_bin"] if kind == "bin" else p["real_tss"]
        w = d[(d >= lo) & (d <= hi)]
        if len(w) > 0:
            x, y = calc_kde(w, lo, hi, scale=1 if kind == "bin" else 1e4)
            ax.plot(x, y, color=REAL_COLORS[i % len(REAL_COLORS)], lw=2, label=f"Real: {p['id']} (n={len(w)})")
    if kind != "bin":
        rng = np.random.default_rng(0)
        for key, ls, lab in (("vert", "--", "Vertical Pooled"), ("hori", ":", "Horizontal Pooled")):
            w = pooled(provs, key, lo, hi)
            if len(w) > 0:
                n = len(w)
                if n > MAX_NOISE_POINTS: w = rng.choice(w, MAX_NOISE_POINTS, replace=False)
                x, y = calc_kde(w, lo, hi)
                ax.plot(x, y, color="gray", ls=ls, lw=1.5, alpha=0.8, label=f"{lab} (n={n})")
    if kind == "up":
        ax.set_title(f"{cname} (-{window} bp)"); ax.set_xlim(-window, 0)
        ax.set_xlabel("Distância ao TSS (bp)")
    elif kind == "flank":
        ax.set_title(f"{cname} (±{window} bp)"); ax.set_xlim(-window, window)
        ax.axvline(0, color="black", lw=0.8, alpha=0.6)
        ax.set_xlabel("Distância ao TSS (bp)")
    else:
        ax.set_title(f"{cname} (Corpo do Gene)"); ax.set_xlim(0, 100)
        ax.set_xlabel("Progresso no Gene (%)")
    if show_ylabel:
        ax.set_ylabel("Densidade (%)" if kind == "bin" else "Densidade prob. (×10⁻⁴)")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, frameon=False, fontsize=7)

def slug(s): return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")

def make_plots():
    classes = load_all()
    if not classes:
        print("\n✗ Nenhum motivo real encontrado. Pastas procuradas:")
        for d in (REF_DIR, VERT_DIR, HORI_DIR): print("  ", os.path.join(d, "non_B_structures"))
        sys.exit(1)
    names = ordered_class_names(classes)
    panels = [("up", w) for w in WINDOWS] + [("flank", FLANK_WINDOW), ("bin", None)]

    fig = plt.figure(figsize=(22, 4.5 * len(names)))
    gs = gridspec.GridSpec(len(names), len(panels), wspace=0.3, hspace=0.7)
    for row, cname in enumerate(names):
        provs = classes[cname]
        for col, (kind, w) in enumerate(panels):
            ax = fig.add_subplot(gs[row, col])
            render_panel(ax, cname, provs, kind, w, show_ylabel=(col == 0))
            # PNG individual do painel
            tag = f"{kind}{w}" if w else "corpo_gene"
            f2, ax2 = plt.subplots(figsize=(5.5, 4.6))
            render_panel(ax2, cname, provs, kind, w)
            f2.savefig(os.path.join(PNG_DIR, f"{slug(cname)}_{tag}.png"), dpi=300, bbox_inches="tight")
            plt.close(f2)

    fig.suptitle("Distribuição posicional (ferramentas separadas + ±500 bp do TSS + corpo do gene)", fontsize=14, y=0.98)
    out_pdf = os.path.join(OUT_DIR, "todas_classes_ferramentas_separadas.pdf")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n✓ Figura completa: {out_pdf}")
    print(f"✓ Painéis individuais (PNG): {PNG_DIR}")

if __name__ == "__main__":
    make_plots()
