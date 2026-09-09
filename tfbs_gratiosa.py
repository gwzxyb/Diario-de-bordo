#!/usr/bin/env python3
"""
Análises de TFBs (ChIP-seq) e estruturas de DNA não-B em Halobacterium salinarum.

Análises incluídas (após modificações):
  1.  Cobertura genômica por replicon
  2.  Promotores COM vs SEM estrutura → recrutamento de TFB (qui-quadrado)
  3.  Cobertura de TFB no promotor × estrutura (Mann‑Whitney, normalizado por kb)
  4.  Expressão gênica × estrutura + TFB no promotor (Mann‑Whitney)
  5.  Orientação de genes vizinhos e estruturas não‑B (tandem/divergente/convergente)
  6.  Proximidade de estruturas não‑B a elementos de inserção (IS)
  7.  Estruturas não‑B a até 50 bp de TFBs (mapeamento fino)
  8.  asRNAs e estruturas de fita simples / R‑loops
  9.  Contatos Hi‑C (lollipop razão rico/pobre)
 10.  Perfil posicional no TSS (±500 bp, contagem absoluta; individual + por classe)
 11.  Distribuição de distâncias mínimas estrutura–TFB (real vs aleatório; KS)
 12.  Enriquecimento de k‑mers em promotores COM vs SEM estrutura (Fisher + FDR)
 13.  Perfil sobreposto de todas as estruturas no TSS (figura única, absoluto + normalizado)
 14.  Genes com/sem 5’UTR (donut)
 15.  Distribuição do tamanho das 5’UTRs
 16.  Distribuição posicional non‑B por região gênica – detalhada (Promotor→5’UTR→decis)
 17.  Distribuição posicional non‑B – agrupada (Promotor/5’UTR/Corpo)
 18.  Contagem de predições por ferramenta (log)
 19.  Conteúdo GC por classe
 20.  Tamanho dos intervalos por classe
"""

import os
import csv
import warnings
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from statsmodels.stats.multitest import multipletests
from Bio import SeqIO
from GRATIOSA import Genome

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURAÇÕES
# =============================================================================

BASE      = "/home/mcavalcante/codigo/Diario-de-bordo/GRATIOSA_DB/data/Hsalinarum/"
CHROMS    = ["NC_002607.1", "NC_001869.1", "NC_002608.1"]
OUT_DIR   = "resultados_tfbs"
SITES_DIR = os.path.join(BASE, "sites/")
ANNOT_DIR = os.path.join(BASE, "annotation/")
CHIPSEQ   = os.path.join(BASE, "chipseq/TFBs.bed")
HIC_FILE  = os.path.join(BASE, "hic/hi-c.bedpe")

REPLICON_LABELS = {
    "NC_002607.1": "Cromossomo principal",
    "NC_001869.1": "Plasmídeo pNRC100",
    "NC_002608.1": "Plasmídeo pNRC200",
}

SS_STRUCTURES = {"R-loop_qmrlfs", "G-Quadruplex_gboost", "G-Quadruplex_nonb"}

PROMOTER_WINDOW = 250
TSS_FLANK       = 500
TFB_PROX        = 50
IS_WINDOW       = 500
HIC_BIN         = 5000
KMER_SIZE       = 6
REAL_FASTA      = os.path.join(OUT_DIR, "real_sequences.fasta")
GENE_IDS_FILE   = os.path.join(OUT_DIR, "gene_ids.txt")

UTR_CUSTOM_SIZE = 100   # tamanho da UTR customizada quando não há UTR real

plt.rcParams.update({
    "font.family"      : "DejaVu Sans",
    "font.size"        : 10,
    "axes.titlesize"   : 11,
    "axes.titleweight" : "bold",
    "axes.spines.top"  : False,
    "axes.spines.right": False,
})

# =============================================================================
# UTILITÁRIOS
# =============================================================================

def _make_dirs():
    for sub in ("tabelas", "figuras"):
        os.makedirs(os.path.join(OUT_DIR, sub), exist_ok=True)

def _iter_tsv(path):
    if not os.path.exists(path):
        return
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip().split("\t")

def _bed_to_mask(path, chrom, length):
    mask = np.zeros(length, dtype=np.int8)
    for cols in _iter_tsv(path):
        if cols[0].startswith("track"):
            continue
        if cols[0] == chrom:
            mask[int(cols[1]):int(cols[2])] = 1
    return mask

def _bed_to_intervals(path, chrom):
    intervals = []
    for cols in _iter_tsv(path):
        if cols[0].startswith("track"):
            continue
        if cols[0] == chrom:
            intervals.append((int(cols[1]), int(cols[2])))
    return intervals

def _gff_to_intervals(path, chrom, feature_type=None):
    intervals = []
    for cols in _iter_tsv(path):
        if len(cols) < 7:
            continue
        if cols[0] != chrom:
            continue
        if feature_type and cols[2] != feature_type:
            continue
        s = int(cols[3]) - 1
        e = int(cols[4])
        strand = cols[6]
        attrs = cols[8] if len(cols) > 8 else ""
        intervals.append((s, e, strand, attrs))
    return intervals

def _write_csv(filename, header, rows):
    path = os.path.join(OUT_DIR, "tabelas", filename)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  [CSV] {path}")
    return path

def list_non_b_files():
    if not os.path.exists(SITES_DIR):
        return []
    return sorted(f for f in os.listdir(SITES_DIR) if f.endswith(".bed"))

def _structure_class(feature):
    f = feature.lower()
    if any(x in f for x in ["z-dna", "zdna", "zseek"]):       return "Z-DNA"
    if any(x in f for x in ["g-quadruplex", "g4", "gboost"]): return "G-Quadruplex"
    if "triplex" in f:                                         return "Triplex"
    if "r-loop" in f:                                          return "R-loop"
    if any(x in f for x in ["short_tandem", "tandem"]):       return "Short_Tandem"
    if "a-phased" in f:                                        return "A-Phased"
    if "cruciform" in f or "inverted" in f:                   return "Cruciform"
    return "Outros"

def _promoter_coords(gene, L):
    strand = "+" if gene.strand else "-"
    tss = gene.start if strand == "+" else gene.end
    if strand == "+":
        return max(0, tss - PROMOTER_WINDOW), tss
    else:
        return tss, min(L, tss + PROMOTER_WINDOW)

def _gene_expr(gene):
    val = getattr(gene, "expression", None)
    if val is None:
        return None
    if hasattr(val, "__iter__"):
        val = float(np.nanmean(list(val)))
    else:
        val = float(val)
    return None if np.isnan(val) else val

def _tss_window(gene, mask, L, flank=TSS_FLANK):
    tss = gene.start if gene.strand else gene.end
    s, e = tss - flank, tss + flank
    if s < 0 or e > L:
        return None
    janela = mask[s:e]
    if not gene.strand:
        janela = janela[::-1]
    if len(janela) != 2 * flank:
        return None
    return janela

def _bar_double(feats, vals_a, vals_b, label_a, label_b, ylabel, title, outpath,
                sigs=None, color_a="#d62728", color_b="#1f77b4"):
    x = np.arange(len(feats))
    w = 0.4
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.bar(x - w/2, vals_a, w, label=label_a, color=color_a, alpha=0.85)
    ax.bar(x + w/2, vals_b, w, label=label_b, color=color_b, alpha=0.85)
    if sigs:
        for i, sig in enumerate(sigs):
            if sig == "*":
                topo = max(vals_a[i], vals_b[i])
                ax.text(i, topo + max(topo * 0.03, 0.5), "*", ha="center", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(feats, fontsize=7)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {outpath}")

CLASS_COLORS = {
    "Z-DNA"        : "#1f77b4",
    "G-Quadruplex" : "#d62728",
    "Triplex"      : "#2ca02c",
    "R-loop"       : "#ff7f0e",
    "Short_Tandem" : "#9467bd",
    "A-Phased"     : "#8c564b",
    "Cruciform"    : "#e377c2",
    "Outros"       : "#7f7f7f",
}

def _bar_by_class(feats, valores, ylabel, title, outpath,
                  sigs=None, hline=None, hline_label=None):
    classes_feat = [_structure_class(f) for f in feats]
    cores = [CLASS_COLORS.get(c, "#7f7f7f") for c in classes_feat]
    fig, ax = plt.subplots(figsize=(13, 6))
    bars = ax.bar(range(len(feats)), valores, color=cores, alpha=0.85)
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", lw=1, label=hline_label)
    if sigs:
        for bar, sig in zip(bars, sigs):
            if sig == "*":
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + abs(bar.get_height()) * 0.02 + 0.01,
                        "*", ha="center", fontsize=14)
    ax.set_xticks(range(len(feats)))
    ax.set_xticklabels([f.replace("_", "\n") for f in feats], fontsize=7)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    classes_presentes = sorted(set(classes_feat))
    handles = [Patch(color=CLASS_COLORS.get(c, "#7f7f7f"), alpha=0.85, label=c)
               for c in classes_presentes]
    if hline is not None:
        handles.append(plt.Line2D([0], [0], color="black", linestyle="--",
                                  label=hline_label))
    ax.legend(handles=handles, fontsize=8, title="Classe estrutural")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {outpath}")

# =============================================================================
# ANÁLISE 1 — COBERTURA GENÔMICA (original, não alterada)
# =============================================================================

def analise_cobertura_genomica(genome, non_b_files):
    print("\n[1] Cobertura genômica das estruturas não-B...")
    rows = []
    cobertura = {}
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        cobertura[feature] = {}
        for chrom in CHROMS:
            chr_obj = genome.genome_dict.get(chrom)
            if chr_obj is None:
                continue
            mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, chr_obj.length)
            n_bp = int(np.sum(mask))
            pct  = n_bp / chr_obj.length * 100
            cobertura[feature][chrom] = pct
            rows.append([feature, chrom, chr_obj.length, n_bp, round(pct, 4)])
    _write_csv("cobertura_estruturas.csv",
               ["feature", "replicon", "comprimento_bp", "bp_cobertos", "pct_replicon"],
               rows)
    features = sorted(cobertura.keys())
    x = np.arange(len(features))
    largura = 0.25
    cores = {"NC_002607.1": "#1f77b4", "NC_001869.1": "#ff7f0e", "NC_002608.1": "#2ca02c"}
    fig, ax = plt.subplots(figsize=(13, 6))
    for i, chrom in enumerate(CHROMS):
        valores = [cobertura[f].get(chrom, 0) for f in features]
        ax.bar(x + i * largura, valores, largura,
               label=REPLICON_LABELS[chrom], color=cores[chrom], alpha=0.85)
    ax.set_xticks(x + largura)
    ax.set_xticklabels([f.replace("_", "\n") for f in features], rotation=0, fontsize=7)
    ax.set_ylabel("% do replicon coberto")
    ax.set_title("Cobertura genômica das estruturas de DNA não-B por replicon")
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "cobertura_genomica.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# =============================================================================
# ANÁLISE 2 — PROMOTORES COM vs SEM ESTRUTURA (original)
# =============================================================================

def analise_promotores_tfb(genome, genes_main, non_b_files):
    print("\n[2] TFBs em promotores com vs sem estrutura...")
    rows = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        L = chr_obj.length
        tfb_mask = _bed_to_mask(CHIPSEQ, chrom, L)
        genes_chrom = [g for g in genes_main.values() if g.chromosome == chrom]
        if not genes_chrom or np.sum(tfb_mask) == 0:
            continue
        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
            if np.sum(nb_mask) == 0:
                continue
            a = b = c = d = 0
            for gene in genes_chrom:
                s, e = _promoter_coords(gene, L)
                tem_nb  = bool(np.any(nb_mask[s:e]))
                tem_tfb = bool(np.any(tfb_mask[s:e]))
                if tem_nb and tem_tfb:   a += 1
                elif tem_nb:             b += 1
                elif tem_tfb:            c += 1
                else:                    d += 1
            n_com = a + b
            n_sem = c + d
            pct_com = a / n_com * 100 if n_com else 0
            pct_sem = c / n_sem * 100 if n_sem else 0
            try:
                chi2, p, _, _ = stats.chi2_contingency(np.array([[a, b], [c, d]]))
            except ValueError:
                chi2, p = np.nan, 1.0
            rows.append([feature, chrom, n_com, a, round(pct_com, 2),
                         n_sem, c, round(pct_sem, 2),
                         round(chi2, 3) if not np.isnan(chi2) else "NA",
                         round(p, 5), "*" if p < 0.05 else "ns"])
    _write_csv("promotores_TFB.csv",
               ["feature", "replicon",
                "n_prom_COM", "n_COM_com_TFB", "pct_COM_com_TFB",
                "n_prom_SEM", "n_SEM_com_TFB", "pct_SEM_com_TFB",
                "chi2", "p_valor", "sig"], rows)
    rows_main = [r for r in rows if r[1] == "NC_002607.1"]
    if rows_main:
        feats = [r[0].replace("_", "\n") for r in rows_main]
        _bar_double(
            feats, [r[4] for r in rows_main], [r[7] for r in rows_main],
            "Promotor COM estrutura não-B", "Promotor SEM estrutura não-B",
            "% de promotores com pelo menos 1 TFB",
            "Recrutamento de TFB em promotores com vs sem estrutura não-B\n"
            "Cromossomo principal  |  * = p < 0,05 (qui-quadrado)",
            os.path.join(OUT_DIR, "figuras", "promotores_TFB.png"),
            sigs=[r[10] for r in rows_main])

# =============================================================================
# ANÁLISE 3 (antiga coocorrência) – REMOVIDA
# =============================================================================
# (removida completamente)

# =============================================================================
# ANÁLISE 4 — COBERTURA DE TFB NO PROMOTOR (NORMALIZADA POR kb)
# =============================================================================

def analise_nTFB_por_gene(genome, genes_main, non_b_files):
    print("\n[4] Cobertura de TFB no promotor (TFB/kb) × estrutura...")
    rows = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        L = chr_obj.length
        tfb_mask = _bed_to_mask(CHIPSEQ, chrom, L)
        genes_chrom = [g for g in genes_main.values() if g.chromosome == chrom]
        if not genes_chrom or np.sum(tfb_mask) == 0:
            continue
        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
            if np.sum(nb_mask) == 0:
                continue
            com_nb, sem_nb = [], []
            for gene in genes_chrom:
                s, e = _promoter_coords(gene, L)
                n_tfb_bp = int(np.sum(tfb_mask[s:e]))
                # Normalizar por kb: bp de TFB / tamanho_janela (em kb)
                tfb_per_kb = n_tfb_bp / ((e - s) / 1000) if (e - s) > 0 else 0
                (com_nb if np.any(nb_mask[s:e]) else sem_nb).append(tfb_per_kb)
            if len(com_nb) < 3 or len(sem_nb) < 3:
                continue
            u, p = stats.mannwhitneyu(com_nb, sem_nb, alternative="two-sided")
            rows.append([feature, chrom,
                         len(com_nb), round(float(np.mean(com_nb)), 2),
                         round(float(np.median(com_nb)), 2),
                         len(sem_nb), round(float(np.mean(sem_nb)), 2),
                         round(float(np.median(sem_nb)), 2),
                         round(u, 1), round(p, 5), "*" if p < 0.05 else "ns"])
    _write_csv("nTFB_por_gene.csv",
               ["feature", "replicon",
                "n_genes_COM", "media_TFB_kb_COM", "mediana_TFB_kb_COM",
                "n_genes_SEM", "media_TFB_kb_SEM", "mediana_TFB_kb_SEM",
                "U", "p_valor", "sig"], rows)
    rows_main = [r for r in rows if r[1] == "NC_002607.1"]
    if rows_main:
        feats = [r[0].replace("_", "\n") for r in rows_main]
        _bar_double(
            feats, [r[4] for r in rows_main], [r[7] for r in rows_main],
            "Promotor COM estrutura não-B", "Promotor SEM estrutura não-B",
            "Mediana de TFB por kb na janela promotora",
            "Cobertura de TFB no promotor (por kb): genes com vs sem estrutura não-B\n"
            "Cromossomo principal  |  * = p < 0,05 (Mann-Whitney)",
            os.path.join(OUT_DIR, "figuras", "nTFB_por_gene.png"),
            sigs=[r[10] for r in rows_main])

# =============================================================================
# ANÁLISE 5 — EXPRESSÃO GÊNICA (original)
# =============================================================================

def analise_expressao_genica(genome, genes_main, non_b_files):
    print("\n[5] Expressão gênica × estrutura + TFB no promotor...")
    try:
        genome.load_expression()
    except Exception as exc:
        print(f"  [AVISO] Não foi possível carregar expressão: {exc}")
        return
    rows = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        L = chr_obj.length
        tfb_mask = _bed_to_mask(CHIPSEQ, chrom, L)
        genes_chrom = [g for g in genes_main.values()
                       if g.chromosome == chrom and _gene_expr(g) is not None]
        if not genes_chrom or np.sum(tfb_mask) == 0:
            continue
        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
            if np.sum(nb_mask) == 0:
                continue
            expr_A, expr_B, expr_C = [], [], []
            for gene in genes_chrom:
                s, e = _promoter_coords(gene, L)
                tem_tfb = bool(np.any(tfb_mask[s:e]))
                tem_nb  = bool(np.any(nb_mask[s:e]))
                val     = _gene_expr(gene)
                if val is None:
                    continue
                if tem_nb and tem_tfb:   expr_A.append(val)
                elif tem_tfb:            expr_B.append(val)
                else:                    expr_C.append(val)
            def _mw(a, b):
                if len(a) >= 3 and len(b) >= 3:
                    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                    return round(float(u), 1), round(float(p), 5)
                return "NA", 1.0
            u_AB, p_AB = _mw(expr_A, expr_B)
            u_AC, p_AC = _mw(expr_A, expr_C)
            rows.append([feature, chrom,
                         len(expr_A),
                         round(float(np.median(expr_A)), 4) if expr_A else "NA",
                         len(expr_B),
                         round(float(np.median(expr_B)), 4) if expr_B else "NA",
                         len(expr_C),
                         round(float(np.median(expr_C)), 4) if expr_C else "NA",
                         u_AB, p_AB,
                         "*" if isinstance(p_AB, float) and p_AB < 0.05 else "ns",
                         u_AC, p_AC,
                         "*" if isinstance(p_AC, float) and p_AC < 0.05 else "ns"])
    if not rows:
        return
    _write_csv("expressao_genica.csv",
               ["feature", "replicon",
                "n_A_estrutura_TFB", "mediana_expr_A",
                "n_B_TFB_only", "mediana_expr_B",
                "n_C_nenhum", "mediana_expr_C",
                "U_A_vs_B", "p_A_vs_B", "sig_A_vs_B",
                "U_A_vs_C", "p_A_vs_C", "sig_A_vs_C"], rows)
    rows_main = [r for r in rows if r[1] == "NC_002607.1"]
    if rows_main:
        feats = [r[0].replace("_", "\n") for r in rows_main]
        med_A = [r[3] if r[3] != "NA" else 0 for r in rows_main]
        med_B = [r[5] if r[5] != "NA" else 0 for r in rows_main]
        med_C = [r[7] if r[7] != "NA" else 0 for r in rows_main]
        x = np.arange(len(feats)); w = 0.27
        fig, ax = plt.subplots(figsize=(14, 6))
        ax.bar(x - w, med_A, w, label="A: estrutura + TFB", color="#d62728", alpha=0.85)
        ax.bar(x,     med_B, w, label="B: TFB sem estrutura", color="#ff7f0e", alpha=0.85)
        ax.bar(x + w, med_C, w, label="C: sem TFB",           color="#1f77b4", alpha=0.85)
        for i, r in enumerate(rows_main):
            topo = max(med_A[i], med_B[i], med_C[i])
            if r[10] == "*":
                ax.text(i - w/2, topo * 1.05, "†", ha="center", fontsize=12, color="#d62728")
            if r[13] == "*":
                ax.text(i + w/2, topo * 1.05, "‡", ha="center", fontsize=12, color="#333")
        ax.set_xticks(x)
        ax.set_xticklabels(feats, fontsize=7)
        ax.set_ylabel("Mediana de expressão gênica")
        ax.set_title("Expressão gênica por grupo (estrutura não-B × TFB no promotor)\n"
                     "Cromossomo principal  |  † A vs B  ‡ A vs C  (p < 0,05, Mann-Whitney)")
        ax.legend(fontsize=9)
        fig.tight_layout()
        out = os.path.join(OUT_DIR, "figuras", "expressao_genica.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  [FIG] {out}")

# =============================================================================
# ANÁLISE 6 — ORIENTAÇÃO (renomeado co‑direcional → tandem)
# =============================================================================

def analise_orientacao_genes(genome, genes_main, non_b_files):
    """
    Classifica o espaço intergênico pelo contexto de orientação dos genes
    flanqueadores:
      - divergente  (← →): supercoiling negativo local
      - convergente (→ ←): supercoiling positivo
      - tandem (co‑direcional): mesmo sentido
    """
    print("\n[6] Orientação de genes vizinhos e estruturas não-B...")
    chrom   = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if chr_obj is None:
        return
    L = chr_obj.length
    genes_sorted = sorted(
        [g for g in genes_main.values() if g.chromosome == chrom],
        key=lambda g: g.start)
    def _strand(g):
        return "+" if g.strand else "-"
    def _context(gl, gr):
        sl, sr = _strand(gl), _strand(gr)
        if sl == "+" and sr == "-":  return "convergente"
        if sl == "-" and sr == "+":  return "divergente"
        return "tandem (co‑direcional)"
    ctx_mask = np.full(L, "sem_vizinho", dtype=object)
    for i in range(len(genes_sorted) - 1):
        g1, g2 = genes_sorted[i], genes_sorted[i + 1]
        s, e = g1.end, g2.start
        if s >= e or s < 0 or e > L:
            continue
        ctx_mask[s:e] = _context(g1, g2)
    ctx_types = ["divergente", "convergente", "tandem (co‑direcional)", "sem_vizinho"]
    rows = []
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
        if np.sum(nb_mask) == 0:
            continue
        nb_pos     = np.where(nb_mask == 1)[0]
        counts_nb  = {c: int(np.sum(ctx_mask[nb_pos] == c)) for c in ctx_types}
        counts_all = {c: int(np.sum(ctx_mask == c))          for c in ctx_types}
        tab = np.array([
            [counts_nb.get(c, 0) for c in ctx_types],
            [counts_all.get(c, 0) - counts_nb.get(c, 0) for c in ctx_types]
        ])
        try:
            chi2, p, _, _ = stats.chi2_contingency(tab)
        except ValueError:
            chi2, p = np.nan, 1.0
        for c in ctx_types:
            pct_nb  = counts_nb[c]  / len(nb_pos) * 100 if len(nb_pos) else 0
            pct_all = counts_all[c] / L           * 100
            rows.append([feature, chrom, c,
                         counts_nb[c], round(pct_nb, 2),
                         counts_all[c], round(pct_all, 2),
                         round(chi2, 3) if not np.isnan(chi2) else "NA",
                         round(p, 5), "*" if p < 0.05 else "ns"])
    _write_csv("orientacao_genes_naoB.csv",
               ["feature", "replicon", "contexto",
                "bp_NB_no_contexto", "pct_NB_no_contexto",
                "bp_total_no_contexto", "pct_genoma_no_contexto",
                "chi2", "p_valor", "sig"], rows)
    features_uniq = sorted(set(r[0] for r in rows))
    cores_ctx = {"divergente": "#d62728", "convergente": "#1f77b4",
                 "tandem (co‑direcional)": "#2ca02c", "sem_vizinho": "#aec7e8"}
    x = np.arange(len(features_uniq))
    bottoms = np.zeros(len(features_uniq))
    fig, ax = plt.subplots(figsize=(13, 6))
    for ctx in ctx_types:
        vals = []
        for feat in features_uniq:
            sub = [r for r in rows if r[0] == feat and r[2] == ctx]
            vals.append(sub[0][4] if sub else 0)
        ax.bar(x, vals, bottom=bottoms, label=ctx, color=cores_ctx[ctx], alpha=0.85)
        bottoms += np.array(vals)
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("_", "\n") for f in features_uniq], fontsize=7)
    ax.set_ylabel("% dos bp da estrutura por contexto de orientação")
    ax.set_title("Contexto de orientação gênica das estruturas não-B\n"
                 "Cromossomo principal  |  divergente = supercoiling negativo local")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "orientacao_genes_naoB.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# =============================================================================
# ANÁLISE 7 — PROXIMIDADE IS (original)
# =============================================================================

def analise_proximidade_IS(genome, non_b_files):
    print(f"\n[7] Proximidade de estruturas não-B a elementos IS (±{IS_WINDOW} bp)...")
    is_path = os.path.join(ANNOT_DIR, "IS_fixed.gff3")
    if not os.path.exists(is_path):
        is_path = os.path.join(ANNOT_DIR, "IS.gff3")
    rows = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        L = chr_obj.length
        is_intervals = _gff_to_intervals(is_path, chrom)
        if not is_intervals:
            continue
        is_win_mask = np.zeros(L, dtype=np.int8)
        for s, e, *_ in is_intervals:
            is_win_mask[max(0, s - IS_WINDOW): min(L, e + IS_WINDOW)] = 1
        n_is_bp    = int(np.sum(is_win_mask))
        n_nonIS_bp = L - n_is_bp
        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
            if np.sum(nb_mask) == 0:
                continue
            nb_em_IS   = int(np.sum((nb_mask == 1) & (is_win_mask == 1)))
            nb_fora_IS = int(np.sum((nb_mask == 1) & (is_win_mask == 0)))
            dens_IS    = nb_em_IS   / n_is_bp    * 1000 if n_is_bp    else 0
            dens_out   = nb_fora_IS / n_nonIS_bp * 1000 if n_nonIS_bp else 0
            tab = np.array([
                [nb_em_IS,   n_is_bp    - nb_em_IS],
                [nb_fora_IS, n_nonIS_bp - nb_fora_IS]
            ])
            try:
                chi2, p, _, _ = stats.chi2_contingency(tab)
            except ValueError:
                chi2, p = np.nan, 1.0
            razao = round(dens_IS / dens_out, 3) if dens_out else "NA"
            rows.append([feature, chrom, len(is_intervals), n_is_bp,
                         nb_em_IS, round(dens_IS, 4),
                         nb_fora_IS, round(dens_out, 4), razao,
                         round(chi2, 3) if not np.isnan(chi2) else "NA",
                         round(p, 5), "*" if p < 0.05 else "ns"])
    _write_csv("proximidade_IS.csv",
               ["feature", "replicon", "n_IS", "bp_janela_IS",
                "bp_NB_em_IS", "densidade_NB_IS_por_kb",
                "bp_NB_fora_IS", "densidade_NB_fora_IS_por_kb",
                "razao_densidade", "chi2", "p_valor", "sig"], rows)
    rows_main = [r for r in rows if r[1] == "NC_002607.1" and r[8] != "NA"]
    if rows_main:
        feats  = [r[0] for r in rows_main]
        razoes = [float(r[8]) for r in rows_main]
        sigs   = [r[11] for r in rows_main]
        _bar_by_class(
            feats, razoes,
            "Razão de densidade (dentro IS / fora IS)",
            f"Enriquecimento de estruturas não-B próximas a elementos IS\n"
            f"Cromossomo principal  |  janela ±{IS_WINDOW} bp  |  * = p < 0,05",
            os.path.join(OUT_DIR, "figuras", "proximidade_IS.png"),
            sigs=sigs, hline=1.0, hline_label="Densidade igual dentro/fora de IS")

# =============================================================================
# ANÁLISE 8 — NÃO‑B PERTO DE TFB (original)
# =============================================================================

def analise_naoB_perto_TFB(genome, non_b_files):
    print(f"\n[8] Estruturas não-B a até {TFB_PROX} bp de TFBs...")
    rows = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        L   = chr_obj.length
        tfb = _bed_to_intervals(CHIPSEQ, chrom)
        if not tfb:
            continue
        n_tfb = len(tfb)
        tfb_prox_mask = np.zeros(L, dtype=np.int8)
        for s, e in tfb:
            tfb_prox_mask[max(0, s - TFB_PROX): min(L, e + TFB_PROX)] = 1
        tfb_exact_mask = _bed_to_mask(CHIPSEQ, chrom, L)
        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
            if np.sum(nb_mask) == 0:
                continue
            viz_mask = (tfb_prox_mask == 1) & (tfb_exact_mask == 0)
            n_sob  = int(np.sum((nb_mask == 1) & (tfb_exact_mask == 1)))
            n_viz  = int(np.sum((nb_mask == 1) & viz_mask))
            n_dist = int(np.sum((nb_mask == 1) & (tfb_prox_mask == 0)))
            n_tot  = int(np.sum(nb_mask))
            pct_sob  = n_sob  / n_tot * 100 if n_tot else 0
            pct_viz  = n_viz  / n_tot * 100 if n_tot else 0
            pct_dist = n_dist / n_tot * 100 if n_tot else 0
            n_tfb_com_viz = sum(
                1 for s, e in tfb
                if np.any(nb_mask[max(0, s - TFB_PROX): min(L, e + TFB_PROX)])
            )
            rows.append([feature, chrom, n_tfb,
                         n_sob, round(pct_sob, 2),
                         n_viz, round(pct_viz, 2),
                         n_dist, round(pct_dist, 2),
                         n_tfb_com_viz,
                         round(n_tfb_com_viz / n_tfb * 100, 2) if n_tfb else 0])
    _write_csv("naoB_perto_TFB.csv",
               ["feature", "replicon", "n_TFB_total",
                "n_NB_sobreposto", "pct_NB_sobreposto",
                "n_NB_vizinho_50bp", "pct_NB_vizinho_50bp",
                "n_NB_distante", "pct_NB_distante",
                "n_TFB_com_NB_vizinho", "pct_TFB_com_NB_vizinho"], rows)
    rows_main = [r for r in rows if r[1] == "NC_002607.1"]
    if rows_main:
        feats   = [r[0].replace("_", "\n") for r in rows_main]
        bp_sob  = [r[3] for r in rows_main]
        bp_viz  = [r[5] for r in rows_main]
        bp_dis  = [r[7] for r in rows_main]
        x = np.arange(len(feats))
        fig, ax = plt.subplots(figsize=(13, 6))
        ax.bar(x, bp_sob, label="Sobreposto ao TFB", color="#d62728", alpha=0.9)
        ax.bar(x, bp_viz, bottom=bp_sob,
               label=f"Vizinho (1–{TFB_PROX} bp)", color="#ff7f0e", alpha=0.9)
        bot2 = [a + b for a, b in zip(bp_sob, bp_viz)]
        ax.bar(x, bp_dis, bottom=bot2,
               label=f"Distante (>{TFB_PROX} bp)", color="#aec7e8", alpha=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(feats, fontsize=7)
        ax.set_ylabel("Nº de bp da estrutura")
        ax.set_title(f"Localização das estruturas não-B em relação aos TFBs\n"
                     f"Cromossomo principal  |  vizinhança = ±{TFB_PROX} bp")
        ax.legend(fontsize=9)
        fig.tight_layout()
        out = os.path.join(OUT_DIR, "figuras", "naoB_perto_TFB.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  [FIG] {out}")

# =============================================================================
# ANÁLISE 9 — asRNA (original)
# =============================================================================

def analise_asRNA_estruturas(genome, non_b_files):
    print("\n[9] asRNAs e estruturas de fita simples / R-loops...")
    asrna_path = os.path.join(ANNOT_DIR, "asRNAs.gff3")
    if not os.path.exists(asrna_path):
        print("  [AVISO] asRNAs.gff3 não encontrado.")
        return
    rows = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        L = chr_obj.length
        asrna_intervals = _gff_to_intervals(asrna_path, chrom)
        if not asrna_intervals:
            continue
        asrna_mask = np.zeros(L, dtype=np.int8)
        for s, e, *_ in asrna_intervals:
            asrna_mask[s:e] = 1
        n_as  = int(np.sum(asrna_mask))
        n_out = L - n_as
        if n_as == 0:
            continue
        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
            if np.sum(nb_mask) == 0:
                continue
            nb_as  = int(np.sum((nb_mask == 1) & (asrna_mask == 1)))
            nb_out = int(np.sum((nb_mask == 1) & (asrna_mask == 0)))
            d_as   = nb_as  / n_as  * 1000 if n_as  else 0
            d_out  = nb_out / n_out * 1000 if n_out else 0
            tab = np.array([
                [nb_as,  n_as  - nb_as],
                [nb_out, n_out - nb_out]
            ])
            try:
                chi2, p, _, _ = stats.chi2_contingency(tab)
            except ValueError:
                chi2, p = np.nan, 1.0
            razao = round(d_as / d_out, 3) if d_out else "NA"
            is_ss = "sim" if feature in SS_STRUCTURES else "nao"
            rows.append([feature, chrom, is_ss,
                         n_as, nb_as, round(d_as, 4),
                         n_out, nb_out, round(d_out, 4), razao,
                         round(chi2, 3) if not np.isnan(chi2) else "NA",
                         round(p, 5), "*" if p < 0.05 else "ns"])
    _write_csv("asRNA_estruturas.csv",
               ["feature", "replicon", "estrutura_ss",
                "bp_asRNA", "bp_NB_em_asRNA", "dens_NB_asRNA_por_kb",
                "bp_nao_asRNA", "bp_NB_fora_asRNA", "dens_NB_fora_asRNA_por_kb",
                "razao_densidade", "chi2", "p_valor", "sig"], rows)
    rows_main = [r for r in rows if r[1] == "NC_002607.1" and r[9] != "NA"]
    if rows_main:
        feats  = [r[0] for r in rows_main]
        razoes = [float(r[9]) for r in rows_main]
        sigs   = [r[12] for r in rows_main]
        _bar_by_class(
            feats, razoes,
            "Razão de densidade (dentro asRNA / fora asRNA)",
            "Estruturas não-B em regiões com transcrição antisense\n"
            "Cromossomo principal  |  * = p < 0,05 (qui-quadrado)",
            os.path.join(OUT_DIR, "figuras", "asRNA_estruturas.png"),
            sigs=sigs, hline=1.0, hline_label="Razão = 1 (densidade igual)")

# =============================================================================
# ANÁLISE 10 — Hi‑C (LOLLIPOP da razão Grupo1/Grupo3)
# =============================================================================

def analise_hic_naoB(genome, non_b_files):
    print(f"\n[10] Contatos Hi-C × riqueza em estruturas não-B (bin {HIC_BIN} bp)...")
    if not os.path.exists(HIC_FILE):
        print(f"  [AVISO] {HIC_FILE} não encontrado.")
        return
    chrom   = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if chr_obj is None:
        return
    L      = chr_obj.length
    n_bins = L // HIC_BIN + 1
    contacts = []
    for cols in _iter_tsv(HIC_FILE):
        if len(cols) < 6:
            continue
        if cols[0] != chrom or cols[3] != chrom:
            continue
        try:
            m1    = (int(cols[1]) + int(cols[2])) // 2
            m2    = (int(cols[4]) + int(cols[5])) // 2
            score = float(cols[6]) if len(cols) > 6 else 1.0
        except ValueError:
            continue
        b1 = min(m1 // HIC_BIN, n_bins - 1)
        b2 = min(m2 // HIC_BIN, n_bins - 1)
        if b1 != b2:
            contacts.append((b1, b2, score))
    if not contacts:
        print("  [AVISO] Nenhum contato Hi-C intra-cromossômico encontrado.")
        return
    rows = []
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
        if np.sum(nb_mask) == 0:
            continue
        frac_bin = np.array([
            np.mean(nb_mask[b * HIC_BIN: min((b + 1) * HIC_BIN, L)])
            for b in range(n_bins)
        ])
        bins_com = frac_bin[frac_bin > 0]
        med  = float(np.median(bins_com)) if len(bins_com) > 0 else 0
        rico = frac_bin >= med
        scores_rr, scores_rp, scores_pp = [], [], []
        for b1, b2, sc in contacts:
            r1, r2 = rico[b1], rico[b2]
            if r1 and r2:    scores_rr.append(sc)
            elif r1 or r2:   scores_rp.append(sc)
            else:            scores_pp.append(sc)
        def _mw(a, b):
            if len(a) >= 3 and len(b) >= 3:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                return round(float(u), 1), round(float(p), 5)
            return "NA", 1.0
        u_rr_pp, p_rr_pp = _mw(scores_rr, scores_pp)
        u_rr_rp, p_rr_rp = _mw(scores_rr, scores_rp)
        media_rr = float(np.mean(scores_rr)) if scores_rr else np.nan
        media_pp = float(np.mean(scores_pp)) if scores_pp else np.nan
        razao = round(media_rr / media_pp, 3) if media_pp and not np.isnan(media_pp) else "NA"
        rows.append([feature, chrom,
                     round(med * 100, 3),
                     len(scores_rr), round(media_rr, 4) if not np.isnan(media_rr) else "NA",
                     len(scores_rp), round(float(np.mean(scores_rp)), 4) if scores_rp else "NA",
                     len(scores_pp), round(media_pp, 4) if not np.isnan(media_pp) else "NA",
                     razao,
                     u_rr_pp, p_rr_pp,
                     "*" if isinstance(p_rr_pp, float) and p_rr_pp < 0.05 else "ns",
                     u_rr_rp, p_rr_rp,
                     "*" if isinstance(p_rr_rp, float) and p_rr_rp < 0.05 else "ns"])
    if not rows:
        return
    _write_csv("hic_naoB.csv",
               ["feature", "replicon", "limiar_rico_pct_cobertura",
                "n_contatos_RR", "media_score_RR",
                "n_contatos_RP", "media_score_RP",
                "n_contatos_PP", "media_score_PP",
                "razao_RR_PP",
                "U_RR_vs_PP", "p_RR_vs_PP", "sig_RR_vs_PP",
                "U_RR_vs_RP", "p_RR_vs_RP", "sig_RR_vs_RP"], rows)
    # Lollipop horizontal da razão Grupo1/Grupo3
    feats   = [r[0] for r in rows]
    razoes  = [float(r[8]) if r[8] != "NA" else np.nan for r in rows]
    sigs    = [r[11] for r in rows]
    # Ordenar por razão
    order = np.argsort([x if not np.isnan(x) else 0 for x in razoes])
    feats_sorted = [feats[i] for i in order]
    razoes_sorted = [razoes[i] for i in order]
    sigs_sorted = [sigs[i] for i in order]
    y_pos = range(len(feats_sorted))
    fig, ax = plt.subplots(figsize=(10, 0.4 * len(feats_sorted) + 2))
    ax.hlines(y_pos, xmin=1, xmax=razoes_sorted, color="gray", alpha=0.6)
    colors = ["#d62728" if s == "*" else "#1f77b4" for s in sigs_sorted]
    ax.scatter(razoes_sorted, y_pos, c=colors, s=80, zorder=3)
    ax.axvline(1.0, color="black", linestyle="--", lw=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f.replace("_", "\n") for f in feats_sorted], fontsize=7)
    ax.set_xlabel("Razão média de contatos (Grupo 1 / Grupo 3)")
    ax.set_title("Enriquecimento de contatos Hi‑C em bins ricos vs pobres em estrutura\n"
                 f"Cromossomo principal  |  bin {HIC_BIN} bp  |  * = p<0.05 (RR vs PP)")
    ax.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728',
                                  markersize=8, label='p<0.05'),
                       plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4',
                                  markersize=8, label='ns')],
              fontsize=9, loc="lower right")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "hic_naoB_lollipop.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# =============================================================================
# ANÁLISE 11 — PERFIL POSICIONAL NO TSS (original)
# =============================================================================

def analise_perfil_posicional(genome, genes_main, non_b_files):
    print(f"\n[11] Perfil posicional no TSS (±{TSS_FLANK} bp, contagem absoluta)...")
    chrom   = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if chr_obj is None:
        return {}
    L = chr_obj.length
    genes_chrom = sorted(
        [g for g in genes_main.values() if g.chromosome == chrom],
        key=lambda g: g.start)
    x = np.arange(-TSS_FLANK, TSS_FLANK)
    perfis   = {}
    csv_rows = []
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
        if np.sum(nb_mask) == 0:
            continue
        contagem = np.zeros(2 * TSS_FLANK, dtype=np.int64)
        n_validos = 0
        for gene in genes_chrom:
            janela = _tss_window(gene, nb_mask, L)
            if janela is None:
                continue
            contagem += janela.astype(np.int64)
            n_validos += 1
        if n_validos < 10:
            continue
        perfis[feature] = contagem
        for pos, val in zip(x, contagem):
            csv_rows.append([feature, int(pos), int(val)])
        # Figura individual
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(x, contagem, color="#d62728", lw=1.5)
        ax.fill_between(x, contagem, color="#d62728", alpha=0.2)
        ax.axvline(0, color="gray", ls=":", lw=1.2, label="TSS")
        ax.set_xlabel("Posição relativa ao TSS (bp)  |  negativo = upstream, positivo = downstream")
        ax.set_ylabel("Nº de genes com a estrutura nessa posição")
        ax.set_title(f"Perfil posicional no TSS: {feature}\n"
                     f"Cromossomo principal  |  janela ±{TSS_FLANK} bp  |  n = {n_validos} genes")
        ax.legend(fontsize=8)
        fig.tight_layout()
        out = os.path.join(OUT_DIR, "figuras", f"perfil_posicional_{feature}.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  [FIG] {out}")
    _write_csv("perfil_posicional.csv",
               ["feature", "posicao_rel_TSS_bp", "n_genes"], csv_rows)
    # Figuras agrupadas por classe
    classes = defaultdict(list)
    for feature in perfis:
        classes[_structure_class(feature)].append(feature)
    paleta = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
    for classe, feats in classes.items():
        if len(feats) < 2:
            continue
        fig, ax = plt.subplots(figsize=(12, 5))
        for i, feature in enumerate(sorted(feats)):
            label = feature.replace(f"{classe}_", "").replace(classe, "").strip("_") or feature
            ax.plot(x, perfis[feature], lw=1.5, alpha=0.85,
                    color=paleta[i % len(paleta)], label=label)
        ax.axvline(0, color="gray", ls=":", lw=1.2, label="TSS")
        ax.set_xlabel("Posição relativa ao TSS (bp)  |  negativo = upstream, positivo = downstream")
        ax.set_ylabel("Nº de genes com a estrutura nessa posição")
        ax.set_title(f"Perfil posicional no TSS — classe {classe}\n"
                     f"Cromossomo principal  |  janela ±{TSS_FLANK} bp  |  {len(feats)} ferramentas")
        ax.legend(fontsize=8, title="Ferramenta")
        fig.tight_layout()
        out = os.path.join(OUT_DIR, "figuras", f"perfil_posicional_CLASSE_{classe}.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  [FIG] {out}")
    return perfis

# =============================================================================
# ANÁLISE 12 — DISTÂNCIA TFB (original)
# =============================================================================

def analise_distancia_TFB(genome, non_b_files):
    print("\n[12] Distribuição de distâncias mínimas estrutura–TFB...")
    chrom   = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if chr_obj is None:
        return
    L = chr_obj.length
    tfb = _bed_to_intervals(CHIPSEQ, chrom)
    if not tfb:
        return
    tfb_centers = np.array([(s + e) // 2 for s, e in tfb])
    rows = []
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        nb = _bed_to_intervals(os.path.join(SITES_DIR, bed_file), chrom)
        if not nb:
            continue
        nb_centers = np.array([(s + e) // 2 for s, e in nb])
        dists_real = np.array([np.min(np.abs(c - tfb_centers)) for c in nb_centers])
        mediana = float(np.median(dists_real))
        media   = float(np.mean(dists_real))
        rows.append([feature, chrom, len(nb_centers),
                     round(mediana, 1), round(media, 1),
                     int(dists_real.min()), int(dists_real.max())])
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(dists_real, bins=60, color="#d62728", alpha=0.8, edgecolor="black", lw=0.3)
        ax.axvline(mediana, color="black", ls="--", lw=1.5,
                   label=f"Mediana = {mediana:.0f} bp")
        ax.set_xlabel("Distância ao TFB mais próximo (bp)")
        ax.set_ylabel("Nº de sítios da estrutura")
        ax.set_title(f"Distância mínima ao TFB: {feature}\n"
                     f"Cromossomo principal  |  n = {len(nb_centers)} sítios")
        ax.legend(fontsize=8)
        fig.tight_layout()
        out = os.path.join(OUT_DIR, "figuras", f"distancia_TFB_{feature}.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  [FIG] {out}")
    _write_csv("distancia_minima_TFB.csv",
               ["feature", "replicon", "n_sitios",
                "mediana_dist_bp", "media_dist_bp",
                "dist_min_bp", "dist_max_bp"], rows)

# =============================================================================
# ANÁLISE 13 — k‑MERS (original)
# =============================================================================

def analise_enriquecimento_kmer(genome, genes_main, non_b_files):
    print(f"\n[13] Enriquecimento de {KMER_SIZE}-mers em promotores COM vs SEM estrutura...")
    chrom   = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if chr_obj is None:
        return
    L = chr_obj.length
    genes_chrom = sorted(
        [g for g in genes_main.values() if g.chromosome == chrom],
        key=lambda g: g.start)
    seqs_por_gene = {}
    if os.path.exists(REAL_FASTA) and os.path.exists(GENE_IDS_FILE):
        with open(GENE_IDS_FILE) as fh:
            gene_ids = [l.strip() for l in fh]
        with open(REAL_FASTA) as fh:
            for rec, gid in zip(SeqIO.parse(fh, "fasta"), gene_ids):
                seqs_por_gene[gid] = str(rec.seq).upper()
        print(f"  Sequências carregadas do FASTA ({len(seqs_por_gene)} genes)")
    else:
        try:
            seq_genoma = np.array(list(chr_obj.sequence.upper()))
        except AttributeError:
            print("  [AVISO] Sequência do genoma não disponível. Pulando análise 13.")
            return
        for gene in genes_chrom:
            s, e = _promoter_coords(gene, L)
            if (e - s) == PROMOTER_WINDOW:
                gid = getattr(gene, "name", str(gene.start))
                seqs_por_gene[gid] = "".join(seq_genoma[s:e])
        print(f"  Sequências extraídas do genoma ({len(seqs_por_gene)} genes)")
    if not seqs_por_gene:
        print("  [AVISO] Nenhuma sequência promotora disponível. Pulando análise 13.")
        return
    def _contar_kmers(lista_seqs, k):
        counts = Counter()
        for seq in lista_seqs:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if "N" not in kmer:
                    counts[kmer] += 1
        return counts
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
        if np.sum(nb_mask) == 0:
            continue
        seqs_com, seqs_sem = [], []
        for gene in genes_chrom:
            s, e = _promoter_coords(gene, L)
            if (e - s) != PROMOTER_WINDOW:
                continue
            gid = getattr(gene, "name", str(gene.start))
            seq = seqs_por_gene.get(gid)
            if seq is None:
                continue
            (seqs_com if np.any(nb_mask[s:e]) else seqs_sem).append(seq)
        if len(seqs_com) < 10 or len(seqs_sem) < 10:
            continue
        counts_com = _contar_kmers(seqs_com, KMER_SIZE)
        counts_sem = _contar_kmers(seqs_sem, KMER_SIZE)
        total_com  = sum(counts_com.values()) or 1
        total_sem  = sum(counts_sem.values()) or 1
        todos_kmers = sorted(set(counts_com) | set(counts_sem))
        resultados  = []
        for kmer in todos_kmers:
            a = counts_com.get(kmer, 0)
            b = counts_sem.get(kmer, 0)
            c = total_com - a
            d = total_sem - b
            _, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            freq_com = a / total_com
            freq_sem = b / total_sem
            fold     = freq_com / freq_sem if freq_sem > 0 else np.inf
            resultados.append((kmer, a, b, round(freq_com, 6),
                               round(freq_sem, 6), round(fold, 4), round(p, 6)))
        pvals = np.array([r[6] for r in resultados])
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        csv_rows = []
        for (kmer, a, b, fc, fs, fold, p), pa in zip(resultados, padj):
            csv_rows.append([feature, kmer, a, b, fc, fs,
                             fold, round(p, 6), round(pa, 6),
                             "*" if pa < 0.05 else "ns"])
        _write_csv(f"kmer_{KMER_SIZE}mer_{feature}.csv",
                   ["feature", "kmer", "n_com", "n_sem",
                    "freq_com", "freq_sem", "fold_enrichment",
                    "p_Fisher", "p_adj_FDR", "sig"], csv_rows)
        sig_rows = sorted(
            [(r[1], r[6], r[8]) for r in csv_rows if r[6] > 1 and r[9] == "*"],
            key=lambda x: -x[1])[:20]
        if sig_rows:
            kmers_plot = [r[0] for r in sig_rows]
            folds_plot = [r[1] for r in sig_rows]
            fig, ax = plt.subplots(figsize=(10, max(5, len(kmers_plot) * 0.35)))
            ax.barh(kmers_plot[::-1], folds_plot[::-1], color="#d62728", alpha=0.85)
            ax.axvline(1.0, color="black", ls="--", lw=1)
            ax.set_xlabel("Fold enrichment (COM / SEM estrutura)")
            ax.set_title(f"K-mers ({KMER_SIZE}-mer) enriquecidos em promotores COM {feature}\n"
                         f"Cromossomo principal  |  p_adj < 0,05 (Fisher + FDR BH)")
            fig.tight_layout()
            out = os.path.join(OUT_DIR, "figuras", f"kmer_{feature}.png")
            fig.savefig(out, dpi=150)
            plt.close(fig)
            print(f"  [FIG] {out}")
        else:
            print(f"  [INFO] Nenhum k-mer significativo após FDR para {feature}")

# =============================================================================
# ANÁLISE 14 — PERFIL SOBREPOSTO (original)
# =============================================================================

def analise_perfil_sobreposto(perfis):
    print(f"\n[14] Perfil sobreposto de todas as estruturas no TSS (±{TSS_FLANK} bp)...")
    if not perfis:
        print("  [AVISO] Nenhum perfil disponível (análise 11 não retornou dados).")
        return
    x = np.arange(-TSS_FLANK, TSS_FLANK)
    features = sorted(perfis.keys())
    cmap = plt.cm.tab20(np.linspace(0, 1, len(features)))
    # Absoluto
    fig, ax = plt.subplots(figsize=(14, 7))
    for feat, cor in zip(features, cmap):
        ax.plot(x, perfis[feat], lw=1.3, alpha=0.8, color=cor,
                label=feat.replace("_", " "))
    ax.axvline(0, color="black", ls=":", lw=1.2, label="TSS")
    ax.set_xlabel("Posição relativa ao TSS (bp)  |  negativo = upstream, positivo = downstream")
    ax.set_ylabel("Nº de genes com a estrutura nessa posição")
    ax.set_title("Perfil posicional sobreposto — todas as estruturas não-B\n"
                 f"Cromossomo principal  |  janela ±{TSS_FLANK} bp  |  contagem absoluta")
    ax.legend(fontsize=6, ncol=2, loc="upper right")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "perfil_sobreposto_absoluto.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")
    # Normalizado
    fig, ax = plt.subplots(figsize=(14, 7))
    for feat, cor in zip(features, cmap):
        perfil = perfis[feat].astype(float)
        norm = perfil / perfil.max() if perfil.max() > 0 else perfil
        ax.plot(x, norm, lw=1.3, alpha=0.8, color=cor, label=feat.replace("_", " "))
    ax.axvline(0, color="black", ls=":", lw=1.2, label="TSS")
    ax.set_xlabel("Posição relativa ao TSS (bp)  |  negativo = upstream, positivo = downstream")
    ax.set_ylabel("Densidade normalizada (0–1, relativa ao pico de cada estrutura)")
    ax.set_title("Perfil posicional sobreposto e normalizado — todas as estruturas não-B\n"
                 f"Cromossomo principal  |  janela ±{TSS_FLANK} bp  |  cada curva escalada ao próprio máximo")
    ax.legend(fontsize=6, ncol=2, loc="upper right")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "perfil_sobreposto_normalizado.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# =============================================================================
# NOVAS ANÁLISES (derivadas do pipeline R)
# =============================================================================

def load_cds_intervals(annot_path, chrom):
    """Retorna dicionário: gene_id -> (cds_start, cds_end, strand) usando o atributo Parent."""
    cds_dict = {}
    for s, e, strand, attrs in _gff_to_intervals(annot_path, chrom, feature_type="CDS"):
        # Extrai gene_id de atributos como "Parent=geneX"
        for part in attrs.split(";"):
            if part.strip().startswith("Parent="):
                gene_id = part.strip().split("=", 1)[1]
                cds_dict[gene_id] = (s, e, strand)
                break
    return cds_dict

def compute_gene_utr_info(genes_main, cds_dict, chrom):
    """
    Para cada gene no cromossomo, define:
      - utr_start, utr_end (5'UTR real ou customizada)
      - utr_length
      - is_real (se CDS existe e a distância TSS->CDS > 0)
    Se não houver CDS ou TSS coincidir com CDS start, usa UTR_CUSTOM_SIZE.
    Retorna dict gene_id -> info.
    """
    info = {}
    genes_chrom = [g for g in genes_main.values() if g.chromosome == chrom]
    for gene in genes_chrom:
        strand = "+" if gene.strand else "-"
        tss = gene.start if strand == "+" else gene.end
        gene_id = getattr(gene, "name", str(gene.start))
        cds = cds_dict.get(gene_id)
        if cds:
            cds_start, cds_end, _ = cds
            if strand == "+":
                utr_len = cds_start - tss
                if utr_len > 0:
                    utr_start, utr_end = tss, cds_start
                    is_real = True
                else:
                    is_real = False
            else:  # strand "-"
                utr_len = tss - cds_end
                if utr_len > 0:
                    utr_start, utr_end = cds_end, tss
                    is_real = True
                else:
                    is_real = False
        else:
            is_real = False

        if not is_real:
            if strand == "+":
                utr_start = tss
                utr_end   = tss + UTR_CUSTOM_SIZE
            else:
                utr_start = tss - UTR_CUSTOM_SIZE
                utr_end   = tss
            utr_len = UTR_CUSTOM_SIZE

        info[gene_id] = {
            "utr_start" : max(0, utr_start),
            "utr_end"   : utr_end,
            "length"    : utr_len,
            "is_real"   : is_real,
            "strand"    : strand,
            "tss"       : tss,
            "gene_start": gene.start,
            "gene_end"  : gene.end,
        }
    return info

# --- Análise 15 (R1): Donut genes com/sem 5'UTR real ---
def analise_genes_sem_utr(genes_utr_info):
    print("\n[15] Genes com vs sem 5'UTR (donut)...")
    real = sum(1 for v in genes_utr_info.values() if v["is_real"])
    custom = len(genes_utr_info) - real
    sizes = [real, custom]
    labels = [f"Com 5'UTR real ({real})", f"Sem 5'UTR real ({custom})"]
    colors = ["#2ca02c", "#d62728"]
    fig, ax = plt.subplots()
    wedges, texts = ax.pie(sizes, labels=None, startangle=90, colors=colors,
                           wedgeprops=dict(width=0.4))
    ax.legend(wedges, labels, title="5'UTR", loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title("Proporção de genes com 5'UTR anotada")
    out = os.path.join(OUT_DIR, "figuras", "donut_genes_utr.png")
    fig.savefig(out, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# --- Análise 16 (R2): Histograma do tamanho das 5'UTRs reais ---
def analise_utr_size_hist(genes_utr_info):
    print("\n[16] Distribuição do tamanho das 5'UTRs (histograma)...")
    real_lengths = [v["length"] for v in genes_utr_info.values() if v["is_real"]]
    if not real_lengths:
        print("  [AVISO] Nenhuma 5'UTR real encontrada.")
        return
    med = np.median(real_lengths)
    media = np.mean(real_lengths)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(real_lengths, bins=30, color="#9467bd", alpha=0.85, edgecolor="black", lw=0.3)
    ax.axvline(media, color="#d62728", ls="--", lw=2, label=f"Média = {media:.0f} bp")
    ax.axvline(med, color="#1f77b4", ls="--", lw=2, label=f"Mediana = {med:.0f} bp")
    ax.set_xlabel("Tamanho da 5'UTR (bp)")
    ax.set_ylabel("Nº de genes")
    ax.set_title("Distribuição do tamanho das 5'UTRs (anotação real)")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "hist_tamanho_utr.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# --- Análise 17 (R4): Distribuição posicional detalhada (empilhado por classe) ---
def analise_posicional_gene_regions_detailed(genome, genes_main, genes_utr_info, non_b_files):
    print("\n[17] Distribuição posicional non-B por região gênica (detalhada)...")
    chrom = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if not chr_obj:
        return
    L = chr_obj.length
    # Para cada gene, construir a lista de bins: Promotor, 5'UTR, 10 deciles do corpo gênico
    bin_names = ["Promotor"] + ["5'UTR"] + [f"Decil_{i+1}" for i in range(10)]
    # Para cada ferramenta, para cada classe, somar a densidade (bp/kb) sobre os genes
    # Resultado: {classe: {bin_name: soma_densidade}}  e depois dividir pelo número de genes com aquele bin válido
    # Armazenaremos contagens por ferramenta, depois agregaremos por classe.
    class_bin_sum = defaultdict(lambda: np.zeros(len(bin_names)))
    class_bin_counts = defaultdict(lambda: np.zeros(len(bin_names)))  # número de genes que contribuíram
    # Iremos processar cada ferramenta
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        classe = _structure_class(feature)
        nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
        if np.sum(nb_mask) == 0:
            continue
        for gene_id, info in genes_utr_info.items():
            strand = info["strand"]
            tss = info["tss"]
            # Promotor
            prom_s, prom_e = _promoter_coords(
                type('Gene', (object,), {'start': info['gene_start'], 'end': info['gene_end'], 'strand': strand=='+'})(),
                L)
            # corrigir para usar o objeto gene diretamente, mas podemos calcular
            # Vamos calcular com base no gene
            if strand == "+":
                prom_s = max(0, tss - PROMOTER_WINDOW)
                prom_e = tss
            else:
                prom_s = tss
                prom_e = min(L, tss + PROMOTER_WINDOW)
            # 5'UTR
            utr_s, utr_e = info["utr_start"], info["utr_end"]
            # Corpo gênico: após UTR até gene end (+), ou gene start até antes da UTR (-)
            if strand == "+":
                corpo_s = utr_e
                corpo_e = info["gene_end"]
            else:
                corpo_s = info["gene_start"]
                corpo_e = utr_s
            # Se corpo não positivo, pula
            if corpo_s >= corpo_e:
                continue
            corpo_len = corpo_e - corpo_s
            decil_len = corpo_len / 10.0
            # Coletar overlap em bp
            bp_prom = np.sum(nb_mask[prom_s:prom_e])
            bp_utr  = np.sum(nb_mask[utr_s:utr_e])
            # Densidades por kb
            dens_prom = bp_prom / ((prom_e - prom_s) / 1000) if (prom_e - prom_s) > 0 else 0
            dens_utr  = bp_utr  / ((utr_e - utr_s) / 1000) if (utr_e - utr_s) > 0 else 0
            # Acumular promotor e UTR
            class_bin_sum[classe][0] += dens_prom
            class_bin_counts[classe][0] += 1
            class_bin_sum[classe][1] += dens_utr
            class_bin_counts[classe][1] += 1
            # Deciles
            for i in range(10):
                d_start = int(corpo_s + i * decil_len)
                d_end   = int(corpo_s + (i+1) * decil_len) if i < 9 else corpo_e
                if d_end <= d_start:
                    continue
                bp_decil = np.sum(nb_mask[d_start:d_end])
                dens_decil = bp_decil / ((d_end - d_start) / 1000) if (d_end - d_start) > 0 else 0
                class_bin_sum[classe][2 + i] += dens_decil
                class_bin_counts[classe][2 + i] += 1
    # Média por classe
    class_dens_media = {}
    for classe in class_bin_sum:
        media = np.divide(class_bin_sum[classe], class_bin_counts[classe],
                          out=np.zeros_like(class_bin_sum[classe]),
                          where=class_bin_counts[classe]!=0)
        class_dens_media[classe] = media
    # Plot stacked bar detalhado (uma barra empilhada por classe)
    if not class_dens_media:
        print("  [AVISO] Sem dados para distribuição posicional detalhada.")
        return
    classes_order = sorted(class_dens_media.keys())
    x = np.arange(len(bin_names))
    width = 0.6
    bottoms = np.zeros(len(bin_names))
    fig, ax = plt.subplots(figsize=(14, 7))
    for classe in classes_order:
        vals = class_dens_media[classe]
        ax.bar(x, vals, bottom=bottoms, label=classe,
               color=CLASS_COLORS.get(classe, "#7f7f7f"), alpha=0.85)
        bottoms += vals
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Densidade média non‑B (bp/kb)")
    ax.set_title("Distribuição posicional das estruturas não‑B por região gênica\n"
                 "Média entre ferramentas de cada classe  |  Cromossomo principal")
    ax.legend(title="Classe", fontsize=8)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "posicional_detalhado.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# --- Análise 18 (R5): Distribuição posicional agrupada (Promotor / 5'UTR / Corpo) ---
def analise_posicional_gene_regions_aggregated(genome, genes_main, genes_utr_info, non_b_files):
    print("\n[18] Distribuição posicional non‑B (agrupada: Promotor/5'UTR/Corpo)...")
    chrom = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if not chr_obj:
        return
    L = chr_obj.length
    bin_names = ["Promotor", "5'UTR", "Corpo gênico"]
    class_bin_sum = defaultdict(lambda: np.zeros(3))
    class_bin_counts = defaultdict(lambda: np.zeros(3))
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        classe = _structure_class(feature)
        nb_mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, L)
        if np.sum(nb_mask) == 0:
            continue
        for gene_id, info in genes_utr_info.items():
            strand = info["strand"]
            tss = info["tss"]
            if strand == "+":
                prom_s = max(0, tss - PROMOTER_WINDOW)
                prom_e = tss
            else:
                prom_s = tss
                prom_e = min(L, tss + PROMOTER_WINDOW)
            utr_s, utr_e = info["utr_start"], info["utr_end"]
            if strand == "+":
                corpo_s, corpo_e = utr_e, info["gene_end"]
            else:
                corpo_s, corpo_e = info["gene_start"], utr_s
            if corpo_s >= corpo_e:
                continue
            # Promotor
            bp_prom = np.sum(nb_mask[prom_s:prom_e])
            dens_prom = bp_prom / ((prom_e - prom_s) / 1000) if (prom_e - prom_s) > 0 else 0
            class_bin_sum[classe][0] += dens_prom
            class_bin_counts[classe][0] += 1
            # UTR
            bp_utr = np.sum(nb_mask[utr_s:utr_e])
            dens_utr = bp_utr / ((utr_e - utr_s) / 1000) if (utr_e - utr_s) > 0 else 0
            class_bin_sum[classe][1] += dens_utr
            class_bin_counts[classe][1] += 1
            # Corpo (soma de todos os deciles)
            corpo_len = corpo_e - corpo_s
            bp_corpo = np.sum(nb_mask[corpo_s:corpo_e])
            dens_corpo = bp_corpo / (corpo_len / 1000) if corpo_len > 0 else 0
            class_bin_sum[classe][2] += dens_corpo
            class_bin_counts[classe][2] += 1
    class_dens_media = {}
    for classe in class_bin_sum:
        media = np.divide(class_bin_sum[classe], class_bin_counts[classe],
                          out=np.zeros_like(class_bin_sum[classe]),
                          where=class_bin_counts[classe]!=0)
        class_dens_media[classe] = media
    if not class_dens_media:
        return
    classes_order = sorted(class_dens_media.keys())
    x = np.arange(3)
    bottoms = np.zeros(3)
    fig, ax = plt.subplots(figsize=(8, 6))
    for classe in classes_order:
        vals = class_dens_media[classe]
        ax.bar(x, vals, bottom=bottoms, label=classe,
               color=CLASS_COLORS.get(classe, "#7f7f7f"), alpha=0.85)
        bottoms += vals
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names)
    ax.set_ylabel("Densidade média non‑B (bp/kb)")
    ax.set_title("Distribuição posicional agregada das estruturas não‑B\n"
                 "Média entre ferramentas por classe  |  Cromossomo principal")
    ax.legend(title="Classe", fontsize=8)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "posicional_agregado.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# --- Análise 19 (R6): Contagem de predições por ferramenta (log) ---
def analise_prediction_counts(non_b_files):
    print("\n[19] Contagem de predições por ferramenta (barras log)...")
    features = []
    n_intervals = []
    for bed_file in non_b_files:
        fname = bed_file.replace(".bed", "")
        # Contar intervalos em todos os cromossomos
        total = 0
        for chrom in CHROMS:
            intervals = _bed_to_intervals(os.path.join(SITES_DIR, bed_file), chrom)
            total += len(intervals)
        features.append(fname)
        n_intervals.append(total)
    x = np.arange(len(features))
    classes_feat = [_structure_class(f) for f in features]
    cores = [CLASS_COLORS.get(c, "#7f7f7f") for c in classes_feat]
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.bar(x, n_intervals, color=cores, alpha=0.85)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("_", "\n") for f in features], fontsize=7)
    ax.set_ylabel("Número de intervalos preditos (log)")
    ax.set_title("Contagem de predições por ferramenta")
    classes_presentes = sorted(set(classes_feat))
    handles = [Patch(color=CLASS_COLORS.get(c, "#7f7f7f"), alpha=0.85, label=c)
               for c in classes_presentes]
    ax.legend(handles=handles, fontsize=8, title="Classe")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "prediction_counts.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# --- Análise 20 (R9): Conteúdo GC por classe ---
def analise_gc_content(genome, non_b_files):
    print("\n[20] Conteúdo GC dos intervalos non‑B por classe...")
    chrom = "NC_002607.1"
    chr_obj = genome.genome_dict.get(chrom)
    if not chr_obj:
        return

    # Obtém a sequência do cromossomo (tenta .sequence, depois .seq)
    try:
        seq = str(chr_obj.sequence).upper()
    except AttributeError:
        try:
            seq = str(chr_obj.seq).upper()
        except AttributeError:
            print("  [AVISO] Não foi possível obter a sequência do cromossomo. Pulando análise 20.")
            return

    L = len(seq)
    gc_by_class = defaultdict(list)
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        classe = _structure_class(feature)
        intervals = _bed_to_intervals(os.path.join(SITES_DIR, bed_file), chrom)
        for s, e in intervals:
            if s < 0 or e > L:
                continue
            subseq = seq[s:e]
            if len(subseq) == 0:
                continue
            gc = (subseq.count("G") + subseq.count("C")) / len(subseq) * 100
            gc_by_class[classe].append(gc)
    if not gc_by_class:
        return
    classes_order = sorted(gc_by_class.keys())
    data = [gc_by_class[c] for c in classes_order]
    fig, ax = plt.subplots(figsize=(10, 6))
    vp = ax.violinplot(data, positions=range(len(classes_order)),
                       showmeans=True, showmedians=True)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(CLASS_COLORS.get(classes_order[i], "#7f7f7f"))
        body.set_alpha(0.7)
    ax.boxplot(data, positions=range(len(classes_order)),
               widths=0.2, showfliers=False, manage_ticks=False)
    ax.set_xticks(range(len(classes_order)))
    ax.set_xticklabels(classes_order, fontsize=9)
    ax.set_ylabel("Conteúdo GC (%)")
    ax.set_title("Conteúdo GC dos intervalos non‑B por classe estrutural\nCromossomo principal")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "gc_por_classe.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# --- Análise 21 (R10): Tamanho dos intervalos por classe ---
def analise_interval_sizes(non_b_files):
    print("\n[21] Tamanho dos intervalos non‑B por classe (log)...")
    chrom = "NC_002607.1"
    sizes_by_class = defaultdict(list)
    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        classe = _structure_class(feature)
        intervals = _bed_to_intervals(os.path.join(SITES_DIR, bed_file), chrom)
        for s, e in intervals:
            sizes_by_class[classe].append(e - s)
    if not sizes_by_class:
        return
    classes_order = sorted(sizes_by_class.keys())
    data = [np.log10(np.array(sizes_by_class[c]) + 1) for c in classes_order]  # log10
    fig, ax = plt.subplots(figsize=(10, 6))
    vp = ax.violinplot(data, positions=range(len(classes_order)),
                       showmeans=True, showmedians=True)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(CLASS_COLORS.get(classes_order[i], "#7f7f7f"))
        body.set_alpha(0.7)
    ax.boxplot(data, positions=range(len(classes_order)),
               widths=0.2, showfliers=False, manage_ticks=False)
    ax.set_xticks(range(len(classes_order)))
    ax.set_xticklabels(classes_order, fontsize=9)
    ax.set_ylabel("log10(Tamanho do intervalo + 1)")
    ax.set_title("Distribuição dos tamanhos dos intervalos non‑B por classe\nCromossomo principal")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "figuras", "tamanho_intervalos.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  [FIG] {out}")

# =============================================================================
# MAIN (atualizado)
# =============================================================================

def main():
    _make_dirs()
    print("=" * 70)
    print("ANÁLISES DE TFB × ESTRUTURAS NÃO‑B — Halobacterium salinarum NRC-1")
    print("=" * 70)

    print("\nCarregando genoma e anotação...")
    g = Genome.Genome("Hsalinarum")
    g.load_seq()
    g.load_annotation(annot_file="sequence.gff3", features=["gene", "CDS"])  # CDS incluído
    genes_main = dict(g.genes)
    non_b_files = list_non_b_files()
    print(f"  Genes carregados : {len(genes_main)}")
    print(f"  Features não‑B   : {len(non_b_files)}")

    # Carregar CDS para análise de 5'UTR
    chrom_principal = "NC_002607.1"
    cds_dict = load_cds_intervals(os.path.join(BASE, "annotation/sequence.gff3"), chrom_principal)
    genes_utr_info = compute_gene_utr_info(genes_main, cds_dict, chrom_principal)

    # Análises originais (com modificações)
    analise_cobertura_genomica(g, non_b_files)            # 1
    analise_promotores_tfb(g, genes_main, non_b_files)    # 2
    # (análise 3 removida)
    analise_nTFB_por_gene(g, genes_main, non_b_files)     # 4 (normalizada)
    analise_expressao_genica(g, genes_main, non_b_files)  # 5
    analise_orientacao_genes(g, genes_main, non_b_files)  # 6 (tandem)
    analise_proximidade_IS(g, non_b_files)                # 7
    analise_naoB_perto_TFB(g, non_b_files)                # 8
    analise_asRNA_estruturas(g, non_b_files)              # 9
    analise_hic_naoB(g, non_b_files)                      # 10 lollipop
    perfis = analise_perfil_posicional(g, genes_main, non_b_files)  # 11
    analise_distancia_TFB(g, non_b_files)                 # 12
    analise_enriquecimento_kmer(g, genes_main, non_b_files)  # 13
    analise_perfil_sobreposto(perfis)                     # 14

    # Novas análises (R)
    analise_genes_sem_utr(genes_utr_info)                 # 15
    analise_utr_size_hist(genes_utr_info)                 # 16
    analise_posicional_gene_regions_detailed(g, genes_main, genes_utr_info, non_b_files)  # 17
    analise_posicional_gene_regions_aggregated(g, genes_main, genes_utr_info, non_b_files) # 18
    analise_prediction_counts(non_b_files)                # 19
    analise_gc_content(g, non_b_files)                    # 20
    analise_interval_sizes(non_b_files)                   # 21

    print("\n" + "=" * 70)
    print("CONCLUÍDO")
    print(f"  Tabelas → {OUT_DIR}/tabelas/")
    print(f"  Figuras → {OUT_DIR}/figuras/")
    print("=" * 70)

if __name__ == "__main__":
    main()
