"""
Pipeline integrado de análise genômica — Halobacterium salinarum NRC-1 [v3]
============================================================================

Experimentos disponíveis:
  A) Cinética temporal : total-RNA TP1 → TP4  (condição padrão 4.3 M NaCl)
  B) Gradiente de sal  : 2.6 M / 4.3 M / 5.1 M NaCl  (rep A + rep B, médias)

Análises realizadas:
  1. RPKM por gene em cada condição/timepoint
  2. Classificação por decis de expressão (top 40% HETUs / bottom 60% LETUs)
     seguindo o critério validado no paper do GRATIOSA (Pineau et al. 2025)
  3. Metagene estratificado por nível de expressão ao redor do TSS
     — perfil de cada estrutura de DNA não-B
     — IC 95% por bin via bootstrap
     — separado por condição de referência, por replicon e por feature
  4. Testes estatísticos formais
     — Proporção (z-score): % promotores com estrutura em HETUs vs LETUs
     — Hipergeométrico (FDR): enriquecimento de estruturas em promotores de genes UP
     — Wilcoxon: comparação de densidade de estruturas entre classes de expressão
  5. Orientação gênica (divergente / convergente / tandem)
     — densidade de cada estrutura não-B por classe de orientação
  6. Cinética temporal (TP1→TP4)
     — metagene estratificado por direção de mudança (UP / DOWN / estável)
  7. Gradiente de sal (2.6 → 4.3 → 5.1 M)
     — metagene estratificado por resposta ao estresse osmótico
  8. RNA ribossomal vs total (análise separada dos genes de RNA estável)

Saídas:
  resultados_v3/
  ├── tabelas/          → CSVs com RPKM, classificações, testes estatísticos
  ├── metagene/         → PNGs por feature × condição × replicon
  ├── orientacao/       → PNGs por feature × classe de orientação
  ├── temporal/         → PNGs cinética TP1→TP4
  ├── salinidade/       → PNGs gradiente de sal
  └── resumo/           → TXTs por replicon

Referência metodológica:
  Pineau et al. (2025) Nucleic Acids Research — GRATIOSA package
  doi: 10.1093/nar/gkaf452
"""

import os
import csv
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import concurrent.futures
from GRATIOSA import Genome

warnings.filterwarnings("ignore")

# =============================================================================
# Paleta e estilo visual
# — cores distintas por decil, do azul escuro (mais expresso) ao vermelho
# — IC 95% como faixa sombreada
# — HETUs (top 40%) em tons quentes, LETUs (bottom 60%) em cinza-azulado
# =============================================================================

DECIL_COLORS = [
    "#d62728",  # decil 10 — mais expresso
    "#e87c1e",
    "#f0a500",
    "#8fbc44",  # decil 7
    "#2ca02c",
    "#17a0a0",
    "#1f77b4",
    "#5a6fa8",
    "#7b5ea7",
    "#9467bd",  # decil 1 — menos expresso
]
HETU_COLOR  = "#d62728"   # vermelho — HETUs
LETU_COLOR  = "#5a7ab5"   # azul — LETUs
UP_COLOR    = "#d62728"
DOWN_COLOR  = "#1f77b4"
STABLE_COLOR= "#888888"
RRNA_COLOR  = "#2ca02c"

REPLICON_LABELS = {
    "NC_002607.1": "Cromossomo principal (2.0 Mb)",
    "NC_001869.1": "Plasmídeo pNRC100 (191 kb)",
    "NC_002608.1": "Plasmídeo pNRC200 (365 kb)",
}

plt.rcParams.update({
    "font.family"      : "DejaVu Sans",
    "font.size"        : 10,
    "axes.titlesize"   : 11,
    "axes.titleweight" : "bold",
    "axes.labelsize"   : 10,
    "legend.fontsize"  : 8,
    "figure.dpi"       : 150,
    "axes.spines.top"  : False,
    "axes.spines.right": False,
})

# =============================================================================
# Configuração
# =============================================================================

BASE      = "/home/mcavalcante/codigo/Diario-de-bordo/GRATIOSA_DB/data/Hsalinarum/"
CHROMS    = ["NC_002607.1", "NC_001869.1", "NC_002608.1"]
OUT_DIR   = "resultados_v3"
SITES_DIR = os.path.join(BASE, "sites/")

# Janelas
PROMOTER_WINDOW = 250   # bp upstream do TSS para análise de promotores
METAGENE_WINDOW = 1500  # bp em cada lado do TSS para metagene
BIN_SIZE        = 50    # resolução do metagene em bp
N_BOOTSTRAP     = 200   # iterações para IC 95%

# Limiar de expressão (paper Pineau et al. 2025)
HETU_QUANTILE = 0.60    # top 40% = HETUs, bottom 60% = LETUs

# Condições experimentais
# — Cinética temporal
TIMEPOINTS = ["total-RNA-TP1", "total-RNA-TP2", "total-RNA-TP3", "total-RNA-TP4"]
TP_LABELS  = {"total-RNA-TP1": "TP1", "total-RNA-TP2": "TP2",
              "total-RNA-TP3": "TP3", "total-RNA-TP4": "TP4"}

# — Gradiente de sal: SRR IDs mapeados para concentração
SAL_CONDITIONS = {
    "2.6M": {"rep_A": "SRR12971714", "rep_B": "SRR12971717"},
    "4.3M": {"rep_A": "SRR12971715", "rep_B": "SRR12971718"},
    "5.1M": {"rep_A": "SRR12971716", "rep_B": "SRR12971719"},
}
SAL_LABELS = {"2.6M": "2,6 M NaCl", "4.3M": "4,3 M NaCl", "5.1M": "5,1 M NaCl"}

# — RNA ribossomal
RIBO_TIMEPOINTS = ["ribosomal_RNA_TP1", "ribosomal_RNA_TP2",
                   "ribosomal_RNA_TP3",  "ribosomal_RNA_TP4"]

RNA_DIR = os.path.join(BASE, "rnaseq_cov/")

# Thresholds FC para cinética temporal
FC_UP   = 2.0
FC_DOWN = 0.5

# Paralelização
MAX_WORKERS = 12


# =============================================================================
# Utilitários de I/O
# =============================================================================

def _make_dirs():
    for sub in ("tabelas", "metagene", "orientacao", "temporal", "salinidade", "resumo"):
        os.makedirs(os.path.join(OUT_DIR, sub), exist_ok=True)


def _iter_tsv(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip().split("\t")


def _bedgraph_to_array(path, chrom, length):
    """Lê bedgraph e retorna array numpy de cobertura por posição."""
    arr = np.zeros(length)
    if not os.path.exists(path):
        return arr
    for cols in _iter_tsv(path):
        if cols[0] == chrom:
            arr[int(cols[1]):int(cols[2])] = float(cols[3])
    return arr


def _bed_to_mask(path, chrom, length):
    """Lê BED e retorna máscara binária numpy."""
    mask = np.zeros(length, dtype=np.int8)
    if not os.path.exists(path):
        return mask
    for cols in _iter_tsv(path):
        if cols[0].startswith("track"):
            continue
        if cols[0] == chrom:
            mask[int(cols[1]):int(cols[2])] = 1
    return mask


def _write_csv(subdir, filename, header, rows):
    path = os.path.join(OUT_DIR, subdir, filename)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  [CSV] {path}")
    return path


def _strand_str(strand):
    if isinstance(strand, bool):
        return "+" if strand else "-"
    return strand


# =============================================================================
# Carregamento do genoma e anotações
# =============================================================================

def load_genome():
    print("Carregando genoma e anotações...")
    g = Genome.Genome("Hsalinarum")
    g.load_seq()
    g.load_annotation(annot_file="sequence.gff3", features=["gene"])
    genes_main = dict(g.genes)
    g.load_annotation(annot_file="asRNAs.gff3", features=["ncRNA"])
    genes_asrna = dict(g.genes)
    print(f"  Genes codificantes : {len(genes_main)}")
    print(f"  asRNAs             : {len(genes_asrna)}")
    return g, genes_main, genes_asrna


def load_tss():
    """Carrega TSSs do GTF por replicon."""
    result = {c: [] for c in CHROMS}
    path = os.path.join(BASE, "annotation/Hsalinarum_TSS.gtf")
    if not os.path.exists(path):
        print("  [AVISO] TSS GTF não encontrado.")
        return result
    for cols in _iter_tsv(path):
        chrom, strand = cols[0], cols[6]
        if chrom not in CHROMS:
            continue
        pos = int(cols[3]) if strand == "+" else int(cols[4])
        result[chrom].append({"pos": pos, "strand": strand})
    return result


def load_hic():
    result = {c: [] for c in CHROMS}
    path = os.path.join(BASE, "hic/hi-c.bedpe")
    if not os.path.exists(path):
        return result
    for cols in _iter_tsv(path):
        chrom = cols[0]
        if chrom not in CHROMS or cols[3] != chrom:
            continue
        result[chrom].append({
            "s1": int(cols[1]), "e1": int(cols[2]),
            "s2": int(cols[4]), "e2": int(cols[5]),
        })
    return result


def load_chipseq():
    result = {c: [] for c in CHROMS}
    path = os.path.join(BASE, "chipseq/TFBs.bed")
    if not os.path.exists(path):
        return result
    for cols in _iter_tsv(path):
        if cols[0] in CHROMS:
            result[cols[0]].append((int(cols[1]), int(cols[2])))
    return result


def list_non_b_files():
    if not os.path.exists(SITES_DIR):
        return []
    return sorted(f for f in os.listdir(SITES_DIR) if f.endswith(".bed"))


# =============================================================================
# Cálculo de RPKM
# — RPKM = (reads no gene) / (tamanho_gene_kb × total_reads_milhões)
# — Para cada condição, soma fwd e rev conforme a fita do gene
# =============================================================================

def calcular_rpkm_por_gene(genes, chrom, chr_length, cond_name,
                            cov_fwd=None, cov_rev=None):
    """
    Calcula RPKM para todos os genes de um replicon em uma condição.
    Se cov_fwd/cov_rev forem fornecidos (arrays pré-carregados), usa-os
    diretamente em vez de ler os bedgraphs — evita I/O repetido.
    Retorna dict {gene_id: rpkm}.
    """
    if cov_fwd is None:
        fwd_path = os.path.join(RNA_DIR, f"{cond_name}-fwd.bedgraph")
        rev_path = os.path.join(RNA_DIR, f"{cond_name}-rev.bedgraph")
        cov_fwd = _bedgraph_to_array(fwd_path, chrom, chr_length)
        cov_rev = _bedgraph_to_array(rev_path, chrom, chr_length)

    total_reads = (np.sum(cov_fwd) + np.sum(cov_rev)) / 1e6
    if total_reads == 0:
        return {}

    rpkm_dict = {}
    for gene_id, gene in genes.items():
        if gene.chromosome != chrom:
            continue
        strand = _strand_str(gene.strand)
        cov = cov_fwd if strand == "+" else cov_rev
        gene_length_kb = max(gene.end - gene.start, 1) / 1000
        reads = np.sum(cov[gene.start:gene.end])
        rpkm_dict[gene_id] = reads / (gene_length_kb * total_reads)

    return rpkm_dict


def calcular_rpkm_media(genes, chrom, chr_length, cond_A, cond_B):
    """
    Calcula RPKM médio de duas replicatas (rep A e rep B) de uma condição de sal.
    Lê cada bedgraph uma única vez.
    """
    # Carrega os 4 arrays de uma vez
    fwd_A = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_A}-fwd.bedgraph"), chrom, chr_length)
    rev_A = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_A}-rev.bedgraph"), chrom, chr_length)
    fwd_B = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_B}-fwd.bedgraph"), chrom, chr_length)
    rev_B = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_B}-rev.bedgraph"), chrom, chr_length)

    # Média dos arrays antes de calcular RPKM (equivalente a média de replicatas)
    cov_fwd = (fwd_A + fwd_B) / 2
    cov_rev = (rev_A + rev_B) / 2

    return calcular_rpkm_por_gene(genes, chrom, chr_length, None,
                                   cov_fwd=cov_fwd, cov_rev=cov_rev)


# =============================================================================
# Classificação por decis e HETUs/LETUs
# =============================================================================

def classificar_por_decis(rpkm_dict):
    """
    Classifica genes em 10 decis de expressão.
    Retorna dict {gene_id: decil} (decil 10 = mais expresso).
    """
    if not rpkm_dict:
        return {}
    ids   = list(rpkm_dict.keys())
    vals  = np.array([rpkm_dict[i] for i in ids])
    decis = pd.qcut(vals, 10, labels=False, duplicates="drop")
    return {ids[i]: int(decis[i]) + 1 for i in range(len(ids))
            if not np.isnan(decis[i])}


def classificar_hetu_letu(rpkm_dict, quantile=HETU_QUANTILE):
    """
    Classifica genes em HETUs (top 40%) e LETUs (bottom 60%)
    conforme critério de Pineau et al. (2025).
    Retorna dict {gene_id: 'HETU' ou 'LETU'}.
    """
    if not rpkm_dict:
        return {}
    limiar = np.quantile(list(rpkm_dict.values()), quantile)
    return {gid: ("HETU" if v >= limiar else "LETU")
            for gid, v in rpkm_dict.items()}


def _worker_rpkm_condicao(args):
    """
    Worker paralelo: calcula RPKM de todos os genes para UMA condição × replicon.
    Lê os bedgraphs uma única vez e processa todos os genes do replicon.
    Retorna lista de (gene_id, cond_nome, rpkm).
    """
    genes_data, chrom, chr_length, cond_nome, cond_val, tipo = args

    if tipo == "single":
        fwd = _bedgraph_to_array(
            os.path.join(RNA_DIR, f"{cond_val}-fwd.bedgraph"), chrom, chr_length)
        rev = _bedgraph_to_array(
            os.path.join(RNA_DIR, f"{cond_val}-rev.bedgraph"), chrom, chr_length)
    else:
        cond_A, cond_B = cond_val
        fwd_A = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_A}-fwd.bedgraph"), chrom, chr_length)
        rev_A = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_A}-rev.bedgraph"), chrom, chr_length)
        fwd_B = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_B}-fwd.bedgraph"), chrom, chr_length)
        rev_B = _bedgraph_to_array(os.path.join(RNA_DIR, f"{cond_B}-rev.bedgraph"), chrom, chr_length)
        fwd = (fwd_A + fwd_B) / 2
        rev = (rev_A + rev_B) / 2

    total_reads = (np.sum(fwd) + np.sum(rev)) / 1e6
    resultados = []
    if total_reads == 0:
        return resultados

    for gene_id, (start, end, strand) in genes_data.items():
        cov = fwd if strand == "+" else rev
        gene_length_kb = max(end - start, 1) / 1000
        reads = np.sum(cov[start:end])
        rpkm = reads / (gene_length_kb * total_reads)
        resultados.append((gene_id, cond_nome, rpkm))

    return resultados


# =============================================================================
# Análise 1 — RPKM e classificação: exporta CSV completo
# =============================================================================

def analise_rpkm_completa(genes, genome, todas_condicoes):
    """
    Calcula RPKM para todos os genes em todas as condições.
    — Cada bedgraph é lido UMA única vez por condição × replicon
    — Paralelizado com MAX_WORKERS processos
    todas_condicoes: dict {nome_condicao: (cond_str_ou_tuple, tipo)}
      tipo = 'single' ou 'media'
    Retorna DataFrame genes × condições com RPKM.
    """
    print(f"\n[1] Calculando RPKM por gene em todas as condições "
          f"({MAX_WORKERS} workers)...")

    # Monta estrutura leve de genes (serializável para multiprocessing)
    # genes_data[chrom] = {gene_id: (start, end, strand)}
    genes_data_por_chrom = {c: {} for c in CHROMS}
    gene_meta = {}  # gene_id → {chrom, start, end, strand, length_bp}
    for gene_id, gene in genes.items():
        chrom = gene.chromosome
        if chrom not in CHROMS:
            continue
        strand = _strand_str(gene.strand)
        genes_data_por_chrom[chrom][gene_id] = (gene.start, gene.end, strand)
        gene_meta[gene_id] = {
            "gene_id"  : gene_id,
            "chrom"    : chrom,
            "start"    : gene.start,
            "end"      : gene.end,
            "strand"   : strand,
            "length_bp": gene.end - gene.start,
        }

    # Monta lista de tarefas: uma por (condição × replicon)
    tarefas = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        gdata = genes_data_por_chrom[chrom]
        for cond_nome, (cond_val, tipo) in todas_condicoes.items():
            tarefas.append((gdata, chrom, chr_obj.length,
                            cond_nome, cond_val, tipo))

    print(f"  Tarefas (condição × replicon) : {len(tarefas)}")

    # Executa em paralelo
    # — usa ProcessPoolExecutor para paralelismo real (CPU-bound + I/O)
    rpkm_acumulado = {}  # {gene_id: {cond_nome: rpkm}}
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(_worker_rpkm_condicao, t): t for t in tarefas}
        concluidos = 0
        for future in concurrent.futures.as_completed(futures):
            resultados = future.result()
            for gene_id, cond_nome, rpkm in resultados:
                if gene_id not in rpkm_acumulado:
                    rpkm_acumulado[gene_id] = {}
                rpkm_acumulado[gene_id][cond_nome] = rpkm
            concluidos += 1
            if concluidos % 10 == 0 or concluidos == len(tarefas):
                print(f"  [{concluidos:3d}/{len(tarefas)}] tarefas concluídas",
                      end="\r")

    print()  # quebra de linha após o \r

    # Monta DataFrame
    cond_names = list(todas_condicoes.keys())
    rows = []
    for gene_id, meta in gene_meta.items():
        row = dict(meta)
        for cond_nome in cond_names:
            row[f"RPKM_{cond_nome}"] = rpkm_acumulado.get(
                gene_id, {}).get(cond_nome, 0.0)
        rows.append(row)

    df = pd.DataFrame(rows)
    _write_csv("tabelas", "rpkm_todos_genes.csv",
               list(df.columns), df.values.tolist())
    print(f"  Genes processados : {len(df)}")
    print(f"  Condições         : {len(cond_names)}")
    return df


# =============================================================================
# Análise 2 — Metagene estratificado por decis e HETU/LETU
# — Perfil de frequência de estruturas não-B ao redor do TSS
# — IC 95% via bootstrap
# — Separado por condição de referência, feature e replicon
# =============================================================================

def _bootstrap_ic95(matrix, n_boot=N_BOOTSTRAP):
    """
    Calcula média e IC 95% por bin via bootstrap.
    matrix: (n_genes × n_bins)
    Retorna (mean, lower, upper).
    """
    if matrix.shape[0] == 0:
        z = np.zeros(matrix.shape[1])
        return z, z, z
    means = []
    n = matrix.shape[0]
    for _ in range(n_boot):
        idx = np.random.randint(0, n, size=n)
        means.append(np.mean(matrix[idx], axis=0))
    means = np.array(means)
    mean  = np.mean(matrix, axis=0)
    lower = np.percentile(means, 2.5,  axis=0)
    upper = np.percentile(means, 97.5, axis=0)
    return mean, lower, upper


def calcular_perfil_metagene(tss_list, mask, chr_length, gene_ids_set=None,
                              tss_gene_map=None):
    """
    Calcula perfis de metagene por bin para um conjunto de TSSs.
    Retorna matrix (n_tss × n_bins) com densidade de estrutura por bin.

    tss_gene_map: dict {(pos, strand): gene_id} para filtrar por conjunto de genes.
    gene_ids_set: conjunto de gene_ids a incluir (None = todos).
    """
    n_bins = (METAGENE_WINDOW * 2) // BIN_SIZE
    perfis = []

    for tss in tss_list:
        if gene_ids_set is not None and tss_gene_map is not None:
            key = (tss["pos"], tss["strand"])
            gid = tss_gene_map.get(key)
            if gid not in gene_ids_set:
                continue

        lo = tss["pos"] - METAGENE_WINDOW
        hi = tss["pos"] + METAGENE_WINDOW
        if lo < 0 or hi >= chr_length:
            continue

        regiao = mask[lo:hi].astype(float)
        if tss["strand"] == "-":
            regiao = regiao[::-1]

        # Binagem por soma
        perfil_bins = regiao.reshape(n_bins, BIN_SIZE).sum(axis=1)
        perfis.append(perfil_bins)

    return np.array(perfis) if perfis else np.zeros((0, n_bins))


def plotar_metagene_decis(perfis_por_decil, feature_name, chrom, cond_label,
                           output_path, n_tss_total):
    """
    Plota metagene estratificado por decis de expressão com IC 95%.
    Segue o estilo da Figura 2 do paper GRATIOSA.
    """
    x = np.arange(-METAGENE_WINDOW + BIN_SIZE // 2,
                   METAGENE_WINDOW,
                   BIN_SIZE)

    fig, ax = plt.subplots(figsize=(10, 5))

    for decil in sorted(perfis_por_decil.keys()):
        matrix = perfis_por_decil[decil]
        if matrix.shape[0] == 0:
            continue
        mean, lower, upper = _bootstrap_ic95(matrix)
        color = DECIL_COLORS[decil - 1]
        label = f"Decil {decil:2d}  (n={matrix.shape[0]})"
        ax.plot(x, mean, color=color, lw=1.8, label=label)
        ax.fill_between(x, lower, upper, color=color, alpha=0.15)

    ax.axvline(0, color="black", linestyle="--", lw=1.2, zorder=5,
               label="TSS (posição 0)")
    ax.axhline(0, color="gray", linestyle=":", lw=0.8)

    # Anotações de contexto
    ax.axvspan(-PROMOTER_WINDOW, 0, alpha=0.04, color="gold",
               label=f"Janela promotor (–{PROMOTER_WINDOW} bp)")

    ax.set_xlabel(f"Distância do TSS (bp)  |  bins de {BIN_SIZE} bp",
                  fontsize=10)
    ax.set_ylabel("Frequência de estrutura por bin\n(soma de bp cobertos)",
                  fontsize=10)
    ax.set_title(
        f"{feature_name}  ×  {REPLICON_LABELS.get(chrom, chrom)}\n"
        f"Condição de referência: {cond_label}  |  n = {n_tss_total} TSSs incluídos",
        fontsize=11, pad=12
    )

    # Legenda em duas colunas (decis + elementos especiais)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, ncol=2, fontsize=7.5,
              title="Decil de expressão (10 = mais expresso)",
              title_fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)

    # Nota metodológica no rodapé
    fig.text(0.01, 0.01,
             "Método: Pineau et al. (2025) NAR — IC 95% via bootstrap "
             f"(n={N_BOOTSTRAP} iterações). "
             f"Janela: ±{METAGENE_WINDOW} bp ao redor do TSS.",
             fontsize=7, color="gray", ha="left")

    fig.tight_layout(rect=[0, 0.04, 0.82, 1])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plotar_metagene_hetu_letu(perfis_hetu, perfis_letu, feature_name, chrom,
                               cond_label, output_path):
    """
    Plota metagene HETU vs LETU com IC 95% e painel inferior mostrando
    a diferença HETU − LETU (enriquecimento específico).
    """
    x = np.arange(-METAGENE_WINDOW + BIN_SIZE // 2,
                   METAGENE_WINDOW,
                   BIN_SIZE)

    fig = plt.figure(figsize=(10, 7))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    for matrix, color, label_prefix in [
        (perfis_hetu, HETU_COLOR, f"HETUs — top 40% mais expressos"),
        (perfis_letu, LETU_COLOR, f"LETUs — bottom 60%"),
    ]:
        if matrix.shape[0] == 0:
            continue
        mean, lower, upper = _bootstrap_ic95(matrix)
        n = matrix.shape[0]
        ax1.plot(x, mean, color=color, lw=2.2,
                 label=f"{label_prefix}  (n={n})")
        ax1.fill_between(x, lower, upper, color=color, alpha=0.18)

    ax1.axvline(0, color="black", linestyle="--", lw=1.2, zorder=5)
    ax1.axvspan(-PROMOTER_WINDOW, 0, alpha=0.05, color="gold")
    ax1.set_ylabel("Frequência de estrutura por bin", fontsize=10)
    if ax1.get_legend_handles_labels()[0]:
        ax1.legend(fontsize=9, loc="upper left")
    ax1.set_title(
        f"{feature_name}  —  HETUs vs LETUs\n"
        f"{REPLICON_LABELS.get(chrom, chrom)}  |  {cond_label}",
        fontsize=11, pad=10
    )

    # Painel inferior: diferença
    if perfis_hetu.shape[0] > 0 and perfis_letu.shape[0] > 0:
        mean_h = np.mean(perfis_hetu, axis=0)
        mean_l = np.mean(perfis_letu, axis=0)
        diff   = mean_h - mean_l
        ax2.bar(x, diff, width=BIN_SIZE * 0.9,
                color=[HETU_COLOR if d > 0 else LETU_COLOR for d in diff],
                alpha=0.7)
        ax2.axhline(0, color="black", lw=0.8)
        ax2.axvline(0, color="black", linestyle="--", lw=1.0)
        ax2.set_ylabel("HETU − LETU", fontsize=9)

    ax2.set_xlabel(f"Distância do TSS (bp)  |  bins de {BIN_SIZE} bp",
                   fontsize=10)
    plt.setp(ax1.get_xticklabels(), visible=False)

    fig.text(0.01, 0.01,
             "Critério HETUs/LETUs: Pineau et al. (2025) NAR — top 40% / bottom 60%"
             " por RPKM na condição de referência. IC 95% via bootstrap.",
             fontsize=7, color="gray")

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _worker_metagene_feature(args):
    """
    Worker paralelo: calcula e salva plots de metagene para UMA feature × replicon.
    Retorna string de status.
    """
    (bed_file, chrom, chr_length, tss_list, gene_tss_map,
     decis_dict, hetu_dict, cond_nome) = args

    feature = bed_file.replace(".bed", "")
    mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file), chrom, chr_length)
    if np.sum(mask) == 0:
        return f"    {feature} × {chrom} — ausente"

    # Metagene por decis
    perfis_decil = {}
    for decil in range(1, 11):
        ids_decil = {gid for gid, d in decis_dict.items() if d == decil}
        perfis_decil[decil] = calcular_perfil_metagene(
            tss_list, mask, chr_length, ids_decil, gene_tss_map)

    n_total = sum(m.shape[0] for m in perfis_decil.values())
    if n_total > 0:
        out_decil = os.path.join(OUT_DIR, "metagene",
                                  f"{feature}_{chrom}_{cond_nome}_decis.png")
        plotar_metagene_decis(perfis_decil, feature, chrom,
                               cond_nome, out_decil, n_total)

    # Metagene HETU vs LETU
    ids_hetu = {gid for gid, cl in hetu_dict.items() if cl == "HETU"}
    ids_letu = {gid for gid, cl in hetu_dict.items() if cl == "LETU"}
    m_hetu = calcular_perfil_metagene(tss_list, mask, chr_length,
                                       ids_hetu, gene_tss_map)
    m_letu = calcular_perfil_metagene(tss_list, mask, chr_length,
                                       ids_letu, gene_tss_map)
    out_hl = os.path.join(OUT_DIR, "metagene",
                           f"{feature}_{chrom}_{cond_nome}_HETU_LETU.png")
    plotar_metagene_hetu_letu(m_hetu, m_letu, feature, chrom, cond_nome, out_hl)

    return f"    {feature} × {chrom} — OK (n={n_total} TSSs)"


def rodar_metagenes(genes, tss_dict, genome, df_rpkm, non_b_files, cond_nome):
    """
    Roda todas as análises de metagene para uma condição de referência.
    Paraleliza por feature × replicon com MAX_WORKERS threads.
    """
    print(f"\n  Metagene — condição de referência: {cond_nome}")
    rpkm_col = f"RPKM_{cond_nome}"
    if rpkm_col not in df_rpkm.columns:
        print(f"  [AVISO] Coluna {rpkm_col} não encontrada no DataFrame.")
        return

    rpkm_dict  = df_rpkm.set_index("gene_id")[rpkm_col].to_dict()
    decis_dict = classificar_por_decis(rpkm_dict)
    hetu_dict  = classificar_hetu_letu(rpkm_dict)

    # Monta tarefas: uma por (feature × replicon)
    tarefas = []
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        tss_list = tss_dict[chrom]
        if not tss_list:
            continue

        gene_tss_map = {}
        for gene_id, gene in genes.items():
            if gene.chromosome != chrom:
                continue
            strand = _strand_str(gene.strand)
            pos = gene.start if strand == "+" else gene.end
            gene_tss_map[(pos, strand)] = gene_id

        for bed_file in non_b_files:
            tarefas.append((bed_file, chrom, chr_obj.length,
                            tss_list, gene_tss_map,
                            decis_dict, hetu_dict, cond_nome))

    print(f"  Tarefas (feature × replicon) : {len(tarefas)}  "
          f"| workers : {MAX_WORKERS}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(_worker_metagene_feature, t) for t in tarefas]
        for i, f in enumerate(concurrent.futures.as_completed(futures), 1):
            msg = f.result()
            print(f"  [{i:3d}/{len(tarefas)}] {msg}", end="\r")
    print()


# =============================================================================
# Análise 3 — Testes estatísticos formais
# =============================================================================

def _proporção_test(n_sucesso_A, n_total_A, n_sucesso_B, n_total_B):
    """
    Teste de proporção z-score (equivalente ao chi-quadrado).
    Retorna (z_stat, p_value).
    Retorna (0.0, 1.0) se qualquer grupo estiver vazio.
    """
    # Guard: grupos vazios → teste não aplicável
    if n_total_A == 0 or n_total_B == 0:
        return 0.0, 1.0
    p_A = n_sucesso_A / n_total_A
    p_B = n_sucesso_B / n_total_B
    p_pool = (n_sucesso_A + n_sucesso_B) / (n_total_A + n_total_B)
    if p_pool in (0.0, 1.0):
        return 0.0, 1.0
    se = np.sqrt(p_pool * (1 - p_pool) * (1/n_total_A + 1/n_total_B))
    if se == 0:
        return 0.0, 1.0
    z = (p_A - p_B) / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return round(z, 4), round(p, 6)


def _hipergeometrico_test(K, M, n, N):
    """
    Teste hipergeométrico.
    K = genes com estrutura no promotor
    M = total de genes com estrutura no genoma
    n = total de genes UP
    N = total de genes no genoma
    Retorna p-value.
    """
    rv = hypergeom(N, M, n)
    return float(rv.sf(K - 1))


def analise_estatistica_formal(genes, genome, tss_dict, df_rpkm,
                                non_b_files, chip_dict, cond_nome):
    """
    Para cada feature não-B e cada replicon:
    1. Teste de proporção: % promotores com estrutura em HETUs vs LETUs
    2. Teste hipergeométrico: enriquecimento em promotores de genes UP (TP1→TP4)
    3. Wilcoxon: RPKM médio em genes com estrutura no promotor vs sem
    Exporta CSV com todos os resultados e FDR aplicado.
    """
    print(f"\n[3] Testes estatísticos — condição: {cond_nome}")
    rpkm_col  = f"RPKM_{cond_nome}"
    if rpkm_col not in df_rpkm.columns:
        return

    rpkm_dict  = df_rpkm.set_index("gene_id")[rpkm_col].to_dict()
    hetu_dict  = classificar_hetu_letu(rpkm_dict)

    rows_prop  = []
    rows_hyper = []
    rows_wilcox= []
    p_vals_prop= []

    for bed_file in non_b_files:
        feature = bed_file.replace(".bed", "")
        for chrom in CHROMS:
            chr_obj = genome.genome_dict.get(chrom)
            if chr_obj is None:
                continue
            mask     = _bed_to_mask(os.path.join(SITES_DIR, bed_file),
                                    chrom, chr_obj.length)
            tss_list = tss_dict[chrom]
            if not tss_list or np.sum(mask) == 0:
                continue

            # Genes deste replicon
            genes_chrom = {gid: g for gid, g in genes.items()
                           if g.chromosome == chrom}

            # Para cada gene: tem estrutura no promotor?
            genes_com_estrutura = set()
            for gid, gene in genes_chrom.items():
                strand = _strand_str(gene.strand)
                tss_pos = gene.start if strand == "+" else gene.end
                s = max(0, tss_pos - PROMOTER_WINDOW) if strand == "+" \
                    else tss_pos
                e = tss_pos if strand == "+" \
                    else min(chr_obj.length, tss_pos + PROMOTER_WINDOW)
                if np.any(mask[s:e]):
                    genes_com_estrutura.add(gid)

            # Classificações
            hetu_ids = {gid for gid, cl in hetu_dict.items()
                        if cl == "HETU" and gid in genes_chrom}
            letu_ids = {gid for gid, cl in hetu_dict.items()
                        if cl == "LETU" and gid in genes_chrom}

            n_hetu_com = len(genes_com_estrutura & hetu_ids)
            n_letu_com = len(genes_com_estrutura & letu_ids)
            n_hetu_tot = len(hetu_ids)
            n_letu_tot = len(letu_ids)

            # Teste de proporção HETUs vs LETUs
            z, p = _proporção_test(n_hetu_com, n_hetu_tot,
                                   n_letu_com, n_letu_tot)
            p_vals_prop.append(p)
            rows_prop.append([
                feature, chrom, cond_nome,
                n_hetu_com, n_hetu_tot,
                round(n_hetu_com/n_hetu_tot*100, 2) if n_hetu_tot > 0 else 0,
                n_letu_com, n_letu_tot,
                round(n_letu_com/n_letu_tot*100, 2) if n_letu_tot > 0 else 0,
                z, p, "pendente_FDR"
            ])

            # Teste de Wilcoxon: RPKM em genes com estrutura vs sem
            rpkm_com = [rpkm_dict[g] for g in genes_com_estrutura
                        if g in rpkm_dict]
            rpkm_sem = [rpkm_dict[g] for g in genes_chrom
                        if g not in genes_com_estrutura and g in rpkm_dict]
            if len(rpkm_com) >= 5 and len(rpkm_sem) >= 5:
                stat_w, p_w = mannwhitneyu(rpkm_com, rpkm_sem,
                                           alternative="two-sided")
                rows_wilcox.append([
                    feature, chrom, cond_nome,
                    len(rpkm_com), round(np.median(rpkm_com), 3),
                    len(rpkm_sem), round(np.median(rpkm_sem), 3),
                    round(stat_w, 1), round(p_w, 6)
                ])

    # Aplica FDR (Benjamini-Hochberg) nos p-values de proporção
    if p_vals_prop:
        _, p_adj, _, _ = multipletests(p_vals_prop, method="fdr_bh")
        for i, row in enumerate(rows_prop):
            row[-1] = round(float(p_adj[i]), 6)
            row.append("*" if p_adj[i] < 0.05 else
                        "**" if p_adj[i] < 0.01 else
                        "ns")

    _write_csv("tabelas",
               f"teste_proporcao_HETU_LETU_{cond_nome}.csv",
               ["feature", "chrom", "condicao",
                "n_HETU_com_estrutura", "n_HETU_total", "pct_HETU",
                "n_LETU_com_estrutura", "n_LETU_total", "pct_LETU",
                "z_stat", "p_value", "p_adj_FDR", "sig"],
               rows_prop)

    _write_csv("tabelas",
               f"teste_wilcoxon_RPKM_{cond_nome}.csv",
               ["feature", "chrom", "condicao",
                "n_com_estrutura", "mediana_RPKM_com",
                "n_sem_estrutura", "mediana_RPKM_sem",
                "U_stat", "p_value"],
               rows_wilcox)


# =============================================================================
# Análise 4 — Orientação gênica
# Classes: divergente, convergente, tandem_fwd, tandem_rev
# =============================================================================

def classificar_orientacao_intergenica(genes, chrom, chr_length):
    """
    Para cada região intergênica, classifica a orientação com base nos
    genes flanqueadores (divergente / convergente / tandem).
    Retorna lista de dicts com {start, end, orientacao}.
    """
    genes_chrom = sorted(
        [g for g in genes.values() if g.chromosome == chrom],
        key=lambda g: g.start
    )
    regioes = []
    for i in range(len(genes_chrom) - 1):
        g1 = genes_chrom[i]
        g2 = genes_chrom[i + 1]
        s1 = _strand_str(g1.strand)
        s2 = _strand_str(g2.strand)
        inicio = g1.end
        fim    = g2.start
        if fim <= inicio:
            continue
        if   s1 == "+" and s2 == "+": orient = "tandem_fwd"
        elif s1 == "-" and s2 == "-": orient = "tandem_rev"
        elif s1 == "+" and s2 == "-": orient = "convergente"
        elif s1 == "-" and s2 == "+": orient = "divergente"
        else: orient = "outro"
        regioes.append({"start": inicio, "end": fim, "orient": orient})
    return regioes


def plotar_orientacao(resultados_orient, feature_name, chrom, output_path):
    """
    Barplot mostrando densidade média de estrutura não-B por classe de orientação.
    Com barras de erro (IC 95%).
    """
    classes = ["divergente", "convergente", "tandem_fwd", "tandem_rev"]
    cores   = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
    labels  = ["Divergente\n(→ ←)", "Convergente\n(← →)",
               "Tandem →", "Tandem ←"]

    medias = []
    erros  = []
    ns     = []
    for cl in classes:
        vals = resultados_orient.get(cl, [])
        if vals:
            medias.append(np.mean(vals))
            erros.append(1.96 * np.std(vals) / np.sqrt(len(vals)))
            ns.append(len(vals))
        else:
            medias.append(0); erros.append(0); ns.append(0)

    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(labels, medias, yerr=erros, capsize=5,
                  color=cores, alpha=0.82, width=0.6,
                  error_kw={"elinewidth": 1.5, "ecolor": "black"})

    for bar, n in zip(bars, ns):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + max(erros) * 0.05,
                f"n={n}", ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("Densidade média de estrutura\n(bp cobertos / bp intergênico)",
                  fontsize=10)
    ax.set_title(
        f"{feature_name}  —  Densidade por orientação gênica\n"
        f"{REPLICON_LABELS.get(chrom, chrom)}",
        fontsize=11, pad=10
    )
    ax.set_ylim(bottom=0)
    fig.text(0.01, 0.01,
             "Barras de erro: IC 95% (1,96 × EP). "
             "Orientação definida pelos genes flanqueadores de cada região intergênica.",
             fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def analise_orientacao(genes, genome, non_b_files):
    """Roda análise de orientação gênica para todas as features e replicons."""
    print("\n[4] Análise de orientação gênica...")
    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue
        regioes = classificar_orientacao_intergenica(genes, chrom, chr_obj.length)
        if not regioes:
            continue

        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            mask    = _bed_to_mask(os.path.join(SITES_DIR, bed_file),
                                   chrom, chr_obj.length)
            if np.sum(mask) == 0:
                continue

            result_orient = {cl: [] for cl in
                             ["divergente", "convergente", "tandem_fwd", "tandem_rev"]}
            for reg in regioes:
                tamanho = reg["end"] - reg["start"]
                if tamanho < 10:
                    continue
                dens = np.sum(mask[reg["start"]:reg["end"]]) / tamanho
                if reg["orient"] in result_orient:
                    result_orient[reg["orient"]].append(dens)

            out = os.path.join(OUT_DIR, "orientacao",
                               f"{feature}_{chrom}_orientacao.png")
            plotar_orientacao(result_orient, feature, chrom, out)

        print(f"    {chrom} — OK")


# =============================================================================
# Análise 5 — Cinética temporal TP1→TP4
# — Metagene estratificado por direção de mudança (UP / DOWN / estável)
# =============================================================================

def analise_temporal(genes, genome, tss_dict, non_b_files):
    """
    Para cada feature não-B, plota metagene separando genes que sobem (UP),
    descem (DOWN) ou se mantêm estáveis entre TP1 e TP4.
    """
    print("\n[5] Cinética temporal TP1 → TP4...")

    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue

        # Calcula FC por gene TP4/TP1
        rpkm_tp1 = calcular_rpkm_por_gene(genes, chrom, chr_obj.length,
                                           "total-RNA-TP1")
        rpkm_tp4 = calcular_rpkm_por_gene(genes, chrom, chr_obj.length,
                                           "total-RNA-TP4")

        direcao = {}
        for gid in rpkm_tp1:
            if gid not in rpkm_tp4:
                continue
            fc = (rpkm_tp4[gid] + 0.1) / (rpkm_tp1[gid] + 0.1)
            if   fc >= FC_UP:   direcao[gid] = "UP"
            elif fc <= FC_DOWN: direcao[gid] = "DOWN"
            else:               direcao[gid] = "Estável"

        # Mapa TSS → gene_id
        gene_tss_map = {}
        for gid, gene in genes.items():
            if gene.chromosome != chrom:
                continue
            strand = _strand_str(gene.strand)
            pos = gene.start if strand == "+" else gene.end
            gene_tss_map[(pos, strand)] = gid

        tss_list = tss_dict[chrom]

        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            mask    = _bed_to_mask(os.path.join(SITES_DIR, bed_file),
                                   chrom, chr_obj.length)
            if np.sum(mask) == 0:
                continue

            perfis = {}
            for classe, cor in [("UP", UP_COLOR), ("DOWN", DOWN_COLOR),
                                 ("Estável", STABLE_COLOR)]:
                ids_cl = {gid for gid, d in direcao.items() if d == classe}
                perfis[classe] = calcular_perfil_metagene(
                    tss_list, mask, chr_obj.length, ids_cl, gene_tss_map
                )

            x = np.arange(-METAGENE_WINDOW + BIN_SIZE//2,
                           METAGENE_WINDOW, BIN_SIZE)
            fig, ax = plt.subplots(figsize=(10, 5))

            for classe, cor in [("UP", UP_COLOR), ("DOWN", DOWN_COLOR),
                                 ("Estável", STABLE_COLOR)]:
                m = perfis[classe]
                if m.shape[0] == 0:
                    continue
                mean, lower, upper = _bootstrap_ic95(m)
                ax.plot(x, mean, color=cor, lw=2,
                        label=f"{classe}  (n={m.shape[0]})")
                ax.fill_between(x, lower, upper, color=cor, alpha=0.15)

            ax.axvline(0, color="black", linestyle="--", lw=1.2)
            ax.axvspan(-PROMOTER_WINDOW, 0, alpha=0.05, color="gold")
            ax.set_xlabel(f"Distância do TSS (bp)  |  bins de {BIN_SIZE} bp")
            ax.set_ylabel("Frequência de estrutura por bin")
            ax.set_title(
                f"{feature}  —  Cinética temporal (TP4 vs TP1)\n"
                f"{REPLICON_LABELS.get(chrom, chrom)}  |  "
                f"FC UP ≥ {FC_UP}  |  FC DOWN ≤ {FC_DOWN}",
                fontsize=11
            )
            if ax.get_legend_handles_labels()[0]:
                ax.legend(title="Direção de expressão TP1→TP4",
                          fontsize=9, title_fontsize=9)
            fig.text(0.01, 0.01,
                     "IC 95% via bootstrap. Estrutura avaliada na janela "
                     f"±{METAGENE_WINDOW} bp ao redor do TSS.",
                     fontsize=7, color="gray")
            fig.tight_layout(rect=[0, 0.05, 1, 1])
            out = os.path.join(OUT_DIR, "temporal",
                               f"{feature}_{chrom}_temporal.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)

        print(f"    {chrom} — OK")


# =============================================================================
# Análise 6 — Gradiente de sal (2.6 → 4.3 → 5.1 M NaCl)
# — Metagene por concentração de sal, com IC 95%
# =============================================================================

def analise_salinidade(genes, genome, tss_dict, df_rpkm, non_b_files):
    """
    Para cada feature não-B, plota metagene por concentração de sal,
    estratificado em HETUs/LETUs definidos pela condição ótima (4.3 M).
    """
    print("\n[6] Gradiente de salinidade (2.6 → 4.3 → 5.1 M NaCl)...")

    # HETUs e LETUs definidos pela condição ótima de crescimento (4.3 M)
    rpkm_ref_col = "RPKM_sal_4.3M"
    if rpkm_ref_col not in df_rpkm.columns:
        print(f"  [AVISO] Coluna {rpkm_ref_col} não encontrada.")
        return

    rpkm_ref   = df_rpkm.set_index("gene_id")[rpkm_ref_col].to_dict()
    hetu_dict  = classificar_hetu_letu(rpkm_ref)
    ids_hetu   = {gid for gid, cl in hetu_dict.items() if cl == "HETU"}
    ids_letu   = {gid for gid, cl in hetu_dict.items() if cl == "LETU"}

    sal_cores  = {"2.6M": "#1f77b4", "4.3M": "#2ca02c", "5.1M": "#d62728"}
    sal_labels = SAL_LABELS

    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue

        gene_tss_map = {}
        for gid, gene in genes.items():
            if gene.chromosome != chrom:
                continue
            strand = _strand_str(gene.strand)
            pos = gene.start if strand == "+" else gene.end
            gene_tss_map[(pos, strand)] = gid

        tss_list = tss_dict[chrom]

        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            mask    = _bed_to_mask(os.path.join(SITES_DIR, bed_file),
                                   chrom, chr_obj.length)
            if np.sum(mask) == 0:
                continue

            # Coleta perfis por concentração × HETU/LETU
            # O perfil não muda com a condição (estrutura é genômica),
            # mas a *classificação* dos genes muda → estratificação muda
            fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
            x = np.arange(-METAGENE_WINDOW + BIN_SIZE//2,
                           METAGENE_WINDOW, BIN_SIZE)

            for ax, (classe, ids_cl, titulo) in zip(
                axes,
                [("HETU", ids_hetu, "HETUs — top 40% (ref: 4,3 M NaCl)"),
                 ("LETU", ids_letu, "LETUs — bottom 60% (ref: 4,3 M NaCl)")]
            ):
                for sal, cor in sal_cores.items():
                    # Reclassifica HETUs/LETUs para ESTA concentração de sal
                    rpkm_sal_col = f"RPKM_sal_{sal}"
                    if rpkm_sal_col in df_rpkm.columns:
                        rpkm_sal = df_rpkm.set_index("gene_id")[rpkm_sal_col].to_dict()
                        hetu_sal = classificar_hetu_letu(rpkm_sal)
                        if classe == "HETU":
                            ids_uso = {g for g, c in hetu_sal.items()
                                       if c == "HETU"}
                        else:
                            ids_uso = {g for g, c in hetu_sal.items()
                                       if c == "LETU"}
                    else:
                        ids_uso = ids_cl

                    m = calcular_perfil_metagene(
                        tss_list, mask, chr_obj.length, ids_uso, gene_tss_map
                    )
                    if m.shape[0] == 0:
                        continue
                    mean, lower, upper = _bootstrap_ic95(m)
                    ax.plot(x, mean, color=cor, lw=2,
                            label=f"{sal_labels[sal]}  (n={m.shape[0]})")
                    ax.fill_between(x, lower, upper, color=cor, alpha=0.15)

                ax.axvline(0, color="black", linestyle="--", lw=1.1)
                ax.axvspan(-PROMOTER_WINDOW, 0, alpha=0.05, color="gold")
                ax.set_xlabel(f"Distância do TSS (bp)")
                ax.set_title(titulo, fontsize=10)
                # Só desenha legenda se houver curvas plotadas
                if ax.get_legend_handles_labels()[0]:
                    ax.legend(fontsize=8)
                else:
                    ax.text(0.5, 0.5, "Sem genes nesta classe",
                            ha="center", va="center", transform=ax.transAxes,
                            fontsize=9, color="gray")

            axes[0].set_ylabel("Frequência de estrutura por bin")
            fig.suptitle(
                f"{feature}  —  Gradiente de salinidade\n"
                f"{REPLICON_LABELS.get(chrom, chrom)}",
                fontsize=12, y=1.01
            )
            fig.text(0.01, -0.02,
                     "HETUs/LETUs reclassificados por concentração de sal. "
                     "IC 95% via bootstrap. Estrutura não-B é propriedade genômica fixa.",
                     fontsize=7, color="gray")
            fig.tight_layout()
            out = os.path.join(OUT_DIR, "salinidade",
                               f"{feature}_{chrom}_salinidade.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)

        print(f"    {chrom} — OK")


# =============================================================================
# Análise 7 — RNA ribossomal vs total
# — Genes de RNA estável têm comportamento distinto (paper GRATIOSA, Fig. 2)
# =============================================================================

def analise_rrna(genes_asrna, genome, tss_dict, non_b_files):
    """
    Perfil metagene para genes de RNA ribossomal/estável separadamente,
    comparado com todos os genes proteicos.
    """
    print("\n[7] RNA ribossomal vs total...")

    for chrom in CHROMS:
        chr_obj = genome.genome_dict.get(chrom)
        if chr_obj is None:
            continue

        ids_rrna = {gid for gid, g in genes_asrna.items()
                    if g.chromosome == chrom}
        if not ids_rrna:
            continue

        gene_tss_map = {}
        for gid, gene in genes_asrna.items():
            if gene.chromosome != chrom:
                continue
            strand = _strand_str(gene.strand)
            pos = gene.start if strand == "+" else gene.end
            gene_tss_map[(pos, strand)] = gid

        tss_list = tss_dict[chrom]
        x = np.arange(-METAGENE_WINDOW + BIN_SIZE//2, METAGENE_WINDOW, BIN_SIZE)

        for bed_file in non_b_files:
            feature = bed_file.replace(".bed", "")
            mask    = _bed_to_mask(os.path.join(SITES_DIR, bed_file),
                                   chrom, chr_obj.length)
            if np.sum(mask) == 0:
                continue

            m_rrna = calcular_perfil_metagene(tss_list, mask, chr_obj.length,
                                              ids_rrna, gene_tss_map)
            m_all  = calcular_perfil_metagene(tss_list, mask, chr_obj.length)

            if m_rrna.shape[0] == 0 and m_all.shape[0] == 0:
                continue

            fig, ax = plt.subplots(figsize=(10, 5))
            for m, cor, lbl in [
                (m_all,  "#888888", f"Todos os genes  (n={m_all.shape[0]})"),
                (m_rrna, RRNA_COLOR, f"RNAs estáveis / asRNAs  (n={m_rrna.shape[0]})"),
            ]:
                if m.shape[0] == 0:
                    continue
                mean, lower, upper = _bootstrap_ic95(m)
                ax.plot(x, mean, color=cor, lw=2, label=lbl)
                ax.fill_between(x, lower, upper, color=cor, alpha=0.18)

            ax.axvline(0, color="black", linestyle="--", lw=1.2)
            ax.axvspan(-PROMOTER_WINDOW, 0, alpha=0.05, color="gold")
            ax.set_xlabel(f"Distância do TSS (bp)  |  bins de {BIN_SIZE} bp")
            ax.set_ylabel("Frequência de estrutura por bin")
            ax.set_title(
                f"{feature}  —  RNA estável vs todos os genes\n"
                f"{REPLICON_LABELS.get(chrom, chrom)}",
                fontsize=11
            )
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=9)
            fig.text(0.01, 0.01,
                     "RNAs estáveis incluem asRNAs e ncRNAs anotados em asRNAs.gff3. "
                     "IC 95% via bootstrap.",
                     fontsize=7, color="gray")
            fig.tight_layout(rect=[0, 0.05, 1, 1])
            out = os.path.join(OUT_DIR, "metagene",
                               f"{feature}_{chrom}_rrna_vs_total.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)

        print(f"    {chrom} — OK")


# =============================================================================
# Painel resumo: todos os replicons × todas as features × uma condição
# =============================================================================

def plotar_painel_resumo(genes, genome, tss_dict, df_rpkm,
                          non_b_files, cond_nome):
    """
    Grid compacto: linha = feature, coluna = replicon.
    Mostra HETU (vermelho) vs LETU (azul) em cada célula.
    Útil para visão panorâmica de um experimento.
    """
    print(f"\n[Painel resumo] condição: {cond_nome}")
    rpkm_col = f"RPKM_{cond_nome}"
    if rpkm_col not in df_rpkm.columns:
        return

    rpkm_dict = df_rpkm.set_index("gene_id")[rpkm_col].to_dict()
    hetu_dict = classificar_hetu_letu(rpkm_dict)
    ids_hetu  = {gid for gid, cl in hetu_dict.items() if cl == "HETU"}
    ids_letu  = {gid for gid, cl in hetu_dict.items() if cl == "LETU"}

    n_feat = len(non_b_files)
    n_rep  = len(CHROMS)
    x = np.arange(-METAGENE_WINDOW + BIN_SIZE//2, METAGENE_WINDOW, BIN_SIZE)

    fig, axes = plt.subplots(n_feat, n_rep,
                              figsize=(5 * n_rep, 3.2 * n_feat),
                              sharex=True)
    if n_feat == 1:
        axes = [axes]

    for fi, bed_file in enumerate(non_b_files):
        feature = bed_file.replace(".bed", "")
        for ci, chrom in enumerate(CHROMS):
            ax = axes[fi][ci]
            chr_obj = genome.genome_dict.get(chrom)
            if chr_obj is None:
                ax.axis("off"); continue
            mask = _bed_to_mask(os.path.join(SITES_DIR, bed_file),
                                chrom, chr_obj.length)
            if np.sum(mask) == 0:
                ax.text(0.5, 0.5, "Ausente", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="gray")
                ax.axis("off"); continue

            tss_list = tss_dict[chrom]
            gene_tss_map = {}
            for gid, gene in genes.items():
                if gene.chromosome != chrom:
                    continue
                strand = _strand_str(gene.strand)
                pos = gene.start if strand == "+" else gene.end
                gene_tss_map[(pos, strand)] = gid

            for ids_cl, cor, lbl in [
                (ids_hetu, HETU_COLOR, "HETUs"),
                (ids_letu, LETU_COLOR, "LETUs"),
            ]:
                m = calcular_perfil_metagene(tss_list, mask, chr_obj.length,
                                             ids_cl, gene_tss_map)
                if m.shape[0] == 0:
                    continue
                mean, lower, upper = _bootstrap_ic95(m)
                ax.plot(x, mean, color=cor, lw=1.5, label=lbl)
                ax.fill_between(x, lower, upper, color=cor, alpha=0.15)

            ax.axvline(0, color="black", linestyle="--", lw=0.8)
            ax.axhline(0, color="gray", linestyle=":", lw=0.5)

            if fi == 0:
                chrom_label = REPLICON_LABELS.get(chrom, chrom).split("(")[0].strip()
                ax.set_title(chrom_label, fontsize=9, pad=4)
            if ci == 0:
                ax.set_ylabel(feature.replace("_", " "), fontsize=8)
            if fi == n_feat - 1:
                ax.set_xlabel("Distância TSS (bp)", fontsize=8)

    # Legenda global
    legend_elements = [
        Line2D([0], [0], color=HETU_COLOR, lw=2, label="HETUs (top 40%)"),
        Line2D([0], [0], color=LETU_COLOR, lw=2, label="LETUs (bottom 60%)"),
        Line2D([0], [0], color="black", lw=1, linestyle="--", label="TSS"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=9, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(
        f"Painel resumo — Todas as features × todos os replicons\n"
        f"Condição de referência: {cond_nome}  |  "
        f"Janela: ±{METAGENE_WINDOW} bp  |  Bins: {BIN_SIZE} bp",
        fontsize=12, y=1.01
    )
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "metagene",
                       f"PAINEL_RESUMO_{cond_nome}.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Painel salvo: {out}")


# =============================================================================
# Main
# =============================================================================

def main():
    _make_dirs()
    print("\n" + "=" * 70)
    print("PIPELINE v3 — Halobacterium salinarum NRC-1")
    print("Metodologia: Pineau et al. (2025) NAR — GRATIOSA")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Carregamento
    # ------------------------------------------------------------------
    genome, genes_main, genes_asrna = load_genome()
    tss_dict  = load_tss()
    hic_dict  = load_hic()
    chip_dict = load_chipseq()
    non_b_files = list_non_b_files()

    print(f"\nFeatures não-B encontradas : {len(non_b_files)}")
    for f in non_b_files:
        print(f"  {f}")

    # ------------------------------------------------------------------
    # Define todas as condições para cálculo de RPKM
    # ------------------------------------------------------------------
    todas_condicoes = {}

    # Cinética temporal (TP1 → TP4)
    for tp in TIMEPOINTS:
        todas_condicoes[tp] = (tp, "single")

    # Gradiente de sal (médias rep A + B)
    for sal, reps in SAL_CONDITIONS.items():
        todas_condicoes[f"sal_{sal}"] = ((reps["rep_A"], reps["rep_B"]), "media")

    # RNA ribossomal
    for tp in RIBO_TIMEPOINTS:
        todas_condicoes[f"ribo_{tp}"] = (tp, "single")

    # ------------------------------------------------------------------
    # Análise 1 — RPKM
    # ------------------------------------------------------------------
    df_rpkm = analise_rpkm_completa(genes_main, genome, todas_condicoes)

    # ------------------------------------------------------------------
    # Análise 2 — Metagene estratificado
    # Roda para cada condição de referência
    # ------------------------------------------------------------------
    print("\n[2] Metagenes estratificados por decis e HETU/LETU...")

    # Condições prioritárias para metagene
    conds_metagene = (
        list(TIMEPOINTS) +
        [f"sal_{sal}" for sal in SAL_CONDITIONS]
    )
    for cond in conds_metagene:
        rodar_metagenes(genes_main, tss_dict, genome,
                        df_rpkm, non_b_files, cond)

    # ------------------------------------------------------------------
    # Análise 3 — Testes estatísticos formais
    # ------------------------------------------------------------------
    for cond in conds_metagene:
        analise_estatistica_formal(genes_main, genome, tss_dict,
                                   df_rpkm, non_b_files, chip_dict, cond)

    # ------------------------------------------------------------------
    # Análise 4 — Orientação gênica
    # ------------------------------------------------------------------
    analise_orientacao(genes_main, genome, non_b_files)

    # ------------------------------------------------------------------
    # Análise 5 — Cinética temporal
    # ------------------------------------------------------------------
    analise_temporal(genes_main, genome, tss_dict, non_b_files)

    # ------------------------------------------------------------------
    # Análise 6 — Gradiente de salinidade
    # ------------------------------------------------------------------
    analise_salinidade(genes_main, genome, tss_dict, df_rpkm, non_b_files)

    # ------------------------------------------------------------------
    # Análise 7 — RNA ribossomal vs total
    # ------------------------------------------------------------------
    analise_rrna(genes_asrna, genome, tss_dict, non_b_files)

    # ------------------------------------------------------------------
    # Painéis resumo (um por condição principal)
    # ------------------------------------------------------------------
    print("\n[Painéis resumo]")
    for cond in ["total-RNA-TP1", "sal_4.3M"]:
        plotar_painel_resumo(genes_main, genome, tss_dict,
                             df_rpkm, non_b_files, cond)

    # ------------------------------------------------------------------
    # Resumo final
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANÁLISES CONCLUÍDAS")
    print(f"  Features não-B   : {len(non_b_files)}")
    print(f"  Replicons        : {len(CHROMS)}")
    print(f"  Condições RPKM   : {len(todas_condicoes)}")
    print(f"\n  Saídas em: {OUT_DIR}/")
    print(f"    tabelas/     → CSVs com RPKM, classificações, testes")
    print(f"    metagene/    → PNGs por feature × condição (decis + HETU/LETU)")
    print(f"    orientacao/  → PNGs por feature × orientação gênica")
    print(f"    temporal/    → PNGs cinética TP1→TP4")
    print(f"    salinidade/  → PNGs gradiente de sal")
    print(f"    resumo/      → TXTs por replicon")
    print()


if __name__ == "__main__":
    main()
