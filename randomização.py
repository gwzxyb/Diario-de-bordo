# RANDOMIZAÇÃO DE JANELAS PROMOTORAS E CDS — Halobacterium salinarum NRC-1
#   1) VERTICAL   — embaralha por coluna (posição): preserva viés AT/GC por posição
#   2) HORIZONTAL — embaralha por linha (gene):     preserva composição de cada gene
#
#   Suporta: Janela [-500, +500] do TSS (ou Start Codon como fallback)
#            Randomização do CDS inteiro (em pasta separada)

import os
import json
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

BASE_DIR       = "/home/mcavalcante/codigo/Diario-de-bordo/dados"  
GENOME_FILE   = os.path.join(BASE_DIR, "Hsalinarum.fa")
TSS_FILE       = os.path.join(BASE_DIR, "Hsalinarum_TSS.gtf")
UPSTREAM_BP    = 500    # pb upstream do TSS/Start
DOWNSTREAM_BP  = 500    # pb downstream do TSS/Start
WINDOW_SIZE    = UPSTREAM_BP + DOWNSTREAM_BP  # Total: 1000 pb
N_PERMUTATIONS = 100    
SEED           = 42     
OUT_DIR        = "output_randomization_py"
MAIN_CHROM     = "NC_002607.1" 

# ESTRUTURA DE SUBDIRETÓRIOS
DIR_REF        = os.path.join(OUT_DIR, "00_reference")
DIR_VERT       = os.path.join(OUT_DIR, "01_vertical_permutations")
DIR_HORIZ      = os.path.join(OUT_DIR, "02_horizontal_permutations")
DIR_CDS        = os.path.join(OUT_DIR, "04_cds_permutations")
DIR_REPORTS    = os.path.join(OUT_DIR, "03_reports")

vert_dir       = os.path.join(DIR_VERT, "individual_fastas")
horiz_dir      = os.path.join(DIR_HORIZ, "individual_fastas")

cds_vert_dir   = os.path.join(DIR_CDS, "01_vertical_permutations")
cds_horiz_dir  = os.path.join(DIR_CDS, "02_horizontal_permutations")

directories = [
    DIR_REF, DIR_VERT, vert_dir, DIR_HORIZ, horiz_dir, 
    DIR_CDS, cds_vert_dir, cds_horiz_dir, DIR_REPORTS
]
for d in directories:
    os.makedirs(d, exist_ok=True)

np.random.seed(SEED)
rng = np.random.default_rng(SEED)

# 1. CARREGAMENTO DO GENOMA
genome = {}
for rec in SeqIO.parse(GENOME_FILE, "fasta"):
    chrom_id = rec.id.split()[0]
    genome[chrom_id] = str(rec.seq).upper()
    print(f"  {chrom_id}: {len(rec.seq)} bp")

if MAIN_CHROM not in genome:
    raise ValueError(f"Cromossomo '{MAIN_CHROM}' não encontrado no FASTA.")

main_seq = genome[MAIN_CHROM]
main_len = len(main_seq)
print(f"  Cromossomo principal: {main_len} bp\n")

# 2. PARSING DO ARQUIVO GTF (TSS E CDS)
tss_dict = {}       # gene_id -> {pos, strand, type: 'TSS'}
cds_start_dict = {} # gene_id -> {pos, strand, type: 'START'}
cds_seq_dict = {}   # gene_id -> sequência CDS inteira

with open(TSS_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        seqname = fields[0].split()[0]
        if seqname != MAIN_CHROM:
            continue

        feature_type = fields[2]
        start  = int(fields[3])  # 1-based
        end    = int(fields[4])
        strand = fields[6]
        attrs  = fields[8] if len(fields) > 8 else ""

        # Extração do gene_id
        gene_id = None
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                break
            if attr.startswith("tss_location") and gene_id is None:
                val = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                for marker in ("_of_gene_", "_of_genes_"):
                    if marker in val:
                        gene_id = val.split(marker)[-1].rstrip(";").strip()
                        break

        if not gene_id:
            continue

        # Processamento de TSS primário ("P")
        if "tss_class" in attrs:
            tss_class = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("tss_class"):
                    tss_class = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                    break
            if tss_class == "P" and gene_id not in tss_dict:
                tss_pos = start if strand == "+" else end
                tss_dict[gene_id] = {"pos": tss_pos, "strand": strand, "type": "TSS"}

        # Processamento de CDS (Start Codon Fallback + CDS Inteiro)
        if feature_type == "CDS":
            start_pos = start if strand == "+" else end
            if gene_id not in cds_start_dict:
                cds_start_dict[gene_id] = {"pos": start_pos, "strand": strand, "type": "START"}
            
            # Extração da sequência de CDS
            s_idx = start - 1
            e_idx = end
            if 0 <= s_idx < main_len and e_idx <= main_len:
                seq_raw = main_seq[s_idx:e_idx]
                if strand == "-":
                    seq_raw = str(Seq(seq_raw).reverse_complement())
                cds_seq_dict[gene_id] = seq_raw

# Unificação das ancoragens (TSS com fallback para Start Codon)
all_genes = set(list(tss_dict.keys()) + list(cds_start_dict.keys()))
anchors = {}
for gid in all_genes:
    if gid in tss_dict:
        anchors[gid] = tss_dict[gid]
    else:
        anchors[gid] = cds_start_dict[gid]

print(f"  Total de genes mapeados: {len(anchors)}")
print(f"  - Ancorados por TSS primário: {len(tss_dict)}")
print(f"  - Ancorados por Start Codon (fallback): {len(anchors) - len(tss_dict)}\n")

# 3. EXTRAÇÃO DAS JANELAS PROMOTORAS [-500, +500]
promoter_sequences = []
promoter_gene_ids  = []

for gid, info in anchors.items():
    pos    = info["pos"]
    strand = info["strand"]

    if strand == "+":
        s = pos - 1 - UPSTREAM_BP     # 0-based, inclusivo
        e = pos - 1 + DOWNSTREAM_BP   # 0-based, exclusivo
    else:
        s = pos - DOWNSTREAM_BP
        e = pos + UPSTREAM_BP

    if s < 0 or e > main_len:
        continue

    window = main_seq[s:e]
    if strand == "-":
        window = str(Seq(window).reverse_complement())

    if len(window) == WINDOW_SIZE and "N" not in window:
        promoter_sequences.append(window)
        promoter_gene_ids.append(gid)

n_genes_prom = len(promoter_sequences)
matrix_prom = np.array([list(seq) for seq in promoter_sequences], dtype="<U1")
print(f"  Matriz de Promotores [-500, +500]: {matrix_prom.shape[0]} genes × {matrix_prom.shape[1]} posições")

# 4. PREPARAÇÃO DOS CDSs INTEIROS
cds_sequences = []
cds_gene_ids  = []
for gid in promoter_gene_ids:
    if gid in cds_seq_dict and "N" not in cds_seq_dict[gid]:
        cds_sequences.append(cds_seq_dict[gid])
        cds_gene_ids.append(gid)

# Padroniza comprimentos dos CDSs via padding à direita com '-' para formar matriz regular
max_cds_len = max(len(s) for s in cds_sequences)
cds_padded = [list(s.ljust(max_cds_len, "-")) for s in cds_sequences]
matrix_cds = np.array(cds_padded, dtype="<U1")
print(f"  Matriz de CDS Inteiros: {matrix_cds.shape[0]} genes × {matrix_cds.shape[1]} posições (máx)")

# Salva referências reais
np.save(os.path.join(DIR_REF, "matrix_promoter_real.npy"), matrix_prom)
np.save(os.path.join(DIR_REF, "matrix_cds_real.npy"), matrix_cds)
with open(os.path.join(DIR_REF, "gene_ids.txt"), "w") as f:
    f.write("\n".join(promoter_gene_ids))

# 5. FUNÇÕES DE RANDOMIZAÇÃO E FASTA
def randomize_vertical(mat, rng):
    out = mat.copy()
    for j in range(out.shape[1]):
        col = out[:, j].copy()
        rng.shuffle(col)
        out[:, j] = col
    return out

def randomize_horizontal(mat, rng):
    out = mat.copy()
    for i in range(out.shape[0]):
        row = out[i, :].copy()
        rng.shuffle(row)
        out[i, :] = row
    return out

def matrix_to_fasta(mat, ids, path, strip_gaps=False):
    with open(path, "w") as f:
        for gid, row in zip(ids, mat):
            seq_str = "".join(row)
            if strip_gaps:
                seq_str = seq_str.replace("-", "")
            f.write(f">{gid}\n{seq_str}\n")

# 6. PERMUTAÇÕES
print("\n GERANDO PERMUTAÇÕES (PROMOTORES E CDSs) ")

all_prom_vert, all_prom_horiz = [], []
all_cds_vert, all_cds_horiz   = [], []

for p in range(1, N_PERMUTATIONS + 1):
    # --- Promotores [-500, +500] ---
    mat_p_v = randomize_vertical(matrix_prom, rng)
    mat_p_h = randomize_horizontal(matrix_prom, rng)
    all_prom_vert.append(mat_p_v)
    all_prom_horiz.append(mat_p_h)
    
    matrix_to_fasta(mat_p_v, promoter_gene_ids, os.path.join(vert_dir, f"promoter_vertical_perm_{p:03d}.fasta"))
    matrix_to_fasta(mat_p_h, promoter_gene_ids, os.path.join(horiz_dir, f"promoter_horizontal_perm_{p:03d}.fasta"))

    # --- CDS Inteiro ---
    mat_c_v = randomize_vertical(matrix_cds, rng)
    mat_c_h = randomize_horizontal(matrix_cds, rng)
    all_cds_vert.append(mat_c_v)
    all_cds_horiz.append(mat_c_h)
    
    matrix_to_fasta(mat_c_v, cds_gene_ids, os.path.join(cds_vert_dir, f"cds_vertical_perm_{p:03d}.fasta"), strip_gaps=True)
    matrix_to_fasta(mat_c_h, cds_gene_ids, os.path.join(cds_horiz_dir, f"cds_horizontal_perm_{p:03d}.fasta"), strip_gaps=True)

    if p % 20 == 0 or p == N_PERMUTATIONS:
        print(f"  {p}/{N_PERMUTATIONS} permutações concluídas")

# Concatenação dos arquivos finais em FASTA
vert_prom_cat  = np.vstack(all_prom_vert)
horiz_prom_cat = np.vstack(all_prom_horiz)
vert_prom_ids  = [f"prom_vert_perm{p:03d}_{i+1:04d}" for p in range(1, N_PERMUTATIONS+1) for i in range(n_genes_prom)]
horiz_prom_ids = [f"prom_horiz_perm{p:03d}_{i+1:04d}" for p in range(1, N_PERMUTATIONS+1) for i in range(n_genes_prom)]

matrix_to_fasta(vert_prom_cat, vert_prom_ids, os.path.join(DIR_VERT, "promoter_vertical_all.fasta"))
matrix_to_fasta(horiz_prom_cat, horiz_prom_ids, os.path.join(DIR_HORIZ, "promoter_horizontal_all.fasta"))

# FASTA de Referência
matrix_to_fasta(matrix_prom, promoter_gene_ids, os.path.join(DIR_REF, "real_promoter_sequences.fasta"))
matrix_to_fasta(matrix_cds, cds_gene_ids, os.path.join(DIR_REF, "real_cds_sequences.fasta"), strip_gaps=True)

print("\nPROCESSAMENTO CONCLUÍDO.")
print(f"Diretório principal de saída: {OUT_DIR}/")
print("  ├── 00_reference/             (Matrizes e FASTAs originais)")
print("  ├── 01_vertical_permutations/  (Promotores [-500, +500] - embaralhado por coluna)")
print("  ├── 02_horizontal_permutations/(Promotores [-500, +500] - embaralhado por linha)")
print("  └── 04_cds_permutations/      (CDSs inteiros - pastas de permutações verticais e horizontais)")
