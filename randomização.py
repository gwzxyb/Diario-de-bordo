# RANDOMIZAÇÃO DE JANELAS PROMOTORAS — Halobacterium salinarum NRC-1
#   1) VERTICAL   — embaralha por coluna (posição): preserva viés AT/GC por posição
#   2) HORIZONTAL — embaralha por linha (gene):     preserva composição de cada gene
#
#  pip install biopython numpy

import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

BASE_DIR     = "/home/mcavalcante/codigo/Diario-de-bordo/dados"  
GENOME_FILE  = os.path.join(BASE_DIR, "Hsalinarum.fa")
TSS_FILE     = os.path.join(BASE_DIR, "Hsalinarum_TSS.gtf")
UPSTREAM_BP    = 250    # pb upstream do TSS (ten Caten 2015)
N_PERMUTATIONS = 100    # permutações por modalidade
SEED           = 42     # a resposta de todo o universo
OUT_DIR        = "output_randomization_py"
MAIN_CHROM = "NC_002607.1" 
np.random.seed(SEED)
os.makedirs(OUT_DIR, exist_ok=True)
genome = {}
for rec in SeqIO.parse(GENOME_FILE, "fasta"):
    chrom_id = rec.id.split()[0]          # pega só o primeiro token
    genome[chrom_id] = str(rec.seq).upper()
    print(f"  {chrom_id}: {len(rec.seq)} bp")
if MAIN_CHROM not in genome:
    raise ValueError(f"Cromossomo '{MAIN_CHROM}' não encontrado no FASTA. "
                     f"IDs disponíveis: {list(genome.keys())}")

main_seq = genome[MAIN_CHROM]
main_len = len(main_seq)
print(f"  Cromossomo principal: {main_len} bp\n")
tss_list = []
with open(TSS_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        seqname = fields[0].split()[0]
        if seqname != MAIN_CHROM:
            continue   # só cromossomo principal

        start  = int(fields[3])   # 1-based
        end    = int(fields[4])
        strand = fields[6]
        attrs  = fields[8] if len(fields) > 8 else ""

        # FILTRO: apenas TSS primário (tss_class "P")
        tss_class = None
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("tss_class"):
                tss_class = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                break
        if tss_class != "P":         
            continue
        # ─────────────────────────────────────────────────────────────────
        # TSS: strand "+" → posição = start; strand "-" → posição = end
        tss_pos = start if strand == "+" else end
        # extrai gene_id — tenta gene_id direto, senão usa tss_location
        gene_id = None
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                break
            if attr.startswith("tss_location") and gene_id is None:
                # extrai "VNG_RS00010" de "82nt_upstream_of_gene_VNG_RS00010;"
                val = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                for marker in ("_of_gene_", "_of_genes_"):
                    if marker in val:
                        gene_id = val.split(marker)[-1].rstrip(";").strip()
                        break
        if gene_id is None:
            gene_id = f"gene_{len(tss_list)}"

        tss_list.append({"pos": tss_pos, "strand": strand, "gene_id": gene_id})
# 1 TSS por gene (mantém o primeiro encontrado)
seen = set()
tss_unique = []
for t in tss_list:
    if t["gene_id"] not in seen:
        seen.add(t["gene_id"])
        tss_unique.append(t)
print(f"  TSS únicos encontrados: {len(tss_unique)}\n")
# EXTRAÇÃO DAS JANELAS PROMOTORAS 
print(" EXTRAÇÃO DE JANELAS ")
# Janela upstream de 250 pb orientada 5'→3' com TSS na borda direita.
#   strand "+" : main_seq[ tss_pos-1-250 : tss_pos-1 ]
#   strand "-" : rev_comp( main_seq[ tss_pos : tss_pos+250 ] )
#   (tss_pos é 1-based; slices Python são 0-based e exclusivos no fim)
sequences = []
gene_ids  = []
for t in tss_unique:
    pos    = t["pos"]
    strand = t["strand"]

    if strand == "+":
        s = pos - 1 - UPSTREAM_BP   # 0-based, inclusivo
        e = pos - 1                  # 0-based, exclusivo
        if s < 0 or e > main_len:
            continue
        window = main_seq[s:e]
    else:
        # upstream de um gene "-" está à DIREITA do TSS nas coord. genômicas
        s = pos       # pos é 1-based; pos em 0-based = pos-1, mas a base imediata
                      # downstream (= upstream do gene "-") começa em pos (0-based)
        e = pos + UPSTREAM_BP
        if s < 0 or e > main_len:
            continue
        window = str(Seq(main_seq[s:e]).reverse_complement())

    if len(window) == UPSTREAM_BP and "N" not in window:
        sequences.append(window)
        gene_ids.append(t["gene_id"])
n_genes = len(sequences)
print(f"  Janelas válidas ({UPSTREAM_BP} pb, sem N): {n_genes} genes\n")
# MATRIZ [n_genes × UPSTREAM_BP]
matrix = np.array([list(seq) for seq in sequences], dtype="<U1")
print(f"  Matriz: {matrix.shape[0]} genes × {matrix.shape[1]} posições")
# Salva matriz real e IDs
np.save(os.path.join(OUT_DIR, "matrix_real.npy"), matrix)
with open(os.path.join(OUT_DIR, "gene_ids.txt"), "w") as f:
    f.write("\n".join(gene_ids))
# FUNÇÕES DE RANDOMIZAÇÃO 
def randomize_vertical(mat, rng):
    """Embaralha cada coluna independentemente (entre genes).
    Preserva a composição de bases POR POSIÇÃO; destrói a identidade do gene."""
    out = mat.copy()
    for j in range(out.shape[1]):
        col = out[:, j].copy()
        rng.shuffle(col)
        out[:, j] = col
    return out

def randomize_horizontal(mat, rng):
    """Embaralha cada linha independentemente (dentro de cada gene).
    Preserva a composição global de cada gene; destrói o contexto posicional."""
    out = mat.copy()
    for i in range(out.shape[0]):
        row = out[i, :].copy()
        rng.shuffle(row)
        out[i, :] = row
    return out

def matrix_to_fasta(mat, ids, path):
    """Salva uma matriz de bases como FASTA (1 sequência por gene)."""
    with open(path, "w") as f:
        for gid, row in zip(ids, mat):
            f.write(f">{gid}\n{''.join(row)}\n")
#  PERMUTAÇÕES
print("\n GERANDO PERMUTAÇÕES ")
vert_dir   = os.path.join(OUT_DIR, "vertical")
horiz_dir  = os.path.join(OUT_DIR, "horizontal")
os.makedirs(vert_dir,   exist_ok=True)
os.makedirs(horiz_dir,  exist_ok=True)
rng = np.random.default_rng(SEED)

# Lista para concatenar todas as permutações em um único FASTA
all_vertical = []
all_horizontal = []
for p in range(1, N_PERMUTATIONS + 1):
    # --- vertical ---
    mat_v = randomize_vertical(matrix, rng)
    all_vertical.append(mat_v)
    # Salva individual (opcional)
    matrix_to_fasta(mat_v, gene_ids,
                    os.path.join(vert_dir, f"vertical_perm_{p:03d}.fasta"))

    # --- horizontal ---
    mat_h = randomize_horizontal(matrix, rng)
    all_horizontal.append(mat_h)
    matrix_to_fasta(mat_h, gene_ids,
                    os.path.join(horiz_dir, f"horizontal_perm_{p:03d}.fasta"))

    if p % 20 == 0 or p == N_PERMUTATIONS:
        print(f"  {p}/{N_PERMUTATIONS} permutações geradas")
# Concatena todas as permutações em um único arquivo FASTA (para usar com o gfa)
vert_cat = np.vstack(all_vertical)
horiz_cat = np.vstack(all_horizontal)
vert_ids = [f"vert_perm{p:03d}_{i+1:04d}" for p in range(1, N_PERMUTATIONS+1) for i in range(n_genes)]
horiz_ids = [f"horiz_perm{p:03d}_{i+1:04d}" for p in range(1, N_PERMUTATIONS+1) for i in range(n_genes)]
matrix_to_fasta(vert_cat, vert_ids, os.path.join(OUT_DIR, "vertical_all.fasta"))
matrix_to_fasta(horiz_cat, horiz_ids, os.path.join(OUT_DIR, "horizontal_all.fasta"))
# Salva também o FASTA real para referência
matrix_to_fasta(matrix, gene_ids,
                os.path.join(OUT_DIR, "real_sequences.fasta"))
# VERIFICAÇÃO DO MÉTODO 
print("\n VERIFICAÇÃO DO MÉTODO ")
def at_per_position(mat):
    """Fração de A/T em cada coluna (posição)."""
    is_at = np.isin(mat, ["A", "T"])
    return is_at.mean(axis=0)
at_real = at_per_position(matrix)
rng_check = np.random.default_rng(SEED + 1)
at_vert   = at_per_position(randomize_vertical(matrix,   rng_check))
at_horiz  = at_per_position(randomize_horizontal(matrix, rng_check))
corr_vert   = np.corrcoef(at_real, at_vert)[0, 1]
corr_horiz  = np.corrcoef(at_real, at_horiz)[0, 1]
print(f"  Correlação AT-por-posição real vs VERTICAL:   r = {corr_vert:.4f}  (esperado ~1.0)")
print(f"  Correlação AT-por-posição real vs HORIZONTAL: r = {corr_horiz:.4f}  (esperado < vertical)")
print("CONCLUÍDO. Arquivos gerados em:", OUT_DIR + "/")
print(f"  matrix_real.npy                      — matriz {n_genes} × {UPSTREAM_BP}")
print(f"  gene_ids.txt                         — IDs dos genes")
print(f"  real_sequences.fasta                 — sequências reais")
print(f"  vertical_all.fasta                   — TODAS as permutações verticais concatenadas")
print(f"  horizontal_all.fasta                 — TODAS as permutações horizontais concatenadas")
print(f"  vertical/vertical_perm_001..{N_PERMUTATIONS:03d}.fasta     — randomização vertical (individuais)")
print(f"  horizontal/horizontal_perm_001..{N_PERMUTATIONS:03d}.fasta — randomização horizontal (individuais)")
