#!/bin/bash
# ============================================================================
# PREDIÇÃO DE ESTRUTURAS NÃO-B SOBRE PERMUTAÇÕES
# ============================================================================

# ============================================================================
# CONFIGURAÇÕES
# ============================================================================

THREADS=8                      # Número de processos simultâneos
BASE_DIR="."                   # Diretório com vertical/, horizontal/, real_sequences.fasta
OUT_BASE="."                   # Diretório de saída (criará subpastas)

# Caminhos para as ferramentas
G4BOOST=~/codigo/Diario-de-bordo/ferramenta_analises_brutas/G4Boost/G4Boost.py
QMRLFS=~/codigo/Diario-de-bordo/ferramenta_analises_brutas/QmRLFS-finder/QmRLFS-finder.py
ZSEEKER="python3 ~/codigo/Diario-de-bordo/ferramenta_analises_brutas/ZSeeker/zseeker.py"

# Parâmetros específicos (incluí --quiet no G4Boost para podermos capturar apenas os dados)
G4BOOST_ARGS="--minloop 1 --maxloop 7 --minG 3 --maxG 4 --quiet"
QMRLFS_ARGS="-m m1,m2 -bed --coordinate_base_zero"
ZSEEKER_ARGS="--threshold 10"

# ----------------------------------------------------------------------------
# FUNÇÃO PARA PROCESSAR UM ARQUIVO COM TODAS AS FERRAMENTAS
# ----------------------------------------------------------------------------

process_file() {
    local fasta="$1"
    local outdir_base="$2"
    local base=$(basename "$fasta" .fasta)

    # 1. G4Boost
    local out_g4="${outdir_base}/g4boost"
    mkdir -p "$out_g4"
    echo "  [G4Boost] $fasta"
    # CORREÇÃO: Usando > para salvar o output ao invés do -o inexistente
    python3 "$G4BOOST" --fasta "$fasta" $G4BOOST_ARGS > "$out_g4/${base}_G4Boost.csv"

    # 2. QmRLFS-finder
    local out_qm="${outdir_base}/qmrlfs"
    mkdir -p "$out_qm"
    echo "  [QmRLFS] $fasta"
    python3 "$QMRLFS" -i "$fasta" -o "$out_qm/${base}_QmRLFS.bed" $QMRLFS_ARGS

    # 3. ZSeeker
    local out_zs="${outdir_base}/zseeker"
    mkdir -p "$out_zs"
    echo "  [ZSeeker] $fasta"
    # CORREÇÃO: Removido o --prefix inexistente
    $ZSEEKER --fasta "$fasta" --output_dir "$out_zs" $ZSEEKER_ARGS
    
    # Busca o arquivo CSV gerado automaticamente pelo ZSeeker para gerar o BED
    ZSEEKER_CSV=$(find "$out_zs" -name "*${base}*.csv" | head -n 1)
    if [ -n "$ZSEEKER_CSV" ] && [ -f "$ZSEEKER_CSV" ]; then
        awk -F',' 'NR>1 && $5 > 0.5 {print $1, $2-1, $3, "ZDNA_"NR-1, $5, "+"}' OFS="\t" \
            "$ZSEEKER_CSV" > "$out_zs/${base}_ZDNA.bed"
    fi
}
export -f process_file
export G4BOOST G4BOOST_ARGS QMRLFS QMRLFS_ARGS ZSEEKER ZSEEKER_ARGS

# ============================================================================
# EXECUÇÃO PRINCIPAL (PARALELIZADA)
# ============================================================================

# --- Vertical ---
echo "=== Processando permutações VERTICAIS ==="
find "$BASE_DIR/vertical" -name "vertical_perm_*.fasta" | sort | \
    while read f; do
        echo "$f $OUT_BASE/vertical"
    done | xargs -P $THREADS -n 2 bash -c 'process_file "$1" "$2"' _

# --- Horizontal ---
echo "=== Processando permutações HORIZONTAIS ==="
find "$BASE_DIR/horizontal" -name "horizontal_perm_*.fasta" | sort | \
    while read f; do
        echo "$f $OUT_BASE/horizontal"
    done | xargs -P $THREADS -n 2 bash -c 'process_file "$1" "$2"' _

# --- Real (opcional) ---
if [ -f "$BASE_DIR/real_sequences.fasta" ]; then
    echo "=== Processando sequências REAIS ==="
    process_file "$BASE_DIR/real_sequences.fasta" "$OUT_BASE/real"
fi

echo "===== PREDIÇÃO CONCLUÍDA ====="
echo "Resultados em: $OUT_BASE/{vertical,horizontal,real}/{g4boost,qmrlfs,zseeker}/"