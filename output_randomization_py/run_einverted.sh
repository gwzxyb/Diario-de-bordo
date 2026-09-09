#!/bin/bash
# ============================================================================
# EINVERTED - Processa cada sequência separadamente (evita Floating point)
# ============================================================================

THREADS=2               # Número de permutações em paralelo
SEQ_THREADS=4           # Número de sequências por permutação em paralelo
BASE_DIR="."
OUT_BASE="."

EINVERTED=$(which einverted 2>/dev/null)
if [ -z "$EINVERTED" ]; then
    for path in /usr/bin/einverted /usr/local/bin/einverted ~/EMBOSS/bin/einverted; do
        [ -x "$path" ] && EINVERTED="$path" && break
    done
fi
if [ -z "$EINVERTED" ]; then
    echo "ERRO: einverted não encontrado."; exit 1
fi
echo "Usando: $EINVERTED"

EINVERTED_ARGS="-gap 3 -threshold 0 -match 1 -mismatch -1 -maxrepeat 2000"

# ----------------------------------------------------------------------------
# Função para processar um arquivo multifasta, sequência por sequência
process_multifasta() {
    local fasta="$1"
    local outdir="$2"
    local base=$(basename "$fasta" .fasta)
    mkdir -p "$outdir"
    local outfile="$outdir/${base}.out"
    local outseq="$outdir/${base}.fasta"

    # Limpa arquivos de saída anteriores
    > "$outfile"
    > "$outseq"

    echo ">>> Processando $fasta (dividindo em sequências)"

    # Extrai cada sequência em um arquivo temporário e processa
    # Usa awk para gerar um FASTA por registro e chama einverted
    # Isso pode ser feito com um loop, mas usamos parallel para acelerar

    # Criar diretório temporário
    tmpdir=$(mktemp -d)
    trap "rm -rf $tmpdir" EXIT

    # Dividir o multifasta em arquivos individuais usando csplit ou awk
    # Opção 1: usar awk para gerar um arquivo por sequência
    awk -v tmp="$tmpdir" '
        /^>/ { if (seq != "") { print seq > tmp "/seq_" (++n) ".fa"; seq="" }
               header=$0; n++; }
        !/^>/ { seq = seq $0 "\n" }
        END { if (seq != "") print seq > tmp "/seq_" (n+1) ".fa" }
    ' "$fasta"

    # Contar quantas sequências foram geradas
    n_seq=$(ls -1 "$tmpdir"/seq_*.fa 2>/dev/null | wc -l)
    echo "  Encontradas $n_seq sequências"

    # Processar cada sequência em paralelo (usando xargs ou parallel)
    process_one() {
        local seqfile="$1"
        local seqbase=$(basename "$seqfile" .fa)
        # Executa einverted para esta sequência
        $EINVERTED -sequence "$seqfile" \
                   -outfile "$tmpdir/${seqbase}.out" \
                   -outseq "$tmpdir/${seqbase}.seq.fa" \
                   $EINVERTED_ARGS 2>/dev/null
        # Se o comando falhar, apenas avança
        return 0
    }
    export -f process_one
    export EINVERTED EINVERTED_ARGS tmpdir

    # Executa em paralelo (limite SEQ_THREADS por vez)
    find "$tmpdir" -name "seq_*.fa" | sort | \
        xargs -P $SEQ_THREADS -I {} bash -c 'process_one "$@"' _ {}

    # Concatena todos os resultados no arquivo de saída principal
    # O .out do einverted tem cabeçalho próprio; vamos juntar todos
    # e adicionar separadores para saber de qual sequência veio
    for f in $(ls -1 "$tmpdir"/seq_*.out 2>/dev/null | sort); do
        echo "===== $(basename $f .out) =====" >> "$outfile"
        cat "$f" >> "$outfile"
        echo "" >> "$outfile"
    done

    # Concatena as sequências de saída (fasta)
    for f in $(ls -1 "$tmpdir"/seq_*.seq.fa 2>/dev/null | sort); do
        cat "$f" >> "$outseq"
    done

    echo ">>> Concluído $fasta (resultados em $outfile)"
}
export -f process_multifasta
export EINVERTED EINVERTED_ARGS SEQ_THREADS

# ----------------------------------------------------------------------------
# Executar sobre as permutações
echo "=== VERTICAL ==="
find "$BASE_DIR/vertical" -name "vertical_perm_*.fasta" | sort | \
    while read f; do echo "$f $OUT_BASE/vertical_einverted"; done | \
    xargs -P $THREADS -n 2 bash -c 'process_multifasta "$1" "$2"' _

echo "=== HORIZONTAL ==="
find "$BASE_DIR/horizontal" -name "horizontal_perm_*.fasta" | sort | \
    while read f; do echo "$f $OUT_BASE/horizontal_einverted"; done | \
    xargs -P $THREADS -n 2 bash -c 'process_multifasta "$1" "$2"' _

if [ -f "$BASE_DIR/real_sequences.fasta" ]; then
    echo "=== REAL ==="
    process_multifasta "$BASE_DIR/real_sequences.fasta" "$OUT_BASE/real_einverted"
fi

echo "===== EINVERTED CONCLUÍDO ====="