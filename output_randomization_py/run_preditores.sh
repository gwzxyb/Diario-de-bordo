#!/bin/bash
# ============================================================================
# PREDIÇÃO DE ESTRUTURAS NÃO-B SOBRE PERMUTAÇÕES E SEQUÊNCIAS REAIS
#   1) G4Boost        — Predição de G-Quadruplexes (G4)
#   2) QmRLFS-finder — Predição de R-loops
#   3) ZSeeker       — Predição de Z-DNA
#   4) einverted     — Predição de Cruciforms / Repetições Invertidas
#   5) GFA2          — Suíte Completa de Motivos Não-B
#
# Uso:
#   bash run_preditores.sh              # retoma de onde parou (pula o que já tem marcador .done)
#   FORCE=1 bash run_preditores.sh      # refaz tudo
#   GFA2_SERIAL=0 bash run_preditores.sh  # roda o GFA2 em paralelo (padrão: 1 por vez)
# ============================================================================
THREADS=10
SEQ_THREADS=4

# Não gerar core dumps (o GFA2 estava abortando e enchendo o disco com cores)
ulimit -c 0

# Caminhos para as ferramentas
G4BOOST="$HOME/codigo/Diario-de-bordo/ferramenta_analises_brutas/G4Boost/G4Boost.py"
QMRLFS="$HOME/codigo/Diario-de-bordo/ferramenta_analises_brutas/QmRLFS-finder/QmRLFS-finder.py"
ZSEEKER="$HOME/codigo/Diario-de-bordo/ferramenta_analises_brutas/ZSeeker/zseeker/zseeker.py"
GFA2="$HOME/codigo/Diario-de-bordo/ferramenta_analises_brutas/non-B_gfa/gfa"

if [ -f "$GFA2" ]; then
    chmod +x "$GFA2"
else
    echo "[AVISO] Binário GFA2 não encontrado em $GFA2"
fi

EINVERTED=$(which einverted 2>/dev/null)
if [ -z "$EINVERTED" ]; then
    for path in /usr/bin/einverted /usr/local/bin/einverted ~/EMBOSS/bin/einverted; do
        [ -x "$path" ] && EINVERTED="$path" && break
    done
fi

G4BOOST_ARGS="--minloop 1 --maxloop 7 --minG 3 --maxG 4 --quiet"
QMRLFS_ARGS="-m m1,m2 -bed --coordinate_base_zero"
ZSEEKER_ARGS="--threshold 10"
EINVERTED_ARGS="-gap 3 -threshold 0 -match 1 -mismatch -1 -maxrepeat 2000"
GFA2_ARGS="-skipWGET"

GFA2_SERIAL="${GFA2_SERIAL:-1}"     # 1 = serializa as chamadas do GFA2 (flock)
GFA2_LOCK="/tmp/gfa2_run.lock"
FORCE="${FORCE:-0}"

if [ -d "output_randomization_py" ]; then
    BASE_ROOT="output_randomization_py"
elif [ -d "00_reference" ]; then
    BASE_ROOT="."
else
    echo "[ERRO] Diretório de permutações não encontrado."
    exit 1
fi

DIR_REF="$BASE_ROOT/00_reference"
DIR_VERT="$BASE_ROOT/01_vertical_permutations"
DIR_HORIZ="$BASE_ROOT/02_horizontal_permutations"
DIR_CDS="$BASE_ROOT/04_cds_permutations"

GFA2_FAIL_LOG="$(readlink -f "$BASE_ROOT")/gfa2_failures.log"
touch "$GFA2_FAIL_LOG"

process_einverted_seqs() {
    local fasta="$1"
    local out_einv="$2"
    local base=$(basename "$fasta" .fasta)
    local outfile="$out_einv/${base}_tool-einverted.out"
    local outseq="$out_einv/${base}_tool-einverted.fasta"

    > "$outfile"
    > "$outseq"

    local tmpdir=$(mktemp -d)

    # Um arquivo por sequência, numeração 1..N (um incremento por cabeçalho)
    awk -v tmp="$tmpdir" '
        /^>/ {
            if (seq != "") { f = tmp "/seq_" n ".fa"; printf "%s", seq > f; close(f) }
            n++; seq = ""; next
        }
        { seq = seq $0 "\n" }
        END { if (seq != "") { f = tmp "/seq_" n ".fa"; printf "%s", seq > f; close(f) } }
    ' "$fasta"

    process_one_einv() {
        local seqfile="$1"
        local seqbase=$(basename "$seqfile" .fa)
        $EINVERTED -sequence "$seqfile" \
                   -outfile "$tmpdir/${seqbase}.out" \
                   -outseq "$tmpdir/${seqbase}.seq.fa" \
                   $EINVERTED_ARGS 2>/dev/null
        return 0
    }
    export -f process_one_einv
    export EINVERTED EINVERTED_ARGS tmpdir

    find "$tmpdir" -name "seq_*.fa" ! -name "*.seq.fa" | sort -V | \
        xargs -P $SEQ_THREADS -I {} bash -c 'process_one_einv "$@"' _ {}

    # Concatena em ordem numérica (sort -V), com o número da sequência no cabeçalho
    for f in $(find "$tmpdir" -name "seq_*.out" | sort -V); do
        echo "===== $(basename $f .out) =====" >> "$outfile"
        cat "$f" >> "$outfile"
        echo "" >> "$outfile"
    done

    for f in $(find "$tmpdir" -name "seq_*.seq.fa" | sort -V); do
        cat "$f" >> "$outseq"
    done

    rm -rf "$tmpdir"
}
export -f process_einverted_seqs

organize_gfa2_outputs() {
    local tmp_gfa_dir="$1"
    local base_out="$2"
    local base="$3"

    # Os arquivos do GFA têm sufixo em MAIÚSCULAS (_GQ, _IR, _Z...): ignora caixa no glob
    shopt -s nocaseglob

    declare -A gfa_map=(
        ["gq"]="05_gfa2_g4"
        ["z"]="06_gfa2_zdna"
        ["cruciform"]="07_gfa2_cruciform"
        ["triplex"]="08_gfa2_triplex"
        ["slipped"]="09_gfa2_slipped"
        ["apr"]="10_gfa2_apr"
        ["str"]="11_gfa2_str"
        ["ir"]="12_gfa2_inverted_repeats"
        ["mr"]="13_gfa2_mirror_repeats"
        ["dr"]="14_gfa2_direct_repeats"
    )

    for key in "${!gfa_map[@]}"; do
        local target_dir="${base_out}/non_B_structures/${gfa_map[$key]}"
        mkdir -p "$target_dir"

        for file in "$tmp_gfa_dir"/*_${key}.*; do
            if [ -f "$file" ]; then
                local ext="${file##*.}"
                mv "$file" "$target_dir/${base}_tool-gfa2_${key}.${ext}" 2>/dev/null
            fi
        done
    done
}
export -f organize_gfa2_outputs

process_file() {
    local fasta="$1"
    local base_out="$2"

    if [ -z "$fasta" ] || [ -z "$base_out" ]; then
        return 1
    fi

    local base=$(basename "$fasta" .fasta)
    local abs_fasta=$(readlink -f "$fasta")
    local abs_base_out=$(readlink -f "$base_out")

    # Retomada: pula o que já terminou com sucesso
    local done_dir="${abs_base_out}/.done"
    local done_marker="${done_dir}/${base}"
    mkdir -p "$done_dir"
    if [ "$FORCE" != "1" ] && [ -f "$done_marker" ]; then
        echo "[PULADO] $base (já concluído)"
        return 0
    fi

    local out_g4="${abs_base_out}/non_B_structures/01_g4boost_g4"
    local out_qm="${abs_base_out}/non_B_structures/02_qmrlfs_rloop"
    local out_zs="${abs_base_out}/non_B_structures/03_zseeker_zdna"
    local out_einv="${abs_base_out}/non_B_structures/04_einverted_cruciform"

    mkdir -p "$out_g4" "$out_qm" "$out_zs" "$out_einv"

    # Cada linha de log leva o nome do arquivo (várias permutações rodam em paralelo)
    echo "[INICIO] Processando arquivo: $base"
    local gfa_ok=0

    # 1. G4Boost
    if [ -f "$G4BOOST" ]; then
        (cd "$(dirname "$G4BOOST")" && python3 "G4Boost.py" --fasta "$abs_fasta" $G4BOOST_ARGS) > "$out_g4/${base}_tool-g4boost.csv" 2>/dev/null
        echo "  [STATUS][$base] G4Boost OK"
    else
        echo "  [STATUS][$base] G4Boost NAO ENCONTRADO (PULADO)"
    fi

    # 2. QmRLFS-finder
    if [ -f "$QMRLFS" ]; then
        python3 "$QMRLFS" -i "$abs_fasta" -o "$out_qm/${base}_tool-qmrlfs.bed" $QMRLFS_ARGS >/dev/null 2>&1
        echo "  [STATUS][$base] QmRLFS OK"
    else
        echo "  [STATUS][$base] QmRLFS NAO ENCONTRADO (PULADO)"
    fi

    # 3. ZSeeker
    if [ -f "$ZSEEKER" ]; then
        python3 "$ZSEEKER" --fasta "$abs_fasta" --output_dir "$out_zs" $ZSEEKER_ARGS >/dev/null 2>&1
        ZSEEKER_OUT=$(find "$out_zs" -type f \( -name "*${base}*.csv" -o -name "*${base}*.tsv" -o -name "*${base}*.bed" \) | head -n 1)
        if [ -n "$ZSEEKER_OUT" ] && [ -f "$ZSEEKER_OUT" ]; then
            if [[ "$ZSEEKER_OUT" == *.csv ]]; then
                awk -F',' 'NR>1 && $5 > 0.5 {print $1, $2-1, $3, "ZDNA_"NR-1, $5, "+"}' OFS="\t" \
                    "$ZSEEKER_OUT" > "$out_zs/${base}_tool-zseeker.bed"
                rm -f "$ZSEEKER_OUT"
            elif [[ "$ZSEEKER_OUT" != *"_tool-zseeker.bed" ]]; then
                mv "$ZSEEKER_OUT" "$out_zs/${base}_tool-zseeker.bed"
            fi
            echo "  [STATUS][$base] ZSeeker OK"
        else
            echo "  [STATUS][$base] ZSeeker OK (SEM MOTIVOS DETECTADOS)"
        fi
    else
        echo "  [STATUS][$base] ZSeeker NAO ENCONTRADO EM $ZSEEKER (PULADO)"
    fi

    # 4. einverted
    if [ -n "$EINVERTED" ] && [ -x "$EINVERTED" ]; then
        process_einverted_seqs "$abs_fasta" "$out_einv"
        echo "  [STATUS][$base] einverted OK"
    else
        echo "  [STATUS][$base] einverted NAO ENCONTRADO (PULADO)"
    fi

    # 5. GFA2 (checa o código de saída; prefixo de saída CURTO "g": o gfa aborta (SIGABRT, exit 134) com caminho/nome longo)
    if [ -f "$GFA2" ] && [ -x "$GFA2" ]; then
        local tmp_gfa=$(mktemp -d)
        local rc
        if [ "$GFA2_SERIAL" = "1" ] && command -v flock >/dev/null 2>&1; then
            flock "$GFA2_LOCK" "$GFA2" -seq "$abs_fasta" -out "$tmp_gfa/g" $GFA2_ARGS >/dev/null 2>&1
            rc=$?
        else
            "$GFA2" -seq "$abs_fasta" -out "$tmp_gfa/g" $GFA2_ARGS >/dev/null 2>&1
            rc=$?
        fi
        if [ $rc -ne 0 ]; then
            echo "  [STATUS][$base] GFA2 FALHOU (exit $rc)"
            echo "$base exit=$rc" >> "$GFA2_FAIL_LOG"
        else
            organize_gfa2_outputs "$tmp_gfa" "$abs_base_out" "$base"
            echo "  [STATUS][$base] GFA2 OK"
            gfa_ok=1
        fi
        rm -rf "$tmp_gfa"
    else
        echo "  [STATUS][$base] GFA2 NAO ENCONTRADO OU SEM PERMISSAO (PULADO)"
    fi

    # Só marca como concluído se o GFA2 deu certo (ou nem existe): falhas são refeitas na próxima execução
    if [ $gfa_ok -eq 1 ] || [ ! -x "$GFA2" ]; then
        touch "$done_marker"
    fi

    echo "[CONCLUIDO] $base"
    echo "--------------------------------------------------"
}

export -f process_file
export G4BOOST G4BOOST_ARGS QMRLFS QMRLFS_ARGS ZSEEKER ZSEEKER_ARGS GFA2 GFA2_ARGS EINVERTED EINVERTED_ARGS SEQ_THREADS
export GFA2_SERIAL GFA2_LOCK GFA2_FAIL_LOG FORCE

concatenate_all_outputs() {
    local target_mode_dir="$1"
    local mode_name="$2"

    local non_b_dir="${target_mode_dir}/non_B_structures"
    if [ ! -d "$non_b_dir" ]; then return; fi

    echo ">>> Concatenando relatórios totais para: $mode_name"

    for subfolder in "$non_b_dir"/*; do
        if [ -d "$subfolder" ]; then
            local folder_name=$(basename "$subfolder")

            local total_fasta="${subfolder}/${folder_name}_all_total.fasta"
            > "$total_fasta"
            find "$subfolder" -name "*_tool-einverted.fasta" -exec cat {} + >> "$total_fasta" 2>/dev/null
            if [ ! -s "$total_fasta" ]; then rm -f "$total_fasta"; fi

            local total_bed="${subfolder}/${folder_name}_all_total.bed"
            > "$total_bed"
            find "$subfolder" -name "*.bed" ! -name "*_all_total.bed" -exec cat {} + >> "$total_bed" 2>/dev/null
            if [ ! -s "$total_bed" ]; then rm -f "$total_bed"; fi

            local total_gff="${subfolder}/${folder_name}_all_total.gff"
            > "$total_gff"
            find "$subfolder" -name "*.gff" ! -name "*_all_total.gff" -exec cat {} + >> "$total_gff" 2>/dev/null
            if [ ! -s "$total_gff" ]; then rm -f "$total_gff"; fi

            local total_csv="${subfolder}/${folder_name}_all_total.csv"
            > "$total_csv"
            local first=1
            for csv_file in $(find "$subfolder" -name "*.csv" ! -name "*_all_total.csv" | sort); do
                if [ $first -eq 1 ]; then
                    head -n 1 "$csv_file" > "$total_csv"
                    first=0
                fi
                tail -n +2 "$csv_file" >> "$total_csv"
            done
            if [ ! -s "$total_csv" ]; then rm -f "$total_csv"; fi
        fi
    done
}

echo "=== Processando sequências REAIS (Referência) ==="
for real_fasta in "$DIR_REF"/*.fasta; do
    if [ -f "$real_fasta" ]; then
        process_file "$real_fasta" "$DIR_REF"
    fi
done
concatenate_all_outputs "$DIR_REF" "REAL"

if [ -d "$DIR_VERT/individual_fastas" ]; then
    echo -e "\n=== Processando permutações VERTICAIS (Promotores) ==="
    find "$DIR_VERT/individual_fastas" -name "*.fasta" | sort | \
        while read f; do
            echo "$f $DIR_VERT"
        done | xargs -P $THREADS -n 2 bash -c 'process_file "$1" "$2"' _
    concatenate_all_outputs "$DIR_VERT" "VERTICAL_PROMOTERS"
fi

if [ -d "$DIR_HORIZ/individual_fastas" ]; then
    echo -e "\n=== Processando permutações HORIZONTAIS (Promotores) ==="
    find "$DIR_HORIZ/individual_fastas" -name "*.fasta" | sort | \
        while read f; do
            echo "$f $DIR_HORIZ"
        done | xargs -P $THREADS -n 2 bash -c 'process_file "$1" "$2"' _
    concatenate_all_outputs "$DIR_HORIZ" "HORIZONTAL_PROMOTERS"
fi

if [ -d "$DIR_CDS/01_vertical_permutations" ]; then
    echo -e "\n=== Processando permutações VERTICAIS (CDSs) ==="
    find "$DIR_CDS/01_vertical_permutations" -name "*.fasta" | sort | \
        while read f; do
            echo "$f $DIR_CDS/01_vertical_permutations"
        done | xargs -P $THREADS -n 2 bash -c 'process_file "$1" "$2"' _
    concatenate_all_outputs "$DIR_CDS/01_vertical_permutations" "VERTICAL_CDS"
fi

if [ -d "$DIR_CDS/02_horizontal_permutations" ]; then
    echo -e "\n=== Processando permutações HORIZONTAIS (CDSs) ==="
    find "$DIR_CDS/02_horizontal_permutations" -name "*.fasta" | sort | \
        while read f; do
            echo "$f $DIR_CDS/02_horizontal_permutations"
        done | xargs -P $THREADS -n 2 bash -c 'process_file "$1" "$2"' _
    concatenate_all_outputs "$DIR_CDS/02_horizontal_permutations" "HORIZONTAL_CDS"
fi

echo -e "\n===== PREDIÇÃO DE ESTRUTURAS NÃO-B CONCLUÍDA ====="
if [ -s "$GFA2_FAIL_LOG" ]; then
    echo "[ATENÇÃO] O GFA2 falhou em $(wc -l < "$GFA2_FAIL_LOG") execução(ões). Veja: $GFA2_FAIL_LOG"
fi