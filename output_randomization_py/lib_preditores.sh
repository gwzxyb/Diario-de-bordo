#!/bin/bash
# lib_preditores.sh — funções comuns aos run_<ferramenta>.sh (carregado por eles; não execute direto)
ulimit -c 0     # sem core dumps

PROJ_DIR="${PROJ_DIR:-$HOME/codigo/Diario-de-bordo}"
FERR="$PROJ_DIR/ferramenta_analises_brutas"
G4BOOST="${G4BOOST:-$FERR/G4Boost/G4Boost.py}"
QMRLFS="${QMRLFS:-$FERR/QmRLFS-finder/QmRLFS-finder.py}"
ZSEEKER="${ZSEEKER:-$FERR/ZSeeker/zseeker/zseeker.py}"
GFA2="${GFA2:-$FERR/non-B_gfa/gfa}"
EINVERTED="${EINVERTED:-$(command -v einverted 2>/dev/null)}"
if [ -z "$EINVERTED" ]; then
    for p in /usr/bin/einverted /usr/local/bin/einverted "$HOME/EMBOSS/bin/einverted"; do
        [ -x "$p" ] && EINVERTED="$p" && break
    done
fi

G4BOOST_ARGS="${G4BOOST_ARGS:---minloop 1 --maxloop 7 --minG 3 --maxG 4 --quiet}"
QMRLFS_ARGS="${QMRLFS_ARGS:--m m1,m2 -bed --coordinate_base_zero}"
ZSEEKER_ARGS="${ZSEEKER_ARGS:---threshold 10}"
EINVERTED_ARGS="${EINVERTED_ARGS:--gap 3 -threshold 0 -match 1 -mismatch -1 -maxrepeat 2000}"
GFA2_ARGS="${GFA2_ARGS:--skipWGET}"

JOBS="${JOBS:-4}"              # arquivos em paralelo
SEQ_THREADS="${SEQ_THREADS:-2}" # (einverted) sequências em paralelo dentro de cada arquivo
FORCE="${FORCE:-0}"            # 1 = refaz mesmo o que já tem marcador
SKIP_CDS="${SKIP_CDS:-1}"      # 1 = ignora as permutações de CDS (o densidade_v2.py não usa)

if   [ -d "output_randomization_py" ]; then BASE_ROOT="$(readlink -f output_randomization_py)"
elif [ -d "00_reference" ];            then BASE_ROOT="$(readlink -f .)"
elif [ -d "$PROJ_DIR/output_randomization_py" ]; then BASE_ROOT="$PROJ_DIR/output_randomization_py"
else echo "[ERRO] output_randomization_py não encontrado (rode a partir de Diario-de-bordo)."; exit 1; fi

DIR_REF="$BASE_ROOT/00_reference"
DIR_VERT="$BASE_ROOT/01_vertical_permutations"
DIR_HORIZ="$BASE_ROOT/02_horizontal_permutations"
DIR_CDS="$BASE_ROOT/04_cds_permutations"

export PROJ_DIR G4BOOST QMRLFS ZSEEKER GFA2 EINVERTED G4BOOST_ARGS QMRLFS_ARGS ZSEEKER_ARGS EINVERTED_ARGS GFA2_ARGS SEQ_THREADS

# Lista "fasta pasta_de_saida" — só FASTAs de ENTRADA (nunca arquivos gerados pelas ferramentas)
build_list() {
    local f pair src out sub d
    for f in "$DIR_REF"/*.fasta; do [ -f "$f" ] && echo "$f $DIR_REF"; done
    for pair in "$DIR_VERT/individual_fastas|$DIR_VERT" "$DIR_HORIZ/individual_fastas|$DIR_HORIZ"; do
        src="${pair%%|*}"; out="${pair##*|}"
        [ -d "$src" ] && find "$src" -name "*.fasta" -not -path "*/non_B_structures/*" ! -name "*_tool-*" | sort | while read -r f; do echo "$f $out"; done
    done
    if [ "$SKIP_CDS" != "1" ]; then
        for sub in 01_vertical_permutations 02_horizontal_permutations; do
            d="$DIR_CDS/$sub"
            [ -d "$d" ] && find "$d" -name "*.fasta" -not -path "*/non_B_structures/*" ! -name "*_tool-*" | sort | while read -r f; do echo "$f $d"; done
        done
    fi
}

tmp_listing() {   # no modo teste, mostra o que a ferramenta gerou na pasta temporária
    if [ "$VERBOSE" = "1" ]; then echo "--- arquivos gerados (pasta temporária $1):"; ls -l "$1"; fi
    return 0
}

run_one() {
    local fasta="$1" base_out="$2"
    local base; base=$(basename "$fasta" .fasta)
    local abs_fasta; abs_fasta=$(readlink -f "$fasta")
    local abs_out; abs_out=$(readlink -f "$base_out")
    local marker="$abs_out/.done_${TOOL}/$base"

    if [ "$FORCE" != "1" ] && [ -f "$marker" ]; then echo "[PULADO][$TOOL][$base]"; return 0; fi
    mkdir -p "$abs_out/.done_${TOOL}"

    local tmp; tmp=$(mktemp -d)
    local log="$tmp/_run.log"
    "$TOOL_FN" "$abs_fasta" "$abs_out" "$base" "$tmp" > "$log" 2>&1
    local rc=$?
    # traceback/usage no log = a ferramenta falhou mesmo que o código de saída seja 0
    if [ $rc -eq 0 ] && grep -qE 'Traceback \(most recent call last\)|^usage:' "$log"; then rc=99; fi
    local n; n=$(grep -o 'NMOTIFS=[0-9]*' "$log" | tail -n1 | cut -d= -f2); n="${n:-NA}"

    if [ $rc -ne 0 ]; then
        echo "[FALHOU][$TOOL][$base] exit=$rc"
        { echo "== $base exit=$rc"; tail -n 10 "$log"; } >> "$FAIL_LOG"
    else
        if [ "$n" = "0" ]; then echo "[OK-SEM-MOTIVOS][$TOOL][$base]"; else echo "[OK][$TOOL][$base] motivos=$n"; fi
        touch "$marker"
    fi
    printf '%s\t%s\t%s\n' "$base" "$rc" "$n" >> "$SUMMARY"
    if [ "$VERBOSE" = "1" ]; then echo "--- log completo da ferramenta:"; cat "$log"; fi
    rm -rf "$tmp"
}
export -f run_one tmp_listing

# tool_main <nome> <função_da_ferramenta> [teste]
tool_main() {
    TOOL="$1"; TOOL_FN="$2"; local MODE="$3"
    FAIL_LOG="$BASE_ROOT/falhas_${TOOL}.log"
    SUMMARY="$BASE_ROOT/resumo_${TOOL}.tsv"
    VERBOSE=0
    if [ "$MODE" = "teste" ]; then VERBOSE=1; JOBS=1; fi
    touch "$FAIL_LOG"; : > "$SUMMARY"
    export TOOL TOOL_FN FORCE VERBOSE FAIL_LOG SUMMARY

    local LIST; LIST=$(mktemp)
    build_list > "$LIST"
    if [ "$MODE" = "teste" ]; then
        grep "real_promoter_sequences.fasta" "$LIST" | head -n 1 > "$LIST.t"; mv "$LIST.t" "$LIST"
    fi
    echo "=== $TOOL: $(wc -l < "$LIST") arquivo(s), $JOBS job(s) ==="
    xargs -a "$LIST" -P "$JOBS" -n 2 bash -c 'run_one "$1" "$2"' _
    rm -f "$LIST"

    echo; echo "===== $TOOL CONCLUÍDO ====="
    awk -F'\t' '{t++; if ($2!=0) f++; else if ($3=="0") z++} END{printf "  processados: %d | falhas: %d | sem motivos: %d\n", t, f, z}' "$SUMMARY"
    [ -s "$FAIL_LOG" ] && echo "  falhas detalhadas em: $FAIL_LOG"
    echo "  resumo por arquivo: $SUMMARY"
    return 0
}
