#!/bin/bash
# GFA2 / non-B gfa (suíte de motivos não-B). Uso: bash run_gfa2.sh teste | bash run_gfa2.sh
# Copia o FASTA para uma pasta temporária com nome curto (in.fa) e usa prefixo de saída "g":
# o gfa aborta (SIGABRT, exit 134) quando o caminho/nome de saída é longo.
source "$(dirname "$(readlink -f "$0")")/lib_preditores.sh"

organize_gfa2_outputs() {
    local tmp_gfa_dir="$1" base_out="$2" base="$3"
    shopt -s nocaseglob    # sufixos do GFA são MAIÚSCULOS (_GQ, _IR, _Z...)
    declare -A gfa_map=(
        ["gq"]="05_gfa2_g4"        ["z"]="06_gfa2_zdna"
        ["cruciform"]="07_gfa2_cruciform" ["triplex"]="08_gfa2_triplex"
        ["slipped"]="09_gfa2_slipped"     ["apr"]="10_gfa2_apr"
        ["str"]="11_gfa2_str"      ["ir"]="12_gfa2_inverted_repeats"
        ["mr"]="13_gfa2_mirror_repeats"   ["dr"]="14_gfa2_direct_repeats")
    local key target_dir file ext
    for key in "${!gfa_map[@]}"; do
        target_dir="${base_out}/non_B_structures/${gfa_map[$key]}"
        mkdir -p "$target_dir"
        for file in "$tmp_gfa_dir"/g_${key}.*; do
            if [ -f "$file" ]; then
                ext="${file##*.}"
                mv "$file" "$target_dir/${base}_tool-gfa2_${key}.${ext}" 2>/dev/null
            fi
        done
    done
}

gfa_run() {
    local fasta="$1" out_root="$2" base="$3" tmp="$4"
    if [ ! -x "$GFA2" ]; then echo "gfa não encontrado ou sem permissão em $GFA2"; return 1; fi
    cp "$fasta" "$tmp/in.fa"
    ( cd "$tmp" && "$GFA2" -seq in.fa -out g $GFA2_ARGS )
    local rc=$?
    tmp_listing "$tmp"
    [ $rc -ne 0 ] && return $rc
    local n
    n=$(ls "$tmp"/g_*.tsv 2>/dev/null | xargs -r awk 'FNR>1' | wc -l)
    organize_gfa2_outputs "$tmp" "$out_root" "$base"
    echo "NMOTIFS=$n"
}
export -f organize_gfa2_outputs gfa_run
tool_main gfa2 gfa_run "$1"
