#!/bin/bash
# G4Boost (G-quadruplex). Uso: bash run_g4boost.sh teste | bash run_g4boost.sh
# O G4Boost grava os resultados AO LADO do FASTA (<fasta>.gff e <fasta>.g4scores.csv) e não no
# stdout; por isso o script antigo (que redirecionava o stdout) gerava CSVs vazios.
source "$(dirname "$(readlink -f "$0")")/lib_preditores.sh"

g4_run() {
    local fasta="$1" out_root="$2" base="$3" tmp="$4"
    local out="$out_root/non_B_structures/01_g4boost_g4"
    mkdir -p "$out"
    if [ ! -f "$G4BOOST" ]; then echo "G4Boost não encontrado em $G4BOOST"; return 1; fi
    cp "$fasta" "$tmp/in.fasta"
    ( cd "$(dirname "$G4BOOST")" && python3 "$(basename "$G4BOOST")" --fasta "$tmp/in.fasta" $G4BOOST_ARGS ) \
        > "$tmp/g4_stdout.txt" 2> "$tmp/g4_stderr.txt"
    local rc=$?
    echo "--- stderr do G4Boost:"; tail -n 15 "$tmp/g4_stderr.txt"
    tmp_listing "$tmp"
    [ $rc -ne 0 ] && return $rc
    local f ext
    for f in "$tmp"/in.fasta.*; do
        [ -f "$f" ] || continue
        ext="${f#"$tmp"/in.fasta.}"
        mv "$f" "$out/${base}_tool-g4boost.${ext}"
    done
    if [ -s "$tmp/g4_stdout.txt" ]; then mv "$tmp/g4_stdout.txt" "$out/${base}_tool-g4boost.csv"; fi
    local n=0
    [ -f "$out/${base}_tool-g4boost.gff" ] && n=$(grep -vc '^#' "$out/${base}_tool-g4boost.gff")
    echo "NMOTIFS=$n"
}
export -f g4_run
tool_main g4boost g4_run "$1"
