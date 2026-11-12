#!/bin/bash
# QmRLFS-finder (R-loops). Uso: bash run_qmrlfs.sh teste | bash run_qmrlfs.sh
source "$(dirname "$(readlink -f "$0")")/lib_preditores.sh"

qm_run() {
    local fasta="$1" out_root="$2" base="$3" tmp="$4"
    local out="$out_root/non_B_structures/02_qmrlfs_rloop"
    mkdir -p "$out"
    if [ ! -f "$QMRLFS" ]; then echo "QmRLFS-finder não encontrado em $QMRLFS"; return 1; fi
    python3 "$QMRLFS" -i "$fasta" -o "$tmp/q.bed" $QMRLFS_ARGS
    local rc=$?
    tmp_listing "$tmp"
    [ $rc -ne 0 ] && return $rc
    local f
    for f in "$tmp"/q.bed*; do
        [ -f "$f" ] || continue
        mv "$f" "$out/${base}_tool-qmrlfs.bed${f#"$tmp"/q.bed}"
    done
    local n=0
    [ -f "$out/${base}_tool-qmrlfs.bed.out.bed" ] && n=$(wc -l < "$out/${base}_tool-qmrlfs.bed.out.bed")
    echo "NMOTIFS=$n"
}
export -f qm_run
tool_main qmrlfs qm_run "$1"
