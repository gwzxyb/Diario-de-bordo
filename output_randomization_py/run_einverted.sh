#!/bin/bash
# einverted / EMBOSS (cruciformes / repetições invertidas). Uso: bash run_einverted.sh teste | bash run_einverted.sh
# Um arquivo por sequência (numeradas 1..N, mesma ordem do FASTA) e o relatório final traz
# "===== seq_N =====" (N = posição da sequência em gene_ids.txt).
source "$(dirname "$(readlink -f "$0")")/lib_preditores.sh"

einv_one() {
    local seqfile="$1" b
    b=$(basename "$seqfile" .fa)
    "$EINVERTED" -sequence "$seqfile" -outfile "$EINV_DIR/$b.out" -outseq "$EINV_DIR/$b.seq.fa" $EINVERTED_ARGS 2>/dev/null
    return 0
}

einv_run() {
    local fasta="$1" out_root="$2" base="$3" tmp="$4"
    local out="$out_root/non_B_structures/04_einverted_cruciform"
    mkdir -p "$out"
    if [ -z "$EINVERTED" ] || [ ! -x "$EINVERTED" ]; then echo "einverted não encontrado (instale o EMBOSS)"; return 1; fi
    local sd="$tmp/seqs"; mkdir -p "$sd"
    awk -v dir="$sd" '
        /^>/ { if (seq != "") { f = dir "/seq_" n ".fa"; printf "%s", seq > f; close(f) } n++; seq = ""; next }
        { seq = seq $0 "\n" }
        END { if (seq != "") { f = dir "/seq_" n ".fa"; printf "%s", seq > f; close(f) } }' "$fasta"
    if [ "$VERBOSE" = "1" ]; then
        echo "--- sondagem do einverted na seq_1 (mensagens de erro, se houver):"
        "$EINVERTED" -sequence "$sd/seq_1.fa" -outfile "$sd/_probe.out" -outseq "$sd/_probe.seq.fa" $EINVERTED_ARGS 2>&1 | head -10
    fi
    export EINV_DIR="$sd"
    find "$sd" -name "seq_*.fa" ! -name "*.seq.fa" | sort -V | xargs -P "$SEQ_THREADS" -I{} bash -c 'einv_one "$1"' _ {}

    local o="$out/${base}_tool-einverted.out" f="$out/${base}_tool-einverted.fasta" x
    : > "$o"; : > "$f"
    for x in $(find "$sd" -name "seq_*.out" | sort -V); do
        echo "===== $(basename "$x" .out) =====" >> "$o"; cat "$x" >> "$o"; echo >> "$o"
    done
    find "$sd" -name "seq_*.seq.fa" | sort -V | xargs -r cat >> "$f"
    echo "NMOTIFS=$(grep -c 'Score' "$o")"
}
export -f einv_one einv_run
tool_main einverted einv_run "$1"
