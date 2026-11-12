#!/bin/bash
# ============================================================================
# RUN ZSEEKER - Leitura Direta de Motivos (Sem Filtro Rígido)
# ============================================================================

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
if [ -f "$SCRIPT_DIR/lib_preditores.sh" ]; then
    source "$SCRIPT_DIR/lib_preditores.sh"
fi

BASE_DIR="${BASE_DIR:-$HOME/codigo/Diario-de-bordo}"
ZSEEKER="${ZSEEKER:-$BASE_DIR/ferramenta_analises_brutas/ZSeeker/zseeker/zseeker.py}"
ZSEEKER_ARGS="${ZSEEKER_ARGS:--threshold 3 --drop_threshold 20}"

zs_run() {
    local fasta="$1" out_root="$2" base="$3" tmp="$4"
    local out="$out_root/non_B_structures/03_zseeker_zdna"
    mkdir -p "$out" "$tmp/zs"

    if [ ! -f "$ZSEEKER" ]; then 
        echo "[ERRO] ZSeeker não encontrado em $ZSEEKER"
        return 1
    fi

    local ZS_PKG_DIR
    ZS_PKG_DIR="$(dirname "$(dirname "$ZSEEKER")")"

    # Executa o ZSeeker com threshold ajustado
    PYTHONPATH="$ZS_PKG_DIR" python3 -m zseeker.zseeker --fasta "$fasta" --output_dir "$tmp/zs" $ZSEEKER_ARGS >/dev/null 2>&1
    local rc=$?

    if [ $rc -ne 0 ]; then
        echo "  [STATUS] ZSeeker FALHOU (exit=$rc)"
        return $rc
    fi

    # Procura qualquer ficheiro de saída (.csv, .tsv ou .bed)
    local ZSEEKER_OUT
    ZSEEKER_OUT=$(find "$tmp/zs" -type f \( -name "*.csv" -o -name "*.tsv" -o -name "*.bed" \) | head -n 1)

    if [ -n "$ZSEEKER_OUT" ] && [ -f "$ZSEEKER_OUT" ]; then
        local dest_csv="$out/${base}_tool-zseeker.csv"
        local dest_bed="$out/${base}_tool-zseeker.bed"

        # Guarda o CSV original bruto
        cp "$ZSEEKER_OUT" "$dest_csv"

        # Converte para BED dinamicamente sem descartar linhas por escore fixo
        awk -F',' 'NR>1 && $1!="" {
            # Converte posições 1-based para 0-based no BED
            s = ($2 ~ /^[0-9]+$/ && $2 > 0) ? $2-1 : $2;
            print $1 "\t" s "\t" $3 "\tZDNA_"NR-1 "\t" $5 "\t+"
        }' "$dest_csv" > "$dest_bed"

        local n_lines
        n_lines=$(awk 'NR>1' "$dest_csv" | wc -l)
        echo "  [STATUS] ZSeeker OK ($n_lines motivos salvos)"
    else
        echo "  [STATUS] ZSeeker OK (SEM MOTIVOS DETECTADOS)"
    fi
}

export -f zs_run

if declare -f tool_main > /dev/null; then
    tool_main zseeker zs_run "$1"
fi