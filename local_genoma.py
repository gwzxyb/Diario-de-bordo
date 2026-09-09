#!/usr/bin/env python3
import sys
import os
import glob
def ler_mapeamento(mapeamento_file):
    """Lê o arquivo BED com 6 colunas: chrom, start, end, name, score, strand.
    Retorna dict {name: (chrom, start, end, strand)}."""
    mapeamento = {}
    with open(mapeamento_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, score, strand = parts[:6]
            start = int(start)
            end = int(end)
            mapeamento[name] = (chrom, start, end, strand)
    return mapeamento

def converter_gff_para_bed(gff_path, mapeamento, out_dir='.'):
    """Lê um GFF local, converte para BED genômico usando o mapeamento."""
    base = os.path.basename(gff_path)
    nome_sem_ext = os.path.splitext(base)[0]
    out_path = os.path.join(out_dir, nome_sem_ext + '.bed')
    
    with open(gff_path) as gff, open(out_path, 'w') as bed_out:
        for line in gff:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            seq_id = parts[0]          # ex: VNG_RS00350
            feature = parts[2]         # ex: Mirror_Repeat
            local_start = int(parts[3])  # 1‑based
            local_end = int(parts[4])    # 1‑based inclusive
            score = parts[5]
            strand = parts[6]          # geralmente '+' no GFF local, mas usamos o da janela
            
            if seq_id not in mapeamento:
                sys.stderr.write(f"Aviso: '{seq_id}' não encontrado no mapeamento\n")
                continue
            
            chrom, win_start, win_end, win_strand = mapeamento[seq_id]
            # win_start e win_end são do BED (0‑based start, 1‑based end)
            
            if win_strand == '+':
                gen_start = win_start + (local_start - 1)
                gen_end   = win_start + local_end
            else:  # fita negativa → sequência é reverso‑complemento
                gen_start = win_end - local_end
                gen_end   = win_end - (local_start - 1)
            
            # Cria um nome descritivo
            nome_bed = f"{seq_id}_{feature}_{local_start}_{local_end}"
            bed_out.write(f"{chrom}\t{gen_start}\t{gen_end}\t{nome_bed}\t{score}\t{win_strand}\n")
    print(f"Gerado: {out_path}")

def main():
    # Caminhos (use absolutos ou relativos a partir da pasta onde o script é executado)
    base_dir = os.path.dirname(os.path.abspath(__file__))  # ~/codigo/Diario-de-bordo
    
    # Arquivo de mapeamento (gerado pelo tabela_prom.py na pasta dados)
    mapeamento_file = os.path.join(base_dir, "dados", "promotor_gene.bed")
    if not os.path.exists(mapeamento_file):
        sys.stderr.write(f"ERRO: Arquivo de mapeamento '{mapeamento_file}' não encontrado.\n")
        sys.exit(1)
    
    mapeamento = ler_mapeamento(mapeamento_file)
    print(f"Mapeamento carregado: {len(mapeamento)} entradas.")
    
    # Pasta onde estão os GFFs locais
    gff_dir = os.path.join(base_dir, "ferramenta_analises_brutas", "non-B_gfa")
    if not os.path.isdir(gff_dir):
        sys.stderr.write(f"ERRO: Pasta '{gff_dir}' não encontrada.\n")
        sys.exit(1)
    
    # Lista apenas os GFFs que começam com "real_" (exclui vertical_*, horizontal_* e analise_completa*)
    padrao = os.path.join(gff_dir, "real_*.gff")
    arquivos = glob.glob(padrao)
    if not arquivos:
        print("Nenhum arquivo real_*.gff encontrado.")
        sys.exit(1)
    
    print(f"Encontrados {len(arquivos)} arquivos GFF reais para processar.")
    
    # Converter cada um, salvando os BEDs na mesma pasta dos GFFs
    for gff in arquivos:
        converter_gff_para_bed(gff, mapeamento, out_dir=gff_dir)
    print("Concluído.")
if __name__ == "__main__":
    main()
