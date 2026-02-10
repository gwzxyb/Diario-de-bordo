#### Non-B DNA — matriz de interseção + GC (script do anexo)

# Objetivo
Classificar arquivos Non-B por tipo, carregar como GRanges e montar uma matriz de interseção entre classes (diagonal = total; fora da diagonal = quantos intervalos de uma classe tocam outra).

Arquivos principais:
- `nonb_classificacao.tsv` (classificação dos arquivos)
- `contagem_real_predicoes.tsv` (número de predições por arquivo/classe)
- `matriz_interseccao_completa.tsv` (matriz classe × classe)

Plots (em `.../plots/`):
- `heatmap_classes_melhorado.png` (heatmap com contagens e % do menor conjunto)
- `violino_gc_por_classe.png` (distribuição do GC por classe, se FASTA disponível)

#### Hi-C vs Non-B DNA (interseções) — sem GC

# Objetivo
Calcular sobreposições entre âncoras de interações Hi-C e regiões Não-B DNA, agrupando por classe (Z-DNA, Triplex, R-loop, Short_tandem, A-phased, Cruciform)

Arquivos principais:
- `nonb_filtrado_classificacao.tsv` (classificação dos arquivos Non-B)
- `contagem_nonb.tsv` e `contagem_hic.tsv` (contagens por arquivo)
- `matriz_hic_nonb_interseccoes.tsv` (interseções por Hi-C, tipo de âncora e classe)
- `relatorio_hic_nonb.txt` (resumo e TOP interseções)
- `hic_data.rds`, `nonb_data.rds`, `interseccoes_data.rds` (objetos salvos)

Plots (em `.../plots/`):
- `heatmap_hic_nonb.png`
- `resumo_classes_nonb.png`

