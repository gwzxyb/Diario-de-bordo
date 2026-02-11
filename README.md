# Non-B DNA — matriz de interseção + GC (script do anexo)

## Objetivo
Classificar arquivos Non-B por tipo, carregar como GRanges e montar uma matriz de interseção entre classes (diagonal = total; fora da diagonal = quantos intervalos de uma classe tocam outra).

Arquivos principais:
- `nonb_classificacao.tsv` (classificação dos arquivos)
- `contagem_real_predicoes.tsv` (número de predições por arquivo/classe)
- `matriz_interseccao_completa.tsv` (matriz classe × classe)

Plots (em `.../plots/`):
- `heatmap_classes_melhorado.png` (heatmap com contagens e % do menor conjunto)
- `violino_gc_por_classe.png` (distribuição do GC por classe, se FASTA disponível)

# Hi-C vs Non-B DNA (interseções)

## Objetivo
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
  
# Análise SMAP × Genes (sobreposição básica)

## Objetivo
- Carregar arquivos de regiões SMAP e anotações de genes como `GRanges`
- Calcular estatísticas básicas (número de regiões, largura, cobertura) por arquivo
- Quantificar sobreposições entre SMAP e genes e gerar algumas visualizações resumidas

### Saídas principais
- Diretório de resultados: `dir_resultados` (no script: `/home/mcavalcante/sobreposição`), com:
  - `estatisticas_gerais.csv`  
    - Tabela com, por arquivo:
      - `Nome`, `Tipo` (SMAP ou Genes), `Regioes`, `Cromossomos`,
      - `Largura_media`, `Largura_mediana`, `Largura_min`, `Largura_max`, `Cobertura_total`.
  - `sobreposicoes_smap_genes_detalhadas.csv`  
    - Para cada par SMAP × Genes:
      - `Regioes_sobrepostas_gr1` (SMAP), `Regioes_sobrepostas_gr2` (Genes),
      - `Percentual_gr1`, `Percentual_gr2`,
      - `Largura_media_sobrepostas_gr1`,
      - `Razao_gr1_sobre_gr2` (relação dos percentuais).
  - `relatorio_analise_smap_genes.rds`  
    - Lista com:
      - número de arquivos SMAP/Genes/Expressão,
      - `Estatisticas_gerais`,
      - `Sobreposicoes_smap_genes`,
      - lista de arquivos de expressão.
        
### Gráficos (em `dir_plots`)
- `estatisticas_gerais.png`  
  - Painel com dois gráficos (escala log10 no eixo Y):
    - Número de regiões por arquivo.
    - Largura média das regiões por arquivo.
- `distribuicao_cromossomos.png`  
  - Distribuição de regiões por cromossomo para:
    - primeiro arquivo SMAP,
    - primeiro arquivo de genes.
