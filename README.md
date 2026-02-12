# Diário de Bordo

**Autor:** Maria Giovana Cavalcante do Nascimento  
**Orientador:** Profa Dra Tie Koide  
**Data:** 12/02/2026  

**Dados brutos (BED/GFF/Hi-C/SMAP):** Clone branch `data-raw` → copie para `data/igv/` → volte `main`.

[![Reprodutível](https://img.shields.io/badge/FAIR-100%25-blue?style=flat-square)] [![R 4.3.0](https://img.shields.io/badge/R-4.3.0-276DC3?style=flat-square&logo=R)](https://cran.r-project.org/)

## Índice Executivo dos Códigos

### 1. DNA Não-B — Matriz de Interseção + GC Content

**Objetivo**  
Classificar Non-B DNA (Z-DNA, G-quadruplex, Triplex, R-loops, etc.) e montar matriz de interseções entre classes (diagonal=total, off-diagonal=overlap).

**Saídas Principais** (`results/nonb_matrix/`)  
- `nonb_classificacao.tsv` — Classificação dos arquivos  
- `contagem_real_predicoes.tsv` — Número de predições por classe  
- `matriz_interseccao_completa.tsv` — Matriz N×N classes  

**Visualizações** (`results/plots/`)  
- `heatmap_classes_melhorado.png` — Heatmap com contagens e percentuais  
- `violino_gc_por_classe.png` — Distribuição de GC por classe  

**Códigos:** `pipeline_sobreposicoes.R` + `Analise_nonb_contGC.r`

---

### 2. DNA Hi-C vs Não-B (Interseções)

**Objetivo**  
Calcular sobreposições entre Hi-C anchors e regiões Non-B DNA, agrupado por classe (Z-DNA, Triplex, R-loop, Short_tandem, A-phased, Cruciform).

**Arquivos Principais** (`results/hic_nonb/`)  
- `nonb_filtrado_classificacao.tsv` — Non-B classificado  
- `contagem_nonb.tsv` + `contagem_hic.tsv` — Contagens por arquivo  
- `matriz_hic_nonb_interseccoes.tsv` — Interseções por Hi-C + âncora + classe  
- `relatorio_hic_nonb.txt` — Resumo e TOP interseções  
- `hic_data.rds`, `nonb_data.rds`, `intersecoes_data.rds` — Objetos serializados  

**Gráficos** (`results/plots/`)  
- `heatmap_hic_nonb.png`  
- `resumo_classes_nonb.png`  

**Códigos:** `Sobreposicao_nonb_e_hic_heatmap` + `pipeline_sobreposicoes.R`

---

### 3. Análise SMAP × Genes (Sobreposição Regulatória)

**Objetivo**  
Análise completa SMAP1 vs genes H. salinarum: estatísticas de regiões, sobreposições (±1kb promotores), cobertura cromossômica.

**Saídas Principais** (`results/sobreposicao/`)  

**estatisticas_gerais.csv
├── Por arquivo SMAP/Genes:
│ ├── Regioes_Totais, Cromossomos_Cobertos
│ ├── Largura_media/min/max/mediana
│ └── Cobertura_Genoma_Total_bp

sobreposicoes_smap_genes_detalhadas.csv
├── SMAP_ID, Gene_ID, Gene_Name
├── SMAP_chr:start-end, Gene_chr:start-end
├── Percentual_Sobreposicao_SMAP, Percentual_Sobreposicao_Gene
├── Distancia_Promotor_bp, Regiao_Regulatoria
└── Razao_SMAP_Gene

relatorio_analise_smap_genes.rds
├── stats_gerais, sobreposicoes_resumo
├── arquivos_expressao_lista
└── top_10_sobreposicoes**


**Visualizações** (`results/plots/`)
estatisticas_gerais.png
├── Número de regiões por arquivo (eixo Y log10)
└── Largura média por arquivo (eixo Y log10)

distribuicao_cromossomos.png
├── 1º arquivo SMAP por cromossomo
└── 1º arquivo Genes por cromossomo

sobreposicoes_promotores.png
└── Histograma distâncias aos promotores


---

# 1. Configuração do Ambiente
```bash
git clone https://github.com/gwzxyb/Diario-de-bordo.git
cd Diario-de-bordo
R -e "renv::restore()"

git checkout data-raw
mkdir -p ../data/igv/
cp *.bed *.gff *.bedpe ../data/igv/
cd .. && git checkout main
mkdir -p results/{sobreposicao,nonb_matrix,hic_nonb,plots}

