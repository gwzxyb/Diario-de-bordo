# --------------Pacotes-----------------------

# ordem - manipulação, genômica e gráficos
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(patchwork)

# ------------- legenda -----------

# smap  - regiões/intervalos (BED/GFF)
# genes - anotações de genes (BED/GFF)
# expr  - sinal (bedGraph) -> fica listado (sem importação por enquanto)

# 1) Diretórios e preparação
dir_dados <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/sobreposição"
dir_plots <- file.path(dir_resultados, "plots")
dir.create(dir_resultados, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# 2) Seleção de arquivos (SEM Hi-C e SEM Non-B)
arquivos <- list.files(dir_dados, full.names = TRUE)
arq_smap <- arquivos[str_detect(arquivos, regex("smap1", ignore_case = TRUE))]
arq_genes <- arquivos[
  str_detect(arquivos, regex("Hsalinarum", ignore_case = TRUE)) &
    str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]
arq_expressao <- arquivos[str_detect(arquivos, regex("bedgraph", ignore_case = TRUE))]
cat("\nARQUIVOS ENCONTRADOS:\n")
cat("SMAP:", length(arq_smap), "\n")
cat("Genes:", length(arq_genes), "\n")
cat("Expressão (bedGraph):", length(arq_expressao), "\n\n")
# 3) Função de importação -> GRanges
# Encapsula a importação; retorna NULL em falhas para o pipeline seguir
carregar_genomico <- function(caminho) {
  tryCatch({
    if (str_detect(caminho, regex("\\.gff3?$", ignore_case = TRUE))) {
      import(caminho, format = "GFF")
    } else if (str_detect(caminho, regex("\\.bed$", ignore_case = TRUE))) {
      import(caminho, format = "BED")
    } else { NULL }
  }, error = function(e) {
    warning("Erro ao carregar ", basename(caminho), ": ", e$message)
    NULL })}

# 4) Carregamento em lista (com nomes)
dados <- list(
  smap = map(set_names(arq_smap, basename(arq_smap)), carregar_genomico) |> discard(is.null),
  genes = map(set_names(arq_genes, basename(arq_genes)), carregar_genomico) |> discard(is.null),
  expressao = set_names(arq_expressao, basename(arq_expressao)))
if (length(dados$smap) == 0) warning("Nenhum arquivo SMAP carregado.")
if (length(dados$genes) == 0) warning("Nenhum arquivo de genes carregado.")

# 5) Estatísticas por GRanges
estatisticas_gr <- function(gr, nome) {
  if (is.null(gr) || length(gr) == 0) return(NULL)
  tibble(
    Nome = nome,
    Regioes = length(gr),
    Cromossomos = length(unique(seqnames(gr))),
    Largura_media = mean(width(gr), na.rm = TRUE),
    Largura_mediana = median(width(gr), na.rm = TRUE),
    Largura_min = suppressWarnings(min(width(gr), na.rm = TRUE)),
    Largura_max = suppressWarnings(max(width(gr), na.rm = TRUE)),
    Cobertura_total = sum(width(gr), na.rm = TRUE))}

estatisticas_gerais <- bind_rows(
  map2(dados$smap, names(dados$smap), estatisticas_gr) %>% bind_rows() %>% mutate(Tipo = "SMAP"),
  map2(dados$genes, names(dados$genes), estatisticas_gr) %>% bind_rows() %>% mutate(Tipo = "Genes"))
write_csv(estatisticas_gerais, file.path(dir_resultados, "estatisticas_gerais.csv"))

# 6) Visualizações principais (log)
if (nrow(estatisticas_gerais) > 0) {
  graf1 <- estatisticas_gerais %>%
    ggplot(aes(x = reorder(Nome, Regioes), y = Regioes, fill = Tipo)) +
    geom_col() +
    coord_flip() +
    labs(title = "Número de regiões por arquivo", x = "Arquivo", y = "Regiões") +
    theme_minimal(base_size = 12) +
    scale_y_log10()
  graf2 <- estatisticas_gerais %>%
    ggplot(aes(x = reorder(Nome, Largura_media), y = Largura_media, fill = Tipo)) +
    geom_col() +
    coord_flip() +
    labs(title = "Largura média por arquivo", x = "Arquivo", y = "Largura média (bp)") +
    theme_minimal(base_size = 12) +
    scale_y_log10()
  graf_comb <- graf1 / graf2
  ggsave(file.path(dir_plots, "estatisticas_gerais.png"),
         graf_comb, width = 12, height = 10, dpi = 300, bg = "white")}

# 7) Sobreposições SMAP x Genes
calc_sobreposicao_detalhe <- function(gr1, gr2, nome1, nome2) {
  if (is.null(gr1) || is.null(gr2) || length(gr1) == 0 || length(gr2) == 0) return(NULL)
  hits <- findOverlaps(gr1, gr2, minoverlap = 1)
  if (length(hits) == 0) return(NULL)
  qhits <- unique(queryHits(hits))
  shits <- unique(subjectHits(hits))
  tibble(
    Conjunto1 = nome1,
    Conjunto2 = nome2,
    Regioes_sobrepostas_gr1 = length(qhits),
    Regioes_sobrepostas_gr2 = length(shits),
    Percentual_gr1 = 100 * length(qhits) / length(gr1),
    Percentual_gr2 = 100 * length(shits) / length(gr2),
    Largura_media_sobrepostas_gr1 = mean(width(gr1[qhits]), na.rm = TRUE),
    Razao_gr1_sobre_gr2 = ifelse(length(shits) > 0, (100 * length(qhits) / length(gr1)) / (100 * length(shits) / length(gr2)), NA_real_)  )}

lista_sobrep_det <- list()
for (nm_smap in names(dados$smap)) {
  for (nm_gene in names(dados$genes)) {
    tmp <- calc_sobreposicao_detalhe(dados$smap[[nm_smap]], dados$genes[[nm_gene]], nm_smap, nm_gene)
    if (!is.null(tmp)) lista_sobrep_det[[paste(nm_smap, nm_gene, sep = "_")]] <- tmp  }}
df_sobrep_det <- bind_rows(lista_sobrep_det)
write_csv(df_sobrep_det, file.path(dir_resultados, "sobreposicoes_smap_genes_detalhadas.csv"))

# 8) Distribuição por cromossomo (primeiro arquivo de cada tipo)
plot_cromossomos <- function(gr, titulo) {
  if (is.null(gr) || length(gr) == 0) return(NULL)
  df <- as.data.frame(gr) %>% count(seqnames, name = "Count")
  ggplot(df, aes(x = seqnames, y = Count, fill = seqnames)) +
    geom_col(alpha = 0.7) +
    labs(title = titulo, x = "Cromossomo", y = "Número de regiões") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))}
plots_cromossomos <- list()

if (length(dados$smap) > 0) {
  plots_cromossomos$smap <- plot_cromossomos(dados$smap[[1]], "Distribuição SMAP por cromossomo")}
if (length(dados$genes) > 0) {
  plots_cromossomos$genes <- plot_cromossomos(dados$genes[[1]], "Distribuição de genes por cromossomo")}
if (length(plots_cromossomos) > 0) {
  p_crom <- wrap_plots(plots_cromossomos, ncol = 1)
  ggsave(file.path(dir_plots, "distribuicao_cromossomos.png"),
         p_crom, width = 10, height = 9, dpi = 300, bg = "white")}

# 9) Relatório (objeto .rds)
relatorio <- list(
  Arquivos = list(
    SMAP = length(dados$smap),
    Genes = length(dados$genes),
    Expressao = length(dados$expressao)
  ),
  Estatisticas_gerais = estatisticas_gerais,
  Sobreposicoes_smap_genes = df_sobrep_det,
  Expressao_arquivos = dados$expressao)

saveRDS(relatorio, file.path(dir_resultados, "relatorio_analise_smap_genes.rds"))
cat("\nAnálise concluída.\n")

