
# ordem - manipulação, genômica, interações e gráficos 
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(InteractionSet)
library(patchwork)
library(ggpubr)
 
dir_dados <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/sobreposição"
 
if (!dir.exists(dir_resultados)) {dir.create(dir_resultados, recursive = TRUE)}
 
arquivos <- list.files(dir_dados, full.names = TRUE)
arq_smap <- arquivos[str_detect(arquivos, regex("smap1", ignore_case = TRUE))]
arq_nonb <- arquivos[str_detect(arquivos, regex("z-dna|gquad|kim|nonb|triplex|r-loop|short_tandem_slipped|A-phased|G-quadruplex|Cruciform", ignore_case = TRUE)) &
                     str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]
arq_genes <- arquivos[str_detect(arquivos, regex("Hsalinarum", ignore_case = TRUE)) &
                      str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]
arq_expressao <- arquivos[str_detect(arquivos, regex("bedgraph", ignore_case = TRUE))]
arq_hic <- arquivos[str_detect(arquivos, regex("bedpe", ignore_case = TRUE))]

carregar_genomico <- function(caminho) {
  tryCatch({
    if (str_detect(caminho, regex("\\.gff3?$", ignore_case = TRUE))) {import.gff(caminho)}
    else if (str_detect(caminho, regex("\\.bed$", ignore_case = TRUE))) {import.bed(caminho)}
    else {NULL}}, error = function(e) {
      warning("Erro ao carregar ", basename(caminho), ": ", e$message)
      NULL})}

carregar_hic <- function(arquivo) {
  tryCatch({
    hic_df <- read_tsv(arquivo, col_names = FALSE)
    gr1 <- GRanges(seqnames = hic_df$X1, ranges = IRanges(hic_df$X2 + 1, hic_df$X3))
    gr2 <- GRanges(seqnames = hic_df$X4, ranges = IRanges(hic_df$X5 + 1, hic_df$X6))
    GInteractions(gr1, gr2)}, error = function(e) {
      message("Erro ao carregar Hi-C ", basename(arquivo), ": ", e$message)
      return(NULL)})}

# =====================================================================
# lista nome
dados <- list(
  smap = map(set_names(arq_smap, basename(arq_smap)), carregar_genomico) |> discard(is.null),
  nonb = map(set_names(arq_nonb, basename(arq_nonb)), carregar_genomico) |> discard(is.null),
  genes = map(set_names(arq_genes, basename(arq_genes)), carregar_genomico) |> discard(is.null),
  expressao = arq_expressao, # expressão fica só listado por enquanto
  hic = map(set_names(arq_hic, basename(arq_hic)), carregar_hic) |> discard(is.null))

# =====================================================================
# GRanges doc
estatisticas_gr <- function(gr, nome) {
  if (is.null(gr)) return(NULL)
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
  map2(dados$nonb, names(dados$nonb), estatisticas_gr) %>% bind_rows() %>% mutate(Tipo = "Non-B"),
  map2(dados$genes, names(dados$genes), estatisticas_gr) %>% bind_rows() %>% mutate(Tipo = "Genes"))

write_csv(estatisticas_gerais, file.path(dir_resultados, "estatisticas_gerais.csv"))

# visualizaçao - log
graf1 <- estatisticas_gerais %>%
  ggplot(aes(x = reorder(Nome, Regioes), y = Regioes, fill = Tipo)) +
  geom_col() +
  coord_flip() +
  labs(title = "Número de regiões por arquivo", x = "Arquivo", y = "Regiões") +
  theme_minimal() +
  scale_y_log10()

graf2 <- estatisticas_gerais %>%
  ggplot(aes(x = reorder(Nome, Largura_media), y = Largura_media, fill = Tipo)) +
  geom_col() +
  coord_flip() +
  labs(title = "Largura media por arquivo", x = "Arquivo", y = "Largura média (bp)") +
  theme_minimal() +
  scale_y_log10()

graf_comb <- graf1 / graf2
ggsave(file.path(dir_resultados, "estatisticas_gerais.png"), graf_comb, width = 12, height = 10)

# nonb - nome do arquivo - categorizar o tipo 
#z dna as vezes vem como outro
analisar_nonb <- function(gr, nome) {
  if (is.null(gr)) return(NULL)
  dados_gc_nonb_classificado <- dados_gc_nonb_df %>%
    mutate(Arquivo_Lower = tolower(Arquivo)) %>%
    mutate(Tipo_Estrutura = case_when(
      str_detect(Arquivo_Lower, "z.dna|zdna|z_dna") ~ "Z-DNA",
      str_detect(Arquivo_Lower, "g.quad|gquad|g_quad|g-quad") ~ "G-Quadruplex",
      str_detect(Arquivo_Lower, "triplex") ~ "Triplex",
      str_detect(Arquivo_Lower, "r.loop|rloop|r_loop") ~ "R-loop",
      str_detect(Arquivo_Lower, "cruciform") ~ "Cruciform",
      str_detect(Arquivo_Lower, "a.phase|aphase|a_phase") ~ "A-phased",
      str_detect(Arquivo_Lower, "slipped") ~ "Slipped",
      str_detect(Arquivo_Lower, "short.tandem|short_tandem|shorttandem") ~ "Short Tandem",
      str_detect(Arquivo_Lower, "inverted.repeat|inverted_repeat") ~ "Inverted Repeat",
      TRUE ~ "Outro"
    )) %>%
  as.data.frame(gr) %>%
    summarise(
      Tipo = tipo,
      Arquivo = nome,
      Total = n(),
      Largura_media = mean(width, na.rm = TRUE),
      Densidade_genomica_Mbp = sum(width, na.rm = TRUE) / 1e6)}

estat_nonb <- map2(dados$nonb, names(dados$nonb), analisar_nonb) %>% bind_rows()
write_csv(estat_nonb, file.path(dir_resultados, "nonb_estatisticas.csv"))

graf_nonb <- estat_nonb %>%
  ggplot(aes(x = Tipo, y = Total, fill = Arquivo)) +
  geom_col(position = "dodge") +
  labs(title = "Distribuição de estruturas não-B", x = "Tipo", y = "Regiões") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(dir_resultados, "nonb_distribuicao.png"), graf_nonb, width = 12, height = 8)

# GR1 e GR2 juntos
calc_sobreposicao <- function(gr1, gr2, nome1, nome2) {
  if (is.null(gr1) || is.null(gr2)) return(NULL)
  hits <- findOverlaps(gr1, gr2)
  perc <- (length(unique(queryHits(hits))) / length(gr1)) * 100
  tibble(Conjunto1 = nome1, Conjunto2 = nome2, Regioes_sobrepostas = length(unique(queryHits(hits))), Percentual_sobreposicao = perc)}

lista_sobrep <- list()
for (nm_smap in names(dados$smap)) {
  for (nm_nonb in names(dados$nonb)) {
    tmp <- calc_sobreposicao(dados$smap[[nm_smap]], dados$nonb[[nm_nonb]], nm_smap, nm_nonb)
    if (!is.null(tmp)) {lista_sobrep[[paste(nm_smap, nm_nonb, sep = "_")]] <- tmp}}}

df_sobrep <- bind_rows(lista_sobrep)
write_csv(df_sobrep, file.path(dir_resultados, "sobreposicoes_smap_nonb.csv"))

# Heatmap
graf_heat <- df_sobrep %>%
  ggplot(aes(x = Conjunto1, y = Conjunto2, fill = Percentual_sobreposicao)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", name = "% Sobreposição") +
  labs(title = "SMAP x Non-B (% com sobreposição)", x = "SMAP", y = "Non-B") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(dir_resultados, "heatmap_sobreposicao.png"), graf_heat, width = 12, height = 10)

# percentuais de ambos lados e medida média de largura
calc_sobreposicao_detalhe <- function(gr1, gr2, nome1, nome2) {
  if (is.null(gr1) || is.null(gr2)) return(NULL)
  hits <- findOverlaps(gr1, gr2)
  if (length(hits) == 0) return(NULL)
  qhits <- unique(queryHits(hits))
  shits <- unique(subjectHits(hits))
  perc_gr1 <- length(qhits) / length(gr1) * 100
  perc_gr2 <- length(shits) / length(gr2) * 100
  # largura média das regiões de gr1 que têm qualquer overlap 
  largura_media_overlap <- mean(width(gr1[qhits]), na.rm = TRUE)
  tibble(
    Conjunto1 = nome1,
    Conjunto2 = nome2,
    Regioes_sobrepostas_gr1 = length(qhits),
    Regioes_sobrepostas_gr2 = length(shits),
    Percentual_gr1 = perc_gr1,
    Percentual_gr2 = perc_gr2,
    Largura_media_sobrepostas_gr1 = largura_media_overlap,
    Razao_gr1_sobre_gr2 = ifelse(perc_gr2 > 0, perc_gr1 / perc_gr2, NA_real_))}

lista_sobrep_det <- list()
for (nm_smap in names(dados$smap)) {
  for (nm_nonb in names(dados$nonb)) {
    tmp <- calc_sobreposicao_detalhe(dados$smap[[nm_smap]], dados$nonb[[nm_nonb]], nm_smap, nm_nonb)
    if (!is.null(tmp)) {lista_sobrep_det[[paste(nm_smap, nm_nonb, sep = "_")]] <- tmp}}}

df_sobrep_det <- bind_rows(lista_sobrep_det)
write_csv(df_sobrep_det, file.path(dir_resultados, "sobreposicoes_detalhadas.csv"))

# hic por arquivo
estat_hic <- function(gi, nome) {
  if (is.null(gi)) return(NULL)
  dist_abs <- abs(start(anchors(gi, type = "second")) - start(anchors(gi, type = "first")))
  tibble(
    Arquivo = nome,
    Interacoes = length(gi),
    Distancia_media = mean(dist_abs, na.rm = TRUE),
    Distancia_mediana = median(dist_abs, na.rm = TRUE),
    Distancia_min = suppressWarnings(min(dist_abs, na.rm = TRUE)),
    Distancia_max = suppressWarnings(max(dist_abs, na.rm = TRUE)),
    Cromossomos_envolvidos = length(unique(c(
      as.character(seqnames(anchors(gi, type = "first"))),
      as.character(seqnames(anchors(gi, type = "second")))))))}

df_hic_stats <- map2(dados$hic, names(dados$hic), estat_hic) %>% bind_rows()
if (nrow(df_hic_stats) > 0) {write_csv(df_hic_stats, file.path(dir_resultados, "hic_estatisticas.csv"))}

# histograma - acho que ta errado
if (length(dados$hic) > 0 && !is.null(dados$hic[[1]])) {
  gi_exemplo <- dados$hic[[1]]
  dist_abs <- abs(start(anchors(gi_exemplo, type = "second")) - start(anchors(gi_exemplo, type = "first")))
  dist_df <- tibble(Distancia = dist_abs) %>% filter(!is.na(Distancia))
  if (nrow(dist_df) > 0) {
    graf_hic <- ggplot(dist_df, aes(x = Distancia)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      scale_x_log10() +
      labs(title = paste("Distância entre âncoras -", names(dados$hic)[1]), x = "Distância (bp, log)", y = "Frequência") +
      theme_minimal()
    ggsave(file.path(dir_resultados, "hic_distancias.png"), graf_hic, width = 10, height = 6)}}

if(exists("hic_stats") && nrow(hic_stats) > 0) {
  # Estatísticas de Hi-C
  p_hic_stats <- hic_stats %>%
    ggplot(aes(x = reorder(Arquivo, Total_Interacoes), y = Total_Interacoes)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Número Total de Interações Hi-C",
         x = "Arquivo Hi-C", y = "Número de Interações") +
    theme_minimal()
  
  # Distribuição de distâncias (apenas para o primeiro arquivo)
  if(length(dados$hic) > 0 && !is.null(dados$hic[[1]])) {
    gi <- dados$hic[[1]]
    distancias <- abs(start(anchors(gi, type = "second")) - start(anchors(gi, type = "first")))
    distancias_df <- tibble(Distancia = distancias) %>% filter(Distancia > 0)
    
    p_hic_dist <- ggplot(distancias_df, aes(x = Distancia)) +
      geom_histogram(fill = "steelblue", alpha = 0.7, bins = 50) +
      scale_x_log10() +
      labs(title = paste("Distribuição de Distâncias de Interação -", names(dados$hic)[1]),
           x = "Distância (bp, escala log)", y = "Frequência") +
      theme_minimal()
    
    ggsave(file.path(dir_resultados, "hic_distancias.png"), p_hic_dist, width = 10, height = 6)
  }
  
  ggsave(file.path(dir_resultados, "hic_estatisticas.png"), p_hic_stats, width = 10, height = 6)
}

#  correlação entre variaveis

p_correlacao <- estatisticas_gerais %>%
  ggplot(aes(x = Numero_Regioes, y = Largura_Media, color = Tipo, size = Cobertura_Total)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = cores_tipo) +
  labs(title = "Correlação entre Número de Regiões e Largura Média",
       x = "Número de Regiões (log)", 
       y = "Largura Média (bp, log)",
       size = "Cobertura Total",
       color = "Tipo") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(dir_resultados, "correlacao_regioes_largura.png"), 
       p_correlacao, width = 10, height = 8)

# cromossomo

# Extrair informação de cromossomos
plot_cromossomos <- function(gr, titulo) {
  if(is.null(gr)) return(NULL)
  
  df <- as.data.frame(gr) %>%
    count(seqnames, name = "Count")
  
  ggplot(df, aes(x = seqnames, y = Count, fill = seqnames)) +
    geom_col(alpha = 0.7) +
    labs(title = titulo, x = "Cromossomo", y = "Número de Regiões") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# plots
plots_cromossomos <- list()

if(length(dados$smap) > 0) {
  plots_cromossomos$smap <- plot_cromossomos(dados$smap[[1]], "Distribuição SMAP por Cromossomo")
}

if(length(dados$nonb) > 0) {
  plots_cromossomos$nonb <- plot_cromossomos(dados$nonb[[1]], "Distribuição não-B por Cromossomo")
}

if(length(dados$genes) > 0) {
  plots_cromossomos$genes <- plot_cromossomos(dados$genes[[1]], "Distribuição de Genes por Cromossomo")
}

# Combinar plots 
if(length(plots_cromossomos) > 0) {
  p_cromossomos <- wrap_plots(plots_cromossomos, ncol = 1)
  ggsave(file.path(dir_resultados, "distribuicao_cromossomos.png"), 
         p_cromossomos, width = 10, height = 12)
}

#data
relatorio <- list(
  Arquivos = list(SMAP = length(dados$smap), NonB = length(dados$nonb), Genes = length(dados$genes), Expressao = length(dados$expressao), HiC = length(dados$hic)),
  Estatisticas_gerais = estatisticas_gerais,
  Sobreposicoes_smap_nonb = df_sobrep,
  NonB_estatisticas = estat_nonb)

saveRDS(relatorio, file.path(dir_resultados, "relatorio_analise.rds"))
