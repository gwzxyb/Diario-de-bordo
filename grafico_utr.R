
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(gggenes)
library(patchwork)

# onde cai no gene
determinar_posicao_relativa <- function(gr_features, gr_genes) {
  if (is.null(gr_features) || is.null(gr_genes) || length(gr_features) == 0) {return(NULL)}
  
  resultados <- list()
  # ver quem encosta em quem 
  overlaps <- findOverlaps(gr_features, gr_genes)  
  if (length(overlaps) == 0) {return(NULL)}
  
  for (i in seq_len(length(overlaps))) {
    feature_idx <- queryHits(overlaps)[i]
    gene_idx <- subjectHits(overlaps)[i]
    feature <- gr_features[feature_idx]
    gene <- gr_genes[gene_idx]
    gene_start <- start(gene); gene_end <- end(gene); gene_length <- width(gene)
    feature_start <- start(feature); feature_end <- end(feature)
    
    # cupstream/downstream ±1kb - regiao promotora, UTR ~10% do gene 
    posicao_relativa <- case_when(
      feature_end < gene_start & feature_start >= (gene_start - 1000) ~ "Promotor",
      feature_start > gene_end & feature_end <= (gene_end + 1000) ~ "Terminator",
      feature_start >= gene_start & feature_end <= (gene_start + gene_length * 0.1) ~ "5' UTR",
      feature_start >= (gene_end - gene_length * 0.1) & feature_end <= gene_end ~ "3' UTR",
      feature_start >= gene_start & feature_end <= (gene_start + gene_length * 0.33) ~ "Início do Gene",
      feature_start >= (gene_start + gene_length * 0.33) & feature_end <= (gene_start + gene_length * 0.66) ~ "Meio do Gene",
      feature_start >= (gene_start + gene_length * 0.66) & feature_end <= gene_end ~ "Fim do Gene",
      feature_start < gene_start & feature_end > gene_end ~ "Spanning_Gene",
      feature_start < gene_start & feature_end <= gene_end ~ "Início_Estendido",
      feature_start >= gene_start & feature_end > gene_end ~ "Fim_Estendido",
      TRUE ~ "Intergênico")
    
    # o que está sobrepondo o gene
    sobreposicao <- intersect(feature, gene)
    fracao_sobreposicao <- ifelse(width(feature) > 0, sum(width(sobreposicao)) / width(feature), 0)
    resultados[[i]] <- tibble(
      Feature_Chrom = as.character(seqnames(feature)),
      Feature_Start = start(feature),
      Feature_End = end(feature),
      Feature_Width = width(feature),
      Gene_Chrom = as.character(seqnames(gene)),
      Gene_Start = start(gene),
      Gene_End = end(gene),
      Gene_Width = width(gene),
      Gene_Strand = as.character(strand(gene)),
      Posicao_Relativa = posicao_relativa,
      Fracao_Sobreposicao = fracao_sobreposicao,
      Distancia_TSS = ifelse(strand(gene) == "+", feature_start - gene_start, gene_end - feature_end),
      Distancia_TES = ifelse(strand(gene) == "+", feature_end - gene_end, gene_start - feature_start) }
  bind_rows(resultados)}

# conjuntos ( tipo - rótulo pra gráficos depois)
analisar_posicoes_genes <- function(dados_genomicos, gr_genes, tipo) {
  res <- list()
  for (nome_arquivo in names(dados_genomicos)) {
    cat("processando:", nome_arquivo, "\n")
    gr_features <- dados_genomicos[[nome_arquivo]]
    if (is.null(gr_features) || length(gr_features) == 0) next
    pos <- determinar_posicao_relativa(gr_features, gr_genes)
    if (!is.null(pos) && nrow(pos) > 0) {
      pos$Arquivo <- nome_arquivo
      pos$Tipo_Feature <- tipo
      res[[nome_arquivo]] <- pos}}
  if (length(res) > 0) {bind_rows(res)} else {NULL}}

# genes: pegando o primeiro da lista
gr_genes <- dados$genes[[1]]
# non-B e smap: mesma lógica
cat("-> Non-B vs genes\n")
posicoes_nonb <- analisar_posicoes_genes(dados$nonb, gr_genes, "Non-B")
cat("-> SMAP vs genes\n")
posicoes_smap <- analisar_posicoes_genes(dados$smap, gr_genes, "SMAP")
posicoes_completas <- bind_rows(posicoes_nonb, posicoes_smap)
write_csv(posicoes_completas, file.path(dir_resultados, "posicoes_relativas_genes.csv"))

#  contagem por posição (percentual dentro de cada tipo)
contagem_posicoes <- posicoes_completas %>%
  group_by(Tipo_Feature, Arquivo, Posicao_Relativa) %>%
  summarise(N = n(), Fracao_Media = mean(Fracao_Sobreposicao, na.rm = TRUE), .groups = "drop") %>%
  group_by(Tipo_Feature, Arquivo) %>%
  mutate(Percentual = N / sum(N) * 100)
write_csv(contagem_posicoes, file.path(dir_resultados, "contagem_posicoes_relativas.csv"))
# grafico 
p_distribuicao <- contagem_posicoes %>%
  ggplot(aes(x = Posicao_Relativa, y = Percentual, fill = Tipo_Feature)) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(title = "Distribuição por posição no gene", x = "Posição", y = "Percentual (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Arquivo, scales = "free_y", ncol = 2)

p_heatmap <- contagem_posicoes %>%
  ggplot(aes(x = Arquivo, y = Posicao_Relativa, fill = Percentual)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Percentual (%)") +
  labs(title = "Heatmap de posições relativas", x = "Arquivo", y = "Posição") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ Tipo_Feature, scales = "free_x")

p_distancia_tss <- posicoes_completas %>%
  filter(abs(Distancia_TSS) < 5000) %>%  # corta extremos 
  ggplot(aes(x = Distancia_TSS, fill = Tipo_Feature)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distância ao TSS (janela ±5kb)", x = "Distância ao TSS (bp)", y = "Densidade") +
  theme_minimal() +
  facet_wrap(~Arquivo, scales = "free_y")

# por tipo de non-B
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
      TRUE ~ "Outro"))
  
  p_estruturas <- posicoes_nonb_classificado %>%
    group_by(Tipo_Estrutura, Posicao_Relativa) %>%
    summarise(N = n(), .groups = "drop") %>%
    group_by(Tipo_Estrutura) %>%
    mutate(Percentual = N / sum(N) * 100) %>%
    ggplot(aes(x = Posicao_Relativa, y = Percentual, fill = Tipo_Estrutura)) +
    geom_col(position = "dodge") +
    labs(title = "Não-B por posição no gene", x = "Posição", y = "Percentual (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(dir_resultados, "posicoes_por_estrutura.png"), p_estruturas, width = 12, height = 8)}
# média global por posição / média das médias
enriquecimento <- contagem_posicoes %>%
  group_by(Posicao_Relativa) %>%
  summarise(Media_Percentual = mean(Percentual, na.rm = TRUE), SD_Percentual = sd(Percentual, na.rm = TRUE), N_Arquivos = n(), .groups = "drop") %>%
  mutate(Enriquecimento_Relativo = Media_Percentual / mean(Media_Percentual))
write_csv(enriquecimento, file.path(dir_resultados, "enriquecimento_posicional.csv"))
ggsave(file.path(dir_resultados, "distribuicao_posicoes.png"), p_distribuicao, width = 14, height = 10)
ggsave(file.path(dir_resultados, "heatmap_posicoes.png"), p_heatmap, width = 12, height = 8)
ggsave(file.path(dir_resultados, "distancia_tss.png"), p_distancia_tss, width = 12, height = 8)

#dá pra plotar setas de genescom gggenes se quiser 
# ver: gggenes::geom_gene_arrow() pra mapa de genes mais visual 
