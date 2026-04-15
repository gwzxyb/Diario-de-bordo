library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(patchwork)
library(ggpubr)
library(InteractionSet)
library(RColorBrewer)
library(VennDiagram)
library(grid)
library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(viridis)
library(Biostrings)

# Configurar paralelização
num_cores <- 12
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Paralelização configurada com", num_cores, "cores\n")

# ------------- CORES PADRÃO -------------
cores_padrao <- c(
  "Z-DNA" = "#F28AA0",
  "G-quadruplex" = "#366899",
  "Triplex" = "#81DACA",
  "R-loop" = "#FC8D62",
  "Short_tandem" = "#66C2A5",
  "A-phased" = "#8DA0CB",
  "Cruciform" = "#377EB8",
  "Outros" = "#BEBADA",
  "Com Hi-C" = "#66C2A5",
  "Sem Hi-C" = "#FC8D62",
  "SMAP em genes" = "#81DACA",
  "SMAP fora de genes" = "#F28AA0")

# Cores específicas para Odds Ratio
cor_significativo <- "#E41A1C"
cor_nao_significativo <- "#377EB8"

# Tema 
theme_artigo <- theme_pubr() + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#2C3E50"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#5D6D7E"),
        axis.text = element_text(size = 11, color = "#2C3E50"),
        axis.title = element_text(size = 12, face = "bold", color = "#2C3E50"),
        panel.background = element_rect(fill = "#F8F9F9"),
        plot.background = element_rect(fill = "#FFFFFF"),
        panel.grid.major = element_line(color = "#EAECEE", size = 0.3),
        panel.grid.minor = element_blank())

# 1) diretorios
data_dir <- "/home/mcavalcante/codigo/Diario-de-bordo/dados/igv/modificados"
dir_base <- "/home/mcavalcante/codigo/Diario-de-bordo/resultados Hsalinarum"
dir_resultados <- dir_base
dir_plots <- file.path(dir_resultados, "plots")
dir_tabelas <- file.path(dir_resultados, "tabelas")  
dir_intersecoes <- file.path(dir_resultados, "intersecoes_nonb")
dir_analise_nonb <- file.path(dir_base, "analise_nonB")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tabelas, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_intersecoes, showWarnings = FALSE, recursive = TRUE)

# 2) Genoma referencia
fa_path <- "/home/mcavalcante/codigo/Diario-de-bordo/dados/Hsalinarum.fa"
if (file.exists(fa_path)) {
  fa <- readDNAStringSet(fa_path)
  genome_size <- sum(width(fa))
  cat("Genoma carregado:", length(fa), "sequências | Tamanho:", format(genome_size, big.mark=","), "bp\n")
} else {
  warning("Arquivo FASTA não encontrado.")
  fa <- NULL
  genome_size <- 2000000 }

# 3) pegar do script ANALISE_NONB
cat("  Diretório:", dir_analise_nonb, "\n")
# 3.1) Non-B por classe (interseção)
arquivo_nonb_intersec <- file.path(dir_analise_nonb, "nonb_por_classe_intersecao.rds")
if(file.exists(arquivo_nonb_intersec)) {
  nonb_por_classe <- readRDS(arquivo_nonb_intersec)
  cat("nonb_por_classe (interseção):", length(nonb_por_classe), "classes\n")
  for(classe in names(nonb_por_classe)) {
    cat("    -", classe, ":", length(nonb_por_classe[[classe]]), "regiões\n")
  }
} else {
  stop("Arquivo nonb_por_classe_intersecao.rds não encontrado!")}

# 3.2) Non-B por classe (todos arquivos combinados)
arquivo_nonb_combinado <- file.path(dir_analise_nonb, "nonb_gr_por_classe.rds")
if(file.exists(arquivo_nonb_combinado)) {
  gr_por_classe <- readRDS(arquivo_nonb_combinado)
  cat("gr_por_classe (combinado):", length(gr_por_classe), "classes\n")
} else {
  gr_por_classe <- nonb_por_classe
  cat("Usando nonb_por_classe como gr_por_classe\n")}

# 3.3) Promotores
arquivo_promotores <- file.path(dir_analise_nonb, "promotores.rds")
if(file.exists(arquivo_promotores)) {
  promotores <- readRDS(arquivo_promotores)
  cat("Promotores:", length(promotores), "regiões\n")
} else {
  promotores <- GRanges()
  cat("Promotores não encontrados\n")}

# 3.4) 5' UTR custom
arquivo_utr5 <- file.path(dir_analise_nonb, "utr5_custom.rds")
if(file.exists(arquivo_utr5)) {
  utr5_custom <- readRDS(arquivo_utr5)
  cat("5' UTR:", length(utr5_custom), "regiões\n")
} else {
  utr5_custom <- GRanges()
  cat("5' UTR não encontrados\n")}

# 4) Carregar arquivos do principaal
# 4.1) SMAP
arquivo_smap <- file.path(dir_resultados, "smap_gr.rds")
if(file.exists(arquivo_smap)) {
  smap_gr <- readRDS(arquivo_smap)
  cat(" SMAP:", length(smap_gr), "regiões\n")
} else {
  # Carregar SMAP do diretório original
  arquivos <- list.files(data_dir, full.names = TRUE)
  arq_smap <- arquivos[str_detect(arquivos, regex("smap1", ignore_case = TRUE))]
  if(length(arq_smap) > 0) {
    smap_gr <- carregar_granges(arq_smap[1])
    cat("SMAP carregado do original:", length(smap_gr), "regiões\n")
  } else {
    stop("Arquivo SMAP não encontrado!")}}

# 4.2) Genes
arquivo_genes <- file.path(dir_resultados, "genes_gr.rds")
if(file.exists(arquivo_genes)) {
  genes_gr <- readRDS(arquivo_genes)
  cat("Genes:", length(genes_gr), "\n")
} else {
  # Carregar GFF
  arq_gff <- arquivos[str_detect(arquivos, regex("Hsalinarum.*gff", ignore_case = TRUE))]
  if(length(arq_gff) > 0) {
    gff_gr <- carregar_granges(arq_gff[1])
    mcols_gff <- as.data.frame(mcols(gff_gr))
    genes_gr <- gff_gr[mcols_gff$type == "gene"]
    cat(" Genes carregado do original:", length(genes_gr), "\n")
  } else {
    stop("Arquivo GFF não encontrado!")}}

# 4.3) TnpB
arquivo_tnpb <- file.path(dir_resultados, "vng_gr.rds")
if(file.exists(arquivo_tnpb)) {
  vng_gr <- readRDS(arquivo_tnpb)
  cat(" TnpB:", length(vng_gr), "loci\n")
} else {
  # Buscar TnpB no GFF
  vng_loci <- c("tnpb")
  vng_gr <- GRanges()
  for(locus in vng_loci) {
    hits <- grepl(locus, mcols_gff$locus_tag, ignore.case=TRUE) |
      grepl(locus, mcols_gff$ID, ignore.case=TRUE) |
      grepl(locus, mcols_gff$product, ignore.case=TRUE)
    if(any(hits)) vng_gr <- c(vng_gr, gff_gr[hits])}
  vng_gr <- unique(vng_gr)
  cat(" TnpB carregado do original:", length(vng_gr), "loci\n")}

# 4.4) Hi-C
arquivo_hic <- file.path(dir_resultados, "hic_anchors.rds")
if(file.exists(arquivo_hic)) {
  hic_anchors_gr <- readRDS(arquivo_hic)
  cat(" Hi-C âncoras:", length(hic_anchors_gr), "\n")
} else {
  arq_hic_bedpe <- arquivos[str_detect(arquivos, regex("\\.bedpe$", ignore_case = TRUE))]
  if(length(arq_hic_bedpe) > 0) {
    hic_bedpe <- read.table(arq_hic_bedpe[1], header=FALSE, stringsAsFactors=FALSE)
    colnames(hic_bedpe)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    anchors1 <- GRanges(hic_bedpe$chr1, IRanges(hic_bedpe$start1, hic_bedpe$end1))
    anchors2 <- GRanges(hic_bedpe$chr2, IRanges(hic_bedpe$start2, hic_bedpe$end2))
    hic_anchors_gr <- unique(c(anchors1, anchors2))
    cat(" Hi-C carregado do original:", length(hic_anchors_gr), "âncoras\n")
  } else {
    hic_anchors_gr <- GRanges()
    cat("Hi-C não encontrado\n")}}

# 5) carregar GRanges 
carregar_granges <- function(arquivo) {
  if (!file.exists(arquivo)) return(GRanges())
  extensao <- tools::file_ext(arquivo)
  tryCatch({
    if (extensao %in% c("bed", "bedgraph")) {
      gr <- import(arquivo, format = "BED")
    } else if (extensao %in% c("gff", "gff3")) {
      gr <- import(arquivo, format = "GFF")
    } else {
      return(GRanges())
    }
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr <- gr[width(gr) > 0]
    return(gr)
  }, error = function(e) return(GRanges()))}

# ============================================================
# BLOCO 1: ANÁLISES COM HI-C
if(length(hic_anchors_gr) > 0) {
  # 1.1) Hi-C vs Non-B (PINGS TOTAIS)]
  #"PING" é uma contagem bruta de eventos de sobreposição (um ping = um par âncora ↔ estrutura não-B)
  if(length(nonb_por_classe) > 0) {
    n_hic_total <- length(hic_anchors_gr)
    tabela_consolidada_nonB <- data.frame()
    for (classe in names(nonb_por_classe)) {
      gr_NB <- nonb_por_classe[[classe]]
      hits <- findOverlaps(hic_anchors_gr, gr_NB)
      n_ancoras_unicas <- length(unique(queryHits(hits)))
      percentual_presenca <- (n_ancoras_unicas / n_hic_total) * 100
      total_pings <- length(subjectHits(hits))
      media_pings <- ifelse(n_ancoras_unicas > 0, round(total_pings / n_ancoras_unicas, 2), 0)
      tabela_consolidada_nonB <- rbind(tabela_consolidada_nonB, data.frame(
        Classe = classe,
        Ancoras_Unicas = n_ancoras_unicas,
        Percentual_Presenca = round(percentual_presenca, 2),
        Total_Pings_Interseccao = total_pings,
        Media_Pings_por_Ancora = media_pings,
        stringsAsFactors = FALSE))
      cat(sprintf("  %-15s | Únicas: %5d | Pings: %6d | Média: %.2f\n", 
                  classe, n_ancoras_unicas, total_pings, media_pings))}
    
    if(nrow(tabela_consolidada_nonB) > 0) {
      tabela_consolidada_nonB <- tabela_consolidada_nonB[order(-tabela_consolidada_nonB$Ancoras_Unicas), ]
      write.csv(tabela_consolidada_nonB, file.path(dir_tabelas, "analise_final_unicas_vs_pings_consenso.csv"), row.names = FALSE)
      tabela_plot_pings <- tabela_consolidada_nonB %>%
        arrange(Total_Pings_Interseccao)
      p_pings_totais <- ggplot(tabela_plot_pings, 
                               aes(x = factor(Classe, levels = Classe), y = Total_Pings_Interseccao, fill = Classe)) +
        geom_col(alpha = 0.9, width = 0.7, color = "black", linewidth = 0.3) +
        geom_text(aes(label = paste0(Total_Pings_Interseccao, " (avg: ", Media_Pings_por_Ancora, ")")), 
                  hjust = -0.1, size = 4, fontface = "bold") +
        coord_flip() +
        scale_fill_manual(values = cores_padrao) +
        labs(title = "Hi-C vs Non-B: Total de Eventos de Sobreposição",
             subtitle = "Volume total de interseções (incluindo múltiplas ocorrências por âncora)",
             x = "Classe Non-B",  y = "Número total de pares (âncora ↔ estrutura não-B)") +
        theme_artigo +
        theme(legend.position = "none") +
        expand_limits(y = max(tabela_plot_pings$Total_Pings_Interseccao) * 1.3)
      
      ggsave(file.path(dir_plots, "barras_hic_pings_totais.png"), 
             p_pings_totais, width = 10, height = 7, dpi = 300)
      cat("  Gráfico de pings totais Hi-C vs Non-B salvo\n")}}
  
  # 1.2) SMAP vs Hi-C (PIZZA)
  smap_hic_int <- findOverlaps(smap_gr, hic_anchors_gr)
  smap_com_hic <- unique(queryHits(smap_hic_int))
  
  if(length(smap_hic_int) > 0) {
    df_smap_hic_coord <- data.frame(
      SMAP_seq = as.character(seqnames(smap_gr[queryHits(smap_hic_int)])),
      SMAP_start = start(smap_gr[queryHits(smap_hic_int)]),
      SMAP_end = end(smap_gr[queryHits(smap_hic_int)]),
      Anchor_seq = as.character(seqnames(hic_anchors_gr[subjectHits(smap_hic_int)])),
      Anchor_start = start(hic_anchors_gr[subjectHits(smap_hic_int)]),
      Anchor_end = end(hic_anchors_gr[subjectHits(smap_hic_int)]),
      stringsAsFactors = FALSE)
    write.csv(df_smap_hic_coord, file.path(dir_tabelas, "SMAP_HiC_interseccao_coordenadas.csv"), row.names = FALSE)}
  
  df_hic_status <- data.frame(
    Status = c("SMAP com interseção Hi-C", "SMAP sem interseção Hi-C"),
    Contagem = c(length(smap_com_hic), length(smap_gr) - length(smap_com_hic)),
    Percentual = c(round(100 * length(smap_com_hic) / length(smap_gr), 1),
                   round(100 * (length(smap_gr) - length(smap_com_hic)) / length(smap_gr), 1)))
  write.csv(df_hic_status, file.path(dir_tabelas, "SMAP_HiC_resumo.csv"), row.names = FALSE)
  
  p_smap_hic <- ggplot(df_hic_status, aes(x = 2, y = Contagem, fill = Status)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) + coord_polar(theta = "y") + xlim(0.5, 2.5) +
    geom_text(aes(label = paste0(Contagem, " (", Percentual, "%)")), 
              position = position_stack(vjust = 0.5), size = 4.5, fontface = "bold", color = "white") +
    scale_fill_manual(values = c(cores_padrao["Com Hi-C"], cores_padrao["Sem Hi-C"])) +
    labs(title = "SMAP vs Hi-C",
         subtitle = paste("Total SMAP:", length(smap_gr), "| Total âncoras Hi-C:", length(hic_anchors_gr)),
         fill = "Categoria") + 
    theme_void() +
    theme(legend.position = "right", 
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#5D6D7E"))
  
  ggsave(file.path(dir_plots, "SMAP_vs_HiC_pizza.png"), p_smap_hic, width = 10, height = 7, dpi = 300)
  cat("Gráfico SMAP vs Hi-C salvo\n")
  
  # 1.3) Distância SMAP-HiC  
  dist_smap_hic <- distanceToNearest(smap_gr, hic_anchors_gr)
  dist_vec <- mcols(dist_smap_hic)$distance / 1000
  df_dist <- data.frame(Distancia_kb = dist_vec)
  media <- mean(dist_vec, na.rm = TRUE)
  mediana <- median(dist_vec, na.rm = TRUE)
  max_y <- max(hist(dist_vec, plot=FALSE, breaks=40)$counts)
  
  p_dist_hic <- ggplot(df_dist, aes(x = Distancia_kb)) +
    geom_histogram(bins = 40, fill = cores_padrao["G-quadruplex"], alpha = 0.7, color = "white", linewidth = 0.3) +
    geom_vline(xintercept = media, color = cores_padrao["R-loop"], linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = mediana, color = cores_padrao["Z-DNA"], linetype = "dotted", linewidth = 1.2) +
    annotate("text", x = media + 3, y = max_y * 0.9, 
             label = paste("Média =", round(media, 2), "kb"), 
             color = cores_padrao["R-loop"], size = 4, hjust = 0, fontface = "bold") +
    annotate("text", x = mediana + 3, y = max_y * 0.8, 
             label = paste("Mediana =", round(mediana, 2), "kb"), 
             color = cores_padrao["Z-DNA"], size = 4, hjust = 0, fontface = "bold") +
    labs(title = "Distância de SMAP à âncora Hi-C mais próxima",
         subtitle = paste("Total SMAP:", length(smap_gr), "| Total âncoras Hi-C:", length(hic_anchors_gr)),
         x = "Distância (kb)",  y = "Frequência (número de SMAP)") +
    theme_artigo + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ggsave(file.path(dir_plots, "distribuicao_distancias_SMAP_HiC.png"), p_dist_hic, width = 12, height = 7, dpi = 300)
  cat("  Gráfico de distâncias SMAP-HiC salvo\n")}
# 1.4)  Âncoras Únicas com Sobreposição 
# Preparar dados para o gráfico
df_plot_ancoras <- tabela_consolidada_nonB %>%
  arrange(desc(Ancoras_Unicas)) %>%
  mutate(Label = paste0(Ancoras_Unicas, " (", Percentual_Presenca, "%)"))
# Cores específicas para este gráfico
cores_ancoras <- c(
  "Cruciform" = "#387EB8",
  "R-loop" = "#FC8D62", 
  "Short_tandem" = "#66C2A5",
  "Z-DNA" = "#F28AA0",
  "G-quadruplex" = "#366899",
  "A-phased" = "#8DA0CB",
  "Triplex" = "#81DBCA")
# Gráfico de barras horizontal
p_ancoras_unicas <- ggplot(df_plot_ancoras, 
                           aes(x = reorder(Classe, Ancoras_Unicas), y = Ancoras_Unicas,  fill = Classe)) +
  geom_col(alpha = 0.9, width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Label), 
            hjust = -0.1, 
            size = 4.5, 
            fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = cores_ancoras) +
  labs(title = "Âncoras Hi-C com sobreposição a estruturas Non-B",
       subtitle = paste("Total de âncoras Hi-C analisadas"),
       x = "Classe Non-B",
       y = "Número de Âncoras Únicas com Sobreposição") +
  theme_artigo +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12, face = "bold")) +
  expand_limits(y = max(df_plot_ancoras$Ancoras_Unicas) * 1.2)
ggsave(file.path(dir_plots, "hic_ancoras_nonb_sobreposicao.png"), 
       p_ancoras_unicas, width = 10, height = 7, dpi = 300)

# 1.5) Heatmap de Enriquecimento (Chance de Non-B estar no Hi-C)    
# Calcular tamanho total do genoma em bp
genome_bp <- genome_size
# Calcular tamanhos totais e cobertura do genoma
tamanhos_classes <- sapply(gr_por_classe, function(x) sum(width(x)))
cobertura_classes <- tamanhos_classes / genome_bp
# Para Hi-C, calcular cobertura total das âncoras
if(length(hic_anchors_gr) > 0) {
  tamanho_hic <- sum(width(hic_anchors_gr))
  cobertura_hic <- tamanho_hic / genome_bp
} else {
  tamanho_hic <- 0
  cobertura_hic <- 0 }
for(classe in names(gr_por_classe)) {
  cat(sprintf("    %-15s: %.4f%% (%d bp)\n", classe, cobertura_classes[classe] * 100, tamanhos_classes[classe]))}
cat(sprintf("    %-15s: %.4f%% (%d bp)\n", 
            "Hi-C", cobertura_hic * 100, tamanho_hic))
# Criar data frame para o gráfico de fração
df_fracao <- data.frame()
# Calcular sobreposição com Hi-C para cada classe
for(classe in names(gr_por_classe)) {
  gr_i <- gr_por_classe[[classe]]
  if(length(hic_anchors_gr) > 0) {
    hits <- findOverlaps(gr_i, hic_anchors_gr)
    if(length(hits) > 0) {
      sobreposicao_bp <- sum(width(pintersect(gr_i[queryHits(hits)], hic_anchors_gr[subjectHits(hits)])))
    } else { sobreposicao_bp <- 0}
    
    fracao_no_hic <- sobreposicao_bp / tamanhos_classes[classe]
    esperado <- tamanhos_classes[classe] * cobertura_hic
    enriquecimento <- ifelse(esperado > 0, sobreposicao_bp / esperado, NA)
  } else {
    sobreposicao_bp <- 0
    fracao_no_hic <- 0
    enriquecimento <- NA}
  df_fracao <- rbind(df_fracao, data.frame(
    Classe = classe,
    Tamanho_bp = tamanhos_classes[classe],
    Cobertura_Genoma = cobertura_classes[classe] * 100,
    Sobreposicao_HiC_bp = sobreposicao_bp,
    Fracao_no_HiC = fracao_no_hic * 100,
    Enriquecimento = enriquecimento,
    stringsAsFactors = FALSE))}
df_fracao <- df_fracao[order(-df_fracao$Fracao_no_HiC), ]
df_fracao$Label_fracao <- paste0(round(df_fracao$Fracao_no_HiC, 1), "%")
df_fracao$Label_enrich <- ifelse(!is.na(df_fracao$Enriquecimento), sprintf("%.2f", df_fracao$Enriquecimento), "NA")

# 1.6) Fração de cada classe no Hi-C
p_fracao <- ggplot(df_fracao, aes(x = reorder(Classe, Fracao_no_HiC), y = Fracao_no_HiC, fill = Classe)) +
  geom_col(alpha = 0.9, width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Label_fracao), hjust = -0.1, size = 4, fontface = "bold") + coord_flip() +
  scale_fill_manual(values = cores_ancoras) +
  labs(title = "Fração de cada estrutura Non-B sobreposta a âncoras Hi-C",
       subtitle = "Percentual do comprimento total da classe que está em âncoras Hi-C",
       x = "Classe Non-B",
       y = "Fração sobreposta (%)") +
  theme_artigo +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  expand_limits(y = max(df_fracao$Fracao_no_HiC) * 1.2)
ggsave(file.path(dir_plots, "fracao_nonb_no_hic.png"), 
       p_fracao, width = 10, height = 7, dpi = 300)

#1.7 Heatmap de Enriquecimento (Fold Enrichment)

# Calcular matriz de enriquecimento entre classes Non-B
n_classes <- length(names(gr_por_classe))
nomes_classes <- names(gr_por_classe)

matriz_enriquecimento <- matrix(1, nrow = n_classes, ncol = n_classes, dimnames = list(nomes_classes, nomes_classes))
for(i in 1:n_classes) {
  for(j in 1:n_classes) {
    if(i != j) {
      gr_i <- gr_por_classe[[nomes_classes[i]]]
      gr_j <- gr_por_classe[[nomes_classes[j]]]
      hits <- findOverlaps(gr_i, gr_j)
      
      if(length(hits) > 0) {
        sobreposicao <- sum(width(pintersect(gr_i[queryHits(hits)], gr_j[subjectHits(hits)])))
        esperado <- (tamanhos_classes[nomes_classes[j]] / genome_bp) * tamanhos_classes[nomes_classes[i]]
        enriquecimento <- ifelse(esperado > 0, sobreposicao / esperado, 1)
        matriz_enriquecimento[i, j] <- enriquecimento
      } else {
        matriz_enriquecimento[i, j] <- 0 }}}}

# Adicionar Hi-C como uma dimensão extra
if(length(hic_anchors_gr) > 0) {
  # Criar matriz expandida
  matriz_expandida <- matrix(1, nrow = n_classes + 1, ncol = n_classes + 1,
                             dimnames = list(c(nomes_classes, "Hi-C"),  c(nomes_classes, "Hi-C")))
  
  # Preencher com valores das classes
  matriz_expandida[1:n_classes, 1:n_classes] <- matriz_enriquecimento
  # Preencher linha/coluna Hi-C
  for(i in 1:n_classes) {
    # Hi-C vs classe i
    enriquecimento_hic <- df_fracao$Enriquecimento[i]
    matriz_expandida["Hi-C", nomes_classes[i]] <- ifelse(is.na(enriquecimento_hic), 1, enriquecimento_hic)
    matriz_expandida[nomes_classes[i], "Hi-C"] <- ifelse(is.na(enriquecimento_hic), 1, enriquecimento_hic)}
  matriz_expandida["Hi-C", "Hi-C"] <- 1
  matriz_plot <- matriz_expandida
} else { matriz_plot <- matriz_enriquecimento}
write.csv(matriz_plot, 
          file.path(dir_tabelas, "matriz_enriquecimento_fold.csv"), row.names = TRUE)

# Preparar dados para heatmap
matriz_melt <- as.data.frame(as.table(matriz_plot))
colnames(matriz_melt) <- c("Classe1", "Classe2", "Fold_Enrichment")
# Remover valores infinitos e NAs
matriz_melt <- matriz_melt[!is.na(matriz_melt$Fold_Enrichment), ]
matriz_melt$Fold_Enrichment[is.infinite(matriz_melt$Fold_Enrichment)] <- 0
matriz_melt$Fold_Enrichment[matriz_melt$Fold_Enrichment > 10] <- 10  # Cap para visualização
# Heatmap
p_heatmap <- ggplot(matriz_melt, 
                    aes(x = Classe2, y = Classe1, fill = log2(Fold_Enrichment + 0.1))) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(Fold_Enrichment > 0 & Fold_Enrichment != 1, sprintf("%.1f", Fold_Enrichment), ifelse(Fold_Enrichment == 1, "1", "0"))), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#E41A1C", midpoint = 0, name = "log2(Fold Enrichment)") +
  labs(title = "Enriquecimento de sobreposição entre estruturas Non-B e Hi-C",
       subtitle = "Fold enrichment = observado / esperado (baseado na cobertura do genoma)", x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  coord_fixed()
ggsave(file.path(dir_plots, "heatmap_enriquecimento_nonb_hic.png"), 
       p_heatmap, width = 10, height = 8, dpi = 300, bg = "white")
# 1.8) Barras de Enriquecimento (Hi-C vs cada classe)
df_enrich <- df_fracao[!is.na(df_fracao$Enriquecimento), ]
df_enrich <- df_enrich[order(-df_enrich$Enriquecimento), ]
df_enrich$Significativo <- ifelse(df_enrich$Enriquecimento > 1.5, "Enriquecido", ifelse(df_enrich$Enriquecimento < 0.67, "Depletado", "Neutro"))

p_enrich_bars <- ggplot(df_enrich, aes(x = reorder(Classe, Enriquecimento), y = Enriquecimento, fill = Significativo)) +
  geom_col(alpha = 0.9, width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Label_enrich), 
            hjust = -0.1, 
            size = 4, 
            fontface = "bold") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Enriquecido" = "#E41A1C", 
                               "Depletado" = "#2166AC", 
                               "Neutro" = "gray70")) +
  labs(title = "Enriquecimento de estruturas Non-B em âncoras Hi-C",
       subtitle = "Fold enrichment > 1 indica enriquecimento; < 1 indica depleção",
       x = "Classe Non-B",
       y = "Fold Enrichment (observado/esperado)") +
  theme_artigo +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  expand_limits(y = max(df_enrich$Enriquecimento) * 1.2)

ggsave(file.path(dir_plots, "enriquecimento_nonb_no_hic.png"), 
       p_enrich_bars, width = 10, height = 7, dpi = 300)
# Resumo
tabela_completa <- df_fracao %>%
  mutate(
    Cobertura_Genoma = round(Cobertura_Genoma, 3),
    Sobreposicao_HiC_bp = Sobreposicao_HiC_bp,
    Sobreposicao_Esperada_bp = round(Tamanho_bp * cobertura_hic, 0),
    Fracao_no_HiC = round(Fracao_no_HiC, 2),
    Fold_Enrichment = round(Enriquecimento, 2),
    Significancia = ifelse(!is.na(Enriquecimento) & Enriquecimento > 1.5, 
                           "Enriquecido", 
                           ifelse(!is.na(Enriquecimento) & Enriquecimento < 0.67, 
                                  "Depletado", "Neutro"))
  ) %>%
  select(Classe, Tamanho_bp, Cobertura_Genoma, Sobreposicao_HiC_bp, 
         Sobreposicao_Esperada_bp, Fracao_no_HiC, Fold_Enrichment, Significancia)

write.csv(tabela_completa, file.path(dir_tabelas, "tabela_enriquecimento_nonb_hic.csv"), row.names = FALSE)
print(tabela_completa)

# ============================================================
# BLOCO 2: ANÁLISES COMBINADAS (VENN E DISTÂNCIAS)

# 2.1) VENN DIAGRAM (SMAP, Hi-C, CRUCIFORM)
if(length(hic_anchors_gr) > 0 && "Cruciform" %in% names(nonb_por_classe)) {  
  cruciform_gr <- nonb_por_classe[["Cruciform"]]
  n_smap <- length(smap_gr)
  n_hic <- length(hic_anchors_gr)
  n_cruc <- length(cruciform_gr)
  n_smap_hic <- length(unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr))))
  n_smap_cruc <- length(unique(queryHits(findOverlaps(smap_gr, cruciform_gr))))
  n_hic_cruc <- length(unique(queryHits(findOverlaps(hic_anchors_gr, cruciform_gr))))
  
  smap_hic_hits <- unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr)))
  n_tripla <- 0
  for(i in smap_hic_hits) {
    if(length(subjectHits(findOverlaps(smap_gr[i], cruciform_gr))) > 0) n_tripla <- n_tripla + 1}
  
  df_venn <- data.frame(
    Conjunto = c("SMAP", "Hi-C", "Cruciform", "SMAP∩Hi-C", "SMAP∩Cruciform", "Hi-C∩Cruciform", "SMAP∩Hi-C∩Cruciform"),
    Contagem = c(n_smap, n_hic, n_cruc, n_smap_hic, n_smap_cruc, n_hic_cruc, n_tripla))
  write.csv(df_venn, file.path(dir_tabelas, "Venn_SMAP_HiC_Cruciform_contagens.csv"), row.names = FALSE)
  png(file.path(dir_plots, "Venn_SMAP_HiC_Cruciform.png"), width = 900, height = 900, res = 150)
  venn.plot <- draw.triple.venn(
    area1 = n_smap, area2 = n_hic, area3 = n_cruc,
    n12 = n_smap_hic, n23 = n_hic_cruc, n13 = n_smap_cruc, n123 = n_tripla,
    category = c("SMAP", "Hi-C", "Cruciform"), 
    fill = c(cores_padrao["Z-DNA"], cores_padrao["G-quadruplex"], cores_padrao["Triplex"]),
    alpha = 0.5, 
    cat.cex = 1.8, 
    cex = 1.8,
    cat.col = c(cores_padrao["Z-DNA"], cores_padrao["G-quadruplex"], cores_padrao["Triplex"]),
    cat.dist = c(0.15, 0.15, 0.1), 
    cat.pos = c(-10, 10, 180),
    main = "SMAP ∩ Hi-C ∩ Cruciform",
    main.cex = 1.8,
    main.fontface = "bold",
    margin = 0.1)
  grid.draw(venn.plot)
  dev.off()
  cat("  Venn Diagram salvo\n")}

# 2.2) DISTÂNCIAS MÉDIAS
if(length(hic_anchors_gr) > 0 && "Cruciform" %in% names(nonb_por_classe)) {
  cruciform_gr <- nonb_por_classe[["Cruciform"]]
  dist_smap_hic <- distanceToNearest(smap_gr, hic_anchors_gr)
  dist_smap_cruc <- distanceToNearest(smap_gr, cruciform_gr)
  dist_hic_cruc <- distanceToNearest(hic_anchors_gr, cruciform_gr)
  
  df_smap_hic_dist <- data.frame(
    SMAP_seq = as.character(seqnames(smap_gr[queryHits(dist_smap_hic)])),
    SMAP_start = start(smap_gr[queryHits(dist_smap_hic)]),
    SMAP_end = end(smap_gr[queryHits(dist_smap_hic)]),
    Anchor_seq = as.character(seqnames(hic_anchors_gr[subjectHits(dist_smap_hic)])),
    Anchor_start = start(hic_anchors_gr[subjectHits(dist_smap_hic)]),
    Anchor_end = end(hic_anchors_gr[subjectHits(dist_smap_hic)]),
    Distancia_kb = mcols(dist_smap_hic)$distance / 1000,
    stringsAsFactors = FALSE)
  write.csv(df_smap_hic_dist, file.path(dir_tabelas, "distancias_SMAP_HiC_coordenadas.csv"), row.names = FALSE)
  df_smap_cruc_dist <- data.frame(
    SMAP_seq = as.character(seqnames(smap_gr[queryHits(dist_smap_cruc)])),
    SMAP_start = start(smap_gr[queryHits(dist_smap_cruc)]),
    SMAP_end = end(smap_gr[queryHits(dist_smap_cruc)]),
    Cruciform_seq = as.character(seqnames(cruciform_gr[subjectHits(dist_smap_cruc)])),
    Cruciform_start = start(cruciform_gr[subjectHits(dist_smap_cruc)]),
    Cruciform_end = end(cruciform_gr[subjectHits(dist_smap_cruc)]),
    Distancia_kb = mcols(dist_smap_cruc)$distance / 1000,
    stringsAsFactors = FALSE)
  write.csv(df_smap_cruc_dist, file.path(dir_tabelas, "distancias_SMAP_Cruciform_coordenadas.csv"), row.names = FALSE)
  df_hic_cruc_dist <- data.frame(
    Anchor_seq = as.character(seqnames(hic_anchors_gr[queryHits(dist_hic_cruc)])),
    Anchor_start = start(hic_anchors_gr[queryHits(dist_hic_cruc)]),
    Anchor_end = end(hic_anchors_gr[queryHits(dist_hic_cruc)]),
    Cruciform_seq = as.character(seqnames(cruciform_gr[subjectHits(dist_hic_cruc)])),
    Cruciform_start = start(cruciform_gr[subjectHits(dist_hic_cruc)]),
    Cruciform_end = end(cruciform_gr[subjectHits(dist_hic_cruc)]),
    Distancia_kb = mcols(dist_hic_cruc)$distance / 1000,
    stringsAsFactors = FALSE)
  write.csv(df_hic_cruc_dist, file.path(dir_tabelas, "distancias_HiC_Cruciform_coordenadas.csv"), row.names = FALSE)
  
  media_smap_hic <- mean(mcols(dist_smap_hic)$distance / 1000, na.rm = TRUE)
  media_smap_cruc <- mean(mcols(dist_smap_cruc)$distance / 1000, na.rm = TRUE)
  media_hic_cruc <- mean(mcols(dist_hic_cruc)$distance / 1000, na.rm = TRUE)
  df_distancias <- data.frame(
    Par = c("SMAP vs Hi-C", "SMAP vs Cruciform", "Hi-C vs Cruciform"),
    Media_kb = c(media_smap_hic, media_smap_cruc, media_hic_cruc))
  write.csv(df_distancias, file.path(dir_tabelas, "distancias_medias_resumo.csv"), row.names = FALSE)
  
  p_distancias <- ggplot(df_distancias, aes(x = Par, y = Media_kb, fill = Par)) +
    geom_col(width = 0.6, alpha = 0.8, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f kb", Media_kb)), vjust = -0.5, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c(cores_padrao["Z-DNA"], cores_padrao["G-quadruplex"], cores_padrao["Triplex"])) +
    labs(title = "Distâncias médias entre conjuntos",
         subtitle = "Distância até o elemento mais próximo da outra categoria", x = "", y = "Distância média (kb)") +
    theme_artigo + 
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.x = element_text(size = 11, angle = 15, hjust = 1, face = "bold"))
  ggsave(file.path(dir_plots, "distancias_medias.png"), p_distancias, width = 9, height = 7, dpi = 300)
  cat("  Gráfico de distâncias médias salvo\n")}

# ============================================================
# SALVAR OBJETOS E FINALIZAR
saveRDS(smap_gr, file.path(dir_resultados, "smap_gr.rds"))
saveRDS(vng_gr, file.path(dir_resultados, "vng_gr.rds"))
saveRDS(genes_gr, file.path(dir_resultados, "genes_gr.rds"))
saveRDS(hic_anchors_gr, file.path(dir_resultados, "hic_anchors.rds"))

stopCluster(cl)