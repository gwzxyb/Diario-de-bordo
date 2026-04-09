# --------------PACOTES-----------------------
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
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Paralelização configurada com", num_cores, "cores\n")

# ------------- CORES PADRÃO -------------
cores_padrao <- c(
  "Z-DNA" = "#F28AA0",        # rosa suave
  "G-quadruplex" = "#366899",  # azul
  "Triplex" = "#81DACA",       # verde água
  "R-loop" = "#FC8D62",        # laranja
  "Short_tandem" = "#66C2A5",  # verde
  "A-phased" = "#8DA0CB",      # roxo
  "Cruciform" = "#377EB8",     # azul escuro
  "Outros" = "#BEBADA",        # cinza
  "Com Hi-C" = "#66C2A5",      # verde para interseção
  "Sem Hi-C" = "#FC8D62",      # laranja para sem interseção
  "SMAP em genes" = "#81DACA", # verde água
  "SMAP fora de genes" = "#F28AA0" # rosa
)

# Cores específicas para Odds Ratio (inseridas manualmente)
cor_significativo <- "#E41A1C"  # vermelho para p < 0.05
cor_nao_significativo <- "#377EB8" # azul para p ≥ 0.05

# Tema consistente
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

# 1) DIRETÓRIOS
data_dir <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/Git limpo/sobreposicao/"
dir_plots <- file.path(dir_resultados, "plots")
dir_tabelas <- file.path(dir_resultados, "tabelas")  
dir_intersecoes <- file.path(dir_resultados, "intersecoes_nonb")

dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tabelas, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_intersecoes, showWarnings = FALSE, recursive = TRUE)

# 2) GENOMA DE REFERÊNCIA
fa_path <- "/home/mcavalcante/dados_brutos/Hsalinarum.fa"
if (file.exists(fa_path)) {
  fa <- readDNAStringSet(fa_path)
  genome_size <- sum(width(fa))
  cat("Genoma carregado:", length(fa), "sequências | Tamanho:", format(genome_size, big.mark=","), "bp\n")
} else {
  warning("Arquivo FASTA não encontrado.")
  fa <- NULL
  genome_size <- 2000000
}

# 3) SELEÇÃO DE ARQUIVOS
arquivos <- list.files(data_dir, full.names = TRUE)

arq_smap <- arquivos[str_detect(arquivos, regex("smap1", ignore_case = TRUE))]
arq_gff <- arquivos[str_detect(arquivos, regex("Hsalinarum.*gff", ignore_case = TRUE))]
arq_hic_bedpe <- arquivos[str_detect(arquivos, regex("\\.bedpe$", ignore_case = TRUE))]

arquivos_nonb <- arquivos[
  str_detect(arquivos, regex("z-dna|triplex|r-loop|short_tandem|A-phased|Cruciform|G-quadruplex", ignore_case = TRUE)) &
    str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]

if(length(arq_smap) == 0) stop("arquivo SMAP nao encontrado")
if(length(arq_gff) == 0) stop("arquivo GFF nao encontrado")

# 4) CLASSIFICAÇÃO NON-B
classificar_nonb <- function(arquivo) {
  nome_lower <- tolower(basename(arquivo))
  if (str_detect(nome_lower, "z-dna|zdna")) return("Z-DNA")
  if (str_detect(nome_lower, "g-quadruplex|g4")) return("G-quadruplex")
  if (str_detect(nome_lower, "triplex")) return("Triplex")
  if (str_detect(nome_lower, "r-loop|rloop")) return("R-loop")
  if (str_detect(nome_lower, "short_tandem|tandem")) return("Short_tandem")
  if (str_detect(nome_lower, "a-phased|aphased")) return("A-phased")
  if (str_detect(nome_lower, "cruciform")) return("Cruciform")
  return("Outros")
}

nonb_classificado <- data.frame(
  arquivo = arquivos_nonb,
  nome_arquivo = basename(arquivos_nonb),
  classe = sapply(arquivos_nonb, classificar_nonb),
  stringsAsFactors = FALSE)

write.table(nonb_classificado, file.path(dir_tabelas, "nonb_classificacao.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nClasses Non-B:\n")
print(table(nonb_classificado$classe))

# 5) FUNÇÃO CARREGAR GRANGES
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
  }, error = function(e) return(GRanges()))
}

# 6) CARREGAR SMAP E GFF
cat("\nCarregando SMAP e GFF...\n")
smap_gr <- carregar_granges(arq_smap[1])
cat("SMAP:", length(smap_gr), "regiões\n")

gff_gr <- carregar_granges(arq_gff[1])
mcols_gff <- as.data.frame(mcols(gff_gr))

genes_gr <- gff_gr[mcols_gff$type == "gene"]
cat("Genes:", length(genes_gr), "\n")

vng_loci <- c("tnpb") 
vng_gr <- GRanges()
for(locus in vng_loci) {
  hits <- grepl(locus, mcols_gff$locus_tag, ignore.case=TRUE) |
    grepl(locus, mcols_gff$ID, ignore.case=TRUE) |
    grepl(locus, mcols_gff$product, ignore.case=TRUE)
  if(any(hits)) vng_gr <- c(vng_gr, gff_gr[hits])
}
vng_gr <- unique(vng_gr)
cat("TnpB:", length(vng_gr), "loci\n")

# 7) CARREGAR Hi-C
hic_anchors_gr <- GRanges()
hic_interacoes <- NULL
if(length(arq_hic_bedpe) > 0) {
  hic_bedpe <- read.table(arq_hic_bedpe[1], header=FALSE, stringsAsFactors=FALSE)
  colnames(hic_bedpe)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  
  anchors1 <- GRanges(hic_bedpe$chr1, IRanges(hic_bedpe$start1, hic_bedpe$end1))
  anchors2 <- GRanges(hic_bedpe$chr2, IRanges(hic_bedpe$start2, hic_bedpe$end2))
  hic_anchors_gr <- unique(c(anchors1, anchors2))
  
  hic_interacoes <- GInteractions(anchors1, anchors2)
  
  cat("Âncoras Hi-C únicas:", length(hic_anchors_gr), "\n")
  cat("Total interações Hi-C:", nrow(hic_bedpe), "\n")
}

# 8) CARREGAR NON-B
cat("\nCarregando Non-B...\n")
gr_list_nonb <- list()
for (i in seq_len(nrow(nonb_classificado))) {
  arquivo <- nonb_classificado$arquivo[i]
  nome_arq <- nonb_classificado$nome_arquivo[i]
  classe <- nonb_classificado$classe[i]
  gr <- carregar_granges(arquivo)
  
  if (length(gr) > 0) {
    mcols(gr)$origem <- nome_arq
    mcols(gr)$classe <- classe
    nome_lista <- gsub("\\.(bed|gff|gff3)$", "", nome_arq)
    gr_list_nonb[[nome_lista]] <- gr
    cat(basename(arquivo), "(", length(gr), "regiões,", classe, ")\n")
  }
}
gr_list_nonb <- gr_list_nonb[sapply(gr_list_nonb, length) > 0]

# 9) CONSOLIDAR NON-B POR CLASSE (INTERSEÇÃO/CONSENSO)
gr_por_classe <- list()
classes_com_dados <- unique(sapply(gr_list_nonb, function(x) unique(mcols(x)$classe)[1]))

for (classe in classes_com_dados) {
  # Filtra a lista para pegar apenas os GRanges desta classe
  grs_classe <- gr_list_nonb[sapply(gr_list_nonb, function(x) unique(mcols(x)$classe)[1] == classe)]
  
  if (length(grs_classe) > 1) {
    # Se houver mais de um arquivo (ex: Triplex), faz a INTERSEÇÃO
    # Isso resulta apenas nas regiões onde todos os arquivos concordam
    gr_consenso <- Reduce(intersect, grs_classe)
    mcols(gr_consenso)$classe <- classe
    gr_por_classe[[classe]] <- gr_consenso
    cat(classe, ": Gerado consenso de", length(gr_consenso), "regiões (Interseção de", length(grs_classe), "arquivos)\n")
  } else if (length(grs_classe) == 1) {
    # Se houver apenas um arquivo, usa ele integralmente
    gr_por_classe[[classe]] <- grs_classe[[1]]
    cat(classe, ":", length(gr_por_classe[[classe]]), "regiões (Arquivo único)\n")
  }
}

# 10) DISTÂNCIA SMAP vs TnpB
if(length(vng_gr) > 0) {
  dist_smap_tnpb <- distanceToNearest(smap_gr, vng_gr)
  dist_vec_tnpb <- mcols(dist_smap_tnpb)$distance / 1000
  
  # TABELA DE COORDENADAS
  df_dist_tnpb <- data.frame(
    SMAP_seq = as.character(seqnames(smap_gr[queryHits(dist_smap_tnpb)])),
    SMAP_start = start(smap_gr[queryHits(dist_smap_tnpb)]),
    SMAP_end = end(smap_gr[queryHits(dist_smap_tnpb)]),
    TnpB_seq = as.character(seqnames(vng_gr[subjectHits(dist_smap_tnpb)])),
    TnpB_start = start(vng_gr[subjectHits(dist_smap_tnpb)]),
    TnpB_end = end(vng_gr[subjectHits(dist_smap_tnpb)]),
    Distancia_kb = dist_vec_tnpb,
    stringsAsFactors = FALSE)
  write.csv(df_dist_tnpb, file.path(dir_tabelas, "distancias_SMAP_TnpB_coordenadas.csv"), row.names = FALSE)
  
  media_tnpb <- mean(dist_vec_tnpb, na.rm = TRUE)
  mediana_tnpb <- median(dist_vec_tnpb, na.rm = TRUE)
  max_y <- max(hist(dist_vec_tnpb, plot=FALSE, breaks=30)$counts)
  
  p_dist_tnpb <- ggplot(data.frame(Distancia_kb = dist_vec_tnpb), aes(x = Distancia_kb)) +
    geom_histogram(bins = 30, fill = cores_padrao["A-phased"], alpha = 0.7, color = "white", linewidth = 0.5) +
    geom_vline(xintercept = media_tnpb, color = cores_padrao["R-loop"], linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = mediana_tnpb, color = cores_padrao["G-quadruplex"], linetype = "dotted", linewidth = 1.2) +
    annotate("text", x = media_tnpb + 2, y = max_y * 0.9, 
             label = paste("Média =", round(media_tnpb, 2), "kb"), 
             color = cores_padrao["R-loop"], size = 4, hjust = 0, fontface = "bold") +
    annotate("text", x = mediana_tnpb + 2, y = max_y * 0.8, 
             label = paste("Mediana =", round(mediana_tnpb, 2), "kb"), 
             color = cores_padrao["G-quadruplex"], size = 4, hjust = 0, fontface = "bold") +
    labs(title = "Distância entre regiões SMAP e loci TnpB",
         subtitle = paste("Total de pares SMAP-TnpB:", length(dist_vec_tnpb)),
         x = "Distância (kb)", 
         y = "Frequência (número de pares)") +
    theme_artigo
  ggsave(file.path(dir_plots, "distancias_SMAP_TnpB.png"), p_dist_tnpb, width = 10, height = 6, dpi = 300)
}

# 11) SMAP vs GENES (COM TABELA DE COORDENADAS)
if(length(genes_gr) > 0) {
  smap_genes_int <- findOverlaps(smap_gr, genes_gr)
  smap_genes_hits <- unique(queryHits(smap_genes_int))
  n_smap_genes <- length(smap_genes_hits)
  perc_smap_genes <- round(100 * n_smap_genes / length(smap_gr), 1)
  
  # TABELA DE COORDENADAS
  if(length(smap_genes_int) > 0) {
    df_smap_genes_coord <- data.frame(
      SMAP_seq = as.character(seqnames(smap_gr[queryHits(smap_genes_int)])),
      SMAP_start = start(smap_gr[queryHits(smap_genes_int)]),
      SMAP_end = end(smap_gr[queryHits(smap_genes_int)]),
      Gene_seq = as.character(seqnames(genes_gr[subjectHits(smap_genes_int)])),
      Gene_start = start(genes_gr[subjectHits(smap_genes_int)]),
      Gene_end = end(genes_gr[subjectHits(smap_genes_int)]),
      Gene_id = mcols(genes_gr[subjectHits(smap_genes_int)])$ID,
      Gene_name = mcols(genes_gr[subjectHits(smap_genes_int)])$Name,
      stringsAsFactors = FALSE)
    write.csv(df_smap_genes_coord, file.path(dir_tabelas, "SMAP_Genes_interseccao_coordenadas.csv"), row.names = FALSE)
  }
  
  df_genes_status <- data.frame(
    Status = c("SMAP em genes", "SMAP fora de genes"),
    Contagem = c(n_smap_genes, length(smap_gr) - n_smap_genes),
    Percentual = c(perc_smap_genes, 100 - perc_smap_genes))
  
  write.csv(df_genes_status, file.path(dir_tabelas, "SMAP_Genes_resumo.csv"), row.names = FALSE)
  
  p_genes <- ggplot(df_genes_status, aes(x = 2, y = Contagem, fill = Status)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) + 
    coord_polar(theta = "y") + 
    xlim(0.5, 2.5) +
    geom_text(aes(label = paste0(Contagem, " (", Percentual, "%)")), 
              position = position_stack(vjust = 0.5), size = 4.5, fontface = "bold", color = "white") +
    scale_fill_manual(values = c(cores_padrao["SMAP em genes"], cores_padrao["SMAP fora de genes"])) +
    labs(title = "SMAP vs Genes",
         subtitle = paste("Total SMAP:", length(smap_gr), "| Total genes:", length(genes_gr)),
         fill = "Categoria") + 
    theme_void() +
    theme(legend.position = "right", 
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#5D6D7E"))
  
  ggsave(file.path(dir_plots, "SMAP_vs_Genes_pizza.png"), p_genes, width = 10, height = 7, dpi = 300)
}


# 12
# 1. Definir o total de âncoras genômicas únicas (N = 20.123)
n_hic_total <- length(hic_anchors_gr)

# 2. Criar data.frame para consolidar todos os resultados
tabela_consolidada_nonB <- data.frame()

cat("\n--- Iniciando Processamento de Interseções (Consenso) ---\n")

for (classe in names(gr_consenso_list)) {
  # Carregar o consensus da classe (DNA Non-B)
  gr_NB <- gr_consenso_list[[classe]]
  
  # Encontrar todas as sobreposições entre Âncoras Hi-C e o DNA Non-B
  # findOverlaps retorna um objeto de hits com pares (ancora, non-b)
  hits <- findOverlaps(hic_anchors_gr, gr_NB)
  
  # --- MÉTRICA 1: ÂNCORAS ÚNICAS (PRESENÇA) ---
  # Quantas âncoras diferentes tocam pelo menos um Non-B desta classe
  n_ancoras_unicas <- length(unique(queryHits(hits)))
  percentual_presenca <- (n_ancoras_unicas / n_hic_total) * 100
  
  # --- MÉTRICA 2: PINGS TOTAIS (VOLUME/INTERSEÇÕES) ---
  # Quantas vezes as estruturas "pingam" nas âncoras no total
  # (Se uma âncora toca 3 estruturas, conta como 3 pings)
  total_pings <- length(subjectHits(hits))
  
  # --- MÉTRICA 3: DENSIDADE (MULTIPLICIDADE) ---
  # Média de estruturas encontradas por ponto de contato atingido
  media_pings <- ifelse(n_ancoras_unicas > 0, round(total_pings / n_ancoras_unicas, 2), 0)
  
  # Adicionar ao data.frame final
  tabela_consolidada_nonB <- rbind(tabela_consolidada_nonB, data.frame(
    Classe = classe,
    Ancoras_Unicas = n_ancoras_unicas,
    Percentual_Presenca = round(percentual_presenca, 2),
    Total_Pings_Interseccao = total_pings,
    Media_Pings_por_Ancora = media_pings,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("Processado: %-15s | Únicas: %5d | Pings: %6d | Média: %.2f\n", 
              classe, n_ancoras_unicas, total_pings, media_pings))
}

# 3. Ordenar a tabela pela presença (Âncoras Únicas) para facilitar a leitura
tabela_consolidada_nonB <- tabela_consolidada_nonB[order(-tabela_consolidada_nonB$Ancoras_Unicas), ]

# 4. Salvar e Exibir Resultados
write.csv(tabela_consolidada_nonB, 
          file.path(dir_tabelas, "analise_final_unicas_vs_pings_consenso.csv"), 
          row.names = FALSE)

cat("\n--- RESULTADOS FINAIS CONSOLIDADOS ---\n")
print(tabela_consolidada_nonB)
# ==============================================================================
# GRÁFICO DE PINGS TOTAIS (INTERSEÇÕES TOTAIS)
# ==============================================================================

# Garantir que a tabela esteja ordenada pelo volume de pings
tabela_plot_pings <- tabela_consolidada_nonB %>%
  arrange(Total_Pings_Interseccao)

# Criar o gráfico
p_pings_totais <- ggplot(tabela_plot_pings, 
                         aes(x = factor(Classe, levels = Classe), 
                             y = Total_Pings_Interseccao, 
                             fill = Classe)) +
  geom_col(alpha = 0.9, width = 0.7, color = "black", linewidth = 0.3) +
  # Adicionar o número total e a média de pings por âncora no rótulo
  geom_text(aes(label = paste0(Total_Pings_Interseccao, " (avg: ", Media_Pings_por_Ancora, ")")), 
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = cores_padrao) +
  labs(title = "Hi-C vs Non-B (Total de Pings)",
       subtitle = "Volume total de interseções (incluindo múltiplas ocorrências por âncora)",
       x = "Classe Non-B", 
       y = "Número Total de Interseções (Pings)") +
  theme_artigo +
  theme(legend.position = "none") +
  # Expandir o limite para o texto não ser cortado
  expand_limits(y = max(tabela_plot_pings$Total_Pings_Interseccao) * 1.3)

# Salvar o gráfico
ggsave(file.path(dir_plots, "barras_hic_pings_totais.png"), 
       p_pings_totais, width = 10, height = 7, dpi = 300)

cat("\nGráfico de Pings Totais gerado com sucesso!\n")

# 13) SMAP vs Non-B (COM TABELA DE COORDENADAS)
if(length(gr_por_classe) > 0) {
  classes <- names(gr_por_classe)
  todas_intersecoes_smap_nonb <- data.frame()
  
  for(classe in classes) {
    inter <- findOverlaps(smap_gr, gr_por_classe[[classe]], minoverlap=1)
    
    if(length(inter) > 0) {
      df_temp <- data.frame(
        classe_nonb = classe,
        SMAP_seq = as.character(seqnames(smap_gr[queryHits(inter)])),
        SMAP_start = start(smap_gr[queryHits(inter)]),
        SMAP_end = end(smap_gr[queryHits(inter)]),
        NonB_seq = as.character(seqnames(gr_por_classe[[classe]][subjectHits(inter)])),
        NonB_start = start(gr_por_classe[[classe]][subjectHits(inter)]),
        NonB_end = end(gr_por_classe[[classe]][subjectHits(inter)]),
        sobreposicao_tamanho = width(pintersect(smap_gr[queryHits(inter)], gr_por_classe[[classe]][subjectHits(inter)])),
        stringsAsFactors = FALSE)
      todas_intersecoes_smap_nonb <- rbind(todas_intersecoes_smap_nonb, df_temp)
    }
  }
  
  if(nrow(todas_intersecoes_smap_nonb) > 0) {
    write.csv(todas_intersecoes_smap_nonb, file.path(dir_tabelas, "SMAP_vs_NonB_interseccoes_coordenadas.csv"), row.names = FALSE)
  }
  
  resultados_overlap <- foreach(classe = classes, .packages = c("GenomicRanges")) %dopar% {
    inter <- findOverlaps(smap_gr, gr_por_classe[[classe]], minoverlap=1)
    n_overlap <- length(unique(queryHits(inter)))
    list(classe = classe, n_overlap = n_overlap)
  }
  
  df_nonb_overlap <- data.frame(
    Classe = classes,
    N_SMAP_sobrepostas = sapply(resultados_overlap, function(x) x$n_overlap))
  df_nonb_overlap$Percentual <- round(100 * df_nonb_overlap$N_SMAP_sobrepostas / length(smap_gr), 1)
  df_nonb_overlap$Label <- paste0(df_nonb_overlap$N_SMAP_sobrepostas, " (", df_nonb_overlap$Percentual, "%)")
  df_nonb_overlap <- df_nonb_overlap[order(-df_nonb_overlap$N_SMAP_sobrepostas), ]
  
  write.csv(df_nonb_overlap, file.path(dir_tabelas, "SMAP_vs_NonB_resumo.csv"), row.names = FALSE)
  
  p1 <- ggplot(df_nonb_overlap, aes(x = reorder(Classe, N_SMAP_sobrepostas), y = N_SMAP_sobrepostas, fill = Classe)) +
    geom_col(alpha = 0.9, width = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = Label), hjust = -0.1, size = 4, fontface = "bold") +
    coord_flip() + 
    scale_fill_manual(values = cores_padrao[intersect(names(cores_padrao), classes)]) +
    labs(title = "SMAP vs Non-B",
         subtitle = paste("Total SMAP:", length(smap_gr), "regiões"),
         x = "Classe Non-B", 
         y = "Número de SMAP em sobreposição") +
    theme_artigo + 
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.y = element_text(size = 11, face = "bold")) +
    expand_limits(y = max(df_nonb_overlap$N_SMAP_sobrepostas) * 1.2)
  ggsave(file.path(dir_plots, "SMAP_vs_NonB_barras.png"), p1, width = 12, height = 8, dpi = 300)
}

# 14) SMAP vs Hi-C (PIZZA PADRONIZADO - COM TABELA DE COORDENADAS)
if(length(hic_anchors_gr) > 0) {
  smap_hic_int <- findOverlaps(smap_gr, hic_anchors_gr)
  smap_com_hic <- unique(queryHits(smap_hic_int))
  
  # TABELA DE COORDENADAS
  if(length(smap_hic_int) > 0) {
    df_smap_hic_coord <- data.frame(
      SMAP_seq = as.character(seqnames(smap_gr[queryHits(smap_hic_int)])),
      SMAP_start = start(smap_gr[queryHits(smap_hic_int)]),
      SMAP_end = end(smap_gr[queryHits(smap_hic_int)]),
      Anchor_seq = as.character(seqnames(hic_anchors_gr[subjectHits(smap_hic_int)])),
      Anchor_start = start(hic_anchors_gr[subjectHits(smap_hic_int)]),
      Anchor_end = end(hic_anchors_gr[subjectHits(smap_hic_int)]),
      stringsAsFactors = FALSE)
    write.csv(df_smap_hic_coord, file.path(dir_tabelas, "SMAP_HiC_interseccao_coordenadas.csv"), row.names = FALSE)
  }
  
  df_hic_status <- data.frame(
    Status = c("SMAP com interseção Hi-C", "SMAP sem interseção Hi-C"),
    Contagem = c(length(smap_com_hic), length(smap_gr) - length(smap_com_hic)),
    Percentual = c(round(100 * length(smap_com_hic) / length(smap_gr), 1),
                   round(100 * (length(smap_gr) - length(smap_com_hic)) / length(smap_gr), 1)))
  
  write.csv(df_hic_status, file.path(dir_tabelas, "SMAP_HiC_resumo.csv"), row.names = FALSE)
  
  p2 <- ggplot(df_hic_status, aes(x = 2, y = Contagem, fill = Status)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) + 
    coord_polar(theta = "y") + 
    xlim(0.5, 2.5) +
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
  
  ggsave(file.path(dir_plots, "SMAP_vs_HiC_pizza.png"), p2, width = 10, height = 7, dpi = 300)
}

# 15) VENN SMAP, Hi-C, CRUCIFORM (COM TABELA)
if("Cruciform" %in% names(gr_por_classe) && length(hic_anchors_gr) > 0) {
  cruciform_gr <- gr_por_classe[["Cruciform"]]
  
  n_smap <- length(smap_gr)
  n_hic <- length(hic_anchors_gr)
  n_cruc <- length(cruciform_gr)
  n_smap_hic <- length(unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr))))
  n_smap_cruc <- length(unique(queryHits(findOverlaps(smap_gr, cruciform_gr))))
  n_hic_cruc <- length(unique(queryHits(findOverlaps(hic_anchors_gr, cruciform_gr))))
  
  smap_hic_hits <- unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr)))
  n_tripla <- 0
  for(i in smap_hic_hits) {
    if(length(subjectHits(findOverlaps(smap_gr[i], cruciform_gr))) > 0) n_tripla <- n_tripla + 1
  }
  
  # TABELA DAS CONTAGENS
  df_venn <- data.frame(
    Conjunto = c("SMAP", "Hi-C", "Cruciform", "SMAP∩Hi-C", "SMAP∩Cruciform", "Hi-C∩Cruciform", "SMAP∩Hi-C∩Cruciform"),
    Contagem = c(n_smap, n_hic, n_cruc, n_smap_hic, n_smap_cruc, n_hic_cruc, n_tripla)
  )
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
}

# 16) DISTÂNCIAS MÉDIAS
if("Cruciform" %in% names(gr_por_classe) && length(hic_anchors_gr) > 0) {
  cruciform_gr <- gr_por_classe[["Cruciform"]]
  
  dist_smap_hic <- distanceToNearest(smap_gr, hic_anchors_gr)
  dist_smap_cruc <- distanceToNearest(smap_gr, cruciform_gr)
  dist_hic_cruc <- distanceToNearest(hic_anchors_gr, cruciform_gr)
  
  # TABELAS DE COORDENADAS
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
  
  p4 <- ggplot(df_distancias, aes(x = Par, y = Media_kb, fill = Par)) +
    geom_col(width = 0.6, alpha = 0.8, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f kb", Media_kb)), vjust = -0.5, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c(cores_padrao["Z-DNA"], cores_padrao["G-quadruplex"], cores_padrao["Triplex"])) +
    labs(title = "Distâncias médias entre conjuntos",
         subtitle = "Distância até o elemento mais próximo da outra categoria",
         x = "", 
         y = "Distância média (kb)") +
    theme_artigo + 
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.x = element_text(size = 11, angle = 15, hjust = 1, face = "bold"))
  ggsave(file.path(dir_plots, "distancias_medias.png"), p4, width = 9, height = 7, dpi = 300)
}

# 17) DISTRIBUIÇÃO DISTÂNCIAS SMAP-HiC
if(length(hic_anchors_gr) > 0) {
  dist_smap_hic <- distanceToNearest(smap_gr, hic_anchors_gr)
  dist_vec <- mcols(dist_smap_hic)$distance / 1000
  df_dist <- data.frame(Distancia_kb = dist_vec)
  
  media <- mean(dist_vec, na.rm = TRUE)
  mediana <- median(dist_vec, na.rm = TRUE)
  max_y <- max(hist(dist_vec, plot=FALSE, breaks=40)$counts)
  
  p5 <- ggplot(df_dist, aes(x = Distancia_kb)) +
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
         x = "Distância (kb)", 
         y = "Frequência (número de SMAP)") +
    theme_artigo + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ggsave(file.path(dir_plots, "distribuicao_distancias_SMAP_HiC.png"), p5, width = 12, height = 7, dpi = 300)
}

# 18) ODDS RATIO FOREST PLOT (COM CORES MANUAIS)
if(exists("fa") && !is.null(fa)) {
  genome_size <- sum(width(fa))
  cat("Tamanho do genoma (do FASTA):", format(genome_size, big.mark=","), "bp\n")
} else {
  genome_size <- 2000000
  cat("Tamanho do genoma (fallback):", format(genome_size, big.mark=","), "bp\n")
}

# Função para calcular odds ratio
calcular_or <- function(a, b, c, d, nome) {
  a <- max(0, a, na.rm = TRUE)
  b <- max(0, b, na.rm = TRUE)
  c <- max(0, c, na.rm = TRUE)
  d <- max(0, d, na.rm = TRUE)
  
  if(a == 0 || b == 0 || c == 0 || d == 0) {
    return(data.frame(
      Comparacao = nome,
      OR = NA,
      OR_lower = NA,
      OR_upper = NA,
      p_value = NA,
      a = a, b = b, c = c, d = d,
      stringsAsFactors = FALSE))
  }
  
  tryCatch({
    mat <- matrix(c(a, b, c, d), nrow = 2)
    test <- fisher.test(mat)
    return(data.frame(
      Comparacao = nome,
      OR = test$estimate,
      OR_lower = test$conf.int[1],
      OR_upper = test$conf.int[2],
      p_value = test$p.value,
      a = a, b = b, c = c, d = d,
      stringsAsFactors = FALSE))
  }, error = function(e) {
    return(data.frame(
      Comparacao = nome,
      OR = NA,
      OR_lower = NA,
      OR_upper = NA,
      p_value = NA,
      a = a, b = b, c = c, d = d,
      stringsAsFactors = FALSE))
  })
}

lista_or <- list()

if(length(hic_anchors_gr) > 0) {
  a1 <- length(unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr))))
  b1 <- length(smap_gr) - a1
  c1 <- length(hic_anchors_gr) - a1
  d1 <- max(0, genome_size - (a1 + b1 + c1))
  lista_or[[1]] <- calcular_or(a1, b1, c1, d1, "SMAP vs Hi-C")
  cat("SMAP vs Hi-C:", a1, "/", length(smap_gr), "\n")
}

if("Cruciform" %in% names(gr_por_classe)) {
  cruciform_gr <- gr_por_classe[["Cruciform"]]
  a2 <- length(unique(queryHits(findOverlaps(smap_gr, cruciform_gr))))
  b2 <- length(smap_gr) - a2
  c2 <- length(cruciform_gr) - a2
  d2 <- max(0, genome_size - (a2 + b2 + c2))
  lista_or[[2]] <- calcular_or(a2, b2, c2, d2, "SMAP vs Cruciform")
  cat("SMAP vs Cruciform:", a2, "/", length(smap_gr), "\n")
}

if(length(hic_anchors_gr) > 0 && "Cruciform" %in% names(gr_por_classe)) {
  a3 <- length(unique(queryHits(findOverlaps(hic_anchors_gr, cruciform_gr))))
  b3 <- length(hic_anchors_gr) - a3
  c3 <- length(cruciform_gr) - a3
  d3 <- max(0, genome_size - (a3 + b3 + c3))
  lista_or[[3]] <- calcular_or(a3, b3, c3, d3, "Hi-C vs Cruciform")
  cat("Hi-C vs Cruciform:", a3, "/", length(hic_anchors_gr), "\n")
}

if(length(vng_gr) > 0) {
  a4 <- length(unique(queryHits(findOverlaps(smap_gr, vng_gr))))
  b4 <- length(smap_gr) - a4
  c4 <- length(vng_gr) - a4
  d4 <- max(0, genome_size - (a4 + b4 + c4))
  lista_or[[4]] <- calcular_or(a4, b4, c4, d4, "SMAP vs TnpB")
  cat("SMAP vs TnpB:", a4, "/", length(smap_gr), "\n")
}

if(length(hic_anchors_gr) > 0 && length(vng_gr) > 0) {
  a5 <- length(unique(queryHits(findOverlaps(hic_anchors_gr, vng_gr))))
  b5 <- length(hic_anchors_gr) - a5
  c5 <- length(vng_gr) - a5
  d5 <- max(0, genome_size - (a5 + b5 + c5))
  lista_or[[5]] <- calcular_or(a5, b5, c5, d5, "Hi-C vs TnpB")
  cat("Hi-C vs TnpB:", a5, "/", length(hic_anchors_gr), "\n")
}

if("Cruciform" %in% names(gr_por_classe) && length(vng_gr) > 0) {
  cruciform_gr <- gr_por_classe[["Cruciform"]]
  a6 <- length(unique(queryHits(findOverlaps(cruciform_gr, vng_gr))))
  b6 <- length(cruciform_gr) - a6
  c6 <- length(vng_gr) - a6
  d6 <- max(0, genome_size - (a6 + b6 + c6))
  lista_or[[6]] <- calcular_or(a6, b6, c6, d6, "Cruciform vs TnpB")
  cat("Cruciform vs TnpB:", a6, "/", length(cruciform_gr), "\n")
}

if(length(hic_anchors_gr) > 0 && "Cruciform" %in% names(gr_por_classe)) {
  smap_em_hic <- unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr)))
  smap_em_cruc <- unique(queryHits(findOverlaps(smap_gr, cruciform_gr)))
  a7 <- length(intersect(smap_em_hic, smap_em_cruc))
  b7 <- length(smap_gr) - a7
  c7 <- length(unique(c(subjectHits(findOverlaps(smap_gr, hic_anchors_gr)), 
                        subjectHits(findOverlaps(smap_gr, cruciform_gr))))) - a7
  d7 <- max(0, genome_size - (a7 + b7 + c7))
  lista_or[[7]] <- calcular_or(a7, b7, c7, d7, "Tripla (SMAP∩Hi-C∩Cruciform)")
  cat("Tripla (SMAP∩Hi-C∩Cruciform):", a7, "/", length(smap_gr), "\n")
}

if(length(hic_anchors_gr) > 0 && length(vng_gr) > 0) {
  smap_em_hic <- unique(queryHits(findOverlaps(smap_gr, hic_anchors_gr)))
  smap_em_tnpb <- unique(queryHits(findOverlaps(smap_gr, vng_gr)))
  a8 <- length(intersect(smap_em_hic, smap_em_tnpb))
  b8 <- length(smap_gr) - a8
  c8 <- length(unique(c(subjectHits(findOverlaps(smap_gr, hic_anchors_gr)), 
                        subjectHits(findOverlaps(smap_gr, vng_gr))))) - a8
  d8 <- max(0, genome_size - (a8 + b8 + c8))
  lista_or[[8]] <- calcular_or(a8, b8, c8, d8, "Tripla (SMAP∩Hi-C∩TnpB)")
  cat("Tripla (SMAP∩Hi-C∩TnpB):", a8, "/", length(smap_gr), "\n")
}

if(length(hic_anchors_gr) > 0 && "Cruciform" %in% names(gr_por_classe) && length(vng_gr) > 0) {
  hic_em_cruc <- unique(queryHits(findOverlaps(hic_anchors_gr, cruciform_gr)))
  hic_em_tnpb <- unique(queryHits(findOverlaps(hic_anchors_gr, vng_gr)))
  a9 <- length(intersect(hic_em_cruc, hic_em_tnpb))
  b9 <- length(hic_anchors_gr) - a9
  c9 <- length(unique(c(subjectHits(findOverlaps(hic_anchors_gr, cruciform_gr)), 
                        subjectHits(findOverlaps(hic_anchors_gr, vng_gr))))) - a9
  d9 <- max(0, genome_size - (a9 + b9 + c9))
  lista_or[[9]] <- calcular_or(a9, b9, c9, d9, "Tripla (Hi-C∩Cruciform∩TnpB)")
  cat("Tripla (Hi-C∩Cruciform∩TnpB):", a9, "/", length(hic_anchors_gr), "\n")
}

df_or <- do.call(rbind, lista_or)
df_or <- df_or[!is.na(df_or$OR), ]

if(nrow(df_or) > 0) {
  df_or$Significativo <- ifelse(df_or$p_value < 0.05, "p < 0.05", "p ≥ 0.05")
  df_or$Contagem <- paste0(df_or$a, "/", df_or$a + df_or$b)
  df_or$Label <- paste0(round(df_or$OR, 2), " [", round(df_or$OR_lower, 2), "-", round(df_or$OR_upper, 2), "]")
  write.csv(df_or, file.path(dir_tabelas, "odds_ratio_resultados.csv"), row.names = FALSE)
  
  p6 <- ggplot(df_or, aes(x = OR, y = reorder(Comparacao, OR))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper, color = Significativo), 
                   height = 0.2, linewidth = 1.2) +
    geom_point(aes(color = Significativo), size = 4) +
    scale_color_manual(values = c("p < 0.05" = cor_significativo, "p ≥ 0.05" = cor_nao_significativo)) +
    scale_x_log10() +
    labs(title = "Odds Ratio de Enriquecimento",
         subtitle = paste("Tamanho do genoma:", format(genome_size, big.mark=","), "bp"),
         x = "Odds Ratio (escala log)", 
         y = "",
         color = "Significância",
         caption = "Linha vertical = OR = 1") +
    theme_artigo +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.y = element_text(size = 11, face = "bold"))
  ggsave(file.path(dir_plots, "odds_ratio_forest.png"), p6, width = 12, height = 9, dpi = 300)
  
  cat("\nResultados Odds Ratio:\n")
  print(df_or[, c("Comparacao", "OR", "p_value", "Significativo")])
} else {
  cat("\nNenhuma comparação com OR válido encontrada\n")
}

# 19) SALVAR OBJETOS
saveRDS(smap_gr, file.path(dir_resultados, "smap_gr.rds"))
saveRDS(vng_gr, file.path(dir_resultados, "vng_gr.rds"))
saveRDS(genes_gr, file.path(dir_resultados, "genes_gr.rds"))
saveRDS(hic_anchors_gr, file.path(dir_resultados, "hic_anchors.rds"))
saveRDS(gr_por_classe, file.path(dir_resultados, "nonb_por_classe.rds"))

# 20) FINALIZAR
stopCluster(cl)