# Pacotes necessários
library(GenomicRanges)
library(data.table)
library(stringr)
library(ggplot2)
library(viridis)
library(dplyr)
library(InteractionSet)
library(rtracklayer)
library(readr)

# 1) Diretórios e preparação
data_dir <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/hic_nonb_intersections"
dir_plots <- file.path(dir_resultados, "plots")
dir.create(dir_resultados, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# 2) Seleção de arquivos
arquivos <- list.files(data_dir, full.names = TRUE)
arquivos_hic <- arquivos[str_detect(arquivos, regex("bedpe$", ignore_case = TRUE))]
arquivos_nonb <- arquivos[
  str_detect(arquivos, regex("z-dna|triplex|r-loop|short_tandem_slipped|A-phased|Cruciform", ignore_case = TRUE)) &
    str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE)) &
    !str_detect(arquivos, regex("_gquad|_kim|gquad|g-quadruplex", ignore_case = TRUE))]

# 3) Classificação de Non-B por nome de arquivo
classificar_nonb <- function(arquivo) {
  nome_lower <- tolower(basename(arquivo))
  if (str_detect(nome_lower, "z-dna|zdna")) return("Z-DNA")
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
write.table(nonb_classificado, file.path(dir_resultados, "nonb_filtrado_classificacao.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("3) CLASSIFICAÇÃO NON-B:\n")
print(table(nonb_classificado$classe))

# 4) Funções utilitárias
carregar_granges <- function(arquivo) {
  if (!file.exists(arquivo)) return(GRanges())
  extensao <- tools::file_ext(arquivo)
  
  tryCatch({
    if (extensao %in% c("bed")) {
      gr <- import(arquivo, format = "BED")
    } else if (extensao %in% c("gff", "gff3")) {
      gr <- import(arquivo, format = "GFF")
    } else {
      return(GRanges())
    }
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr
  }, error = function(e) {
    warning(paste("Erro ao carregar", basename(arquivo), ":", e$message))
    GRanges()
  })
}

carregar_hic <- function(arquivo) {
  tryCatch({
    hic_df <- read_tsv(arquivo, col_names = FALSE, show_col_types = FALSE)
    gr1 <- GRanges(seqnames = hic_df$X1, ranges = IRanges(hic_df$X2 + 1, hic_df$X3))
    gr2 <- GRanges(seqnames = hic_df$X4, ranges = IRanges(hic_df$X5 + 1, hic_df$X6))
    seqlevels(gr1) <- gsub("^chr", "", seqlevels(gr1))
    seqlevels(gr2) <- gsub("^chr", "", seqlevels(gr2))
    GInteractions(gr1, gr2)
  }, error = function(e) {
    message("Erro ", basename(arquivo), ": ", e$message)
    NULL
  })
}

# 5) Carregamento de dados
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
    cat("5.1)", basename(arquivo), "(", length(gr), "regioes)\n")
  }
}
gr_list_nonb <- gr_list_nonb[sapply(gr_list_nonb, length) > 0]

gr_list_hic <- list()
for (hic_file in arquivos_hic) {
  nome_arq <- basename(hic_file)
  cat("5.2) Hi-C:", nome_arq, "...")
  hic <- carregar_hic(hic_file)
  
  if (!is.null(hic) && length(hic) > 0) {
    mcols(hic)$origem <- nome_arq
    nome_lista <- gsub("\\.bedpe$", "", nome_arq)
    gr_list_hic[[nome_lista]] <- hic
    cat("OK (", length(hic), "interacoes)\n")
  } else {
    cat("VAZIO\n")
  }
}

# 6) Contagens e tabelas básicas
contagem_nonb <- data.frame(
  Arquivo = names(gr_list_nonb),
  Classe = sapply(gr_list_nonb, function(x) unique(mcols(x)$classe)[1]),
  Num_Regioes = sapply(gr_list_nonb, length),
  stringsAsFactors = FALSE
)

contagem_hic <- data.frame(
  Arquivo = names(gr_list_hic),
  Num_Interacoes = sapply(gr_list_hic, length),
  stringsAsFactors = FALSE
)

write.table(contagem_nonb, file.path(dir_resultados, "contagem_nonb.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(contagem_hic, file.path(dir_resultados, "contagem_hic.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n6) RESUMO CARREGAMENTO:\n")
cat("Non-B:", sum(contagem_nonb$Num_Regioes), "regioes\n")
cat("Hi-C:", sum(contagem_hic$Num_Interacoes), "interacoes\n")

# 7) Consolidação por classe Non-B
gr_por_classe <- list()
classes_com_dados <- unique(contagem_nonb$Classe)
for (classe in classes_com_dados) {
  arquivos_classe <- contagem_nonb$Arquivo[contagem_nonb$Classe == classe]
  grs_classe <- gr_list_nonb[arquivos_classe]
  if (length(grs_classe) > 0) {
    gr_por_classe[[classe]] <- Reduce(c, grs_classe)
    cat("7)", classe, ":", length(gr_por_classe[[classe]]), "regioes\n")
  }
}

# 8) Interseções Hi-C anchors vs Non-B classes
matriz_hic_nonb <- data.frame()
for (hic_nome in names(gr_list_hic)) {
  hic <- gr_list_hic[[hic_nome]]
  anchors1 <- anchors(hic, type = "first")
  anchors2 <- anchors(hic, type = "second")
  
  for (classe in names(gr_por_classe)) {
    nonb <- gr_por_classe[[classe]]
    ov1 <- findOverlaps(anchors1, nonb, minoverlap = 1)
    n1 <- length(unique(queryHits(ov1)))
    ov2 <- findOverlaps(anchors2, nonb, minoverlap = 1)
    n2 <- length(unique(queryHits(ov2)))
    
    matriz_hic_nonb <- rbind(matriz_hic_nonb,
                             data.frame(hic_arquivo = hic_nome, hic_tipo = "anchor1", nonb_classe = classe,
                                        n_interseccoes = n1, total_hic = length(hic), total_nonb = length(nonb)),
                             data.frame(hic_arquivo = hic_nome, hic_tipo = "anchor2", nonb_classe = classe,
                                        n_interseccoes = n2, total_hic = length(hic), total_nonb = length(nonb))
    )
  }
}

matriz_hic_nonb <- matriz_hic_nonb %>%
  mutate(perc_hic = round(n_interseccoes / total_hic * 100, 2),
         perc_nonb = round(n_interseccoes / total_nonb * 100, 2))

write.table(matriz_hic_nonb, file.path(dir_resultados, "matriz_hic_nonb_interseccoes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("8) Matriz salva (", nrow(matriz_hic_nonb), "linhas)\n")

# 9) Plots principais
heatmap_hic <- ggplot(matriz_hic_nonb, aes(x = hic_tipo, y = nonb_classe, fill = n_interseccoes)) +
  facet_wrap(~hic_arquivo, scales = "free") +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(n_interseccoes > 0, n_interseccoes, "")),
            color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#E3F2FD", mid = "#FFF3E0", high = "#FFCDD2",
                       midpoint = median(matriz_hic_nonb$n_interseccoes), name = "Interseções") +
  labs(title = "Hi-C Âncoras vs Classes Non-B DNA", subtitle = "Halobacterium salinarum | Overlap ≥ 1bp",
       x = "Tipo Âncora", y = "Classe Non-B") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

ggsave(file.path(dir_plots, "9.1_heatmap_hic_nonb.png"), heatmap_hic, width = 14, height = 10, dpi = 300, bg = "white")
print(heatmap_hic)

# Resumo por classe
resumo_classes <- matriz_hic_nonb %>%
  group_by(nonb_classe) %>%
  summarise(total_interseccoes = sum(n_interseccoes), media_perc_hic = round(mean(perc_hic), 2),
            n_anchor1 = sum(n_interseccoes[hic_tipo == "anchor1"]),
            n_anchor2 = sum(n_interseccoes[hic_tipo == "anchor2"]),
            n_hic_arquivos = length(unique(hic_arquivo)), .groups = "drop") %>%
  arrange(desc(total_interseccoes)) %>%
  mutate(label = paste0(total_interseccoes, "(", media_perc_hic, "%)"))

cat("\n9.2) RESUMO POR CLASSE:\n")
print(resumo_classes)

p2 <- ggplot(resumo_classes, aes(x = reorder(nonb_classe, total_interseccoes), y = total_interseccoes, fill = nonb_classe)) +
  geom_bar(stat = "identity", alpha = 0.88, width = 0.75, color = "white", linewidth = 1.2) +
  geom_text(aes(label = label), hjust = -0.12, size = 4.5, fontface = "bold", color = "white") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Z-DNA" = "red", "Cruciform" = "#1F78B4", "R-loop" = "#FF7F00",
                               "Triplex" = "#33A02C", "A-phased" = "#6A3D9A", "Short_tandem" = "#FDBF6F"),
                    name = "Classe Non-B") +
  labs(title = "Total de Sobreposições Hi-C vs Classes Non-B DNA",
       subtitle = "Âncoras Hi-C (anchor1 + anchor2) | Halobacterium salinarum",
       x = "Classe Non-B DNA", y = "Âncoras Hi-C Sobrepostas") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        legend.position = "bottom", legend.direction = "horizontal")

ggsave(file.path(dir_plots, "9.2_resumo_classes_nonb.png"), p2, width = 14, height = 10, dpi = 400, bg = "white")
print(p2)

# 10) Salvar objetos principais
saveRDS(gr_list_hic, file.path(dir_resultados, "10_hic_data.rds"))
saveRDS(gr_list_nonb, file.path(dir_resultados, "10_nonb_data.rds"))
saveRDS(matriz_hic_nonb, file.path(dir_resultados, "10_interseccoes_data.rds"))

# 11) Relatório final
cat("\n11) ARQUIVOS GERADOS:\n")
cat("-", length(list.files(dir_resultados, pattern = "\\.tsv$")), "arquivos TSV\n")
cat("-", length(list.files(dir_plots, pattern = "\\.png$")), "plots PNG\n")
cat("-", length(names(gr_por_classe)), "classes Non-B\n")
cat("-", length(gr_list_hic), "datasets Hi-C\n")

cat("\nANALISE CONCLUIDA\n")
cat("Resultados:", dir_resultados, "\n")
cat("Plots:", dir_plots, "\n")
