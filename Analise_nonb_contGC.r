# --------------Pacotes-----------------------
library(GenomicRanges)
library(data.table)
library(stringr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(Biostrings)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(ggforce)
library(dplyr)
library(grid)
library(rtracklayer)
library(tidyverse)
library(InteractionSet)
library(parallel)
library(foreach)
library(doParallel)

# Configurar paralelização
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Paralelização configurada com", num_cores, "cores\n")

# ------------- legenda -----------
# nonb - non-B-form DNA

# 1) Diretórios e preparação
data_dir <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/Git limpo/sobreposicao/"
dir_plots <- file.path(dir_resultados, "plots")
dir_tabelas <- file.path(dir_resultados, "tabelas") 
dir_intersecoes <- file.path(dir_resultados, "intersecoes_nonb")

dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_resultados, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tabelas, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_intersecoes, showWarnings = FALSE, recursive = TRUE)

# 2) Genoma de referência (para GC)
fa_path <- "/home/mcavalcante/dados_brutos/Hsalinarum.fa"
if (file.exists(fa_path)) {
  fa <- readDNAStringSet(fa_path)
  cat("Genoma carregado:", length(fa), "sequências\n")
} else {
  warning("Arquivo não encontrado.")
  fa <- NULL}

# 3) Seleção de arquivos Non-B
arquivos <- list.files(data_dir, full.names = TRUE)
arquivos_nonb <- arquivos[
  str_detect(arquivos, regex("z-dna|nonb|triplex|r-loop|short_tandem|A-phased|G-quadruplex|Cruciform",
                             ignore_case = TRUE)) &
    str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]
cat("\nArquivos Non-B encontrados:", length(arquivos_nonb), "\n")

# 4) Classificação por classe de DNA alternativo
classificar_nonb <- function(arquivo) {
  nome <- basename(arquivo)
  nome_lower <- tolower(nome)
  if (str_detect(nome_lower, "z-dna|zdna")) return("Z-DNA")
  if (str_detect(nome_lower, "g-quadruplex|g4")) return("G-quadruplex")
  if (str_detect(nome_lower, "triplex")) return("Triplex")
  if (str_detect(nome_lower, "r-loop|rloop")) return("R-loop")
  if (str_detect(nome_lower, "short_tandem|tandem")) return("Short_tandem")
  if (str_detect(nome_lower, "a-phased|aphased")) return("A-phased")
  if (str_detect(nome_lower, "cruciform")) return("Cruciform")
  return("Outros")}

nonb_classificado <- data.frame(
  arquivo = arquivos_nonb,
  nome_arquivo = basename(arquivos_nonb),
  classe = sapply(arquivos_nonb, classificar_nonb),
  stringsAsFactors = FALSE)
# Salva a tabela de classificação
write.table(nonb_classificado,
            file.path(dir_tabelas, "nonb_classificacao.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResumo por classe:\n")
print(table(nonb_classificado$classe))

# 5) Função para carregar GRanges
carregar_granges <- function(arquivo) {
  if (!file.exists(arquivo)) {
    warning(paste("Arquivo não encontrado:", arquivo))
    return(GRanges())}
  extensao <- tools::file_ext(arquivo)
  tryCatch({
    if (extensao %in% c("bed", "bedgraph")) {
      gr <- import(arquivo, format = "BED")
    } else if (extensao %in% c("gff", "gff3")) {
      gr <- import(arquivo, format = "GFF")
    } else {
      warning(paste("Formato não suportado:", arquivo))
      return(GRanges()) }
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr <- gr[width(gr) > 0]
    return(gr)
  }, error = function(e) {
    warning(paste("Erro:", e$message))
    return(GRanges())})}

# 6) Carregamento de dados 
carregar_com_metadados <- function(i, df) {
  arquivo <- df$arquivo[i]
  nome_arq <- df$nome_arquivo[i]
  classe <- df$classe[i]
  gr <- carregar_granges(arquivo)
  if (length(gr) > 0) {
    mcols(gr)$origem <- nome_arq
    mcols(gr)$classe <- classe
    mcols(gr)$arquivo_id <- i
    nome_lista <- gsub("\\.(bed|gff|gff3)$", "", nome_arq)
    return(list(nome = nome_lista, gr = gr, status = "OK", n = length(gr)))
  } else {
    return(list(nome = NULL, gr = NULL, status = "VAZIO", n = 0))}}
resultados_carregamento <- foreach(i = 1:nrow(nonb_classificado), 
                                   .packages = c("GenomicRanges", "rtracklayer")) %dopar% {
                                     carregar_com_metadados(i, nonb_classificado)}

# Processar resultados
gr_list <- list()
for (res in resultados_carregamento) {
  if (!is.null(res$gr)) {
    gr_list[[res$nome]] <- res$gr
    cat("  Carregado:", res$nome, "-", res$n, "elementos\n")}}

cat("\nTotal de arquivos carregados:", length(gr_list), "\n")

# 7) Agrupar arquivos por classe
arquivos_por_classe <- list()
for (nome in names(gr_list)) {
  classe <- unique(mcols(gr_list[[nome]])$classe)[1]
  
  if (is.null(arquivos_por_classe[[classe]])) {
    arquivos_por_classe[[classe]] <- list()}
  arquivos_por_classe[[classe]][[nome]] <- gr_list[[nome]]}

cat("\nClasses encontradas:\n")
for (classe in names(arquivos_por_classe)) {
  n_arquivos <- length(arquivos_por_classe[[classe]])
  cat("  ", classe, ":", n_arquivos, "arquivo(s)\n")}

# 8) Calculando interseções por classe
calcular_interseccao_classe <- function(classe, arquivos_classe) {
  n_arquivos <- length(arquivos_classe)
  if (n_arquivos == 1) {
    gr_intersect <- arquivos_classe[[1]]
    return(list(
      classe = classe,
      gr_intersect = gr_intersect,
      n_intersect = length(gr_intersect),
      status = "único arquivo",
      arquivos = names(arquivos_classe)))}
  gr_intersect <- arquivos_classe[[1]]
  interseccao_valida <- TRUE
  for (i in 2:n_arquivos) {
    hits <- findOverlaps(gr_intersect, arquivos_classe[[i]])
    if (length(hits) == 0) {
      interseccao_valida <- FALSE
      break}
    gr_intersect <- pintersect(gr_intersect[queryHits(hits)], 
                               arquivos_classe[[i]][subjectHits(hits)])
    if (length(gr_intersect) == 0) {
      interseccao_valida <- FALSE
      break}}
  if (interseccao_valida && length(gr_intersect) > 0) {
    gr_intersect <- reduce(gr_intersect)
    return(list(
      classe = classe,
      gr_intersect = gr_intersect,
      n_intersect = length(gr_intersect),
      status = "interseção",
      arquivos = names(arquivos_classe)))
  } else {
    return(list(
      classe = classe,
      gr_intersect = GRanges(),
      n_intersect = 0,
      status = "vazio",
      arquivos = names(arquivos_classe)))}}
classes <- names(arquivos_por_classe)
resultados_interseccao <- foreach(classe = classes, 
                                  .packages = c("GenomicRanges")) %dopar% {
                                    calcular_interseccao_classe(classe, arquivos_por_classe[[classe]])}
nonb_por_classe <- list()
resumo_intersecoes <- data.frame()
for (res in resultados_interseccao) {
  classe <- res$classe
  cat("\n", classe, ":\n")
  cat("  Arquivos:", paste(res$arquivos, collapse = ", "), "\n")
  if (res$n_intersect > 0) {
    nonb_por_classe[[classe]] <- res$gr_intersect
    mcols(res$gr_intersect)$classe <- classe
    mcols(res$gr_intersect)$n_arquivos_origem <- length(res$arquivos)
    mcols(res$gr_intersect)$score <- 0
    mcols(res$gr_intersect)$name <- paste0(gsub("-", "_", classe), "_", 
                                           sprintf("%04d", 1:res$n_intersect))
    nome_bed <- file.path(dir_intersecoes, sprintf("nonB_%s_interseccao.bed", 
                                                   gsub("-", "_", classe)))
    export(res$gr_intersect, nome_bed, format = "BED")
    cat("  BED salvo:", basename(nome_bed), "\n")
    cat("  Regiões:", res$n_intersect, "\n")
    resumo_intersecoes <- rbind(resumo_intersecoes, data.frame(
      Classe = classe,
      Regioes_Interseccao = res$n_intersect,
      Arquivos_Origem = length(res$arquivos),
      Arquivos = paste(res$arquivos, collapse = "; "),
      stringsAsFactors = FALSE))
  } else {
    cat("Interseção vazia\n")}}
if (nrow(resumo_intersecoes) > 0) {
  write.table(resumo_intersecoes, 
              file.path(dir_intersecoes, "resumo_intersecoes_nonB.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  print(resumo_intersecoes)
  cat("\nTotal de classes com interseção:", nrow(resumo_intersecoes), "\n")
  cat("Total de regiões em interseção:", sum(resumo_intersecoes$Regioes_Interseccao), "\n")}

# 9) Contagem por arquivo
contagem_real <- data.frame(
  Arquivo = names(gr_list),
  Classe = sapply(gr_list, function(x) unique(mcols(x)$classe)[1]),
  Num_Predicoes = sapply(gr_list, length),
  stringsAsFactors = FALSE)
write.table(contagem_real,
            file.path(dir_tabelas, "contagem_real_predicoes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
# 10) Gráfico de barras das interseções
if (nrow(resumo_intersecoes) > 0) {
  p_intersec <- ggplot(resumo_intersecoes, 
                       aes(x = reorder(Classe, Regioes_Interseccao), 
                           y = Regioes_Interseccao, 
                           fill = Classe)) +
    geom_col(alpha = 0.8, width = 0.7, color = "white") +
    geom_text(aes(label = Regioes_Interseccao), 
              hjust = -0.1, size = 4) +
    coord_flip() +
    scale_fill_viridis_d(option = "turbo") +
    labs(title = "Regiões em Interseção por Classe Non-B",
         x = "Classe", y = "Número de regiões em interseção") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(dir_plots, "interseccoes_por_classe.png"),
         p_intersec, width = 10, height = 6, dpi = 300)}

# 11) Comparação: total vs interseção
if (exists("contagem_real") && nrow(resumo_intersecoes) > 0) {
  totais_classe <- contagem_real %>% group_by(Classe) %>% summarise(Total_Regioes = sum(Num_Predicoes), .groups = "drop")
  comparacao <- merge(totais_classe, resumo_intersecoes, by = "Classe", all = TRUE)
  comparacao$Regioes_Interseccao[is.na(comparacao$Regioes_Interseccao)] <- 0
  comparacao$Percentual_Interseccao <- round(100 * comparacao$Regioes_Interseccao / comparacao$Total_Regioes, 1)
  comparacao_long <- comparacao %>% pivot_longer(cols = c(Total_Regioes, Regioes_Interseccao), names_to = "Tipo", values_to = "Contagem") %>% mutate(Tipo = ifelse(Tipo == "Total_Regioes", "Total", "Interseção"))
  p_comp <- ggplot(comparacao_long, aes(x = reorder(Classe, -Contagem), y = Contagem, fill = Tipo)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7, color = "white") +
    geom_text(aes(label = Contagem), position = position_dodge(0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Total" = "#366899", "Interseção" = "#F28AA0")) +
    labs(title = "Comparação: Total vs Interseção por Classe", x = "Classe", y = "Número de regiões") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(dir_plots, "comparacao_total_vs_interseccao.png"), p_comp, width = 12, height = 7, dpi = 300)
  write.table(comparacao, file.path(dir_tabelas, "comparacao_total_vs_interseccao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)}

# 12) Consolidar por classe (para matriz)
gr_por_classe <- list()
classes_com_dados <- unique(contagem_real$Classe)
for (classe in classes_com_dados) {
  arquivos_classe <- contagem_real$Arquivo[contagem_real$Classe == classe]
  grs_classe <- gr_list[arquivos_classe]
  if (length(grs_classe) > 0) {
    gr_por_classe[[classe]] <- Reduce(c, grs_classe)
    cat("  Classe", classe, ":", length(gr_por_classe[[classe]]), "elementos\n")}}

# 13) Matriz de interseção completa
classes_todas <- unique(nonb_classificado$classe)
n_todas <- length(classes_todas)
matriz_completa <- matrix(0, nrow = n_todas, ncol = n_todas,
                          dimnames = list(classes_todas, classes_todas))
# Diagonal
for (classe in classes_todas) {
  if (classe %in% names(gr_por_classe)) {
    matriz_completa[classe, classe] <- length(gr_por_classe[[classe]])}}
# Fora da diagonal
for (i in 1:n_todas) {
  for (j in 1:n_todas) {
    if (i != j) {
      classe_i <- classes_todas[i]
      classe_j <- classes_todas[j]
      if (classe_i %in% names(gr_por_classe) && classe_j %in% names(gr_por_classe)) {
        if (length(gr_por_classe[[classe_i]]) > 0 && length(gr_por_classe[[classe_j]]) > 0) {
          inter <- findOverlaps(gr_por_classe[[classe_i]], gr_por_classe[[classe_j]], minoverlap = 1)
          matriz_completa[i, j] <- length(unique(queryHits(inter)))}}}}}
write.table(matriz_completa,
            file.path(dir_tabelas, "matriz_interseccao_completa.tsv"),
            sep = "\t", quote = FALSE)

# 14) Heatmap da matriz de interseção
matriz_melt <- melt(matriz_completa)
colnames(matriz_melt) <- c("Classe1", "Classe2", "Interseccao")
# Calcular porcentagem do menor conjunto
matriz_melt$Porcentagem <- NA
for (i in 1:nrow(matriz_melt)) {
  total1 <- matriz_completa[matriz_melt$Classe1[i], matriz_melt$Classe1[i]]
  total2 <- matriz_completa[matriz_melt$Classe2[i], matriz_melt$Classe2[i]]
  
  if (matriz_melt$Classe1[i] == matriz_melt$Classe2[i]) {
    matriz_melt$Porcentagem[i] <- 100
  } else if (total1 > 0 && total2 > 0) {
    min_total <- min(total1, total2)
    if (min_total > 0) {
      matriz_melt$Porcentagem[i] <- round((matriz_melt$Interseccao[i] / min_total) * 100, 1)}}}
# Heatmap
heatmap_classes <- ggplot(matriz_melt, aes(x = Classe1, y = Classe2, fill = Interseccao)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(Classe1 == Classe2,
                               paste0(Interseccao),
                               ifelse(Interseccao > 0,
                                      paste0(Interseccao, "\n(",
                                             ifelse(is.na(Porcentagem), "0", Porcentagem), "%)"),
                                      ""))),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC",
                       mid = "lightblue",
                       high = "red",
                       midpoint = ifelse(length(matriz_melt$Interseccao[matriz_melt$Interseccao > 0]) > 0,
                                         median(matriz_melt$Interseccao[matriz_melt$Interseccao > 0], na.rm = TRUE),
                                         0),
                       name = "Número de\nInterseções") +
  labs(title = "Matriz de Interseção entre Classes de Estruturas Non-B DNA",
       x = "Classe", y = "Classe") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank()) +
  coord_fixed()
ggsave(file.path(dir_plots, "heatmap_interseccao_classes.png"),
       heatmap_classes,
       width = max(12, n_todas * 0.8),
       height = max(10, n_todas * 0.8),
       dpi = 300, bg = "white")

# 15) Análise de conteúdo GC
calcular_gc_real <- function(gr, fa_reference) {
  if (is.null(fa_reference) || length(gr) == 0) return(NA)
  tryCatch({
    # Harmonizar seqlevels com nomes exatos do FASTA
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    valid_seqs <- intersect(seqlevels(gr), names(fa_reference))
    if(length(valid_seqs) == 0) {
      cat("Nenhuma sequência válida encontrada\n")
      return(NA) }
    
    gr <- gr[seqnames(gr) %in% valid_seqs]
    gr <- trim(gr)
    if(length(gr) == 0) return(NA)
    seqs <- getSeq(fa_reference, gr)
    gc_values <- sapply(as.character(seqs), function(s) {
      if (nchar(s) == 0) return(NA)
      s <- DNAString(s)
      af <- alphabetFrequency(s, baseOnly = TRUE)
      total <- sum(af[c("A", "C", "G", "T")])
      if (total > 0) return((af["C"] + af["G"]) / total)
      return(NA)})
    return(gc_values[!is.na(gc_values)])
  }, error = function(e) {
    cat("    ⚠ Erro ao calcular GC:", e$message, "\n")
    return(NA)})}

# Calcular GC para todas as regiões
gc_violin_data <- data.frame()
if (!is.null(fa)) {
  cat("  Processando", length(gr_list), "arquivos...\n")
  for (nome_arq in names(gr_list)) {
    gr <- gr_list[[nome_arq]]
    classe <- unique(mcols(gr)$classe)[1]
    cat("    ", sprintf("%-25s", nome_arq), "...")
    gc_vals <- calcular_gc_real(gr, fa)
    if (!all(is.na(gc_vals)) && length(gc_vals) > 0) {
      temp_df <- data.frame(
        Classe = classe,
        Arquivo = nome_arq,
        GC = gc_vals,
        stringsAsFactors = FALSE)
      gc_violin_data <- rbind(gc_violin_data, temp_df)
      cat(" ", length(gc_vals), "valores GC\n")
    } else {
      cat(" sem dados GC\n")}}
  gc_violin_data <- gc_violin_data[!is.na(gc_violin_data$GC), ]
  if (nrow(gc_violin_data) > 0) {
    # Estatísticas por classe
    stats_gc <- gc_violin_data %>%
      group_by(Classe) %>%
      summarise(
        Mediana = median(GC, na.rm = TRUE),
        Media = mean(GC, na.rm = TRUE),
        Q1 = quantile(GC, 0.25, na.rm = TRUE),
        Q3 = quantile(GC, 0.75, na.rm = TRUE),
        Min = min(GC, na.rm = TRUE),
        Max = max(GC, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>% arrange(Mediana)
    
    write.csv(stats_gc, file.path(dir_tabelas, "estatisticas_gc_por_classe.csv"), row.names = FALSE)
    gc_violin_data$Classe <- factor(gc_violin_data$Classe, levels = stats_gc$Classe)
    # Violino + boxplot
    plot_violino <- ggplot(gc_violin_data, aes(x = Classe, y = GC, fill = Classe)) +
      geom_violin(alpha = 0.7, trim = TRUE, scale = "width",
                  color = "black", linewidth = 0.5) +
      geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
      stat_summary(fun = median, geom = "point",
                   shape = 18, size = 3, color = "black") +
      geom_text(data = stats_gc,
                aes(x = Classe, y = Mediana, 
                    label = paste0(round(Mediana * 100, 1), "%")),
                vjust = -1, size = 3.5, fontface = "bold") +
      geom_text(data = stats_gc,
                aes(x = Classe, y = max(gc_violin_data$GC, na.rm = TRUE) * 0.95,
                    label = paste0("n=", n)),
                size = 3, color = "gray40") +
      scale_fill_viridis_d(option = "turbo", alpha = 0.8) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         breaks = seq(0, 1, 0.2),
                         limits = c(0, 1.05)) +
      labs(title = "Distribuição do Conteúdo GC por Classe Non-B",
           subtitle = paste("Total de regiões:", format(nrow(gc_violin_data), big.mark=",")),
           x = "Classe", y = "Conteúdo GC (%)") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")
    ggsave(file.path(dir_plots, "violino_gc_por_classe.png"),
           plot_violino,
           width = max(12, length(unique(gc_violin_data$Classe)) * 0.8),
           height = 8, dpi = 300, bg = "white")
    print(stats_gc)
    } else {
    cat("Nenhum dado GC disponível\n")}
} else {
  cat("Genoma não carregado - pulando análise GC\n")}
# 16) Análise de tamanho dos intervalos
tamanhos_data <- data.frame()
for (nome_arq in names(gr_list)) {
  gr <- gr_list[[nome_arq]]
  classe <- unique(mcols(gr)$classe)[1]
  temp_df <- data.frame(
    Classe = classe,
    Arquivo = nome_arq,
    Tamanho = width(gr),
    stringsAsFactors = FALSE)
  tamanhos_data <- rbind(tamanhos_data, temp_df)}
stats_tamanhos <- tamanhos_data %>%
  group_by(Classe) %>%
  summarise(
    Media = mean(Tamanho, na.rm = TRUE),
    Mediana = median(Tamanho, na.rm = TRUE),
    Q1 = quantile(Tamanho, 0.25, na.rm = TRUE),
    Q3 = quantile(Tamanho, 0.75, na.rm = TRUE),
    Min = min(Tamanho, na.rm = TRUE),
    Max = max(Tamanho, na.rm = TRUE),
    n = n(),
    .groups = "drop")
write.csv(stats_tamanhos, file.path(dir_tabelas, "estatisticas_tamanhos.csv"), row.names = FALSE)
p_tamanhos <- ggplot(tamanhos_data, aes(x = Tamanho, fill = Classe)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Distribuição dos Tamanhos dos Intervalos",
       x = "Tamanho (bp) - escala log", y = "Densidade") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(dir_plots, "densidade_tamanhos.png"), p_tamanhos, width = 10, height = 6, dpi = 300)
p_box_tamanhos <- ggplot(tamanhos_data, aes(x = Classe, y = Tamanho, fill = Classe)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Distribuição dos Tamanhos por Classe",
       x = "Classe", y = "Tamanho (bp) - escala log") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(file.path(dir_plots, "boxplot_tamanhos.png"), p_box_tamanhos, width = 10, height = 6, dpi = 300)

# 17) Salvar objetos R para uso futuro
saveRDS(gr_list, file.path(dir_resultados, "nonb_gr_list_completa.rds"))
saveRDS(nonb_por_classe, file.path(dir_resultados, "nonb_por_classe_intersecao.rds"))
saveRDS(gr_por_classe, file.path(dir_resultados, "nonb_gr_por_classe.rds"))
saveRDS(matriz_completa, file.path(dir_resultados, "nonb_matriz_interseccao.rds"))
# 18) Relatório final
cat("\nResumo dos arquivos gerados:\n")
cat("  - Plots:", length(list.files(dir_plots, pattern = "\\.png$")), "arquivos\n")
cat("  - Tabelas:", length(list.files(dir_tabelas, pattern = "\\.(tsv|csv)$")), "arquivos\n")
cat("  - BED interseções:", length(list.files(dir_intersecoes, pattern = "\\.bed$")), "arquivos\n")
cat("  - Objetos R:", length(list.files(dir_resultados, pattern = "\\.rds$")), "arquivos\n")
for (classe in names(gr_por_classe)) {
  n_total <- length(gr_por_classe[[classe]])
  n_intersec <- ifelse(classe %in% names(nonb_por_classe), length(nonb_por_classe[[classe]]), 0)
  n_arquivos <- sum(contagem_real$Classe == classe)
  cat(sprintf("  %-15s: %6d regiões totais | %6d em interseção | %d arquivo(s)\n",
              classe, n_total, n_intersec, n_arquivos))}

cat("\nPasta de resultados:", dir_resultados, "\n")
cat("Pasta de plots:", dir_plots, "\n")
cat("Pasta de tabelas:", dir_tabelas, "\n")
cat("Pasta de interseções:", dir_intersecoes, "\n")
stopCluster(cl)