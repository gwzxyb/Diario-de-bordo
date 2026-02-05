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
# ------------- legenda -----------
# nonb - non-B-form DNA

# 1) Diretórios e preparação

# Defição dos pontos de entrada e saída 
data_dir <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/igv/sobreposição"
dir_plots <- file.path(dir_resultados, "plots")

# Criação de diretórios para os resutados.
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_resultados, showWarnings = FALSE, recursive = TRUE)


# 2) Genoma de referência (para GC)

# Se o genoma não estiver disponível, o pipeline segue sem a etapa de GC

fa_path <- "/home/mcavalcante/dados_brutos/Hsalinarum.fa"
if (file.exists(fa_path)) {
  fa <- readDNAStringSet(fa_path)
  cat("Genoma carregado:", length(fa), "sequências\n")
} else {
  warning("Arquivo não encontrado.")
  fa <- NULL}


# 3) Seleção de arquivos Non-B
# filtra por arquivos por meio de : 
## a)palavras-chave no nome do arquivo (classe/tipo)
## b) extensões suportadas (BED/GFF/GFF3)
arquivos <- list.files(data_dir, full.names = TRUE)

arquivos_nonb <- arquivos[
  str_detect(arquivos, regex("z-dna|nonb|triplex|r-loop|short_tandem|A-phased|G-quadruplex|Cruciform",
                             ignore_case = TRUE)) &
    str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]


# 4) Classificação por classe de DNA alternativo
# feita a partir do nome do arquivo 
# define os agrupamentos e afeta todo o processo 
classificar_nonb <- function(arquivo) {
  nome <- basename(arquivo)
  nome_lower <- tolower(nome) #sem Case sensitive
  
  if (str_detect(nome_lower, "z-dna|zdna")) return("Z-DNA")
  if (str_detect(nome_lower, "g-quadruplex|g4")) return("G-quadruplex")
  if (str_detect(nome_lower, "triplex")) return("Triplex")
  if (str_detect(nome_lower, "r-loop|rloop")) return("R-loop")
  if (str_detect(nome_lower, "short_tandem|tandem")) return("Short_tandem")
  if (str_detect(nome_lower, "a-phased|aphased")) return("A-phased")
  if (str_detect(nome_lower, "cruciform")) return("Cruciform")
  
  # arquivos fora das regras caem em "Outros".
  return("Outros")}

nonb_classificado <- data.frame(
  arquivo = arquivos_nonb,
  nome_arquivo = basename(arquivos_nonb),
  classe = sapply(arquivos_nonb, classificar_nonb),
  stringsAsFactors = FALSE)

# Salva a tabela de classificação
write.table(nonb_classificado,
            file.path(dir_resultados, "nonb_classificacao.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Mostra quantos arquivos existem por classe
print(table(nonb_classificado$classe))
cat("\n")

# 5) Importação -> GRanges 
# Encapsula a importação e padronização de seqnames.
## Retorna o GRanges() vazio em falhas para o pipeline não apresentar erros
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
    
    #  remove "chr" para compatibilizar seqnames entre fontes.
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    return(gr)
  }, error = function(e) {
    warning(paste("Erro", ":", e$message))
    return(GRanges())})}


# 6) Carregamento de dados
#cada GRanges recebe dados p/ rastrear origem e classe 
gr_list <- list()

for (i in seq_len(nrow(nonb_classificado))) {
  arquivo <- nonb_classificado$arquivo[i]
  nome_arq <- nonb_classificado$nome_arquivo[i]
  classe <- nonb_classificado$classe[i]
  
  cat("Carregando:", nome_arq, "...")
  gr <- carregar_granges(arquivo)
  
  if (length(gr) > 0) {
    mcols(gr)$origem <- nome_arq
    mcols(gr)$classe <- classe
    mcols(gr)$arquivo_id <- i
    
    #Nome da lista: sem extensão 
    nome_lista <- gsub("\\.(bed|gff|gff3)$", "", nome_arq)
    gr_list[[nome_lista]] <- gr
    cat(" OK (", length(gr), " elementos)\n")
  } else {
    cat(" VAZIO ou ERRO\n")}}

# Remoção de entradas vazias; se tudo falhar, interrompe aqui
gr_list <- gr_list[sapply(gr_list, length) > 0]
if (length(gr_list) == 0) {
  stop("Nenhum arquivo Non-B foi carregado com sucesso.")}


# 7) Contagem por arquivo -> dado para fazer o grafico de quantidade de predições 

contagem_real <- data.frame(
  Arquivo = names(gr_list),
  Classe = sapply(gr_list, function(x) unique(mcols(x)$classe)[1]),
  Num_Predicoes = sapply(gr_list, length),
  stringsAsFactors = FALSE)

write.table(contagem_real,
            file.path(dir_resultados, "contagem_real_predicoes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 8) Consolidação por classe
# Junção de todos os arquivos da mesma classe em um único GRanges por classe
#-----> Reduce(c, ...) concatena; NÃO faz "reduce()" (não funde intervalos).
gr_por_classe <- list()
classes_com_dados <- unique(contagem_real$Classe)

for (classe in classes_com_dados) {
  arquivos_classe <- contagem_real$Arquivo[contagem_real$Classe == classe]
  grs_classe <- gr_list[arquivos_classe]
  
  if (length(grs_classe) > 0) {
    gr_por_classe[[classe]] <- Reduce(c, grs_classe)
    cat("Classe", classe, ":", length(gr_por_classe[[classe]]), "elementos\n")}}

# 9)Matriz de interseção completa

# Inclui todas as classes (mesmo sem dados) para manter matriz completa
classes_todas <- unique(nonb_classificado$classe)
n_todas <- length(classes_todas)

matriz_completa <- matrix(0, nrow = n_todas, ncol = n_todas,
                          dimnames = list(classes_todas, classes_todas))

## Diagonal -> total de predições por classe (contagem de intervalos)
for (classe in classes_todas) {
  if (classe %in% names(gr_por_classe)) {
    matriz_completa[classe, classe] <- length(gr_por_classe[[classe]])}}

# Fora da diagonal -> interseções segundo definição:
#   conta quantos intervalos da classe_i tocam pelo menos um intervalo da classe_j
for (i in 1:n_todas) {
  for (j in 1:n_todas) {
    if (i != j) {
      classe_i <- classes_todas[i]
      classe_j <- classes_todas[j]
      
      if (classe_i %in% names(gr_por_classe) && classe_j %in% names(gr_por_classe)) {
        if (length(gr_por_classe[[classe_i]]) > 0 && length(gr_por_classe[[classe_j]]) > 0) {
          inter <- findOverlaps(gr_por_classe[[classe_i]],
                                gr_por_classe[[classe_j]],
                                minoverlap = 1)
          
          #-------> unique(queryHits) -> conta intervalos únicos de classe_i que possuem pelo menos um overlap com classe_j.
          matriz_completa[i, j] <- length(unique(queryHits(inter)))}}}}}

matriz_classes <- matriz_completa
n_classes <- n_todas
classes_validas <- classes_todas

write.table(matriz_classes,
            file.path(dir_resultados, "matriz_interseccao_completa.tsv"),
            sep = "\t", quote = FALSE)

# 10) Visualização: heatmap com contagens e % do menor conjunto

# Converte matriz para formato longo para ggplot2.
matriz_melt <- melt(matriz_classes)
colnames(matriz_melt) <- c("Classe1", "Classe2", "Interseccao")

#   Para classes diferentes, usa Interseccao / min(total Classe1, total Classe2).
#   quanto do menor conjunto é coberto
matriz_melt$Porcentagem <- NA
for (i in 1:nrow(matriz_melt)) {
  total1 <- matriz_classes[matriz_melt$Classe1[i], matriz_melt$Classe1[i]]
  total2 <- matriz_classes[matriz_melt$Classe2[i], matriz_melt$Classe2[i]]
  
  if (matriz_melt$Classe1[i] == matriz_melt$Classe2[i]) {
    matriz_melt$Porcentagem[i] <- 100
  } else if (total1 > 0 && total2 > 0) {
    min_total <- min(total1, total2)
    if (min_total > 0) {
      matriz_melt$Porcentagem[i] <- round((matriz_melt$Interseccao[i] / min_total) * 100, 1)}}}

#Heatmap: diagonal mostra total; fora diagonal mostra contagem e %
# midpoint por mediana dos valores >0 melhora contraste em matrizes esparsas
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
                                         quantile(matriz_melt$Interseccao[matriz_melt$Interseccao > 0], 0.5, na.rm = TRUE),
                                         0),
                       na.value = "grey90",
                       name = "Número de\nInterseções") +
  labs(title = "Matriz de Interseção entre Classes de Estruturas Non-B DNA",
       subtitle = "Diagonal: total de predições | Fora diagonal: interseções (porcentagem do menor conjunto)",
       x = "Classe", y = "Classe") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "darkred"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  ) + coord_fixed(ratio = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

ggsave(file.path(dir_plots, "heatmap_classes_melhorado.png"),
       heatmap_classes,
       width = max(10, n_classes * 1.2),
       height = max(8, n_classes * 1.1),
       dpi = 300, bg = "white")

# 11) GC por intervalo:calcula GC por intervalo usando getSeq
## Seqnames de GRanges e FASTA precisam ser iguais
## GC = (G + C) / (A + C + G + T)
calcular_gc_real <- function(gr, fa_reference) {
  if (is.null(fa_reference) || length(gr) == 0) return(NA)
  tryCatch({
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    seqlevels(fa_reference) <- gsub("^chr", "", seqlevels(fa_reference))
    
    seqs <- getSeq(fa_reference, gr)
    gc_values <- sapply(as.character(seqs), function(s) {
      if (nchar(s) == 0) return(NA)
      s <- DNAString(s)
      af <- alphabetFrequency(s, baseOnly = TRUE)
      total <- sum(af[c("A", "C", "G", "T")])
      if (total > 0) return((af["C"] + af["G"]) / total)
      return(NA) })
    
    ## Remove NAs (intervalos fora do FASTA / problemas de seqname).
    return(gc_values[!is.na(gc_values)])
  }, error = function(e) {
    cat("Erro ao calcular GC:", e$message, "\n")
    return(NA)  })}

# Empilha GC em formato “tidy”: uma linha por intervalo 
gc_violin_data <- data.frame()

for (nome_arq in names(gr_list)) {
  gr <- gr_list[[nome_arq]]
  classe <- unique(mcols(gr)$classe)[1]
  
  if (!is.null(fa)) {
    gc_vals <- calcular_gc_real(gr, fa)
    
    if (!all(is.na(gc_vals)) && length(gc_vals) > 0) {
      temp_df <- data.frame(
        Classe = classe,
        Arquivo = nome_arq,
        GC = gc_vals,
        stringsAsFactors = FALSE )
      gc_violin_data <- rbind(gc_violin_data, temp_df)}}}

gc_violin_data <- gc_violin_data[!is.na(gc_violin_data$GC), ]

if (nrow(gc_violin_data) > 0) {
  cat("gerando plot de violino para conteúdo GC...\n")
  
  # Ordenação de classes e anotação de mediana e n.
  stats_gc <- gc_violin_data %>%
    group_by(Classe) %>%
    summarise(
      Mediana = median(GC, na.rm = TRUE),
      Media = mean(GC, na.rm = TRUE),
      Q1 = quantile(GC, 0.25, na.rm = TRUE),
      Q3 = quantile(GC, 0.75, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>% arrange(Mediana)
  
  gc_violin_data$Classe <- factor(gc_violin_data$Classe, levels = stats_gc$Classe)
  
  #Violino + boxplot
  plot_violino <- ggplot(gc_violin_data, aes(x = Classe, y = GC, fill = Classe)) +
    geom_violin(alpha = 0.7, trim = TRUE, scale = "width",
                color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
    stat_summary(fun = median, geom = "point",
                 shape = 18, size = 3, color = "black") +
    geom_text(data = stats_gc,
              aes(x = Classe, y = Mediana, label = paste0(round(Mediana * 100, 1), "%")),
              vjust = -1, size = 3.5, fontface = "bold", color = "black") +
    geom_text(data = stats_gc,
              aes(x = Classe, y = max(gc_violin_data$GC, na.rm = TRUE) * 0.95,
                  label = paste0("n=", n)),
              size = 3, color = "gray40") +
    scale_fill_viridis_d(option = "turbo", alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, 0.2),
                       limits = c(0, 1.05)) +
    labs(
      title = "Distribuição do Conteúdo GC por Classe de Estruturas Non-B",
      subtitle = "Valores calculados a partir do genoma de referência",
      x = "Classe de Estrutura Non-B",
      y = "Conteúdo GC (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "darkblue"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
  
  ggsave(file.path(dir_plots, "violino_gc_por_classe.png"),
         plot_violino,
         width = max(12, length(unique(gc_violin_data$Classe)) * 1.0),
         height = 8, dpi = 300, bg = "white")
} else {
  cat("Nenhum dado GC disponível para plot.\n")}
