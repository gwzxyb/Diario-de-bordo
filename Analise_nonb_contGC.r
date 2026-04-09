# --------------PACOTES-----------------------
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
data_dir <- "/home/mcavalcante/codigo/Diario-de-bordo/dados/igv/modificados"
dir_base <- "/home/mcavalcante/codigo/Diario-de-bordo/resultados Hsalinarum"
dir_resultados <- file.path(dir_base, "analise_nonB")  
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

# 3) Seleção de arquivos
arquivos <- list.files(data_dir, full.names = TRUE)
arq_gff <- arquivos[str_detect(arquivos, regex("Hsalinarum.*gff", ignore_case = TRUE))]
arq_tss <- "/home/mcavalcante/dados_brutos/Hsalinarum_TSS.gtf"

# CORRIGIDO: usar regex() para ignore_case
arquivos_nonb <- arquivos[
  str_detect(arquivos, regex("\\.bed$", ignore_case = TRUE)) &  
    str_detect(arquivos, regex("z-dna|nonb|triplex|R-loop|short_tandem|A-phased|G-quadruplex|Cruciform",
                               ignore_case = TRUE))]

if(length(arquivos_nonb) == 0) {
  warning("Nenhum arquivo BED de estruturas não-B encontrado!")
  cat("\nArquivos BED disponíveis:\n")
  # CORRIGIDO: usar regex() também aqui
  arquivos_bed <- arquivos[str_detect(arquivos, regex("\\.bed$", ignore_case = TRUE))]
  print(arquivos_bed)
} else {
  cat("\nArquivos BED encontrados:\n")
  print(basename(arquivos_nonb))
}

if(length(arq_gff) == 0) warning("arquivo GFF nao encontrado")
if(length(arq_tss) == 0) warning("arquivo TSS nao encontrado, algumas análises serão puladas")


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
    } else if (extensao %in% c("gff", "gff3", "gtf")) {
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

# 6) Carregamento de dados Non-B
cat("\n6) Carregando arquivos Non-B...\n")
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

gr_list <- list()
for (res in resultados_carregamento) {
  if (!is.null(res$gr)) {
    gr_list[[res$nome]] <- res$gr
    cat("  Carregado:", res$nome, "-", res$n, "elementos\n")}}

cat("\nTotal de arquivos carregados:", length(gr_list), "\n")

# 7) Agrupar arquivos por classe
cat("\n7) Agrupando por classe...\n")
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
  
  # Se só tem 1 arquivo, retorna ele mesmo
  if (n_arquivos == 1) {
    gr_intersect <- arquivos_classe[[1]]
    # Simplificar metadados para evitar problemas
    mcols(gr_intersect) <- NULL
    return(list(
      classe = classe, 
      gr_intersect = gr_intersect, 
      n_intersect = length(gr_intersect), 
      status = "único arquivo", 
      arquivos = names(arquivos_classe)
    ))
  }
  
  # Para múltiplos arquivos, calcular interseção
  # Converter todos para GRanges simples (sem metadados)
  gr_simplificados <- list()
  for (i in 1:n_arquivos) {
    gr_temp <- arquivos_classe[[i]]
    gr_simples <- GRanges(
      seqnames = seqnames(gr_temp),
      ranges = ranges(gr_temp),
      strand = strand(gr_temp)
    )
    gr_simplificados[[i]] <- gr_simples
  }
  
  # Calcular interseção progressiva
  gr_intersect <- gr_simplificados[[1]]
  for (i in 2:n_arquivos) {
    hits <- findOverlaps(gr_intersect, gr_simplificados[[i]])
    if (length(hits) == 0) {
      return(list(
        classe = classe, 
        gr_intersect = GRanges(), 
        n_intersect = 0, 
        status = "vazio", 
        arquivos = names(arquivos_classe)
      ))
    }
    gr_intersect <- pintersect(
      gr_intersect[queryHits(hits)], 
      gr_simplificados[[i]][subjectHits(hits)]
    )
    if (length(gr_intersect) == 0) {
      return(list(
        classe = classe, 
        gr_intersect = GRanges(), 
        n_intersect = 0, 
        status = "vazio", 
        arquivos = names(arquivos_classe)
      ))
    }
  }
  
  gr_intersect <- reduce(gr_intersect)
  
  return(list(
    classe = classe, 
    gr_intersect = gr_intersect, 
    n_intersect = length(gr_intersect), 
    status = "interseção", 
    arquivos = names(arquivos_classe)
  ))
}

# Executar em paralelo
classes <- names(arquivos_por_classe)
resultados_interseccao <- foreach(classe = classes, .packages = c("GenomicRanges")) %dopar% {
  calcular_interseccao_classe(classe, arquivos_por_classe[[classe]])
}

# Processar resultados
nonb_por_classe <- list()
resumo_intersecoes <- data.frame()

for (res in resultados_interseccao) {
  classe <- res$classe
  cat("\n", classe, ":\n")
  cat("  Arquivos:", paste(res$arquivos, collapse = ", "), "\n")
  
  if (res$n_intersect > 0) {
    nonb_por_classe[[classe]] <- res$gr_intersect
    
    # Adicionar metadados básicos
    mcols(res$gr_intersect)$classe <- classe
    mcols(res$gr_intersect)$n_arquivos_origem <- length(res$arquivos)
    mcols(res$gr_intersect)$name <- paste0(gsub("-", "_", classe), "_", sprintf("%04d", 1:res$n_intersect))
    
    # Salvar BED
    nome_bed <- file.path(dir_intersecoes, sprintf("nonB_%s_interseccao.bed", gsub("-", "_", classe)))
    export(res$gr_intersect, nome_bed, format = "BED")
    cat("  BED salvo:", basename(nome_bed), "\n")
    cat("  Regiões:", res$n_intersect, "\n")
    
    resumo_intersecoes <- rbind(resumo_intersecoes, data.frame(
      Classe = classe,
      Regioes_Interseccao = res$n_intersect,
      Arquivos_Origem = length(res$arquivos),
      Arquivos = paste(res$arquivos, collapse = "; "),
      stringsAsFactors = FALSE
    ))
  } else {
    cat("  Interseção vazia\n")
  }
}

# Salvar resumo
if (nrow(resumo_intersecoes) > 0) {
  write.table(resumo_intersecoes, 
              file.path(dir_intersecoes, "resumo_intersecoes_nonB.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\n Resumo final:\n")
  print(resumo_intersecoes)
  cat("\nTotal de classes com interseção:", nrow(resumo_intersecoes), "\n")
  cat("Total de regiões em interseção:", sum(resumo_intersecoes$Regioes_Interseccao), "\n")
}
# 9) Carregar GFF e TSS
cat("\n9) Carregando GFF e TSS...\n")
if(length(arq_gff) > 0) {
  genes_all <- carregar_granges(arq_gff[1])
  if(!is.null(genes_all) && length(genes_all) > 0) {
    stopifnot("type" %in% colnames(mcols(genes_all)))
    genes <- genes_all[genes_all$type == "gene"]
    genes$region <- "gene"
    cat("  Genes carregados:", length(genes), "\n") } }

if(length(arq_tss) > 0) {
  tss <- carregar_granges(arq_tss)
  if(!is.null(tss) && length(tss) > 0) {
    stopifnot("type" %in% colnames(mcols(tss)))
    
    # Filtrar apenas TSS (remover outras features se houver)
    tss <- tss[tss$type %in% c("TSS", "transcription_start_site")]
    
    # Filtrar APENAS TSS classe P
    if("tss_class" %in% colnames(mcols(tss))) {
      tss_P <- tss[tss$tss_class == "P"]
      cat("  TSS classe P:", length(tss_P), "\n")
      
      # Mostrar exemplo de um TSS classe P para confirmar
      if(length(tss_P) > 0) {
        cat("  Exemplo de TSS classe P:\n")
        exemplo <- head(tss_P, 1)
        cat("    ", as.character(seqnames(exemplo)), ":", start(exemplo), "-", end(exemplo), 
            "| strand:", as.character(strand(exemplo)), 
            "| tss_class:", mcols(exemplo)$tss_class, "\n")
      }
    } else {
      warning("Coluna 'tss_class' não encontrada. Usando todos os TSS.")
      tss_P <- tss
      cat("  TSS (sem classificação):", length(tss_P), "\n")
    }
    
    # Promotores (250bp upstream do TSS)
    promotores <- promoters(tss_P, upstream = 250, downstream = 0)
    promotores$region <- "promoter"
    cat("  Promotores gerados:", length(promotores), "\n")
    
    # Salvar tabela de promotores
    promotores_df <- as.data.frame(promotores) %>% 
      transmute(seqname = seqnames, start = start, end = end, strand = strand, width = width, tss_id = names(promotores))
    write.csv(promotores_df, file.path(dir_tabelas, "promotores_250bp.csv"), row.names = FALSE)
    
  } } else {
    tss_P <- GRanges()
    promotores <- GRanges() }
# 10) Identificar genes com e sem 5' UTR
cat("\n10) Identificando genes com 5' UTR...\n")

# Função para encontrar o gene associado a cada TSS
encontrar_gene_por_tss <- function(tss_gr, genes_gr) {
  result <- data.frame()
  for (i in seq_along(tss_gr)) {
    tss_pos <- tss_gr[i]
    seqn <- as.character(seqnames(tss_pos))
    strnd <- as.character(strand(tss_pos))
    
    # Genes no mesmo cromossomo e strand
    genes_cand <- genes_gr[seqnames(genes_gr) == seqn & strand(genes_gr) == strnd]
    if (length(genes_cand) == 0) next
    
    if (strnd == "+") {
      # Gene mais próximo à direita do TSS
      genes_dir <- genes_cand[end(genes_cand) > start(tss_pos)]
      if (length(genes_dir) > 0) {
        gene_closest <- genes_dir[which.min(start(genes_dir))]
        result <- rbind(result, data.frame(
          tss_id = i,
          gene_id = mcols(gene_closest)$ID,
          gene_name = ifelse("Name" %in% colnames(mcols(gene_closest)), mcols(gene_closest)$Name, NA),
          seqname = seqn,
          tss_start = start(tss_pos),
          tss_end = end(tss_pos),
          gene_start = start(gene_closest),
          gene_end = end(gene_closest),
          strand = strnd,
          stringsAsFactors = FALSE)) } }
    else if (strnd == "-") {
      # Gene mais próximo à esquerda do TSS (invertido)
      genes_esq <- genes_cand[start(genes_cand) < end(tss_pos)]
      if (length(genes_esq) > 0) {
        gene_closest <- genes_esq[which.max(end(genes_esq))]
        result <- rbind(result, data.frame(
          tss_id = i,
          gene_id = mcols(gene_closest)$ID,
          gene_name = ifelse("Name" %in% colnames(mcols(gene_closest)), mcols(gene_closest)$Name, NA),
          seqname = seqn,
          tss_start = start(tss_pos),
          tss_end = end(tss_pos),
          gene_start = start(gene_closest),
          gene_end = end(gene_closest),
          strand = strnd,
          stringsAsFactors = FALSE)) } } }
  return(result) }

# Encontrar associação TSS-genes
if(exists("tss_P") && length(tss_P) > 0 && exists("genes") && length(genes) > 0) {
  tss_gene_assoc <- encontrar_gene_por_tss(tss_P, genes)
  
  if(nrow(tss_gene_assoc) > 0) {
    # Calcular distância TSS ao start do gene
    tss_gene_assoc <- tss_gene_assoc %>%
      mutate(distancia_tss_gene = ifelse(strand == "+", gene_start - tss_end, tss_start - gene_end))
    
    # Identificar genes com 5' UTR (distância positiva = tem espaço entre TSS e gene)
    genes_com_utr <- tss_gene_assoc %>% filter(distancia_tss_gene > 0)
    genes_sem_utr <- tss_gene_assoc %>% filter(distancia_tss_gene <= 0)
    
    cat("\n  Genes com 5' UTR:", nrow(genes_com_utr), "\n")
    cat("  Genes sem 5' UTR:", nrow(genes_sem_utr), "\n")
    
    # Tabela completa de genes e status 5' UTR
    tabela_genes_utr <- tss_gene_assoc %>%
      mutate(tem_5utr = distancia_tss_gene > 0) %>%
      select(gene_id, gene_name, seqname, strand, tss_start, gene_start, gene_end, distancia_tss_gene, tem_5utr)
    
    write.csv(tabela_genes_utr, file.path(dir_tabelas, "genes_status_5utr.csv"), row.names = FALSE)
    cat("  Tabela de genes com status 5' UTR salva\n")
    
    # Estatísticas
    stats_utr <- data.frame(
      Categoria = c("Genes com 5' UTR", "Genes sem 5' UTR", "Total genes com TSS"),
      Quantidade = c(nrow(genes_com_utr), nrow(genes_sem_utr), nrow(tss_gene_assoc)),
      Percentual = c(round(100 * nrow(genes_com_utr)/nrow(tss_gene_assoc), 1),
                     round(100 * nrow(genes_sem_utr)/nrow(tss_gene_assoc), 1),
                     100)
    )
    write.csv(stats_utr, file.path(dir_tabelas, "estatisticas_5utr.csv"), row.names = FALSE)
    print(stats_utr)
    
    # Gráfico de pizza: genes com e sem 5' UTR
    p_utr_pie <- ggplot(stats_utr[1:2,], aes(x = 2, y = Quantidade, fill = Categoria)) +
      geom_col(width = 1, color = "white", linewidth = 0.5) +
      coord_polar(theta = "y") + xlim(0.5, 2.5) +
      geom_text(aes(label = paste0(Quantidade, " (", Percentual, "%)")), 
                position = position_stack(vjust = 0.5), size = 4.5, fontface = "bold", color = "white") +
      scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
      labs(title = "Genes com e sem 5' UTR", fill = "") +
      theme_void() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    ggsave(file.path(dir_plots, "genes_com_sem_5utr.png"), p_utr_pie, width = 8, height = 6, dpi = 300)
  }
}

# 11) Criar 5' UTRs
cat("\n11) Calculando 5' UTRs...\n")
criar_utr5_custom <- function(genes_gr, tss_gr) {
  if(length(tss_gr) == 0) return(GRanges())
  utr_list <- vector("list", length(tss_gr))
  for (i in seq_along(tss_gr)) {
    tss_pos <- tss_gr[i]; seqn <- as.character(seqnames(tss_pos)); strnd <- as.character(strand(tss_pos))
    genes_cand <- genes_gr[seqnames(genes_gr) == seqn & strand(genes_gr) == strnd]
    if (length(genes_cand) == 0) next
    if (strnd == "+") {
      genes_dir <- genes_cand[end(genes_cand) > start(tss_pos)]
      if (length(genes_dir) > 0) {
        gene_closest <- genes_dir[which.min(start(genes_dir))]
        utr_start <- end(tss_pos); utr_end <- start(gene_closest) - 1
        if (utr_end >= utr_start) {
          utr_list[[i]] <- GRanges(seqn, IRanges(utr_start, utr_end), strand = strnd, 
                                   region = "5_UTR", gene_id = mcols(gene_closest)$ID) } }
    } else if (strnd == "-") {
      genes_esq <- genes_cand[start(genes_cand) < end(tss_pos)]
      if (length(genes_esq) > 0) {
        gene_closest <- genes_esq[which.max(end(genes_esq))]
        utr_start <- end(gene_closest) + 1; utr_end <- start(tss_pos)
        if (utr_end >= utr_start) {
          utr_list[[i]] <- GRanges(seqn, IRanges(utr_start, utr_end), strand = strnd, 
                                   region = "5_UTR", gene_id = mcols(gene_closest)$ID) } } } }
  utr_valid <- utr_list[!sapply(utr_list, is.null)]
  if (length(utr_valid) == 0) return(GRanges())
  utr5_final <- unlist(GRangesList(utr_valid))
  return(trim(utr5_final)) }

if(length(tss_P) > 0 && length(genes) > 0) {
  utr5_custom <- criar_utr5_custom(genes, tss_P)
  cat("  5' UTRs gerados:", length(utr5_custom), "\n")
  
  if(length(utr5_custom) > 0) {
    # Tabela 5' UTR com informações dos genes
    utr5_table <- as.data.frame(utr5_custom) %>% 
      transmute(seqname = seqnames, start = start, end = end, strand = strand, 
                width = width, gene_id = gene_id, region = "5_UTR") %>% 
      arrange(seqname, start)
    write.csv(utr5_table, file.path(dir_resultados, "5_UTR_custom.csv"), row.names = FALSE)
    cat("  Tabela 5' UTR salva\n")
    
    # Estatísticas de tamanho das 5' UTRs
    stats_utr_tamanho <- utr5_table %>%
      summarise(media = mean(width), mediana = median(width), 
                min = min(width), max = max(width), q1 = quantile(width, 0.25), q3 = quantile(width, 0.75))
    write.csv(stats_utr_tamanho, file.path(dir_tabelas, "estatisticas_tamanho_5utr.csv"), row.names = FALSE)
    
    # Histograma do tamanho das 5' UTRs
    p_utr_hist <- ggplot(utr5_table, aes(x = width)) +
      geom_histogram(bins = 30, fill = "#377EB8", alpha = 0.7, color = "white") +
      geom_vline(xintercept = stats_utr_tamanho$media, color = "#E41A1C", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = stats_utr_tamanho$mediana, color = "#4DAF4A", linetype = "dotted", linewidth = 1) +
      annotate("text", x = stats_utr_tamanho$media + 5, y = max(table(cut(utr5_table$width, 30))) * 0.8, 
               label = paste("Média =", round(stats_utr_tamanho$media, 1), "bp"), 
               color = "#E41A1C", size = 3.5, hjust = 0, fontface = "bold") +
      annotate("text", x = stats_utr_tamanho$mediana + 5, y = max(table(cut(utr5_table$width, 30))) * 0.7, 
               label = paste("Mediana =", round(stats_utr_tamanho$mediana, 1), "bp"), 
               color = "#4DAF4A", size = 3.5, hjust = 0, fontface = "bold") +
      labs(title = "Distribuição do tamanho das 5' UTRs", x = "Tamanho (bp)", y = "Frequência") +
      theme_minimal(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(dir_plots, "histograma_tamanho_5utr.png"), p_utr_hist, width = 10, height = 6, dpi = 300)
    
    # Boxplot do tamanho por strand
    p_utr_box <- ggplot(utr5_table, aes(x = strand, y = width, fill = strand)) +
      geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.3) +
      scale_fill_manual(values = c("#F28AA0", "#366899")) +
      labs(title = "Tamanho das 5' UTRs por fita", x = "Fita", y = "Tamanho (bp)") +
      theme_minimal() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(dir_plots, "boxplot_5utr_por_strand.png"), p_utr_box, width = 6, height = 6, dpi = 300)
    
  } } else { utr5_custom <- GRanges() }

# 12) Criar bins dos genes
cat("\n12) Gerando bins dos genes...\n")
bins_por_gene <- function(gr_gene, n_bins = 10) {
  gene_length <- width(gr_gene)
  starts <- start(gr_gene) - 1 + floor(seq(0, gene_length, length.out = n_bins + 1)[-(n_bins + 1)]) + 1
  ends <- start(gr_gene) - 1 + floor(seq(0, gene_length, length.out = n_bins + 1)[-1])
  bins <- GRanges(seqnames = seqnames(gr_gene), ranges = IRanges(start = starts, end = ends), strand = strand(gr_gene))
  bins$region <- paste0("gene_bin_", seq_len(n_bins))
  valid <- start(bins) <= end(bins)
  return(bins[valid]) }

genes_list <- split(genes, seq_along(genes))
gene_bins_list <- foreach(gene = genes_list, .packages = c("GenomicRanges"), .combine = c) %dopar% { bins_por_gene(gene) }
gene_bins <- sort(gene_bins_list)
cat("  Bins gerados:", length(gene_bins), "\n")

# 13) Concatenar todas as regiões
cat("\n13) Consolidando regiões gênicas...\n")
regioes_list <- list()
if(length(promotores) > 0) regioes_list$promotores <- promotores
if(length(utr5_custom) > 0) regioes_list$utr5 <- utr5_custom
if(length(gene_bins) > 0) regioes_list$bins <- gene_bins
todas_regioes <- trim(Reduce(c, regioes_list))
cat("  Total de regiões para análise:", length(todas_regioes), "\n")

# 14) Analisar sobreposições Non-B vs regiões gênicas
cat("\n14) Analisando sobreposições Non-B vs regiões gênicas...\n")
if(length(nonb_por_classe) > 0 && length(todas_regioes) > 0) {
  overlap_list <- list()
  for (classe in names(nonb_por_classe)) {
    nonb_gr <- nonb_por_classe[[classe]]
    ov <- findOverlaps(todas_regioes, nonb_gr, ignore.strand = TRUE)
    if (length(ov) > 0) {
      reg_names <- as.character(todas_regioes$region)[queryHits(ov)]
      df_temp <- tibble(classe = classe, region = reg_names, nb_index = subjectHits(ov))
      overlap_list[[classe]] <- df_temp } }
  overlap_df <- bind_rows(overlap_list)
  if (nrow(overlap_df) > 0) {
    overlap_df <- overlap_df %>%
      mutate(region = recode(region, promoter = "Promotor", `5_UTR` = "5' UTR")) %>%
      mutate(region = ifelse(str_detect(region, "^gene_bin_"), paste0(as.numeric(str_extract(region, "\\d+")) * 10, "%"), region))
    write.csv(overlap_df, file.path(dir_resultados, "sobreposicoes_nonB_regioes.csv"), row.names = FALSE)
    cat("  Tabela de sobreposições salva\n")
    
    ordem_regioes <- c("Promotor", "5' UTR", paste0(seq(10, 100, by = 10), "%"))
    overlap_plot <- overlap_df %>% mutate(region = factor(region, levels = ordem_regioes)) %>% group_by(classe, region) %>% summarise(frequencia = n(), .groups = "drop") %>% filter(!is.na(region))
    
    cores_padrao <- c("Z-DNA" = "#F28AA0", "G-quadruplex" = "#366899", "Triplex" = "#81DACA", "R-loop" = "#FC8D62", "Short_tandem" = "#66C2A5", "A-phased" = "#8DA0CB", "Cruciform" = "#377EB8", "Outros" = "#BEBADA")
    
    for (i in unique(overlap_plot$classe)) {
      df_plot <- overlap_plot %>% filter(classe == i)
      if (nrow(df_plot) == 0) next
      cor_classe <- ifelse(is.na(cores_padrao[i]), "#BEBADA", cores_padrao[i])
      p_categ <- ggplot(df_plot, aes(x = region, y = frequencia)) +
        geom_col(fill = cor_classe, color = "black", width = 0.7, alpha = 0.8) +
        labs(title = paste("Distribuição", i, "por região gênica"), subtitle = paste("Regiões:", sum(df_plot$frequencia)), x = "Região", y = "Non-B sobrepostos") +
        theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
      ggsave(file.path(dir_plots, paste0("intersec_", gsub(" ", "_", i), ".png")), p_categ, width = 12, height = 8, dpi = 300) }
    
    p_stack <- ggplot(overlap_plot, aes(x = region, y = frequencia, fill = classe)) +
      geom_col(position = "stack", color = "white", linewidth = 0.3) +
      scale_fill_manual(values = cores_padrao[intersect(names(cores_padrao), unique(overlap_plot$classe))]) +
      labs(title = "Distribuição Non-B por região gênica", x = "Região", y = "Total", fill = "Classe") +
      theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    ggsave(file.path(dir_plots, "stacked_nonB_regioes.png"), p_stack, width = 14, height = 9, dpi = 300)
    
    resumo <- overlap_plot %>% group_by(classe, region) %>% summarise(total = sum(frequencia), .groups = "drop") %>% pivot_wider(names_from = region, values_from = total, values_fill = 0)
    write.csv(resumo, file.path(dir_tabelas, "resumo_nonB_por_regiao.csv"), row.names = FALSE)
    print(resumo) } }

# 15) Contagem por arquivo
cat("\n15) Contagem por arquivo...\n")
contagem_real <- data.frame(
  Arquivo = names(gr_list),
  Classe = sapply(gr_list, function(x) unique(mcols(x)$classe)[1]),
  Num_Predicoes = sapply(gr_list, length),
  stringsAsFactors = FALSE)
write.table(contagem_real, file.path(dir_tabelas, "contagem_real_predicoes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# 16) Gráfico de barras das interseções
if (nrow(resumo_intersecoes) > 0) {
  p_intersec <- ggplot(resumo_intersecoes, aes(x = reorder(Classe, Regioes_Interseccao), y = Regioes_Interseccao, fill = Classe)) +
    geom_col(alpha = 0.8, width = 0.7, color = "white") + geom_text(aes(label = Regioes_Interseccao), hjust = -0.1, size = 4) +
    coord_flip() + scale_fill_viridis_d(option = "turbo") +
    labs(title = "Regiões em Interseção por Classe Non-B", x = "Classe", y = "Número de regiões") +
    theme_minimal(base_size = 12) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(dir_plots, "interseccoes_por_classe.png"), p_intersec, width = 10, height = 6, dpi = 300) }

# 17) Comparação: total vs interseção
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
    labs(title = "Comparação: Total vs Interseção", x = "Classe", y = "Número de regiões") +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(dir_plots, "comparacao_total_vs_interseccao.png"), p_comp, width = 12, height = 7, dpi = 300)
  write.table(comparacao, file.path(dir_tabelas, "comparacao_total_vs_interseccao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE) }

# 18) Consolidar por classe (para matriz)
cat("\n18) Consolidando por classe para matriz...\n")
gr_por_classe <- list()
classes_com_dados <- unique(contagem_real$Classe)
for (classe in classes_com_dados) {
  arquivos_classe <- contagem_real$Arquivo[contagem_real$Classe == classe]
  grs_classe <- gr_list[arquivos_classe]
  if (length(grs_classe) > 0) {
    gr_por_classe[[classe]] <- Reduce(c, grs_classe)
    cat("  Classe", classe, ":", length(gr_por_classe[[classe]]), "elementos\n") } }

# 19) Matriz de interseção completa
cat("\n19) Calculando matriz de interseção...\n")
classes_todas <- unique(nonb_classificado$classe)
n_todas <- length(classes_todas)
matriz_completa <- matrix(0, nrow = n_todas, ncol = n_todas, dimnames = list(classes_todas, classes_todas))
for (classe in classes_todas) if (classe %in% names(gr_por_classe)) matriz_completa[classe, classe] <- length(gr_por_classe[[classe]])
for (i in 1:n_todas) for (j in 1:n_todas) if (i != j) {
  classe_i <- classes_todas[i]; classe_j <- classes_todas[j]
  if (classe_i %in% names(gr_por_classe) && classe_j %in% names(gr_por_classe)) {
    if (length(gr_por_classe[[classe_i]]) > 0 && length(gr_por_classe[[classe_j]]) > 0) {
      inter <- findOverlaps(gr_por_classe[[classe_i]], gr_por_classe[[classe_j]], minoverlap = 1)
      matriz_completa[i, j] <- length(unique(queryHits(inter))) } } }
write.table(matriz_completa, file.path(dir_tabelas, "matriz_interseccao_completa.tsv"), sep = "\t", quote = FALSE)

# 20) Heatmap da matriz de interseção
cat("\n20) Gerando heatmap...\n")
matriz_melt <- melt(matriz_completa)
colnames(matriz_melt) <- c("Classe1", "Classe2", "Interseccao")
matriz_melt$Porcentagem <- NA
for (i in 1:nrow(matriz_melt)) {
  total1 <- matriz_completa[matriz_melt$Classe1[i], matriz_melt$Classe1[i]]
  total2 <- matriz_completa[matriz_melt$Classe2[i], matriz_melt$Classe2[i]]
  if (matriz_melt$Classe1[i] == matriz_melt$Classe2[i]) matriz_melt$Porcentagem[i] <- 100
  else if (total1 > 0 && total2 > 0) {
    min_total <- min(total1, total2)
    if (min_total > 0) matriz_melt$Porcentagem[i] <- round((matriz_melt$Interseccao[i] / min_total) * 100, 1) } }

heatmap_classes <- ggplot(matriz_melt, aes(x = Classe1, y = Classe2, fill = Interseccao)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(Classe1 == Classe2, paste0(Interseccao), ifelse(Interseccao > 0, paste0(Interseccao, "\n(", ifelse(is.na(Porcentagem), "0", Porcentagem), "%)"), ""))), color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "lightblue", high = "red", midpoint = ifelse(length(matriz_melt$Interseccao[matriz_melt$Interseccao > 0]) > 0, median(matriz_melt$Interseccao[matriz_melt$Interseccao > 0], na.rm = TRUE), 0), name = "Interseções") +
  labs(title = "Matriz de Interseção Non-B", x = "Classe", y = "Classe") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.text.y = element_text(face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank()) + coord_fixed()
ggsave(file.path(dir_plots, "heatmap_interseccao_classes.png"), heatmap_classes, width = max(12, n_todas * 0.8), height = max(10, n_todas * 0.8), dpi = 300, bg = "white")

# 21) Análise de conteúdo GC
cat("\n21) Análise de conteúdo GC...\n")
calcular_gc_real <- function(gr, fa_reference) {
  if (is.null(fa_reference) || length(gr) == 0) return(NA)
  tryCatch({
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    valid_seqs <- intersect(seqlevels(gr), names(fa_reference))
    if(length(valid_seqs) == 0) return(NA)
    gr <- gr[seqnames(gr) %in% valid_seqs]; gr <- trim(gr)
    if(length(gr) == 0) return(NA)
    seqs <- getSeq(fa_reference, gr)
    gc_values <- sapply(as.character(seqs), function(s) {
      if (nchar(s) == 0) return(NA)
      s <- DNAString(s); af <- alphabetFrequency(s, baseOnly = TRUE)
      total <- sum(af[c("A", "C", "G", "T")])
      if (total > 0) return((af["C"] + af["G"]) / total)
      return(NA) })
    return(gc_values[!is.na(gc_values)])
  }, error = function(e) { return(NA) }) }

gc_violin_data <- data.frame()
if (!is.null(fa)) {
  cat("  Processando", length(gr_list), "arquivos...\n")
  for (nome_arq in names(gr_list)) {
    gr <- gr_list[[nome_arq]]; classe <- unique(mcols(gr)$classe)[1]
    cat("    ", sprintf("%-25s", nome_arq), "...")
    gc_vals <- calcular_gc_real(gr, fa)
    if (!all(is.na(gc_vals)) && length(gc_vals) > 0) {
      temp_df <- data.frame(Classe = classe, Arquivo = nome_arq, GC = gc_vals, stringsAsFactors = FALSE)
      gc_violin_data <- rbind(gc_violin_data, temp_df)
      cat(" ", length(gc_vals), "valores GC\n") } else { cat(" sem dados GC\n") } }
  gc_violin_data <- gc_violin_data[!is.na(gc_violin_data$GC), ]
  if (nrow(gc_violin_data) > 0) {
    stats_gc <- gc_violin_data %>% group_by(Classe) %>% summarise(Mediana = median(GC, na.rm = TRUE), Media = mean(GC, na.rm = TRUE), Q1 = quantile(GC, 0.25, na.rm = TRUE), Q3 = quantile(GC, 0.75, na.rm = TRUE), Min = min(GC, na.rm = TRUE), Max = max(GC, na.rm = TRUE), n = n(), .groups = "drop") %>% arrange(Mediana)
    write.csv(stats_gc, file.path(dir_tabelas, "estatisticas_gc_por_classe.csv"), row.names = FALSE)
    gc_violin_data$Classe <- factor(gc_violin_data$Classe, levels = stats_gc$Classe)
    plot_violino <- ggplot(gc_violin_data, aes(x = Classe, y = GC, fill = Classe)) +
      geom_violin(alpha = 0.7, trim = TRUE, scale = "width", color = "black", linewidth = 0.5) +
      geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
      stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black") +
      geom_text(data = stats_gc, aes(x = Classe, y = Mediana, label = paste0(round(Mediana * 100, 1), "%")), vjust = -1, size = 3.5, fontface = "bold") +
      geom_text(data = stats_gc, aes(x = Classe, y = max(gc_violin_data$GC, na.rm = TRUE) * 0.95, label = paste0("n=", n)), size = 3, color = "gray40") +
      scale_fill_viridis_d(option = "turbo", alpha = 0.8) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.2), limits = c(0, 1.05)) +
      labs(title = "Distribuição do GC por Classe", x = "Classe", y = "GC (%)") +
      theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
    ggsave(file.path(dir_plots, "violino_gc_por_classe.png"), plot_violino, width = max(12, length(unique(gc_violin_data$Classe)) * 0.8), height = 8, dpi = 300, bg = "white")
    print(stats_gc) } else { cat("  Nenhum dado GC disponível\n") } } else { cat("  Genoma não carregado - pulando análise GC\n") }

# 22) Análise de tamanho dos intervalos
cat("\n22) Análise de tamanho dos intervalos...\n")
tamanhos_data <- data.frame()
for (nome_arq in names(gr_list)) {
  gr <- gr_list[[nome_arq]]; classe <- unique(mcols(gr)$classe)[1]
  temp_df <- data.frame(Classe = classe, Arquivo = nome_arq, Tamanho = width(gr), stringsAsFactors = FALSE)
  tamanhos_data <- rbind(tamanhos_data, temp_df) }
stats_tamanhos <- tamanhos_data %>% group_by(Classe) %>% summarise(Media = mean(Tamanho, na.rm = TRUE), Mediana = median(Tamanho, na.rm = TRUE), Q1 = quantile(Tamanho, 0.25, na.rm = TRUE), Q3 = quantile(Tamanho, 0.75, na.rm = TRUE), Min = min(Tamanho, na.rm = TRUE), Max = max(Tamanho, na.rm = TRUE), n = n(), .groups = "drop")
write.csv(stats_tamanhos, file.path(dir_tabelas, "estatisticas_tamanhos.csv"), row.names = FALSE)

p_tamanhos <- ggplot(tamanhos_data, aes(x = Tamanho, fill = Classe)) + geom_density(alpha = 0.5) + scale_x_log10() + scale_fill_viridis_d(option = "turbo") + labs(title = "Distribuição dos Tamanhos", x = "Tamanho (bp) - log", y = "Densidade") + theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(dir_plots, "densidade_tamanhos.png"), p_tamanhos, width = 10, height = 6, dpi = 300)

p_box_tamanhos <- ggplot(tamanhos_data, aes(x = Classe, y = Tamanho, fill = Classe)) + geom_boxplot(alpha = 0.7) + scale_y_log10() + scale_fill_viridis_d(option = "turbo") + labs(title = "Tamanhos por Classe", x = "Classe", y = "Tamanho (bp) - log") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(dir_plots, "boxplot_tamanhos.png"), p_box_tamanhos, width = 10, height = 6, dpi = 300)

# 23) Salvar objetos R
cat("\n23) Salvando objetos R...\n")
saveRDS(gr_list, file.path(dir_resultados, "nonb_gr_list_completa.rds"))
saveRDS(nonb_por_classe, file.path(dir_resultados, "nonb_por_classe_intersecao.rds"))
saveRDS(gr_por_classe, file.path(dir_resultados, "nonb_gr_por_classe.rds"))
saveRDS(matriz_completa, file.path(dir_resultados, "nonb_matriz_interseccao.rds"))
if(exists("utr5_custom") && length(utr5_custom) > 0) saveRDS(utr5_custom, file.path(dir_resultados, "utr5_custom.rds"))
if(exists("promotores") && length(promotores) > 0) saveRDS(promotores, file.path(dir_resultados, "promotores.rds"))

# 24) Relatório final
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RELATÓRIO FINAL\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("\nResumo dos arquivos gerados:\n")
cat("  - Plots:", length(list.files(dir_plots, pattern = "\\.png$")), "arquivos\n")
cat("  - Tabelas:", length(list.files(dir_tabelas, pattern = "\\.(tsv|csv)$")), "arquivos\n")
cat("  - BED interseções:", length(list.files(dir_intersecoes, pattern = "\\.bed$")), "arquivos\n")
cat("  - Objetos R:", length(list.files(dir_resultados, pattern = "\\.rds$")), "arquivos\n")

cat("\nEstruturas Non-B processadas:\n")
for (classe in names(gr_por_classe)) {
  n_total <- length(gr_por_classe[[classe]])
  n_intersec <- ifelse(classe %in% names(nonb_por_classe), length(nonb_por_classe[[classe]]), 0)
  n_arquivos <- sum(contagem_real$Classe == classe)
  cat(sprintf("  %-15s: %6d regiões totais | %6d em interseção | %d arquivo(s)\n", classe, n_total, n_intersec, n_arquivos)) }

# Estatísticas de 5' UTR
if(exists("utr5_custom") && length(utr5_custom) > 0) {
  cat("\nEstatísticas 5' UTR:\n")
  cat("  Total 5' UTRs:", length(utr5_custom), "\n")
  cat("  Tamanho médio:", round(mean(width(utr5_custom)), 1), "bp\n")
  cat("  Tamanho mediano:", median(width(utr5_custom)), "bp\n")
  cat("  Mínimo:", min(width(utr5_custom)), "bp | Máximo:", max(width(utr5_custom)), "bp\n") }

cat("\nPastas de resultados:\n")
cat("  Resultados principais:", dir_resultados, "\n")
cat("  Plots:", dir_plots, "\n")
cat("  Tabelas:", dir_tabelas, "\n")
cat("  Interseções BED:", dir_intersecoes, "\n")

# 25) Finalizar
stopCluster(cl)
cat("\nANÁLISE CONCLUÍDA\n")