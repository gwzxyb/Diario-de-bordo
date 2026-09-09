# ANÁLISE DE ESTRUTURAS DE DNA NÃO-B EM H. salinarum
# Descrição : Identifica, classifica e analisa regiões de DNA não-canônico
#             (non-B DNA) no genoma de H. salinarum, cruzando com anotações
#             gênicas (GFF), sítios de início de transcrição (TSS) e regiões
#             regulatórias (promotores e 5' UTRs) 

# Pacotes
suppressPackageStartupMessages({
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
  library(scales)})

#namespace
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by  <- dplyr::group_by

# Configuração
## Paralelização
NUM_CORES <- min(12L, detectCores(logical = FALSE))
cl <- makeCluster(NUM_CORES)
registerDoParallel(cl)

## Paleta de cores non-B
# Usar e em TODOS os gráficos
CORES_NONB <- c(
  "Z-DNA"        = "#E63946",
  "G-quadruplex" = "#1D3557",
  "Triplex"      = "#2A9D8F",
  "R-loop"       = "#F4A261",
  "Short_tandem" = "#57CC99",
  "A-phased"     = "#6A4C93",
  "Cruciform"    = "#264653",
  "Outros"       = "#ADB5BD")

## Tema ggplot2
tema_nonb <- theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    plot.caption  = element_text(hjust = 0, color = "grey50", size = 9,
                                 margin = margin(t = 8)),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    strip.text       = element_text(face = "bold"))
theme_set(tema_nonb)

# Arquivos
data_dir        <- "./Diario-de-bordo/dados/igv/modificados"
dir_base        <- "./Diario-de-bordo/resultados Hsalinarum"
dir_resultados  <- file.path(dir_base, "analise_nonB")
dir_plots       <- file.path(dir_resultados, "plots")
dir_tabelas     <- file.path(dir_resultados, "tabelas")
dir_intersecoes <- file.path(dir_resultados, "intersecoes_nonb")
for (d in c(dir_resultados, dir_plots, dir_tabelas, dir_intersecoes)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)}

# Genoma
fa_path <- "./Diario-de-bordo/dados/Hsalinarum.fa"
if (file.exists(fa_path)) {
  fa <- readDNAStringSet(fa_path)
  cat(sprintf("Genoma carregado: %d sequências\n", length(fa)))
} else {
  warning("Arquivo FASTA não encontrado — análise de GC será pulada.")
  fa <- NULL}


# Seleção e classificaão de arquivos
arquivos    <- list.files(data_dir, full.names = TRUE)
arq_gff     <- arquivos[str_detect(arquivos, regex("Hsalinarum.*gff", ignore_case = TRUE))]
arq_tss     <- "./Diario-de-bordo/dados/Hsalinarum_TSS.gtf"
PATTERN_NONB <- "z-dna|nonb|triplex|R-loop|short_tandem|A-phased|G-quadruplex|Cruciform"
arquivos_nonb <- arquivos[
  str_detect(arquivos, regex("\\.bed$",    ignore_case = TRUE)) &
    str_detect(arquivos, regex(PATTERN_NONB, ignore_case = TRUE))]
if (length(arquivos_nonb) == 0) {
  warning("Nenhum arquivo BED non-B encontrado")
  print(arquivos[str_detect(arquivos, regex("\\.bed$", ignore_case = TRUE))])
} else {
  cat("\nArquivos BED encontrados:\n"); print(basename(arquivos_nonb))}

if (length(arq_gff) == 0) warning("Arquivo GFF não encontrado.")
if (!file.exists(arq_tss))  warning("Arquivo TSS não encontrado — algumas análises serão puladas.")


# Classificar por tipo de estrutura non-B
classificar_nonb <- function(arquivo) {
  n <- tolower(basename(arquivo))
  if (str_detect(n, "z-dna|zdna"))                return("Z-DNA")
  if (str_detect(n, "g-quadruplex|g4"))           return("G-quadruplex")
  if (str_detect(n, "triplex"))                   return("Triplex")
  if (str_detect(n, "r-loop|rloop"))              return("R-loop")
  if (str_detect(n, "short_tandem|tandem"))       return("Short_tandem")
  if (str_detect(n, "a-phased|aphased"))          return("A-phased")
  if (str_detect(n, "cruciform"))                 return("Cruciform")
  return("Outros")}

nonb_classificado <- data.frame(
  arquivo      = arquivos_nonb,
  nome_arquivo = basename(arquivos_nonb),
  classe       = sapply(arquivos_nonb, classificar_nonb),
  stringsAsFactors = FALSE)

write.table(nonb_classificado,
            file.path(dir_tabelas, "nonb_classificacao.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResumo por classe:\n"); print(table(nonb_classificado$classe))


# FUNÇÕES 
## Carrega qualquer arquivo BED/GFF como GRanges
carregar_granges <- function(arquivo) {
  if (!file.exists(arquivo)) {
    warning(paste("Arquivo não encontrado:", arquivo))
    return(GRanges())}
  ext <- tools::file_ext(arquivo)
  tryCatch({
    gr <- if (ext %in% c("bed", "bedgraph")) {
      import(arquivo, format = "BED")
    } else if (ext %in% c("gff", "gff3", "gtf")) {
      import(arquivo, format = "GFF")
    } else {
      warning(paste("Formato não suportado:", ext))
      return(GRanges())}
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr[width(gr) > 0]          # remove intervalos de largura zero
  }, error = function(e) {
    warning(paste("Erro ao carregar", arquivo, ":", e$message))
    GRanges()})}

## Salva um gráfico ggplot com configurações padrão
salvar_plot <- function(p, nome, w = 12, h = 7) {
  ggsave(
    filename = file.path(dir_plots, nome),
    plot     = p,
    width    = w, height = h, dpi = 300, bg = "white")
  cat("  Gráfico salvo:", nome, "\n")}


#  CARREGAMENTO DOS NON-B 
carregar_com_metadados <- function(i, df) {
  gr    <- carregar_granges(df$arquivo[i])
  nome  <- gsub("\\.(bed|gff|gff3)$", "", df$nome_arquivo[i])
  if (length(gr) > 0) {
    mcols(gr)$origem     <- df$nome_arquivo[i]
    mcols(gr)$classe     <- df$classe[i]
    mcols(gr)$arquivo_id <- i
    return(list(nome = nome, gr = gr, status = "OK", n = length(gr)))}
  list(nome = NULL, gr = NULL, status = "VAZIO", n = 0)}

resultados_carregamento <- foreach(
  i = seq_len(nrow(nonb_classificado)),
  .packages = c("GenomicRanges", "rtracklayer")
) %dopar% {
  carregar_com_metadados(i, nonb_classificado)}

gr_list <- list()
for (res in resultados_carregamento) {
  if (!is.null(res$gr)) {
    gr_list[[res$nome]] <- res$gr
    cat(sprintf("  Carregado: %-50s — %d regiões\n", res$nome, res$n))}}
cat(sprintf("\nTotal de arquivos carregados: %d\n", length(gr_list)))


#  AGRUPAR E CALCULAR INTERSEÇÕES classe

## Organiza os GRanges por classe
arquivos_por_classe <- list()
for (nome in names(gr_list)) {
  classe <- unique(mcols(gr_list[[nome]])$classe)[1]
  arquivos_por_classe[[classe]][[nome]] <- gr_list[[nome]]}

cat("\nClasses encontradas:\n")
for (cl_nome in names(arquivos_por_classe)) {
  cat(sprintf("  %-15s: %d arquivo(s)\n", cl_nome,
              length(arquivos_por_classe[[cl_nome]])))}

## interseção progressiva entre réplicas/ferramentas da mesma classe
calcular_interseccao_classe <- function(classe, arquivos_classe) {
  n <- length(arquivos_classe)
  
  # Apenas um arquivo → não há o que interseccionar
  if (n == 1) {
    gr <- arquivos_classe[[1]]
    mcols(gr) <- NULL
    return(list(classe = classe, gr_intersect = gr,
                n_intersect = length(gr), status = "único arquivo",
                arquivos = names(arquivos_classe)))}
  
  # Simplificando
  simplificar <- function(g) GRanges(seqnames(g), ranges(g), strand(g))
  grs <- lapply(arquivos_classe, simplificar)
  
  resultado <- grs[[1]]
  for (i in 2:n) {
    hits <- findOverlaps(resultado, grs[[i]])
    if (length(hits) == 0) {
      return(list(classe = classe, gr_intersect = GRanges(),
                  n_intersect = 0, status = "vazio",
                  arquivos = names(arquivos_classe)))}
    resultado <- pintersect(resultado[queryHits(hits)], grs[[i]][subjectHits(hits)])
    if (length(resultado) == 0) break}
  
  resultado <- reduce(resultado)
  list(classe = classe, gr_intersect = resultado,
       n_intersect = length(resultado), status = "interseção",
       arquivos = names(arquivos_classe))}

resultados_interseccao <- foreach(
  classe = names(arquivos_por_classe),
  .packages = "GenomicRanges"
) %dopar% {
  calcular_interseccao_classe(classe, arquivos_por_classe[[classe]])}

## resultados das interseções
nonb_por_classe    <- list()
resumo_intersecoes <- data.frame()
for (res in resultados_interseccao) {
  cat(sprintf("\n[%s] Arquivos: %s\n", res$classe,
              paste(res$arquivos, collapse = ", ")))
  
  if (res$n_intersect > 0) {
    gr_out <- res$gr_intersect
    mcols(gr_out)$classe            <- res$classe
    mcols(gr_out)$n_arquivos_origem <- length(res$arquivos)
    mcols(gr_out)$name <- paste0(
      gsub("-", "_", res$classe), "_",
      sprintf("%04d", seq_len(res$n_intersect)))
    nonb_por_classe[[res$classe]] <- res$gr_intersect
    
    nome_bed <- file.path(dir_intersecoes,
                          sprintf("nonB_%s_interseccao.bed",
                                  gsub("-", "_", res$classe)))
    export(gr_out, nome_bed, format = "BED")
    cat(sprintf("  BED salvo: %s (%d regiões)\n",
                basename(nome_bed), res$n_intersect))
    
    resumo_intersecoes <- rbind(resumo_intersecoes, data.frame(
      Classe             = res$classe,
      Regioes_Interseccao = res$n_intersect,
      Arquivos_Origem    = length(res$arquivos),
      Arquivos           = paste(res$arquivos, collapse = "; "),
      stringsAsFactors   = FALSE))
  } else {
    cat("  Interseção resultou em conjunto vazio.\n")}}

if (nrow(resumo_intersecoes) > 0) {
  write.table(resumo_intersecoes,
              file.path(dir_intersecoes, "resumo_intersecoes_nonB.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nResumo das interseções:\n"); print(resumo_intersecoes)
  cat(sprintf("Total de classes com interseção: %d\n", nrow(resumo_intersecoes)))
  cat(sprintf("Total de regiões em interseção : %d\n",
              sum(resumo_intersecoes$Regioes_Interseccao)))}


# GFF E TSS
genes    <- GRanges()
tss_P    <- GRanges()
promotores <- GRanges()
if (length(arq_gff) > 0) {
  genes_all <- carregar_granges(arq_gff[1])
  if (length(genes_all) > 0 && "type" %in% colnames(mcols(genes_all))) {
    genes <- genes_all[genes_all$type == "gene"]
    genes$region <- "gene"
    cat(sprintf("Genes carregados: %d\n", length(genes)))}}

if (file.exists(arq_tss)) {
  tss_raw <- carregar_granges(arq_tss)
  if (length(tss_raw) > 0 && "type" %in% colnames(mcols(tss_raw))) {
    tss_raw <- tss_raw[tss_raw$type %in% c("TSS", "transcription_start_site")]
    if ("tss_class" %in% colnames(mcols(tss_raw))) {
      tss_P <- tss_raw[tss_raw$tss_class == "P"]
      cat(sprintf("TSS classe P: %d\n", length(tss_P)))
    } else {
      warning("Coluna 'tss_class' não encontrada — usando todos os TSS.")
      tss_P <- tss_raw}
    
    # Promotores: 250 bp upstream do TSS (conforme Caten et al.)
    promotores         <- promoters(tss_P, upstream = 250, downstream = 0)
    promotores$region  <- "promoter"
    cat(sprintf("Promotores gerados: %d\n", length(promotores)))
    promotores_df <- as.data.frame(promotores) %>%
      transmute(seqname = seqnames, start, end, strand, width,
                tss_id  = names(promotores))
    write.csv(promotores_df, file.path(dir_tabelas, "promotores_250bp.csv"), row.names = FALSE)}}


# GENES COM E SEM 5' UTR 
## Para cada TSS, encontra o gene mais próximo na mesma fita
encontrar_gene_por_tss <- function(tss_gr, genes_gr) {
  result <- vector("list", length(tss_gr))
  for (i in seq_along(tss_gr)) {
    tp     <- tss_gr[i]
    seqn   <- as.character(seqnames(tp))
    strnd  <- as.character(strand(tp))
    cand   <- genes_gr[seqnames(genes_gr) == seqn & strand(genes_gr) == strnd]
    if (length(cand) == 0) next
    
    gene_closest <- if (strnd == "+") {
      cand_dir <- cand[end(cand) > start(tp)]
      if (length(cand_dir) == 0) next
      cand_dir[which.min(start(cand_dir))]
    } else {
      cand_esq <- cand[start(cand) < end(tp)]
      if (length(cand_esq) == 0) next
      cand_esq[which.max(end(cand_esq))]}
    
    result[[i]] <- data.frame(
      tss_id    = i,
      gene_id   = mcols(gene_closest)$ID,
      gene_name = ifelse("Name" %in% colnames(mcols(gene_closest)),
                         mcols(gene_closest)$Name, NA_character_),
      seqname   = seqn,
      tss_start = start(tp),  tss_end = end(tp),
      gene_start = start(gene_closest), gene_end = end(gene_closest),
      strand    = strnd,
      stringsAsFactors = FALSE)}
  bind_rows(result)}

if (length(tss_P) > 0 && length(genes) > 0) {
  tss_gene_assoc <- encontrar_gene_por_tss(tss_P, genes) %>%
    mutate(distancia_tss_gene = ifelse(
      strand == "+",
      gene_start - tss_end,
      tss_start  - gene_end))
  
  genes_com_utr <- tss_gene_assoc %>% filter(distancia_tss_gene >  0)
  genes_sem_utr <- tss_gene_assoc %>% filter(distancia_tss_gene <= 0)
  
  cat(sprintf("\nGenes com 5' UTR: %d | sem 5' UTR: %d\n",
              nrow(genes_com_utr), nrow(genes_sem_utr)))
  
  tabela_genes_utr <- tss_gene_assoc %>%
    mutate(tem_5utr = distancia_tss_gene > 0) %>%
    select(gene_id, gene_name, seqname, strand,
           tss_start, gene_start, gene_end, distancia_tss_gene, tem_5utr)
  
  write.csv(tabela_genes_utr,file.path(dir_tabelas, "genes_status_5utr.csv"),row.names = FALSE)
  
  stats_utr <- data.frame(
    Categoria  = c("Com 5' UTR", "Sem 5' UTR"),
    Quantidade = c(nrow(genes_com_utr), nrow(genes_sem_utr))
  ) %>% mutate(Percentual = round(100 * Quantidade / sum(Quantidade), 1),
               Label      = paste0(Quantidade, "\n(", Percentual, "%)"))
  
  write.csv( stats_utr %>% select(-Label), file.path(dir_tabelas, "estatisticas_5utr.csv"), row.names = FALSE)
  
  # genes com/sem 5' UTR 
  # Donut
  p_utr_pie <- ggplot(stats_utr,
                      aes(x = 2, y = Quantidade, fill = Categoria)) +
    geom_col(width = 1, color = "white", linewidth = 0.8) +
    geom_text(aes(label = Label),
              position = position_stack(vjust = 0.5),
              size = 4.5, fontface = "bold", color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = c("Com 5' UTR" = "#1D3557",
                                 "Sem 5' UTR"  = "#E63946")) +
    labs(
      title   = "Presença de 5' UTR nos genes de H. salinarum",
      caption = "Baseado em TSS de classe P; distância TSS → início do gene"
    ) + theme_void() + theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 13,
                                  margin = margin(b = 8)),
      plot.caption = element_text(hjust = 0.5, color = "grey50", size = 9),
      legend.title = element_blank(),
      legend.position = "bottom")
  salvar_plot(p_utr_pie, "genes_com_sem_5utr.png", w = 7, h = 7)}


# 5' UTRs CUSTOMIZADAS 
criar_utr5_custom <- function(genes_gr, tss_gr) {
  if (length(tss_gr) == 0) return(GRanges())
  utr_list <- vector("list", length(tss_gr))
  
  for (i in seq_along(tss_gr)) {
    tp    <- tss_gr[i]
    seqn  <- as.character(seqnames(tp))
    strnd <- as.character(strand(tp))
    cand  <- genes_gr[seqnames(genes_gr) == seqn & strand(genes_gr) == strnd]
    if (length(cand) == 0) next
    
    if (strnd == "+") {
      cand_dir <- cand[end(cand) > start(tp)]
      if (length(cand_dir) == 0) next
      gc     <- cand_dir[which.min(start(cand_dir))]
      u_s    <- end(tp);           u_e <- start(gc) - 1
      gene_id <- mcols(gc)$ID
    } else {
      cand_esq <- cand[start(cand) < end(tp)]
      if (length(cand_esq) == 0) next
      gc     <- cand_esq[which.max(end(cand_esq))]
      u_s    <- end(gc) + 1;      u_e <- start(tp)
      gene_id <- mcols(gc)$ID}
    
    if (u_e >= u_s) {
      utr_list[[i]] <- GRanges(seqn, IRanges(u_s, u_e),
                               strand  = strnd,
                               region  = "5_UTR",
                               gene_id = gene_id)}}
  
  valid <- Filter(Negate(is.null), utr_list)
  if (length(valid) == 0) return(GRanges())
  trim(unlist(GRangesList(valid)))}

utr5_custom <- GRanges()
if (length(tss_P) > 0 && length(genes) > 0) {
  utr5_custom <- criar_utr5_custom(genes, tss_P)
  cat(sprintf("5' UTRs geradas: %d\n", length(utr5_custom)))
  
  if (length(utr5_custom) > 0) {
    utr5_table <- as.data.frame(utr5_custom) %>%
      transmute(seqname = seqnames, start, end, strand, width,
                gene_id, region = "5_UTR") %>%
      arrange(seqname, start)
    
    write.csv(utr5_table,file.path(dir_resultados, "5_UTR_custom.csv"),row.names = FALSE)
    
    stats_utr_tam <- utr5_table %>%
      summarise(media   = mean(width),  mediana = median(width),
                min     = min(width),   max     = max(width),
                q1      = quantile(width, 0.25),
                q3      = quantile(width, 0.75))
    
    write.csv(stats_utr_tam, file.path(dir_tabelas, "estatisticas_tamanho_5utr.csv"),row.names = FALSE)
    
    # distribuição do tamanho das 5' UTRs ─────────────────────────
    p_utr_hist <- ggplot(utr5_table, aes(x = width)) +
      geom_histogram(bins = 35, fill = "#1D3557", alpha = 0.8, color = "white") +
      geom_vline(xintercept = stats_utr_tam$media,   color = "#E63946",
                 linetype = "dashed",  linewidth = 1) +
      geom_vline(xintercept = stats_utr_tam$mediana, color = "#2A9D8F",
                 linetype = "dotted",  linewidth = 1) +
      annotate("label",
               x     = stats_utr_tam$media * 1.05,
               y     = Inf, vjust = 1.5,
               label = sprintf("Média = %.1f bp", stats_utr_tam$media),
               color = "#E63946", fill = "white", size = 3.5,
               label.size = 0.3, fontface = "bold") +
      annotate("label",
               x     = stats_utr_tam$mediana * 1.05,
               y     = Inf, vjust = 3.2,
               label = sprintf("Mediana = %.1f bp", stats_utr_tam$mediana),
               color = "#2A9D8F", fill = "white", size = 3.5,
               label.size = 0.3, fontface = "bold") +
      labs(
        title   = "Distribuição do tamanho das 5' UTRs",
        subtitle = "Estimadas a partir de TSS de classe P",
        x       = "Tamanho (bp)",
        y       = "Frequência",
        caption = "Linhas: média (vermelho tracejado) e mediana (verde pontilhado)")
    salvar_plot(p_utr_hist, "histograma_tamanho_5utr.png", w = 10, h = 6)
    
    # tamanho por fita ─────────────────────────────────────────────
    p_utr_box <- ggplot(utr5_table, aes(x = strand, y = width, fill = strand)) +
      geom_violin(alpha = 0.6, trim = TRUE, color = "grey40") +
      geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA) +
      geom_jitter(width = 0.08, alpha = 0.25, size = 1.2, color = "grey30") +
      scale_fill_manual(values = c("+" = "#1D3557", "-" = "#E63946")) +
      labs(
        title   = "Tamanho das 5' UTRs por fita",
        x       = "Fita",
        y       = "Tamanho (bp)",
        caption = "Violino + boxplot; pontos individuais em cinza"
      ) + theme(legend.position = "none")
    salvar_plot(p_utr_box, "boxplot_5utr_por_strand.png", w = 6, h = 6)}}


# BINS DO CORPO GÊNICO ──────────────────────────────────────────────────
# Divide cada gene em N bins de tamanho proporcional para análise posicional
bins_por_gene <- function(gr_gene, n_bins = 10) {
  gl    <- width(gr_gene)
  pts   <- floor(seq(0, gl, length.out = n_bins + 1))
  s_pos <- start(gr_gene) - 1 + pts[-(n_bins + 1)] + 1
  e_pos <- start(gr_gene) - 1 + pts[-1]
  bins  <- GRanges(
    seqnames = seqnames(gr_gene),
    ranges   = IRanges(start = s_pos, end = e_pos),
    strand   = strand(gr_gene))
  bins$region <- paste0("gene_bin_", seq_len(n_bins))
  bins[start(bins) <= end(bins)]}

gene_bins <- GRanges()
if (length(genes) > 0) {
  gene_bins_list <- foreach(
    gene = split(genes, seq_along(genes)),
    .packages = "GenomicRanges",
    .combine  = c
  ) %dopar% { bins_por_gene(gene) }
  gene_bins <- sort(gene_bins_list)
  cat(sprintf("Bins gerados: %d\n", length(gene_bins)))}


# CONCATENAR TODAS AS REGIÕES 
regioes_list <- list()
if (length(promotores)  > 0) regioes_list$promotores <- promotores
if (length(utr5_custom) > 0) regioes_list$utr5       <- utr5_custom
if (length(gene_bins)   > 0) regioes_list$bins       <- gene_bins
todas_regioes <- if (length(regioes_list) > 0) {
  trim(Reduce(c, regioes_list))
} else { GRanges()}
cat(sprintf("Total de regiões para análise: %d\n", length(todas_regioes)))

# SOBREPOSIÇÕES NON-B vs REGIÕES GÊNICAS 
##  Distribuição posicional detalhada (bins individuais)
if (length(nonb_por_classe) > 0 && length(todas_regioes) > 0) {
  overlap_list <- list()
  for (classe in names(nonb_por_classe)) {
    ov <- findOverlaps(todas_regioes, nonb_por_classe[[classe]],
                       ignore.strand = TRUE)
    if (length(ov) > 0) {
      overlap_list[[classe]] <- tibble(
        classe   = classe,
        region   = as.character(todas_regioes$region)[queryHits(ov)],
        nb_index = subjectHits(ov))}}
  
  overlap_df <- bind_rows(overlap_list) %>%
    mutate(
      region = recode(region, promoter = "Promotor", `5_UTR` = "5' UTR"),
      region = ifelse(str_detect(region, "^gene_bin_"),
                      paste0(as.numeric(str_extract(region, "\\d+")) * 10, "%"),
                      region))
  
  write.csv(overlap_df,file.path(dir_resultados, "sobreposicoes_nonB_regioes.csv"),row.names = FALSE)
  
  ordem_det <- c("Promotor", "5' UTR", paste0(seq(10, 100, 10), "%"))
  overlap_plot_det <- overlap_df %>%
    mutate(region = factor(region, levels = ordem_det)) %>%
    group_by(classe, region) %>%
    summarise(frequencia = n(), .groups = "drop") %>%
    filter(!is.na(region))
  
  ## distribuição ao longo do gene
  for (cl_nome in unique(overlap_plot_det$classe)) {
    df_cl <- overlap_plot_det %>% filter(classe == cl_nome)
    if (nrow(df_cl) == 0) next
    cor_cl <- CORES_NONB[cl_nome] %||% "#ADB5BD"
    p_cl <- ggplot(df_cl, aes(x = region, y = frequencia)) +
      geom_col(fill = cor_cl, color = "white", width = 0.75, alpha = 0.9) +
      geom_text(aes(label = frequencia),
                vjust = -0.4, size = 3.2, fontface = "bold",
                color = "grey20") +
      labs(
        title    = sprintf("Distribuição de %s por região gênica", cl_nome),
        subtitle = sprintf("Total de sobreposições: %d", sum(df_cl$frequencia)),
        x        = "Região (Promotor → 3' do gene)",
        y        = "Sobreposições non-B",
        caption  = "Promotor = 250 bp upstream do TSS (Caten et al.); corpo gênico dividido em decis"
      ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    salvar_plot(p_cl,sprintf("intersec_%s_detalhado.png", gsub("[ -]", "_", cl_nome)),w = 13, h = 7)}
  
  ## todas as classes juntas
  p_stack <- ggplot(overlap_plot_det,
                    aes(x = region, y = frequencia, fill = classe)) +
    geom_col(position = "stack", color = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = CORES_NONB[intersect(names(CORES_NONB),
                                    unique(overlap_plot_det$classe))],
      name   = "Classe"
    ) + labs(
      title    = "Distribuição de estruturas non-B por região gênica",
      subtitle = "Todas as classes sobrepostas",
      x        = "Região",
      y        = "Total de sobreposições",
      caption  = "Promotor = 250 bp upstream; corpo gênico em decis de 10%"
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  salvar_plot(p_stack, "stacked_nonB_regioes_detalhado.png", w = 14, h = 8)
  resumo_det <- overlap_plot_det %>%
    pivot_wider(names_from = region, values_from = frequencia,
                values_fill = 0)
  write.csv(resumo_det,file.path(dir_tabelas, "resumo_nonB_por_regiao.csv"), row.names = FALSE)}


## Distribuição agrupada (Promotor / 5' UTR / Corpo do gene)
if (length(nonb_por_classe) > 0 && length(todas_regioes) > 0) {
  overlap_agrup <- overlap_df %>%
    mutate(regiao_agrupada = case_when(
      region == "Promotor" ~ "Promotor",
      region == "5' UTR"  ~ "5' UTR",
      str_detect(region, "%") ~ "Corpo do gene",
      TRUE ~ region))
  
  write.csv(overlap_agrup, file.path(dir_resultados, "sobreposicoes_nonB_regioes_agrupado.csv"),row.names = FALSE)
  
  ordem_agrup <- c("Promotor", "5' UTR", "Corpo do gene")
  overlap_plot_agrup <- overlap_agrup %>%
    group_by(classe, regiao_agrupada) %>%
    summarise(frequencia = n(), .groups = "drop") %>%
    filter(!is.na(regiao_agrupada)) %>%
    mutate(regiao_agrupada = factor(regiao_agrupada, levels = ordem_agrup))
  
  ## Gráfico por classe (agrupado)
  for (cl_nome in unique(overlap_plot_agrup$classe)) {
    df_cl <- overlap_plot_agrup %>% filter(classe == cl_nome)
    if (nrow(df_cl) == 0) next
    cor_cl <- CORES_NONB[cl_nome] %||% "#ADB5BD"
    p_ag <- ggplot(df_cl, aes(x = regiao_agrupada, y = frequencia)) +
      geom_col(fill = cor_cl, color = "white", width = 0.6, alpha = 0.9) +
      geom_text(aes(label = frequencia),
                vjust = -0.4, size = 4, fontface = "bold") +
      labs(
        title    = sprintf("Distribuição de %s por região gênica", cl_nome),
        subtitle = sprintf("Total: %d sobreposições", sum(df_cl$frequencia)),
        x        = NULL,
        y        = "Sobreposições non-B",
        caption  = "Corpo do gene agrupa todos os decis (10%–100%)" )
    
    salvar_plot(p_ag, sprintf("intersec_%s_agrupado.png", gsub("[ -]", "_", cl_nome)),w = 7, h = 6)}
  
  ## Gráfico empilhado agrupado
  p_stack_ag <- ggplot(overlap_plot_agrup,
                       aes(x = regiao_agrupada, y = frequencia, fill = classe)) +
    geom_col(position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = CORES_NONB[intersect(names(CORES_NONB),unique(overlap_plot_agrup$classe))],
      name   = "Classe"
    ) +
    labs(
      title    = "Distribuição de estruturas non-B por região gênica (agrupado)",
      subtitle = "Promotor · 5' UTR · Corpo do gene",
      x        = NULL,
      y        = "Total de sobreposições",
      caption  = "Corpo do gene = junção dos decis 10%–100% do CDS")
  
  salvar_plot(p_stack_ag, "stacked_nonB_regioes_agrupado.png", w = 9, h = 7)
  
  resumo_agrup <- overlap_plot_agrup %>%
    pivot_wider(names_from = regiao_agrupada, values_from = frequencia,
                values_fill = 0)
  write.csv(resumo_agrup, file.path(dir_tabelas, "resumo_nonB_por_regiao_agrupado.csv"),row.names = FALSE)}


# CONTAGEM REAL POR ARQUIVO/FERRAMENTA ──────────────────────────────────

contagem_real <- data.frame(
  Arquivo       = names(gr_list),
  Classe        = sapply(gr_list, function(x) unique(mcols(x)$classe)[1]),
  Num_Predicoes = sapply(gr_list, length),
  stringsAsFactors = FALSE
) %>% mutate(
  ferramenta_exibicao = case_when(
    grepl("nonb",     Arquivo, ignore.case = TRUE) ~ "Non-B DB v2.0",
    grepl("einverted",Arquivo, ignore.case = TRUE) ~ "einverted",
    grepl("gboost",   Arquivo, ignore.case = TRUE) ~ "G4Boost",
    grepl("qmrlfs",   Arquivo, ignore.case = TRUE) ~ "QmRLFS-finder",
    grepl("triplex",  Arquivo, ignore.case = TRUE) ~ "triplex (Bioconductor)",
    grepl("zseek",    Arquivo, ignore.case = TRUE) ~ "ZSeeker",
    TRUE ~ Arquivo
  ), Classe = factor(Classe,
                     levels = c("A-phased", "Short_tandem", "Cruciform",
                                "G-quadruplex", "R-loop", "Triplex", "Z-DNA")))

write.table(contagem_real,file.path(dir_tabelas, "contagem_real_predicoes.tsv"),sep = "\t", row.names = FALSE, quote = FALSE)

#  INTERSEÇÃO E COMPARAÇÃO
if (nrow(resumo_intersecoes) > 0) {
  
  ## 15a. Regiões em interseção por classe
  p_intersec <- ggplot(
    resumo_intersecoes,
    aes(x = reorder(Classe, Regioes_Interseccao),
        y = Regioes_Interseccao, fill = Classe)
  ) +
    geom_col(alpha = 0.85, width = 0.65, color = "white") +
    geom_text(aes(label = format(Regioes_Interseccao, big.mark = ",")),
              hjust = -0.15, size = 3.8, fontface = "bold", color = "grey20") +
    coord_flip() +
    scale_fill_manual(values = CORES_NONB) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = "Regiões consenso por classe non-B",
      subtitle = "Interseção entre todas as ferramentas da mesma classe",
      x        = NULL,
      y        = "Regiões em interseção",
      caption  = "Interseção progressiva; apenas regiões presentes em TODAS as ferramentas"
    ) + theme(legend.position = "none")
  salvar_plot(p_intersec, "interseccoes_por_classe.png", w = 10, h = 6)
  
  ## Total vs. interseção
  totais_classe <- contagem_real %>%
    group_by(Classe) %>%
    summarise(Total_Regioes = sum(Num_Predicoes), .groups = "drop") %>%
    mutate(Classe = as.character(Classe))
  
  comparacao <- merge(totais_classe, resumo_intersecoes, by = "Classe",
                      all = TRUE) %>%
    replace_na(list(Regioes_Interseccao = 0)) %>%
    mutate(Percentual_Interseccao =
             round(100 * Regioes_Interseccao / Total_Regioes, 1))
  
  write.table(comparacao,file.path(dir_tabelas, "comparacao_total_vs_interseccao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  comp_long <- comparacao %>%
    pivot_longer(c(Total_Regioes, Regioes_Interseccao),
                 names_to = "Tipo", values_to = "Contagem") %>%
    mutate(Tipo = recode(Tipo,
                         Total_Regioes      = "Total predito",
                         Regioes_Interseccao = "Consenso (interseção)"))
  
  p_comp <- ggplot(comp_long,
                   aes(x = reorder(Classe, -Contagem),
                       y = Contagem, fill = Tipo)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.72,
             color = "white") +
    geom_text(aes(label = format(Contagem, big.mark = ",")),
              position = position_dodge(0.72),
              vjust = -0.4, size = 3, color = "grey20") +
    scale_fill_manual(
      values = c("Total predito" = "#1D3557",
                 "Consenso (interseção)" = "#E63946"),
      name   = NULL
    ) + scale_y_continuous(labels = label_number(big.mark = ","),expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = "Total de predições vs. regiões consenso",
      subtitle = "Por classe de DNA não-B",
      x        = NULL,
      y        = "Número de regiões",
      caption  = "Consenso = regiões presentes em todas as ferramentas para a classe"
    ) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  salvar_plot(p_comp, "comparacao_total_vs_interseccao.png", w = 12, h = 7)}


# PREDIÇÕES POR FERRAMENTA ─────────────────────────────────────────────
p_ferramentas <- ggplot(
  contagem_real,
  aes(x = Classe, y = Num_Predicoes, fill = ferramenta_exibicao)
) +
  geom_col(
    position = position_dodge2(preserve = "single"),
    width = 0.8, color = "white", alpha = 0.88
  ) +
  geom_text(
    aes(label = format(Num_Predicoes, big.mark = ",")),
    position = position_dodge2(width = 0.8, preserve = "single"),
    vjust = -0.45, size = 2.6, color = "grey25"
  ) +
  scale_y_log10(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0.02, 0.15))
  ) + scale_fill_brewer(palette = "Set1", na.value = "grey70", name = "Ferramenta") +
  labs(
    title    = "Predições de DNA não-B por ferramenta e classe",
    subtitle = "Escala logarítmica — cada barra representa uma ferramenta",
    x        = "Classe non-B",
    y        = "Regiões preditas (escala log₁₀)",
    caption  = "Ferramentas usadas: Non-B DB v2.0, einverted, G4Boost, QmRLFS-finder, triplex, ZSeeker"
  ) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1),
    panel.grid.major.x = element_blank()
  ) salvar_plot(p_ferramentas, "predicoes_por_ferramenta.png", w = 13, h = 7)

# GRanges POR CLASSE 
limpar_granges <- function(gr) {
  cols_manter <- intersect(c("classe", "origem", "name", "score"),
                           colnames(mcols(gr)))
  if (length(cols_manter) > 0) {
    mcols(gr) <- mcols(gr)[, cols_manter, drop = FALSE]
  } else {
    mcols(gr) <- NULL}
  gr}
gr_por_classe <- list()
for (classe in unique(as.character(contagem_real$Classe))) {
  arqs   <- contagem_real$Arquivo[contagem_real$Classe == classe]
  grs    <- lapply(gr_list[arqs], limpar_granges)
  tryCatch({
    gr_por_classe[[classe]] <- unlist(GRangesList(grs))
  }, error = function(e) {
    warning(paste("Erro ao combinar classe", classe, ":", e$message,
                  "— usando método sequencial."))
    gr_por_classe[[classe]] <<- Reduce(c, grs)})
  cat(sprintf("  %-15s: %d elementos consolidados\n",
              classe, length(gr_por_classe[[classe]])))}


# CONTEÚDO GC
calcular_gc_real <- function(gr, fa_ref, max_amostras = 10000L) {
  if (is.null(fa_ref) || length(gr) == 0) return(NULL)
  tryCatch({
    seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr <- gr[seqnames(gr) %in% intersect(seqlevels(gr), names(fa_ref))]
    gr <- trim(gr)
    if (length(gr) == 0) return(NULL)
    if (length(gr) > max_amostras) gr <- gr[sample(length(gr), max_amostras)]
    
    seqs <- getSeq(fa_ref, gr)
    gc_vals <- sapply(as.character(seqs), function(s) {
      if (nchar(s) == 0) return(NA_real_)
      af <- alphabetFrequency(DNAString(s), baseOnly = TRUE)
      tot <- sum(af[c("A", "C", "G", "T")])
      if (tot > 0) (af["C"] + af["G"]) / tot else NA_real_})
    gc_vals[!is.na(gc_vals)]
  }, error = function(e) {
    warning(paste("Erro no cálculo de GC:", e$message))
    NULL })}

if (!is.null(fa) && length(gr_por_classe) > 0) {
  gc_violin_data <- bind_rows(lapply(names(gr_por_classe), function(cl_nome) {
    cat(sprintf("    GC — %-15s: ", cl_nome))
    vals <- calcular_gc_real(gr_por_classe[[cl_nome]], fa)
    if (!is.null(vals) && length(vals) > 0) {
      cat(sprintf("%d valores\n", length(vals)))
      data.frame(Classe = cl_nome, GC = vals, stringsAsFactors = FALSE)
    } else {
      cat("sem dados\n"); NULL }}))
  
  if (nrow(gc_violin_data) > 0) {
    stats_gc <- gc_violin_data %>%
      group_by(Classe) %>%
      summarise(
        Mediana = median(GC, na.rm = TRUE),
        Media   = mean(GC,   na.rm = TRUE),
        Q1      = quantile(GC, 0.25, na.rm = TRUE),
        Q3      = quantile(GC, 0.75, na.rm = TRUE),
        n       = n(), .groups = "drop"
      ) %>% arrange(Mediana)
    
    write.csv(stats_gc, file.path(dir_tabelas, "estatisticas_gc_por_classe.csv"), row.names = FALSE)
    
    gc_violin_data$Classe <- factor(gc_violin_data$Classe, levels = stats_gc$Classe)
    
    # violino GC por classe 
    p_gc <- ggplot(gc_violin_data, aes(x = Classe, y = GC, fill = Classe)) +
      geom_violin(alpha = 0.7, trim = TRUE, scale = "width", color = "grey40", linewidth = 0.4) +
      geom_boxplot(width = 0.12, alpha = 0.9, outlier.shape = NA, color = "grey20") +
      stat_summary(fun = median, geom = "point", shape = 18, size = 3.5, color = "white") +
      geom_text(data = stats_gc,
                aes(x = Classe, y = Q3 + 0.05,
                    label = sprintf("%.1f%%", Mediana * 100)),
                size = 3.5, fontface = "bold", color = "grey20",
                inherit.aes = FALSE) +
      scale_fill_manual(values = CORES_NONB) +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        breaks = seq(0, 1, 0.2),
        limits = c(0, 1.1)
      ) + labs(
        title    = "Conteúdo de GC por classe de DNA não-B",
        subtitle = sprintf("Genoma: %s", basename(fa_path)),
        x        = NULL,
        y        = "Conteúdo GC (%)",
        caption  = "Violino: distribuição; caixa: IQR; losango: mediana. Texto = mediana (%)."
      ) + theme(
        axis.text.x  = element_text(angle = 35, hjust = 1, face = "bold"),
        legend.position = "none")
    
    salvar_plot(p_gc, "violino_gc_por_classe.png", w = max(12, length(unique(gc_violin_data$Classe)) * 1.2), h = 7)
    print(stats_gc)
  } else {
    cat("Nenhum dado de GC disponível.\n")}}


# TAMANHO DOS INTERVALOS
tamanhos_data <- bind_rows(lapply(names(gr_list), function(nm) {
  data.frame(
    Classe  = unique(mcols(gr_list[[nm]])$classe)[1],
    Arquivo = nm,
    Tamanho = width(gr_list[[nm]]),
    stringsAsFactors = FALSE )}))

stats_tamanhos <- tamanhos_data %>%
  group_by(Classe) %>%
  summarise(
    Media   = mean(Tamanho),
    Mediana = median(Tamanho),
    Q1      = quantile(Tamanho, 0.25),
    Q3      = quantile(Tamanho, 0.75),
    Min     = min(Tamanho),
    Max     = max(Tamanho),
    n       = n(), .groups = "drop")
write.csv(stats_tamanhos,file.path(dir_tabelas, "estatisticas_tamanhos.csv"),row.names = FALSE)

## Densidade (escala log)
p_dens <- ggplot(tamanhos_data, aes(x = Tamanho, fill = Classe, color = Classe)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_x_log10(labels = label_number(big.mark = ",")) +
  scale_fill_manual(values  = CORES_NONB, name = "Classe") +
  scale_color_manual(values = CORES_NONB, name = "Classe") +
  labs(
    title    = "Distribuição do tamanho das estruturas non-B",
    subtitle = "Escala logarítmica no eixo x",
    x        = "Tamanho (bp) — log₁₀",
    y        = "Densidade",
    caption  = "Cada curva representa todas as predições de uma classe (antes da interseção)"
  ) + facet_wrap(~Classe, scales = "free_y", ncol = 3) +
  theme(legend.position = "none",strip.background = element_rect(fill = "grey95", color = NA))
salvar_plot(p_dens, "densidade_tamanhos.png", w = 13, h = 9)

## Boxplot (escala log)
p_box_tam <- ggplot(tamanhos_data,
                    aes(x = Classe, y = Tamanho, fill = Classe)) +
  geom_violin(alpha = 0.5, trim = TRUE, color = "grey50") +
  geom_boxplot(width = 0.14, alpha = 0.9, outlier.shape = NA) +
  scale_y_log10(labels = label_number(big.mark = ",")) +
  scale_fill_manual(values = CORES_NONB) +
  labs(
    title    = "Tamanho das estruturas non-B por classe",
    subtitle = "Violino + boxplot; escala log₁₀",
    x        = NULL,
    y        = "Tamanho (bp) — log₁₀",
    caption  = "Inclui todas as predições (antes da interseção entre ferramentas)"
  ) + theme(
    axis.text.x     = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "none" )
salvar_plot(p_box_tam, "boxplot_tamanhos.png", w = 11, h = 7)

#HEATMAP DE DENSIDADE (INTERSEÇÕES)
if (length(nonb_por_classe) > 0) {
  calc_densidade_intersecao <- function(regioes, nonb_gr) {
    if (length(regioes) == 0 || length(nonb_gr) == 0)
      return(list(total = 0, kb = 0, densidade = 0))
    total <- 0L
    kb    <- 0
    for (i in seq_along(regioes)) {
      hits  <- findOverlaps(regioes[i], nonb_gr, ignore.strand = TRUE)
      total <- total + length(unique(subjectHits(hits)))
      kb    <- kb + width(regioes[i]) / 1000 }
    list(total = total, kb = round(kb, 2),
         densidade = round(total / max(kb, 1e-9), 4))}
  
  df_dens <- bind_rows(lapply(names(nonb_por_classe), function(cl_nome) {
    cat(sprintf("  Densidade — %s\n", cl_nome))
    nonb_gr <- nonb_por_classe[[cl_nome]]
    data.frame(
      Classe            = cl_nome,
      Promotor_dens     = calc_densidade_intersecao(promotores,  nonb_gr)$densidade,
      UTR5_dens         = calc_densidade_intersecao(utr5_custom, nonb_gr)$densidade,
      Gene_dens         = calc_densidade_intersecao(genes,       nonb_gr)$densidade,
      stringsAsFactors  = FALSE )}))
  
  write.csv(df_dens, file.path(dir_tabelas, "densidade_nonB_por_regiao.csv"), row.names = FALSE)
  mat_heat <- as.matrix(df_dens[, -1])
  rownames(mat_heat) <- df_dens$Classe
  colnames(mat_heat) <- c("Promotor", "5' UTR", "Corpo Gênico")
  
  pheatmap(
    mat_heat,
    main            = "Densidade de estruturas non-B (consenso) por região\n(regiões / kb)",
    color           = viridis(100, option = "D"),
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    display_numbers = TRUE,
    number_format   = "%.4f",
    number_color    = "white",
    fontsize         = 11,
    fontsize_number  = 9,
    border_color     = "white",
    filename        = file.path(dir_plots, "heatmap_densidade_intersecoes.png"),
    width           = 8,
    height          = max(5, nrow(mat_heat) * 0.8),
    bg              = "white"
  ) cat("  Heatmap salvo.\n") }

# SALVAR OBJETOS R 
cat("\nSalvando objetos RDS...\n")
saveRDS(gr_list,      file.path(dir_resultados, "nonb_gr_list_completa.rds"))
saveRDS(nonb_por_classe, file.path(dir_resultados, "nonb_por_classe_intersecao.rds"))
saveRDS(gr_por_classe, file.path(dir_resultados, "nonb_gr_por_classe.rds"))
if (length(utr5_custom) > 0)
  saveRDS(utr5_custom, file.path(dir_resultados, "utr5_custom.rds"))
if (length(promotores) > 0)
  saveRDS(promotores,  file.path(dir_resultados, "promotores.rds"))

# RELATÓRIO
stopCluster(cl)
cat("\n", strrep("=", 60), "\n", sep = "")
cat("  ANÁLISE CONCLUÍDA\n")
cat(strrep("=", 60), "\n\n", sep = "")
cat(sprintf("  Plots     : %d arquivos\n",
            length(list.files(dir_plots,       pattern = "\\.png$"))))
cat(sprintf("  Tabelas   : %d arquivos\n",
            length(list.files(dir_tabelas,     pattern = "\\.(tsv|csv)$"))))
cat(sprintf("  BED (inter): %d arquivos\n",
            length(list.files(dir_intersecoes, pattern = "\\.bed$"))))
cat(sprintf("  Objetos R  : %d arquivos\n",
            length(list.files(dir_resultados,  pattern = "\\.rds$"))))
cat("\nEstrutura non-B:\n")
for (cl_nome in names(gr_por_classe)) {
  n_tot    <- length(gr_por_classe[[cl_nome]])
  n_inter  <- ifelse(cl_nome %in% names(nonb_por_classe),
                     length(nonb_por_classe[[cl_nome]]), 0L)
  n_arqs   <- sum(contagem_real$Classe == cl_nome)
  cat(sprintf("  %-15s: %7s regiões totais | %6s consenso | %d ferramenta(s)\n",
              cl_nome,
              format(n_tot,   big.mark = ","),
              format(n_inter, big.mark = ","),
              n_arqs))}

if (length(utr5_custom) > 0) {
  cat(sprintf(
    "\n5' UTRs: %d | média %.1f bp | mediana %d bp | min %d | max %d bp\n",
    length(utr5_custom), mean(width(utr5_custom)),
    median(width(utr5_custom)),
    min(width(utr5_custom)), max(width(utr5_custom))))}

cat(sprintf("\nResultados em: %s\n", dir_resultados))