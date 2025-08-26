# tidyverse -> tabelas e listas
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(InteractionSet)
# BiocManager::install("tidyverse")
data_dir <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/sobreposição"
arquivos <- list.files(data_dir, full.names = TRUE)

# separar arq por categoria usando regex
arquivos_smap <- arquivos[str_detect(arquivos, regex("smap1", ignore_case = TRUE))]
arquivos_nonb <- arquivos[str_detect(arquivos, regex("z-dna|gquad|kim|nonb|triplex|r-loop|short_tandem_slipped|A-phased|G-quadruplex|Cruciform", ignore_case = TRUE)) & str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]
arquivos_genes <- arquivos[str_detect(arquivos, regex("Hsalinarum", ignore_case = TRUE)) & str_detect(arquivos, regex("\\.(bed|gff3?)$", ignore_case = TRUE))]
arquivos_expressao <- arquivos[str_detect(arquivos, regex("bedgraph", ignore_case = TRUE))]
arquivos_hic <- arquivos[str_detect(arquivos, regex("bedpe", ignore_case = TRUE))]

# ===================================================================== n mexer - ta funcionando
#tryCatch - evita erros e ajuda a importar o tipo
carregar_genomic <- function(caminho) {
  tryCatch({
    if (str_detect(caminho, regex("\\.gff3?$", ignore_case = TRUE))) {import.gff(caminho)}
    else if (str_detect(caminho, regex("\\.bed$", ignore_case = TRUE))) {import.bed(caminho)}
    else {NULL}}, error = function(e) {
      warning("Problema ", basename(caminho), ": ", e$message)
      NULL})}

# para Hi-C -> separar manualmente as âncoras (não achei jeito mais direto)
carregar_hic <- function(arquivo) {
  tryCatch({
    hic_df <- read_tsv(arquivo, col_names = FALSE)
    gr1 <- GRanges(seqnames = hic_df$X1, ranges = IRanges(hic_df$X2 + 1, hic_df$X3))
    gr2 <- GRanges(seqnames = hic_df$X4, ranges = IRanges(hic_df$X5 + 1, hic_df$X6))
    GInteractions(gr1, gr2)}, error = function(e) {
      message("Erro", basename(arquivo), ": ", e$message)
      return(NULL)})}

# =====================================================================
#  lista nomeada
dados <- list(
  smap = map(set_names(arquivos_smap, basename(arquivos_smap)), carregar_genomic) |> discard(is.null),
  nonb = map(set_names(arquivos_nonb, basename(arquivos_nonb)), carregar_genomic) |> discard(is.null),
  genes = map(set_names(arquivos_genes, basename(arquivos_genes)), carregar_genomic) |> discard(is.null),
  expressao = arquivos_expressao, # ainda não estou usando expressao
  hic = map(set_names(arquivos_hic, basename(arquivos_hic)), carregar_hic) |> discard(is.null))

# cria pasta se não existir
if (!dir.exists(dir_resultados)) {dir.create(dir_resultados)}

# =====================================================================
# Comparação SmAP1 vs NonB
for (nome_smap in names(dados$smap)) {
  smap <- dados$smap[[nome_smap]]
  for (nome_nonb in names(dados$nonb)) {
    nonb <- dados$nonb[[nome_nonb]]
    sobreposicoes <- findOverlaps(smap, nonb)
    if (length(sobreposicoes) > 0) {
      coord_smap <- smap[queryHits(sobreposicoes)]
      coord_nonb <- nonb[subjectHits(sobreposicoes)]
      resultado <- tibble(
        SmAP1_ID = if (!is.null(mcols(coord_smap)$ID)) mcols(coord_smap)$ID else paste0("SmAP1_", seq_along(coord_smap)),
        NonB_ID = if (!is.null(mcols(coord_nonb)$ID)) mcols(coord_nonb)$ID else paste0("NonB_", seq_along(coord_nonb)),
        SmAP1_chr = as.character(seqnames(coord_smap)),
        SmAP1_inicio = start(coord_smap),
        SmAP1_fim = end(coord_smap),
        SmAP1_tam = width(coord_smap),
        NonB_chr = as.character(seqnames(coord_nonb)),
        NonB_inicio = start(coord_nonb),
        NonB_fim = end(coord_nonb),
        NonB_tam = width(coord_nonb),
        Sobreposicao_inicio = pmax(start(coord_smap), start(coord_nonb)),
        Sobreposicao_fim = pmin(end(coord_smap), end(coord_nonb)),
        Sobreposicao_tam = Sobreposicao_fim - Sobreposicao_inicio + 1,
        Arquivo_SmAP1 = nome_smap,
        Arquivo_NonB = nome_nonb) %>%
        mutate(Porc_SmAP1 = round(Sobreposicao_tam / SmAP1_tam * 100, 2),
               Porc_NonB = round(Sobreposicao_tam / NonB_tam * 100, 2))
      arq_saida <- file.path(dir_resultados, sprintf("sobreposicoes_%s_vs_%s.csv", tools::file_path_sans_ext(nome_smap), tools::file_path_sans_ext(nome_nonb)))
      write_csv(resultado, arq_saida)
      message("Salvo: ", arq_saida)}}}

# =====================================================================
# Genes vs NonB -> 1000bp de promotor
for (gene_nome in names(dados$genes)) {
  genes <- dados$genes[[gene_nome]]
  reg_regulatorias <- promoters(genes, upstream=1000, downstream=1000)
  for (nonb_nome in names(dados$nonb)) {
    nonb <- dados$nonb[[nonb_nome]]
    overlaps <- findOverlaps(reg_regulatorias, nonb)
    if (length(overlaps) > 0) {
      reg_hit <- reg_regulatorias[queryHits(overlaps)]
      nonb_hit <- nonb[subjectHits(overlaps)]
      gene_hit <- genes[queryHits(overlaps)]
      gene_ids <- if (!is.null(mcols(gene_hit)$ID)) mcols(gene_hit)$ID else gene_hit$ID
      nonb_ids <- if (!is.null(mcols(nonb_hit)$ID)) mcols(nonb_hit)$ID else nonb_hit$ID
      resultado <- tibble(
        Gene_ID = gene_ids,
        NonB_ID = nonb_ids,
        Reg_chr = as.character(seqnames(reg_hit)),
        Reg_inicio = start(reg_hit),
        Reg_fim = end(reg_hit),
        Reg_tam = width(reg_hit),
        NonB_chr = as.character(seqnames(nonb_hit)),
        NonB_inicio = start(nonb_hit),
        NonB_fim = end(nonb_hit),
        NonB_tam = width(nonb_hit),
        Gene_chr = as.character(seqnames(gene_hit)),
        Gene_inicio = start(gene_hit),
        Gene_fim = end(gene_hit),
        Gene_tam = width(gene_hit),
        Gene_strand = as.character(strand(gene_hit)),
        Distancia = distance(gene_hit, nonb_hit),
        Arquivo_Gene = gene_nome,
        Arquivo_NonB = nonb_nome)
      arq_saida <- file.path(dir_resultados, sprintf("genes_regulatorio_%s_vs_%s.csv", tools::file_path_sans_ext(gene_nome), tools::file_path_sans_ext(nonb_nome)))
      write_csv(resultado, arq_saida)}}}

# =====================================================================
# Hi-C vs NonB -> aqui juntei first e second na mesma tabela de saída
# cada linha tem Tipo_Anchor = "first" ou "second"
if (length(dados$hic) > 0) {
  for (hic_nome in names(dados$hic)) {
    hic <- dados$hic[[hic_nome]]
    anchors_first <- anchors(hic, type = "first")
    anchors_second <- anchors(hic, type = "second")
    for (nonb_nome in names(dados$nonb)) {
      nonb <- dados$nonb[[nonb_nome]]
      overlaps_first <- findOverlaps(anchors_first, nonb)
      overlaps_second <- findOverlaps(anchors_second, nonb)
      resultados <- list()
      if (length(overlaps_first) > 0) {
        qh_f <- queryHits(overlaps_first); sh_f <- subjectHits(overlaps_first)
        segunda_f <- anchors_second[qh_f]
        resultados[["first"]] <- tibble(
          Interacao_ID = hic$id[qh_f],
          Tipo_Anchor = "first",
          HiC_chr = as.character(seqnames(anchors_first[qh_f])),
          HiC_inicio = start(anchors_first[qh_f]),
          HiC_fim = end(anchors_first[qh_f]),
          NonB_ID = nonb[sh_f]$ID,
          NonB_chr = as.character(seqnames(nonb[sh_f])),
          NonB_inicio = start(nonb[sh_f]),
          NonB_fim = end(nonb[sh_f]),
          Segunda_chr = as.character(seqnames(segunda_f)),
          Segunda_inicio = start(segunda_f),
          Segunda_fim = end(segunda_f),
          Arquivo_HiC = hic_nome,
          Arquivo_NonB = nonb_nome)}
      if (length(overlaps_second) > 0) {
        qh_s <- queryHits(overlaps_second); sh_s <- subjectHits(overlaps_second)
        segunda_s <- anchors_first[qh_s]
        resultados[["second"]] <- tibble(
          Interacao_ID = hic$id[qh_s],
          Tipo_Anchor = "second",
          HiC_chr = as.character(seqnames(anchors_second[qh_s])),
          HiC_inicio = start(anchors_second[qh_s]),
          HiC_fim = end(anchors_second[qh_s]),
          NonB_ID = nonb[sh_s]$ID,
          NonB_chr = as.character(seqnames(nonb[sh_s])),
          NonB_inicio = start(nonb[sh_s]),
          NonB_fim = end(nonb[sh_s]),
          Segunda_chr = as.character(seqnames(segunda_s)),
          Segunda_inicio = start(segunda_s),
          Segunda_fim = end(segunda_s),
          Arquivo_HiC = hic_nome,
          Arquivo_NonB = nonb_nome)}
      if (length(resultados) > 0) {
        resultado_final <- bind_rows(resultados)
        arq_saida <- file.path(dir_resultados, sprintf("hic_%s_vs_%s.csv", tools::file_path_sans_ext(hic_nome), tools::file_path_sans_ext(nonb_nome)))
        write_csv(resultado_final, arq_saida)}}}}
