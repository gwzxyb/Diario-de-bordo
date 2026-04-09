library(Biostrings)
library(triplex)
library(rtracklayer)
library(parallel)

if (!dir.exists("triplex")) {
  dir.create("triplex")
}

fa <- readDNAStringSet("/home/mcavalcante/dados_brutos/Hsalinarum.fa")

n_cores <- 8
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(Biostrings)
  library(triplex)
})

clusterExport(cl, c("fa"))

processar_sequencia <- function(i) {
  
  seqname <- names(fa)[i]
  short_name <- substr(seqname, 1, 11)
  
  t <- triplex.search(fa[[i]], 
                      min_score = 2.0,
                      p_value = 1.0,
                      min_len = 5)
  
  if (length(t) == 0) { return(NULL) }
  
  gr <- as(t, "GRanges")
  
  seqlevels(gr) <- short_name
  seqnames(gr) <- short_name
  
  mcols(gr)$sequencia_original <- seqname
  
  return(gr)
}

gr_list <- parLapply(cl, seq_along(fa), processar_sequencia)

stopCluster(cl)

# Remover elementos NULL da lista
gr_list <- gr_list[!sapply(gr_list, is.null)]

if (length(gr_list) == 0) {
  cat("NENHUM TRIPLEX ENCONTRADO!\n")
  
  # Segunda tentativa com parâmetros mínimos (sem paralelismo)
  gr_list <- list()
  
  for (i in seq_along(fa)) {
    
    seqname <- names(fa)[i]
    short_name <- substr(seqname, 1, 11)
    
    cat("Tentando", short_name, "com parâmetros mínimos...\n")
    
    # Parâmetros mínimos
    t <- triplex.search(fa[[i]], 
                        min_score = 1.0,
                        p_value = 1.0,
                        min_len = 3)
    
    if (length(t) > 0) {
      cat("  Encontrado!", length(t), "triplex\n")
      gr <- as(t, "GRanges")
      seqlevels(gr) <- short_name
      seqnames(gr) <- short_name
      mcols(gr)$sequencia_original <- seqname
      gr_list[[i]] <- gr
    } else {
      cat("  Nenhum triplex encontrado\n")
    }
  }
  
  gr_list <- gr_list[!sapply(gr_list, is.null)]
}

if (length(gr_list) == 0) {
  stop("Execução interrompida: nenhum triplex encontrado")
}

gr_all <- do.call(c, gr_list)

cat("Total de triplex encontrados:", length(gr_all), "\n\n")

# FILTROS

# Verificar se há resultados
if (length(gr_all) == 0) {
  stop()
}

# Filtro 1: Comprimento do triplex (stem) entre 10 e 100 pb
cat("Filtro 1: Comprimento entre 10 e 100 pb\n")

widths <- width(gr_all)
idx_len <- which(widths >= 10 & widths <= 100)

if (length(idx_len) > 0) {
  gr_len <- gr_all[idx_len]
} else {
  gr_len <- GRanges()
}

cat("  Antes:", length(gr_all), "\n")
cat("  Depois:", length(gr_len), "\n")
cat("  Removidos:", length(gr_all) - length(gr_len), "\n\n")

# Filtro 2: Espaçador (loop) ≤ 100 pb
cat("Filtro 2: Espaçador (loop) ≤ 100 pb\n")

if (length(gr_len) > 0) {
  spacer_size <- mcols(gr_len)$lend - mcols(gr_len)$lstart
  idx_spacer <- which(spacer_size <= 100)
  gr_spacer <- gr_len[idx_spacer]
} else {
  gr_spacer <- GRanges()
}

cat("  Antes:", length(gr_len), "\n")
cat("  Depois:", length(gr_spacer), "\n")
cat("  Removidos:", length(gr_len) - length(gr_spacer), "\n\n")

# Filtro 3: Score ≥ 6
cat("Filtro 3: Score ≥ 6\n")

if (length(gr_spacer) > 0) {
  scores <- mcols(gr_spacer)$score
  idx_score <- which(scores >= 6)
  gr_score <- gr_spacer[idx_score]
} else {
  gr_score <- GRanges()
}

cat("  Antes:", length(gr_spacer), "\n")
cat("  Depois:", length(gr_score), "\n")
cat("  Removidos:", length(gr_spacer) - length(gr_score), "\n\n")

# Filtro 4: Composição ≥ 90% purinas OU ≥ 90% pirimidinas
cat("Filtro 4: Composição ≥ 90% purinas OU ≥ 90% pirimidinas\n")

if (length(gr_score) > 0) {
  
  seqs_list <- list()
  
  for (i in seq_along(gr_score)) {
    seqname <- as.character(seqnames(gr_score)[i])
    
    fa_idx <- which(names(fa) == seqname)
    
    if (length(fa_idx) == 1) {
      start_pos <- start(gr_score)[i]
      end_pos <- end(gr_score)[i]
      
      seq_sub <- subseq(fa[[fa_idx]], start = start_pos, end = end_pos)
      seqs_list[[i]] <- seq_sub
    } else {
      seqs_list[[i]] <- DNAString("")
    }
  }
  
  seqs <- DNAStringSet(seqs_list)
  
  af <- alphabetFrequency(seqs, baseOnly = TRUE)
  widths_seq <- width(seqs)
  
  # Calcular porcentagem de purinas e pirimidinas
  purine_pct <- (af[, "A"] + af[, "G"]) / widths_seq
  pyrimidine_pct <- (af[, "C"] + af[, "T"]) / widths_seq
  
  # Identificar regiões que atendem ao critério
  idx_comp <- which(purine_pct >= 0.90 | pyrimidine_pct >= 0.90)
  gr_final <- gr_score[idx_comp]
  
} else {
  gr_final <- GRanges()
}

cat("  Antes:", length(gr_score), "\n")
cat("  Depois:", length(gr_final), "\n")
cat("  Removidos:", length(gr_score) - length(gr_final), "\n\n")

# Resultados finais
if (length(gr_final) > 0) {
  
  # Exportar para BED
  export(gr_final, "triplex/Hsalinarum_triplex.bed")
  cat("Arquivo salvo: triplex/Hsalinarum_triplex.bed\n")
  
  # Salvar também como RDS
  saveRDS(gr_final, "triplex/Hsalinarum_triplex.rds")
  cat("Arquivo salvo: triplex/Hsalinarum_triplex.rds\n")
  
  # Salvar tabela com estatísticas
  df <- data.frame(
    seqnames = as.character(seqnames(gr_final)),
    start = start(gr_final),
    end = end(gr_final),
    width = width(gr_final),
    spacer = mcols(gr_final)$lend - mcols(gr_final)$lstart,
    score = mcols(gr_final)$score,
    purine_pct = purine_pct[idx_comp],
    pyrimidine_pct = pyrimidine_pct[idx_comp],
    sequencia_original = mcols(gr_final)$sequencia_original
  )
  write.csv(df, "triplex/Hsalinarum_triplex.csv", row.names = FALSE)
  cat("Arquivo salvo: triplex/Hsalinarum_triplex.csv\n\n")
  
} else {
  cat("Nenhum resultado para salvar após os filtros\n\n")
}

# Estatísticas finais
if (length(gr_final) > 0) {
  
  cat("Total de triplex preditos:", length(gr_final), "\n\n")
  
  cat("Por sequência (plasmídeo):\n")
  seqnames_gr <- as.character(seqnames(gr_final))
  seq_counts <- table(seqnames_gr)
  print(seq_counts)
  
  cat("\n")
  
  # Estatísticas de comprimento
  cat("Comprimento dos triplex (pb):\n")
  cat("  Mínimo:", min(width(gr_final)), "\n")
  cat("  Médio:", round(mean(width(gr_final)), 1), "\n")
  cat("  Máximo:", max(width(gr_final)), "\n")
  cat("  Mediana:", median(width(gr_final)), "\n\n")
  
  # Estatísticas do espaçador (loop)
  cat("Tamanho do espaçador/loop (pb):\n")
  spacer_sizes <- mcols(gr_final)$lend - mcols(gr_final)$lstart
  cat("  Mínimo:", min(spacer_sizes), "\n")
  cat("  Médio:", round(mean(spacer_sizes), 1), "\n")
  cat("  Máximo:", max(spacer_sizes), "\n")
  cat("  Mediana:", median(spacer_sizes), "\n\n")
  
  # Estatísticas de score
  cat("Score dos triplex:\n")
  cat("  Mínimo:", min(mcols(gr_final)$score), "\n")
  cat("  Médio:", round(mean(mcols(gr_final)$score), 1), "\n")
  cat("  Máximo:", max(mcols(gr_final)$score), "\n")
  cat("  Mediana:", median(mcols(gr_final)$score), "\n\n")
  
  # Estatísticas de composição
  cat("Composição de purinas/pirimidinas:\n")
  cat("  % Purinas - Mínimo:", round(min(purine_pct[idx_comp]) * 100, 1), "%\n")
  cat("  % Purinas - Médio:", round(mean(purine_pct[idx_comp]) * 100, 1), "%\n")
  cat("  % Purinas - Máximo:", round(max(purine_pct[idx_comp]) * 100, 1), "%\n")
  cat("  % Pirimidinas - Mínimo:", round(min(pyrimidine_pct[idx_comp]) * 100, 1), "%\n")
  cat("  % Pirimidinas - Médio:", round(mean(pyrimidine_pct[idx_comp]) * 100, 1), "%\n")
  cat("  % Pirimidinas - Máximo:", round(max(pyrimidine_pct[idx_comp]) * 100, 1), "%\n\n")
  
  print(head(gr_final, 10))
  
} else {
  cat("Nenhum triplex encontrado após os filtros!\n")
  cat("Tente rodar com parâmetros ainda mais brandos.\n")
}