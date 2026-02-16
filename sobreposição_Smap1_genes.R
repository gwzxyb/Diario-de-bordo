

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(patchwork)
library(ggpubr)
library(InteractionSet)
library(RColorBrewer)

# CONFIGURAÇÃO INICIAL
# Diretórios de trabalho
dir_dados <- "/home/mcavalcante/igv/modificados"
dir_resultados <- "/home/mcavalcante/Git limpo/sobreposicao/"
dir_plots <- file.path(dir_resultados, "plots")

# Criar diretórios se necessário
dir.create(dir_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

# Tema consistente para publicações com cores suaves
theme_artigo <- theme_pubr() + 
  theme(
    legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#2C3E50"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#5D6D7E"),
    axis.text = element_text(size = 11, color = "#2C3E50"),
    axis.title = element_text(size = 12, face = "bold", color = "#2C3E50"),
    legend.text = element_text(size = 10, color = "#2C3E50"),
    legend.title = element_text(size = 11, face = "bold", color = "#2C3E50"),
    panel.background = element_rect(fill = "#F8F9F9"),
    plot.background = element_rect(fill = "#FFFFFF"),
    panel.grid.major = element_line(color = "#EAECEE", size = 0.3)
  )


# 1. LOCALIZAÇÃO DE ARQUIVOS

cat("\n1. LOCALIZAÇÃO DE ARQUIVOS:\n")
arquivos <- list.files(dir_dados, full.names = TRUE)

# Arquivos principais
arq_smap <- arquivos[str_detect(arquivos, regex("smap1", ignore_case = TRUE))]
arq_gff <- arquivos[str_detect(arquivos, regex("Hsalinarum.*gff", ignore_case = TRUE))]

# asRNA e IS elements
arq_asrna <- arquivos[str_detect(arquivos, regex("asRNA|asrna", ignore_case = TRUE))]
arq_is <- arquivos[str_detect(arquivos, regex("IS|insertion.*sequence", ignore_case = TRUE))]

# Arquivos Non-B (todos que parecem ser Non-B DNA)
arquivos_nonb <- arquivos[str_detect(arquivos, regex("z-dna|zdna|g-quadruplex|g4|triplex|r-loop|cruciform|tandem|a-phased", ignore_case = TRUE))]

cat("SMAP:", length(arq_smap), "| GFF:", length(arq_gff), "\n")
cat("asRNA:", length(arq_asrna), "| IS:", length(arq_is), "\n")
cat("Non-B:", length(arquivos_nonb), "arquivos encontrados\n")

# Verificar arquivos essenciais
if(length(arq_smap) == 0) stop("ERRO: Arquivo SMAP não encontrado!")
if(length(arq_gff) == 0) stop("ERRO: Arquivo GFF não encontrado!")

# 2. CLASSIFICAÇÃO DE ARQUIVOS Non-B

cat("\n2. CLASSIFICAÇÃO DE ARQUIVOS Non-B:\n")

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
  
  return("Outros")
}

if(length(arquivos_nonb) > 0) {
  nonb_classificado <- data.frame(
    arquivo = arquivos_nonb,
    nome_arquivo = basename(arquivos_nonb),
    classe = sapply(arquivos_nonb, classificar_nonb),
    stringsAsFactors = FALSE
  )
  
  # Salva a tabela de classificação
  write.table(nonb_classificado,
              file.path(dir_resultados, "nonb_classificacao.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Mostra quantos arquivos existem por classe
  cat("\n  Distribuição de arquivos Non-B por classe:\n")
  print(table(nonb_classificado$classe))
  
  # Gráfico da distribuição de classes Non-B
  df_nonb_files <- as.data.frame(table(nonb_classificado$classe))
  names(df_nonb_files) <- c("Classe", "Frequencia")
  
  p_nonb_files <- df_nonb_files %>%
    ggplot(aes(x = reorder(Classe, Frequencia), y = Frequencia, fill = Classe)) +
    geom_col(alpha = 0.8, width = 0.7, color = "white") +
    geom_text(aes(label = Frequencia), hjust = -0.1, size = 5, fontface = "bold", color = "#2C3E50") +
    coord_flip() +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set3"))(nrow(df_nonb_files))) +
    labs(title = "Distribuição de arquivos Non-B por classe",
         x = NULL, y = "Número de arquivos") +
    theme_artigo + theme(legend.position = "none")
  
  ggsave(file.path(dir_plots, "00_Distribuicao_NonB_Arquivos.png"), 
         p_nonb_files, width = 10, height = 6, dpi = 400)
}


# 3. CARREGAMENTO DE DADOS

cat("\n3. CARREGAMENTO DE DADOS:\n")

# SMAP GFF3 com validação
cat("  Carregando SMAP... ")
smap_gr <- import(arq_smap[1], format = "GFF")
seqlevels(smap_gr) <- gsub("^chr", "", seqlevels(smap_gr))
smap_gr <- smap_gr[width(smap_gr) > 0]
cat(sprintf("%d regiões válidas\n", length(smap_gr)))

# GFF principal
cat("  Carregando GFF... ")
gff_gr <- import(arq_gff[1], format = "GFF")
seqlevels(gff_gr) <- gsub("^chr", "", seqlevels(gff_gr))
mcols_gff <- as.data.frame(mcols(gff_gr))
cat(sprintf("%d features\n", length(gff_gr)))

# Features GFF - filtros precisos
genes_gr <- gff_gr[mcols_gff$type == "gene"]

# LISTA ESPECÍFICA DE LOCI VNG (TnpB) - Corrigido para os loci reais
vng_loci <- c(
 "tnpb"
)

cat("\n  Buscando loci VNG (TnpB) no GFF...\n")

# Buscar em todos os campos relevantes do GFF
tnpbs_gr_list <- list()

for(locus in vng_loci) {
  locus_gr <- gff_gr[
    grepl(locus, mcols_gff$locus_tag, ignore.case = TRUE) |
      grepl(locus, mcols_gff$ID, ignore.case = TRUE) |
      grepl(locus, mcols_gff$Name, ignore.case = TRUE) |
      grepl(locus, mcols_gff$product, ignore.case = TRUE) |
      grepl(locus, mcols_gff$gene, ignore.case = TRUE)
  ]
  
  if(length(locus_gr) > 0) {
    cat(sprintf("    ✓ %s encontrado (%d features)\n", locus, length(locus_gr)))
    tnpbs_gr_list[[locus]] <- locus_gr
  } else {
    cat(sprintf("    ✗ %s NÃO encontrado\n", locus))
  }
}

# Combinar todos os loci encontrados
if(length(tnpbs_gr_list) > 0) {
  tnpbs_gr <- Reduce(c, tnpbs_gr_list)
  tnpbs_gr <- tnpbs_gr[!duplicated(tnpbs_gr)]
  cat(sprintf("\n  Total de loci VNG encontrados: %d/%d\n", length(tnpbs_gr_list), length(vng_loci)))
  cat(sprintf("  Total de features VNG: %d\n", length(tnpbs_gr)))
} else {
  tnpbs_gr <- GRanges()
  cat("\n  ⚠ NENHUM loci VNG encontrado! Verifique o arquivo GFF.\n")
}

# Transposases (outras, excluindo VNG loci)
transp_gr <- gff_gr[
  grepl("transposase|transposon", mcols_gff$product, ignore.case = TRUE) &
    !grepl(paste(vng_loci, collapse = "|"), mcols_gff$locus_tag, ignore.case = TRUE) &
    !grepl(paste(vng_loci, collapse = "|"), mcols_gff$ID, ignore.case = TRUE) &
    !grepl(paste(vng_loci, collapse = "|"), mcols_gff$Name, ignore.case = TRUE)
]

# asRNA e IS elements
asrna_gr <- GRanges()
if(length(arq_asrna) > 0) {
  asrna_gr <- tryCatch(
    import(arq_asrna[1], format = "BED"),
    error = function(e) import(arq_asrna[1], format = "GFF")
  )
  seqlevels(asrna_gr) <- gsub("^chr", "", seqlevels(asrna_gr))
}

is_gr <- GRanges()
if(length(arq_is) > 0) {
  is_gr <- tryCatch(
    import(arq_is[1], format = "BED"),
    error = function(e) import(arq_is[1], format = "GFF")
  )
  seqlevels(is_gr) <- gsub("^chr", "", seqlevels(is_gr))
}


# 4. CARREGAMENTO E CLASSIFICAÇÃO DE DADOS Non-B

cat("\n4. CARREGAMENTO DE DADOS Non-B:\n")

nonb_gr_list <- list()
nonb_all <- GRanges()

if(length(arquivos_nonb) > 0) {
  for(i in 1:nrow(nonb_classificado)) {
    arquivo <- nonb_classificado$arquivo[i]
    classe <- nonb_classificado$classe[i]
    
    gr <- tryCatch(
      import(arquivo, format = "BED"),
      error = function(e) {
        tryCatch(
          import(arquivo, format = "GFF"),
          error = function(e2) GRanges()
        )
      }
    )
    
    if(length(gr) > 0) {
      seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
      
      if(!classe %in% names(nonb_gr_list)) {
        nonb_gr_list[[classe]] <- gr
      } else {
        nonb_gr_list[[classe]] <- c(nonb_gr_list[[classe]], gr)
      }
      
      cat(sprintf("  ✓ %s: %d regiões\n", classe, length(gr)))
    }
  }
  
  # Combinar todas as regiões Non-B
  if(length(nonb_gr_list) > 0) {
    nonb_all <- Reduce(c, nonb_gr_list)
    cat(sprintf("\n  Total Non-B: %d regiões (%d classes)\n", 
                length(nonb_all), length(nonb_gr_list)))
  }
}


# 5. ANÁLISE DE DISTÂNCIA SMAP - VNG (TnpB)

cat("\n5. ANÁLISE DE DISTÂNCIA SMAP - VNG (TnpB):\n")
cat("============================================\n")

if(length(tnpbs_gr) > 0 && length(smap_gr) > 0) {
  # Features expandidas para análise de proximidade
  tnpbs_near_80bp <- resize(tnpbs_gr, width(tnpbs_gr) + 160, fix = "center")
  tnpbs_near_100bp <- resize(tnpbs_gr, width(tnpbs_gr) + 200, fix = "center")
  tnpbs_near_200bp <- resize(tnpbs_gr, width(tnpbs_gr) + 400, fix = "center")
  
  dists <- distanceToNearest(smap_gr, tnpbs_gr)
  dist_vec <- mcols(dists)$distance
  
  media_dist_bp <- mean(dist_vec, na.rm = TRUE)
  mediana_dist_bp <- median(dist_vec, na.rm = TRUE)
  media_dist_kb <- round(media_dist_bp / 1000, 2)
  mediana_dist_kb <- round(mediana_dist_bp / 1000, 2)
  min_dist_kb <- round(min(dist_vec) / 1000, 2)
  max_dist_kb <- round(max(dist_vec) / 1000, 2)
  
  cat(sprintf("  SMAP analisados: %d | VNG loci: %d\n", length(smap_gr), length(tnpbs_gr)))
  cat(sprintf("  DISTÂNCIA MÉDIA: %.2f kb (%.0f bp)\n", media_dist_kb, media_dist_bp))
  cat(sprintf("  DISTÂNCIA MEDIANA: %.2f kb (%.0f bp)\n", mediana_dist_kb, mediana_dist_bp))
  cat(sprintf("  Mínimo: %.2f kb | Máximo: %.2f kb\n", min_dist_kb, max_dist_kb))
  
  df_distancias <- tibble(
    smap_index = queryHits(dists),
    vng_index = subjectHits(dists),
    distancia_bp = dist_vec,
    distancia_kb = dist_vec / 1000
  )
  
  # Estatísticas por faixa de distância
  df_dist_stats <- df_distancias %>%
    mutate(faixa = cut(distancia_kb, 
                       breaks = c(0, 1, 5, 10, 50, 100, 200, Inf),
                       labels = c("0-1kb", "1-5kb", "5-10kb", "10-50kb", "50-100kb", "100-200kb", ">200kb"))) %>%
    group_by(faixa) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(percentual = round(100 * n / sum(n), 1))
  
  cat("\n  Distribuição por faixa:\n")
  print(df_dist_stats)
  
  # Criar histograma
  p_dist <- ggplot(df_distancias, aes(x = distancia_kb)) +
    geom_histogram(bins = 30, fill = "#8DA0CB", alpha = 0.8, color = "white") +
    geom_vline(xintercept = media_dist_kb, color = "#FC8D62", linetype = "dashed", size = 1.2) +
    geom_vline(xintercept = mediana_dist_kb, color = "#66C2A5", linetype = "dotted", size = 1) +
    annotate("text", x = media_dist_kb + 10, y = 20, 
             label = paste("Média =", media_dist_kb, "kb"), 
             color = "#FC8D62", size = 4, hjust = 0, fontface = "bold") +
    annotate("text", x = mediana_dist_kb + 10, y = 15, 
             label = paste("Mediana =", mediana_dist_kb, "kb"), 
             color = "#66C2A5", size = 4, hjust = 0, fontface = "bold") +
    labs(title = "Distância SMAP - VNG loci (TnpB)",
         subtitle = paste("n =", length(dist_vec), "SMAP | Média:", media_dist_kb, "kb | Mediana:", mediana_dist_kb, "kb"),
         x = "Distância (kb)", y = "Frequência") +
    theme_artigo
  
  ggsave(file.path(dir_plots, "01_Distancia_SMAP_VNG.png"), 
         p_dist, width = 12, height = 7, dpi = 400)
  
  write_csv(df_distancias, file.path(dir_resultados, "Distancias_SMAP_VNG.csv"))
  write_csv(df_dist_stats, file.path(dir_resultados, "Distribuicao_Distancias.csv"))
  
} else {
  cat("  ⚠ Dados insuficientes para análise de distância\n")
  media_dist_kb <- NA
  tnpbs_near_80bp <- tnpbs_near_100bp <- tnpbs_near_200bp <- GRanges()
}


# 6. MATRIZ DE SOBREPOSIÇÕES

cat("\n6. MATRIZ DE SOBREPOSIÇÕES:\n")

features <- list(
  Genes = genes_gr, 
  VNG_loci = tnpbs_gr, 
  Transposases = transp_gr, 
  asRNA = asrna_gr, 
  IS_elements = is_gr
)

# Adicionar features expandidas de VNG
if(length(tnpbs_gr) > 0) {
  features$VNG_80bp <- tnpbs_near_80bp
  features$VNG_100bp <- tnpbs_near_100bp
  features$VNG_200bp <- tnpbs_near_200bp
}

# Adicionar classes Non-B
if(length(nonb_gr_list) > 0) {
  for(classe in names(nonb_gr_list)) {
    if(length(nonb_gr_list[[classe]]) > 0) {
      features[[paste0("NonB_", classe)]] <- nonb_gr_list[[classe]]
    }
  }
}

matriz_sobre <- map_dfr(names(features), function(nome) {
  if(length(features[[nome]]) == 0) {
    return(tibble(Feature = nome, N_sobre = 0, Perc_SM = 0))
  }
  
  hits <- findOverlaps(smap_gr, features[[nome]])
  n_hits <- length(unique(queryHits(hits)))
  
  tibble(
    Feature = nome,
    N_sobre = n_hits,
    Perc_SM = round(100 * n_hits / length(smap_gr), 2)
  )
}) %>% arrange(desc(N_sobre))

write_csv(matriz_sobre, file.path(dir_resultados, "Matriz_Sobreposicoes.csv"))
print(matriz_sobre)

# Gráfico da matriz (apenas features com sobreposição)
matriz_plot <- matriz_sobre %>% 
  filter(N_sobre > 0) %>% 
  mutate(Feature = str_trunc(Feature, 25))

if(nrow(matriz_plot) > 0) {
  n_features <- nrow(matriz_plot)
  
  p_matriz <- matriz_plot %>%
    ggplot(aes(x = reorder(Feature, N_sobre), y = N_sobre, fill = Feature)) +
    geom_col(alpha = 0.8, width = 0.7, color = "white") +
    geom_text(aes(label = paste0(N_sobre, " (", Perc_SM, "%)")), 
              hjust = -0.1, size = 4, fontface = "bold", color = "#2C3E50") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(n_features)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Sobreposições SMAP x Features Genômicas",
         subtitle = paste("Total SMAP:", length(smap_gr), "| Features:", n_features),
         x = NULL, y = "Número de SMAP") +
    theme_artigo + theme(legend.position = "none")
  
  ggsave(file.path(dir_plots, "02_Matriz_Sobreposicoes.png"), 
         p_matriz, width = 14, height = 10, dpi = 400)
}


# 7. ANÁLISE VNG (TnpB) EM DIFERENTES JANELAS

cat("\n7. ANÁLISE VNG (TnpB) EM DIFERENTES JANELAS:\n")

if(length(tnpbs_gr) > 0) {
  # Sobreposição direta
  vng_hits <- findOverlaps(smap_gr, tnpbs_gr)
  vng_smap <- length(unique(queryHits(vng_hits)))
  vng_pct <- round(100 * vng_smap / length(smap_gr), 1)
  
  # Janelas de proximidade
  vng_hits_80 <- findOverlaps(smap_gr, tnpbs_near_80bp)
  vng_hits_100 <- findOverlaps(smap_gr, tnpbs_near_100bp)
  vng_hits_200 <- findOverlaps(smap_gr, tnpbs_near_200bp)
  
  vng_smap_80 <- length(unique(queryHits(vng_hits_80)))
  vng_smap_100 <- length(unique(queryHits(vng_hits_100)))
  vng_smap_200 <- length(unique(queryHits(vng_hits_200)))
  
  df_vng_analysis <- tibble(
    Feature = c("VNG direto", "VNG ≤80bp", "VNG ≤100bp", "VNG ≤200bp"),
    N_smap = c(vng_smap, vng_smap_80, vng_smap_100, vng_smap_200),
    Perc_SM = round(100 * N_smap / length(smap_gr), 1)
  )
  
  print(df_vng_analysis)
  write_csv(df_vng_analysis, file.path(dir_resultados, "Analise_VNG_Janelas.csv"))
  
  p_vng <- df_vng_analysis %>%
    ggplot(aes(x = reorder(Feature, N_smap), y = N_smap, fill = Feature)) +
    geom_col(alpha = 0.8, width = 0.7, color = "white") +
    geom_text(aes(label = paste0(N_smap, " (", Perc_SM, "%)")), 
              hjust = -0.1, size = 5, fontface = "bold", color = "#2C3E50") +
    coord_flip() +
    scale_fill_manual(values = brewer.pal(4, "Set2")) +
    labs(title = "VNG loci (TnpB) próximos a SMAP",
         subtitle = paste("Total SMAP:", length(smap_gr)),
         y = "Número de SMAP", x = NULL) +
    theme_artigo + theme(legend.position = "none")
  
  ggsave(file.path(dir_plots, "03_VNG_Janelas.png"), 
         p_vng, width = 10, height = 6, dpi = 400)
}


# 9. ANÁLISE DE CLASSES Non-B

if(length(nonb_gr_list) > 0) {
  cat("\n9. ANÁLISE DE CLASSES Non-B:\n")
  
  df_nonb_classes <- map_dfr(names(nonb_gr_list), function(classe) {
    gr <- nonb_gr_list[[classe]]
    if(length(gr) == 0) return(NULL)
    
    hits <- findOverlaps(smap_gr, gr)
    n_smap <- length(unique(queryHits(hits)))
    n_nonb <- length(unique(subjectHits(hits)))
    
    tibble(
      Classe = classe,
      N_smap = n_smap,
      Perc_SM = round(100 * n_smap / length(smap_gr), 1),
      N_nonb = n_nonb,
      Perc_nonb = round(100 * n_nonb / length(gr), 1)
    )
  }) %>% 
    filter(N_smap > 0) %>% 
    arrange(desc(N_smap))
  
  if(nrow(df_nonb_classes) > 0) {
    print(df_nonb_classes)
    write_csv(df_nonb_classes, file.path(dir_resultados, "Analise_NonB_Classes.csv"))
    
    p_nonb <- df_nonb_classes %>%
      ggplot(aes(x = reorder(Classe, N_smap), y = N_smap, fill = Classe)) +
      geom_col(alpha = 0.8, width = 0.7, color = "white") +
      geom_text(aes(label = paste0(N_smap, " (", Perc_SM, "%)")), 
                hjust = -0.1, size = 5, fontface = "bold", color = "#2C3E50") +
      coord_flip() +
      scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(nrow(df_nonb_classes))) +
      labs(title = "SMAP sobrepostos a Non-B DNA por classe",
           subtitle = paste("Total SMAP:", length(smap_gr)),
           y = "Número de SMAP", x = NULL) +
      theme_artigo + theme(legend.position = "none")
    
    ggsave(file.path(dir_plots, "05_NonB_Classes.png"), 
           p_nonb, width = 12, height = 8, dpi = 400)
  }
}


# 10. ANÁLISE asRNA E IS ELEMENTS

cat("\n10. ANÁLISE asRNA E IS ELEMENTS:\n")

df_rna_is <- bind_rows(
  if(length(asrna_gr) > 0) {
    hits <- findOverlaps(smap_gr, asrna_gr)
    tibble(
      Feature = "asRNA",
      N_smap = length(unique(queryHits(hits))),
      Perc_SM = round(100 * length(unique(queryHits(hits))) / length(smap_gr), 1)
    )
  } else NULL,
  if(length(is_gr) > 0) {
    hits <- findOverlaps(smap_gr, is_gr)
    tibble(
      Feature = "IS elements",
      N_smap = length(unique(queryHits(hits))),
      Perc_SM = round(100 * length(unique(queryHits(hits))) / length(smap_gr), 1)
    )
  } else NULL
)

if(nrow(df_rna_is) > 0) {
  print(df_rna_is)
  write_csv(df_rna_is, file.path(dir_resultados, "Analise_RNA_IS.csv"))
  
  p_rna_is <- df_rna_is %>%
    ggplot(aes(x = Feature, y = N_smap, fill = Feature)) +
    geom_col(alpha = 0.8, width = 0.6, color = "white") +
    geom_text(aes(label = paste0(N_smap, " (", Perc_SM, "%)")), 
              vjust = -0.5, size = 6, fontface = "bold", color = "#2C3E50") +
    scale_fill_manual(values = brewer.pal(nrow(df_rna_is), "Set3")) +
    labs(title = "SMAP sobrepostos a asRNA e IS elements",
         subtitle = paste("Total SMAP:", length(smap_gr)),
         y = "Número de SMAP", x = NULL) +
    theme_artigo + theme(legend.position = "none")
  
  ggsave(file.path(dir_plots, "06_asRNA_IS_elements.png"), 
         p_rna_is, width = 8, height = 6, dpi = 400)
}

