library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
library(parallel)
# Janelas -> região próxima ao gene (classificação de vizinhança)
janela_250 <- 250L
janela_500 <- 500L

# paralelismo -> g4 é grande -> acelera execução
n_threads <- 8
diretorio_dados <- "/home/mcavalcante/igv/modificados"
diretorio_resultados <- "/home/mcavalcante/sobreposição"
if (!dir.exists(diretorio_resultados)) dir.create(diretorio_resultados, recursive = TRUE)
to_GRanges <- function(x){
  if (is.null(x)) return(NULL)
  if (is(x,"GRanges")) return(x)
  if (is(x,"GRangesList")) { y <- unlist(x,use.names=FALSE); if (length(y)==0) return(NULL); return(y) }
  if (is.data.frame(x) && all(c("seqnames","start","end") %in% names(x)))
    return(GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns=TRUE))
  return(NULL)}
# tem lugar que está "chr1" ou "1" -> mapear para acessões NC_*
mapear_cromossomos <- function(granges_obj) {
  cromossomos <- as.character(seqnames(granges_obj))
  cromossomos_mapeados <- dplyr::case_when(
    str_detect(cromossomos, "NC_002607|chr1|^1$") ~ "NC_002607.1",
    str_detect(cromossomos, "NC_001869|chr2|^2$") ~ "NC_001869.1",
    str_detect(cromossomos, "NC_002608|chr3|^3$") ~ "NC_002608.1",
    TRUE ~ cromossomos  )
  novo <- granges_obj
  seqnames(novo) <- factor(cromossomos_mapeados, levels = unique(cromossomos_mapeados))
  novo }
carregar_genomic <- function(caminho) {
  tryCatch({
    if (str_detect(caminho, regex("\\.gff3?$", ignore_case = TRUE))) import.gff(caminho)
    else if (str_detect(caminho, regex("\\.bed$",   ignore_case = TRUE))) import.bed(caminho)
    else NULL }, error = function(e) NULL)}
# separar por tipo (usa nome do arquivo)
extrair_tipo_nonb <- function(caminho_arquivo) {
  nome <- basename(caminho_arquivo)
  if (str_detect(nome, regex("A-phased", ignore_case = TRUE))) return("A-Phased")
  if (str_detect(nome, regex("G-quadruplex|gquad|g_quad|g-quad", ignore_case = TRUE))) return("G-Quadruplex")
  if (str_detect(nome, regex("Z-DNA|zdna|z_dna", ignore_case = TRUE))) return("Z-DNA")
  if (str_detect(nome, regex("triplex", ignore_case = TRUE))) return("Triplex")
  if (str_detect(nome, regex("r-?loop|r_loop", ignore_case = TRUE))) return("R-Loop")
  if (str_detect(nome, regex("short_tandem|slipped", ignore_case = TRUE))) return("Short_Tandem")
  if (str_detect(nome, regex("Cruciform", ignore_case = TRUE))) return("Cruciform")
  "Non-B"}


classificar_regiao <- function(nonb, gene, win250 = 250L, win500 = 500L) {
  # pega fita e coordenadas do gene
  gstr <- as.character(strand(gene))         # "+" ou "-"
  gstart <- start(gene)                      # inicio do gene
  gend   <- end(gene)                        # fim
  glen   <- width(gene)                      # tamanho
  
  # coordenadas da nonb
  nstart <- start(nonb)
  nend   <- end(nonb)
  meio_gene <- gstart + floor(glen/2)
  mid <- floor((nstart + nend) / 2) #medio nb
  # todo no gene?
  inside <- (nstart >= gstart & nend <= gend)
  #se ta no gene
  reg_corpo <- dplyr::case_when(
    inside & nend   <= meio_gene ~ "Inicio",
    inside & nstart >  meio_gene ~ "Final",
    TRUE ~ NA_character_  )
  #sentido do gene
  if (gstr == "+") {
    # antes = mid < gstart depois = mid > gend
    dist_antes  <- gstart - mid   # >0 => antes
    dist_depois <- mid - gend  }    # >0 => depois
  else {
    # - ->  antes = mid > gend ; depois = mid < gstart
    dist_antes  <- mid - gend     # >0 => antes
    dist_depois <- gstart - mid   }  # >0 => depoi
  
  #fora 
  reg_antes <- dplyr::case_when(
    !inside & dist_antes  > 0 & dist_antes  <= win250 ~ "antes_250",
    !inside & dist_antes  > win250 & dist_antes  <= win500 ~ "antes_500",
    TRUE ~ NA_character_ )
  reg_depois <- dplyr::case_when(
    !inside & dist_depois > 0 & dist_depois <= win250 ~ "depois_250",
    !inside & dist_depois > win250 & dist_depois <= win500 ~ "depois_500",
    TRUE ~ NA_character_)
  # corpo > antes > depois > outro
  if (!is.na(reg_corpo))  return(reg_corpo)
  if (!is.na(reg_antes))  return(reg_antes)
  if (!is.na(reg_depois)) return(reg_depois)
  "Outro"}
# gff
arquivo_genes <- file.path(diretorio_dados, "Hsalinarum.gff")
genes <- import(arquivo_genes, format = "gff3")
genes <- to_GRanges(genes)
if (!is.null(genes) && length(genes) > 0 && "type" %in% colnames(mcols(genes))) {
  sel <- which(mcols(genes)$type %in% c("gene","Gene"))
  if (length(sel) > 0) genes <- genes[sel]
}
if (is.null(genes) || length(genes) == 0) stop("Arquivo GFF sem entradas 'gene' após filtro.")
# genes <- mapear_cromossomos(genes)
#---------------------------------------- n mexer
# arquivos de non-B 
arquivos <- list.files(diretorio_dados, full.names = TRUE)
padroes_nonb <- c("z-dna","gquad","nonb","triplex","r-loop","short_tandem","A-phased","G-quadruplex","Cruciform")
arquivos_nonb <- arquivos[
  sapply(arquivos, function(x) any(sapply(padroes_nonb, function(p) str_detect(basename(x), regex(p, ignore_case = TRUE))))) &
    str_detect(arquivos, regex("\\.(bed|gff3?|gff)$", ignore_case = TRUE))]

# overlaps -> classifica -> fração 
processar_arquivo <- function(arquivo_nonb){
  nonb_data <- tryCatch({ if (str_detect(arquivo_nonb,regex("\\.bed$",TRUE))) import.bed(arquivo_nonb) else import.gff(arquivo_nonb) }, error=function(e) NULL)
  nonb_data <- to_GRanges(nonb_data); if (is.null(nonb_data) || length(nonb_data)==0) return(tibble())
  
  out <- tibble()
  
  crom_genes <- unique(as.character(seqnames(genes))); crom_nonb <- unique(as.character(seqnames(nonb_data))); comuns <- intersect(crom_genes,crom_nonb)
  if (length(comuns)==0){
    nonb_map <- mapear_cromossomos(nonb_data)
    crom_nonb2 <- unique(as.character(seqnames(nonb_map))); comuns2 <- intersect(crom_genes,crom_nonb2)
    if (length(comuns2)==0) return(out)
    nonb_filtrado <- nonb_map[seqnames(nonb_map) %in% comuns2]; genes_filtro <- genes[seqnames(genes) %in% comuns2] } 
  else {
    nonb_filtrado <- nonb_data[seqnames(nonb_data) %in% comuns]; genes_filtro <- genes[seqnames(genes) %in% comuns] }
  
  if (length(nonb_filtrado)==0 || length(genes_filtro)==0) return(out)
  
  ov <- findOverlaps(nonb_filtrado, genes_filtro, ignore.strand=FALSE); if (length(ov)==0) return(out)
  
  tipo <- extrair_tipo_nonb(arquivo_nonb); nome_arq <- basename(arquivo_nonb)
  qh <- queryHits(ov); sh <- subjectHits(ov); res <- vector("list", length(qh))
  for (i in seq_along(qh)){
    nf <- nonb_filtrado[qh[i]]; gf <- genes_filtro[sh[i]]
    regiao <- classificar_regiao(nf,gf,win250=janela_250,win500=janela_500)
    ovseg <- GenomicRanges::intersect(nf,gf); frac <- if (width(nf)>0) sum(width(ovseg))/width(nf) else 0
    res <- vector("list", length(qh))
    for (i in seq_along(qh)){
      nf <- nonb_filtrado[qh[i]]; gf <- genes_filtro[sh[i]]
      regiao <- classificar_regiao(nf,gf,win250=janela_250,win500=janela_500)
      ovseg <- GenomicRanges::intersect(nf,gf); frac <- if (width(nf)>0) sum(width(ovseg))/width(nf) else 0
      res[[i]] <- tibble(Tipo_NonB=tipo, Arquivo=nome_arq, Regiao=regiao, Fracao_Sobreposicao=frac,
                         Tamanho_NonB=width(nf), Cromossomo=as.character(seqnames(nf)),
                         Start_NonB=start(nf), End_NonB=end(nf), Start_Gene=start(gf), End_Gene=end(gf),
                         Strand_Gene=as.character(strand(gf)))
    }
    out <- dplyr::bind_rows(res)  # out sempre existe a partir daqui
    
    if (nrow(out)>0){
      nome_clean <- gsub("[^A-Za-z0-9]","_",nome_arq)
      nome_clean <- gsub("\\.(bed|gff3?|gff)$","",nome_clean, ignore.case=TRUE)
      write.csv(out, file.path(diretorio_resultados, paste0("resultado_",nome_clean,".csv")), row.names=FALSE)
    }
    out
 
    
# grafico ---------------------------------
    plot_tipo_nonb <- function(dados_tipo, tipo_rotulo){
      ordem <- c("antes_500","antes_250","Inicio","Final","depois_250","depois_500")
      dados_tipo$Regiao <- factor(dados_tipo$Regiao, levels=ordem)
      
      # tabela base com todos os níveis e n=0
      todos <- tibble::tibble(Regiao=factor(ordem, levels=ordem), n_base=0L)
      obs <- dados_tipo %>% dplyr::count(Regiao, name="n_obs")
      rotulos_tab <- dplyr::left_join(todos, obs, by="Regiao") %>%
        dplyr::mutate(n = dplyr::coalesce(n_obs, n_base),
                      lbl = paste0(as.character(Regiao), " (n=", n, ")")) %>%
        dplyr::select(Regiao, lbl)
      rotulos_map <- setNames(rotulos_tab$lbl, rotulos_tab$Regiao)  # cobre todos os níveis, sem NULL
      
      ggplot2::ggplot(dados_tipo, ggplot2::aes(x=Regiao, y=Fracao_Sobreposicao, fill=Regiao)) +
        ggplot2::geom_boxplot(width=0.6, size=0.7, color="black", outlier.shape=NA) +
        ggplot2::scale_x_discrete(drop=FALSE, labels=rotulos_map) +  # não gera NULL
        ggplot2::labs(title=paste("Distribuição de", tipo_rotulo), x="Região", y="Fração de Sobreposição") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1), legend.position="none") +
        ggplot2::scale_fill_brewer(palette="Set2")
    }

    for (tipo in tipos_unicos){
      dt <- dplyr::filter(dados_completos, Tipo_NonB==tipo)
      if (NROW(dt)==0) next  # nada a plotar
      p <- plot_tipo_nonb(dt, tipo)
      ggsave(file.path(diretorio_resultados, paste0("boxplot_", gsub("[^A-Za-z0-9]","_", tipo), ".png")),
             p, width=10, height=6, dpi=300)
    }
