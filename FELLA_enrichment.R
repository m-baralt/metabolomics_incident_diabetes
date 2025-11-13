library(FELLA)
library(ggplot2)
library(dplyr)
library(viridis)
library(forcats)

graph <- buildGraphFromKEGGREST(
  organism = "hsa",
  filter.path = c("01100", "01200", "01210", "01212", "01230"))


tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(
  keggdata.graph = graph,
  databaseDir = tmpdir,
  internalDir = FALSE,
  matrices = c("diffusion", "hypergeom"),
  normality = "diffusion",
  niter = 50)

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = c("diffusion", "hypergeom")
)

cat(getInfo(fella.data))

fella_enrichment <- function(sign_cmps, fella_data){
  
  sign_cmps <- unique(sign_cmps)
  hgeom_matrix <- fella.data@hypergeom@matrix
  
  analysis <- defineCompounds(
    compounds = sign_cmps,
    data = fella_data)
  
  analysis_diff <- runDiffusion(
    object = analysis,
    data = fella_data,
    approx = "normality")
  
  tab_diff <- generateResultsTable(
    method = "diffusion",
    nlimit = 1000,
    object = analysis_diff,
    data = fella_data)
  
  tab_diff <- tab_diff[tab_diff$Entry.type %in% "pathway",]
  
  tab_diff <- plyr::ldply(tab_diff$KEGG.id, function(hsa){
    cmps_hsa <- names(which(hgeom_matrix[,hsa]))
    #sign_cmps <- unique(sign_cmps)
    sign_cmps_hsa <- sum(sign_cmps %in% cmps_hsa)
    compl_name <- fella_data@keggdata@id2name[[hsa]]
    tab_diff <- data.frame(tab_diff[tab_diff$KEGG.id %in% hsa,], Hits = sign_cmps_hsa, Total_cmps = length(cmps_hsa))
    tab_diff$KEGG.name <- compl_name
    tab_diff
  })
  
  analysis_hgeom <- runHypergeom(
    object = analysis,
    data = fella_data,
    p.adjust = "fdr")
  
  tab_hgeom <- generateResultsTable(
    method = "hypergeom",
    nlimit = 100,
    object = analysis_hgeom,
    data = fella_data)
  
  return(list(table_diff = tab_diff,
              table_hgeom = tab_hgeom))
}

sign_cmps_analysis <- readRDS("processed_files/list_sign_cmps.rds")

enrich_r1 <- fella_enrichment(sign_cmps = names(sign_cmps_analysis$R1), fella_data = fella.data)

enrich_r3 <- fella_enrichment(sign_cmps = names(sign_cmps_analysis$R3), fella_data = fella.data)

enrich_r4 <- fella_enrichment(sign_cmps = names(sign_cmps_analysis$R4), fella_data = fella.data)

enrich_r5 <- fella_enrichment(sign_cmps = names(sign_cmps_analysis$R5), fella_data = fella.data)

enrich_ir <- fella_enrichment(sign_cmps = names(sign_cmps_analysis$IR), fella_data = fella.data)

diff_r1 <- enrich_r1$table_diff
diff_r1$analysis <- "R1"
diff_r3 <- enrich_r3$table_diff
diff_r3$analysis <- "R3"
diff_r4 <- enrich_r4$table_diff
diff_r4$analysis <- "R4"
diff_r5 <- enrich_r5$table_diff
diff_r5$analysis <- "R5"
diff_ir <- enrich_ir$table_diff
diff_ir$analysis <- "IR"

hgeom_r3 <- enrich_r3$table_hgeom
hgeom_r4 <- enrich_r4$table_hgeom
hgeom_r3$KEGG.name <- gsub(hgeom_r3$KEGG.name, pattern = " [-] Homo sapiens [(]human[)]", replacement = "")
hgeom_r4$KEGG.name <- gsub(hgeom_r4$KEGG.name, pattern = " [-] Homo sapiens [(]human[)]", replacement = "")

diff_all <- rbind(diff_r1, diff_r3, diff_r4, diff_ir)
#diff_all <- rbind(diff_r1, diff_r3, diff_r4)
diff_all$KEGG.name <- gsub(diff_all$KEGG.name, pattern = " [-] Homo sapiens [(]human[)]", replacement = "")

analysis_names <- c(
  `R1` = "T2D development - All subjects",
  `R3` = "Glycemia status",
  `R4` = "T2D development - No prediabetes tf",
  `IR` = "Insulin resistance"
)

diff_all$signif <- diff_all$p.score
diff_all$signif[diff_all$p.score<=0.05] <- "*"
diff_all$signif[diff_all$p.score<=0.01] <- "**"
diff_all$signif[diff_all$p.score<=0.001] <- "***"
diff_all$signif[diff_all$p.score<=0.0001] <- "****"

diff_all$approach <- "Diffusion"
diff_all$approach[(diff_all$KEGG.id %in% hgeom_r3$KEGG.id) & (diff_all$analysis %in% "R3")] <- "Diffusion | Hypergeometric"
diff_all$approach[(diff_all$KEGG.id %in% hgeom_r4$KEGG.id) & (diff_all$analysis %in% "R4")] <- "Diffusion | Hypergeometric"

hgeom_r3$signif <- "*"
hgeom_r3$signif[hgeom_r3$p.value<=0.01] <- "**"
hgeom_r3$signif[hgeom_r3$p.value<=0.001] <- "***"
hgeom_r3$signif[hgeom_r3$p.value<=0.0001] <- "****"

hgeom_r4$signif <- "*"
hgeom_r4$signif[hgeom_r4$p.value<=0.01] <- "**"
hgeom_r4$signif[hgeom_r4$p.value<=0.001] <- "***"
hgeom_r4$signif[hgeom_r4$p.value<=0.0001] <- "****"


total_r3 <- plyr::ldply(hgeom_r3$KEGG.id, function(pway_id){
  dd <- diff_all[(diff_all$KEGG.id %in% pway_id) & (diff_all$analysis %in% "R3"),]
  dd$signif <- paste(dd$signif, hgeom_r3$signif[hgeom_r3$KEGG.id %in% pway_id], sep = "|")
  dd
})

total_r4 <- plyr::ldply(hgeom_r4$KEGG.id, function(pway_id){
  dd <- diff_all[(diff_all$KEGG.id %in% pway_id) & (diff_all$analysis %in% "R4"),]
  if (nrow(dd)>0){
    dd$signif <- paste(dd$signif, hgeom_r4$signif[hgeom_r4$KEGG.id %in% pway_id], sep = "|")
    dd
  }
})


hgeom_res <- data.frame(KEGG.id = hgeom_r4$KEGG.id[3], Entry.type = "pathway", 
                        KEGG.name = hgeom_r4$KEGG.name[3], p.score = hgeom_r4$p.value[3],
                        Hits = hgeom_r4$CompoundHits[3], Total_cmps = hgeom_r4$CompoundsInPathway[3], 
                        analysis = "R4", signif = "*", approach = "Hypergeometric")

diff_all <- rbind(diff_all[diff_all$approach %in% "Diffusion",], total_r3, total_r4, hgeom_res)

diff_all <- diff_all[diff_all$Hits>0,]

hsa_minhits <- unlist(plyr::llply(unique(diff_all$KEGG.id), function(i){
  dx <- diff_all[diff_all$KEGG.id %in% i,]
  if (sum(dx$Hits>2)>0) return(i)
}))

p_bar <- ggplot(diff_all[diff_all$KEGG.id %in% hsa_minhits,], aes(x = KEGG.name, y = Hits, fill = approach)) +
  geom_col(width = 0.8, alpha = 0.7) + 
  coord_flip() + facet_grid(~analysis, labeller = labeller(analysis = analysis_names))+
  theme_bw() +
  geom_text(aes(label = Hits), hjust = 1.8, size = 5.5, 
            color = "white", fontface = "bold") + 
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_text(size = 17),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 13), 
        axis.title.x = element_text(size = 16), legend.position = "bottom", 
        legend.text = element_text(size = 15), legend.title = element_blank()) +
  labs(y = "Number of significant compounds per pathway and analysis", 
       fill = "Hypergeometric significance") +
  geom_text(aes(x= KEGG.name, y = Hits+(0.4*nchar(signif)), label = signif), 
            vjust = 0.7,size = 5, family = "Courier") + ylim(0, max(diff_all$Hits)+2) + 
  scale_fill_manual(values =c("#5E3C99", "goldenrod1", "#E66100"))


ggsave(plot = p_bar, filename = "figures/Figure1.jpg", 
       width = 26, height = 12, dpi = 900)

library(pathview)
setwd("figures/pathview/")
for (i in names(sign_cmps_analysis)){
  metab.data <- sign_cmps_analysis[[i]]
  hsa_ids <- unique(diff_all$KEGG.id[diff_all$analysis %in% i])
  pathview(
    cpd.data = metab.data,        # metabolite data
    pathway.id = hsa_ids,      # pyruvate metabolism pathway
    species = "hsa",              # human
    out.suffix = paste0("metab_", i),    # output file suffix
    #kegg.native = TRUE,           # use KEGG native pathway diagram
    low = list(cpd = "red"),    # color for downregulated
    mid = list(cpd = "gray"),     # color for unchanged (optional)
    high = list(cpd = "green"),      # color for upregulated
    limit = list(cpd = c(min(metab.data), max(metab.data)))
  ) 
}

