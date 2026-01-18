library(ggfortify)
library(magrittr)
library(plyr)
library(kableExtra)
library(RColorBrewer)
library(ggplot2)
library(qvalue)

path <- "/home/maria/metabolomics_incident_diabetes/"
setwd(path)
source("Functions_discovery.R")

################################################################################
############################## POSITIVE IONIZATION #############################
################################################################################

# data loading - positive ionization
MS_pos <- read.csv("data/DM_FIS2018_Hilic_pos_results2023_filled.csv")
rownames(MS_pos) <- paste0("SOI", MS_pos$Qidx)
MS_pos_t <- as.data.frame(t(MS_pos[,-c(1:9)]))
MS_pos_t <- MS_pos_t[!(grepl(pattern = "Bi", rownames(MS_pos_t)) | grepl(pattern = "Bf", rownames(MS_pos_t))),]

# data loading - injection order
inj_order_pos <- read.csv("data/DidacMauricio_hilic_pos_injectionorder.csv")
rownames(inj_order_pos) <- inj_order_pos$SampleName
inj_order_pos <- inj_order_pos[,-1]

# data loading - data base
db <- read.csv("data/rebuilt_db_29052023.csv")
db <- db[,-1]
rownames(db) <- db$ID.sample

# list of SOIs
feat.cols <- grep(pattern = "SOI", colnames(MS_pos_t), value = TRUE)
feats.before80 <- length(feat.cols)

# 80% rule filter and imputation with minimum/2
proc_object <- filtering_imputing(df = MS_pos_t, features.names = feat.cols, imp.method = "min/2")

# number of features after filtering
feats.after80 <- length(grep(pattern = "SOI", x = colnames(proc_object$Imputed), value = TRUE))

proc_object$log_transf <- proc_object$Imputed

for (j in 1:ncol(proc_object$log_transf)){
  proc_object$log_transf[,j] <- log10(proc_object$log_transf[,j])
}

res_prep_viz <- data_prep_vis(proc_object = proc_object, inj_order = inj_order_pos)
pca.initial.pos <- ggpubr::ggarrange(res_prep_viz$pca_data, 
                                     res_prep_viz$pca_inj, 
                                     ncol = 2)

pca_pos_none <- res_prep_viz$pca_data + labs(title = "No correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))

# outliers removals
data_log <- proc_object$log_transf
outlier_list <- plyr::llply(1:ncol(data_log), function(j){
  idx.rows <- which(data_log[,j] %in% boxplot.stats(data_log[,j])$out)
  out.samps <- rownames(data_log)[idx.rows]
  return(out.samps)
}) %>% unlist()

freq_out_samps <- as.data.frame(table(outlier_list))
freq_out_samps$prop <- freq_out_samps$Freq/feats.before80
out_samps <- as.character(freq_out_samps$outlier_list[freq_out_samps$prop>0.2])

data_log <- data_log[!(rownames(data_log) %in% out_samps),]
proc_object$log_transf <- data_log

res_prep_viz <- data_prep_vis(proc_object = proc_object, inj_order = inj_order_pos)
pca.initial.pos.nooutliers <- ggpubr::ggarrange(res_prep_viz$pca_data, 
                                                res_prep_viz$pca_inj, 
                                                ncol = 2)

proc_object <- res_prep_viz$proc_object

common_proc_object <- proc_object

###################################### CPCA ####################################

euclidean_dist <- function(p1, p2){
  p1 <- unname(p1) 
  p2 <- unname(p2) 
  sqrt(((p2[2]-p1[2])^2)+((p2[1]-p1[1])^2))
}

normalization_ns_pcas <- plyr::llply(seq_len(4), function(n){
  normalization_res <- data_normalization(proc_object = proc_object, ncomps = n)
  pca.data <- prcomp(normalization_res$norm_data[,-c(1,2)], center = TRUE, scale = TRUE)
  pcs <- pca.data$x[grepl("QC", rownames(pca.data$x)),c(1,2)]
  centroid_coords <- rearrr::centroid(pcs[,1], pcs[,2])
  dist.centroid <- unlist(plyr::llply(1:nrow(pcs), function(i) euclidean_dist(p1 = pcs[i,], p2 = centroid_coords)))
  
  df.dispersion <- data.frame(ncomps = n, dist.centroid = mean(dist.centroid))
  p <- normalization_res$pca_norm + 
    labs(title = paste0("Number of CPCs = ",n, " (dist=",round(mean(dist.centroid),2),")")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          legend.position = "bottom", legend.title = element_blank(), 
          legend.text = element_text(size =12))
  list(df.dispersion = df.dispersion, pca = p)
})

df.dispersion.all <- ldply(normalization_ns_pcas, function(n) n$df.dispersion)
pcas_norms <- llply(normalization_ns_pcas, function(n) n$pca)

pcas_ncomps <- ggpubr::ggarrange(plotlist = pcas_norms, nrow = 2, 
                                 ncol = 2, common.legend = TRUE, legend = "bottom") +
  theme(plot.title = element_text(size = 10))

ggsave(plot=pcas_ncomps, filename = "figures/PCA_cpca_ns_pos.jpeg", 
       width = 10, height = 8, dpi = 500)

# looking at the PCAs, we select the first n cpcs where technical bias is 
# suposedly removed in the first two PCs (n = 3)
normalization_res <- data_normalization(proc_object = proc_object, ncomps = 3)

pca_pos_cpca <- normalization_res$pca_norm + labs(title = "CPCA correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))


proc_object$log_transf <- normalization_res$norm_data
proc_object$scaled <- proc_object$log_transf
proc_object$scaled[,-c(1,2)] <- scale(proc_object$scaled[,-c(1,2)], 
                                      center = TRUE, scale = TRUE)

feats.to.keep <- CV_filtering(proc_object = proc_object)
feats.to.keep <- names(feats.to.keep)[feats.to.keep]

feats.afterCV <- length(feats.to.keep)

feats.df <- data.frame(Feats0 = feats.before80, feats80 = feats.after80, featsCV = feats.afterCV)

proc_object$scaled <- proc_object$scaled[,c("sampleID", "class", feats.to.keep)]

scaled_dat <- proc_object$scaled

db_original <- db

db[,c("Age","BMI","HT","HOMA_IR","Smoking","Family_history_DM", "Glucose","seguimiento")] <- 
  scale(db[,c("Age","BMI","HT","HOMA_IR","Smoking","Family_history_DM", "Glucose","seguimiento")], 
        center = T, scale = T)

merge.pos.cpca <- merge(db,proc_object$scaled[,-c(1,2)], by = 0)
rownames(merge.pos.cpca) <- merge.pos.cpca$Row.names
merge.pos.cpca <- merge.pos.cpca[,-c(1,2)]

data.frame(feats.before80, feats.after80, feats.afterCV)

merge_pos_original <- merge(db_original, proc_object$scaled[,-c(1,2)], by = 0)
merge_pos_original$HT <- as.factor(merge_pos_original$HT)
merge_pos_original$Smoking <- as.factor(merge_pos_original$Smoking)

merge_pos_original <- merge_pos_original[!is.na(merge_pos_original$DM_seguimiento),]

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

merge_pos_original <- merge_pos_original[,c("DM_INCIDENTE",covars)]
merge_pos_original <- na.omit(merge_pos_original) %>% as.data.frame()
merge_pos_original$Family_history_DM <- as.factor(merge_pos_original$Family_history_DM)
res1 <- compareGroups::compareGroups(DM_INCIDENTE ~ Sex+Age+BMI+HT+HOMA_IR+Smoking+
                                       seguimiento+Glucose+Family_history_DM, 
                                     data = merge_pos_original)
compTab <- compareGroups::createTable(res1, show.ratio = FALSE)

compTab

compareGroups::export2csv(compTab, "figures/comparegroups_discovery.csv")

############################# Statistical analysis #############################

merge.pos.cpca <- merge.pos.cpca[!is.na(merge.pos.cpca$DM_seguimiento),]

write.csv(merge.pos.cpca, "processed_files/cpca_processed_pos.csv")

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.pos.cpca[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos.cpca), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

#plot_perf <- perf_ncomp(x = x, y = merge.pos.cpca$DM_INCIDENTE, ncomp = 10) #ncomp = 3

# incident diabetis using all subjects
pls_result <- pls_stats(x = x, y = merge.pos.cpca$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

pls_result$plot

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos.cpca,
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos.cpca), value = TRUE),
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos.cpca, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15, psign = 0.05)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r1.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R1_pos_cpca.csv")

# incident diabetis with only prediabetics at t=0
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos.pred <- merge.pos.cpca[merge.pos.cpca$Prediabetes==1,]
xpred <- merge.pos.pred[,c(covars,grep(pattern = "SOI", x = colnames(merge.pos.pred), value = TRUE))]
xpred$Sex <- as.numeric(as.factor(xpred$Sex))
xpred$HT <- as.numeric(as.factor(xpred$HT))

#plot_perf <- perf_ncomp(x = xpred, y = merge.pos.pred$DM_INCIDENTE, ncomp = 10) #ncomp = 2

pls_result <- pls_stats(x = xpred, y = merge.pos.pred$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos.pred, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos.pred), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos.pred, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)


t(stats_results$summary.stats)

r2.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

write.csv(stats_results$Pvals.df, "results/discovery/R2_pos_cpca.csv")

# Incident diabetes status
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos2 <- merge.pos.cpca[!is.na(merge.pos.cpca$DM_seguimiento),]
x2 <- merge.pos2[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

#plot_perf <- perf_ncomp(x = x2, y = merge.pos2$DM_INCIDENTE, ncomp = 10) #ncomp = 5

pls_result <- pls_stats(x = x2, y = merge.pos2$DM_seguimiento, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)
volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)
t(stats_results$summary.stats)

r3.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

write.csv(stats_results$Pvals.df, "results/discovery/R3_pos_cpca.csv")

# incident diabetes without prediabetes at the end
merge.pos2 <- merge.pos.cpca[merge.pos.cpca$DM_seguimiento %in% c(0,2),
                             c(covars,"DM_seguimiento", "ID.sample", 
                               grep(pattern = "SOI", x = colnames(merge.pos.cpca), 
                                    value = TRUE))]

x2 <- merge.pos2[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

merge.pos2$DM_seguimiento[merge.pos2$DM_seguimiento==2] <- 1

# tune = tune.splsda(x2, merge.pos2$DM_seguimiento, ncomp = 10, nrepeat = 50, 
#                    logratio = "none", 
#                    folds = 10, dist = "max.dist",
#                    progressBar = TRUE)
# tune$choice.ncomp
# plot_perf <- perf_ncomp(x = x2, y = merge.pos2$DM_seguimiento, ncomp = 10) #ncomp = 6

pls_result <- pls_stats(x = x2, 
                        y = merge.pos2$DM_seguimiento, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos2, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos2, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

r4.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

write.csv(stats_results$Pvals.df, "results/discovery/R4_pos_cpca.csv")

# Glycemia change
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")

# glycemia change
# normoglycemia at t=0
merge.pos2 <- merge.pos.cpca[!is.na(merge.pos.cpca$DM_seguimiento),]
merge.pos2$Glyc_change <- paste(merge.pos2$Prediabetes, merge.pos2$DM_seguimiento, sep = "->")
merge.pos2$Glyc_change <- as.factor(merge.pos2$Glyc_change)
table(merge.pos2$Glyc_change)

merge.pos20 <- merge.pos2[merge.pos2$Glyc_change %in% c("0->0", "0->1", "0->2"),]
merge.pos20$Glyc_change <- as.factor(as.character(merge.pos20$Glyc_change))
levels(merge.pos20$Glyc_change) <- c("0", "1", "2")
merge.pos20$Glyc_change <- as.numeric(as.character(merge.pos20$Glyc_change))

x2 <- merge.pos20[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos20), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.pos20$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos20, 
                                         feats.names = grep(pattern = "SOI", 
                                                            x = colnames(merge.pos20), 
                                                            value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos20, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r5.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R5_pos_cpca.csv")

# prediabetes at t=0
merge.pos21 <- merge.pos2[merge.pos2$Glyc_change %in% c("1->0", "1->1", "1->2"),]
merge.pos21$Glyc_change <- as.factor(as.character(merge.pos21$Glyc_change))
levels(merge.pos21$Glyc_change) <- c("0", "1", "2")
merge.pos21$Glyc_change <- as.numeric(as.character(merge.pos21$Glyc_change))

x2 <- merge.pos21[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos21), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.pos21$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos21, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos21), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos21, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r6.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R6_pos_cpca.csv")

# Glycemia at follow-up
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

stats_results <- statistical_analysis_lm(df = merge.pos.cpca, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos.cpca), value = TRUE), 
                                         covars = covars, outcome = "glucosa_FU",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

r7.pos.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R7_pos_cpca.csv")

###################################### LOESS ###################################

proc_object <- common_proc_object

res <- proc_object$log_transf
res <- merge(inj_order_pos, res, by = 0)
rownames(res) <- res$Row.names
res <- res[,-1]
res2 <- res[,-c(1:4)]
res2 <- t(res2)
batch <- rep(1, times=ncol(res2))

sum(is.na(res2))

corrected_data <- pmp::QCRSC(df=res2, order=res$injectionOrder, 
                             batch=batch,
                             log = F, # log = F log-transformation previously applied
                             classes=res$class,
                             spar_lim = c(-2,2),
                             minQC=10)

corrdata <- as.data.frame(t(corrected_data))
sum(is.na(corrdata))
corrdata <- data.frame(sampleID = rownames(corrdata), class = "sample", corrdata)
corrdata$class[grepl(pattern = "QC", x = corrdata$sampleID)] <- "QC"

proc_object$scaled <- corrdata

proc_object$scaled[,-c(1,2)] <- scale(proc_object$scaled[,-c(1,2)], 
                                      center = TRUE, scale = TRUE)

p.loess <- autoplot(prcomp(proc_object$scaled[,-c(1,2)], center = FALSE, scale. = FALSE), 
                    data = proc_object$scaled, colour = "class", size = 2) 

pca_pos_loess <- p.loess +
  theme_bw() + scale_color_brewer(palette = "Set1") + labs(title = "QC-RSC correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))

feats.to.keep <- CV_filtering(proc_object = proc_object)
feats.to.keep <- names(feats.to.keep)[feats.to.keep]

feats.afterCV <- length(feats.to.keep)

feats.df <- data.frame(Feats0 = feats.before80, feats80 = feats.after80, featsCV = feats.afterCV)

proc_object$scaled <- proc_object$scaled[,c("sampleID", "class", feats.to.keep)]

rownames(db) <- db$ID.sample
merge.pos.loess <- merge(db,proc_object$scaled[,-c(1,2)], by = 0)
rownames(merge.pos.loess) <- merge.pos.loess$Row.names
merge.pos.loess <- merge.pos.loess[,-c(1,2)]

######################## Statistical analysis ##################################

merge.pos.loess <- merge.pos.loess[!is.na(merge.pos.loess$DM_seguimiento),]

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.pos.loess[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos.loess), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

# incident diabetis using all subjects
pls_result <- pls_stats(x = x, y = merge.pos.loess$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

pls_result$plot

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos.loess, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos.loess), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats) 

volcano_plot(stats_results = stats_results, pls_result = pls_result)

#write.csv(stats_results$Pvals.df, "~/fis2018-s1/Results_discovery/R11_pos.csv")

stats_results <- statistical_analysis(df = merge.pos.loess, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15, psign = 0.05)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r1.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R1_pos_loess.csv")

# incident diabetis with only prediabetics at t=0
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos.pred <- merge.pos.loess[merge.pos.loess$Prediabetes==1,]
xpred <- merge.pos.pred[,c(covars,grep(pattern = "SOI", x = colnames(merge.pos.pred), value = TRUE))]
xpred$Sex <- as.numeric(as.factor(xpred$Sex))
xpred$HT <- as.numeric(as.factor(xpred$HT))

pls_result <- pls_stats(x = xpred, y = merge.pos.pred$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos.pred, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos.pred), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos.pred, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)


t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r2.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R2_pos_loess.csv")

# Incident diabetes status
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos2 <- merge.pos.loess[!is.na(merge.pos.loess$DM_seguimiento),]
x2 <- merge.pos2[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

pls_result <- pls_stats(x = x2, y = merge.pos2$DM_seguimiento, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r3.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R3_pos_loess.csv")

# incident diabetes without prediabetes at the end
merge.pos2 <- merge.pos.loess[merge.pos.loess$DM_seguimiento %in% c(0,2),
                              c(covars,"DM_seguimiento", "ID.sample", grep(pattern = "SOI", x = colnames(merge.pos.loess), value = TRUE))]

x2 <- merge.pos2[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

merge.pos2$DM_seguimiento[merge.pos2$DM_seguimiento==2] <- 1

pls_result <- pls_stats(x = x2, 
                        y = merge.pos2$DM_seguimiento, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos2, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos2, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r4.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R4_pos_loess.csv")

# Glycemia change

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")

# glycemia change
# normoglycemia at t=0
merge.pos2 <- merge.pos.loess[!is.na(merge.pos.loess$DM_seguimiento),]
merge.pos2$Glyc_change <- paste(merge.pos2$Prediabetes, merge.pos2$DM_seguimiento, sep = "->")
merge.pos2$Glyc_change <- as.factor(merge.pos2$Glyc_change)
table(merge.pos2$Glyc_change)

merge.pos20 <- merge.pos2[merge.pos2$Glyc_change %in% c("0->0", "0->1", "0->2"),]
merge.pos20$Glyc_change <- as.factor(as.character(merge.pos20$Glyc_change))
levels(merge.pos20$Glyc_change) <- c("0", "1", "2")
merge.pos20$Glyc_change <- as.numeric(as.character(merge.pos20$Glyc_change))

x2 <- merge.pos20[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos20), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.pos20$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos20, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos20), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos20, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r5.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R5_pos_loess.csv")

# prediabetes at t=0
merge.pos21 <- merge.pos2[merge.pos2$Glyc_change %in% c("1->0", "1->1", "1->2"),]
merge.pos21$Glyc_change <- as.factor(as.character(merge.pos21$Glyc_change))
levels(merge.pos21$Glyc_change) <- c("0", "1", "2")
merge.pos21$Glyc_change <- as.numeric(as.character(merge.pos21$Glyc_change))

x2 <- merge.pos21[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos21), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.pos21$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos21, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos21), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos21, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r6.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R6_pos_loess.csv")

# Glycemia at follow-up
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

stats_results <- statistical_analysis_lm(df = merge.pos.loess, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos.loess), value = TRUE), 
                                         covars = covars, outcome = "glucosa_FU",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

r7.pos.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R7_pos_loess.csv")


###################################### NONE ####################################

proc_object <- common_proc_object

proc_object$scaled <- proc_object$log_transf
proc_object$scaled[,-c(1,2)] <- scale(proc_object$scaled[,-c(1,2)], 
                                      center = TRUE, scale = TRUE)

feats.to.keep <- CV_filtering(proc_object = proc_object)
feats.to.keep <- names(feats.to.keep)[feats.to.keep]

feats.afterCV <- length(feats.to.keep)

feats.df <- data.frame(Feats0 = feats.before80, feats80 = feats.after80, featsCV = feats.afterCV)

proc_object$scaled <- proc_object$scaled[,c("sampleID", "class", feats.to.keep)]

merge.pos.none <- merge(db,proc_object$scaled[,-c(1,2)], by = 0)
rownames(merge.pos.none) <- merge.pos.none$Row.names
merge.pos.none <- merge.pos.none[,-c(1,2)]

######################## Statistical analysis ##################################

merge.pos.none <- merge.pos.none[!is.na(merge.pos.none$DM_seguimiento),]

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.pos.none[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos.none), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

# incident diabetis using all subjects
pls_result <- pls_stats(x = x, y = merge.pos.none$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

pls_result$plot

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos.none, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos.none), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats) 

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos.none, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15, psign = 0.05)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r1.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R1_pos_none.csv")

# incident diabetis with only prediabetics at t=0
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos.pred <- merge.pos.none[merge.pos.none$Prediabetes==1,]
xpred <- merge.pos.pred[,c(covars,grep(pattern = "SOI", x = colnames(merge.pos.pred), value = TRUE))]
xpred$Sex <- as.numeric(as.factor(xpred$Sex))
xpred$HT <- as.numeric(as.factor(xpred$HT))

pls_result <- pls_stats(x = xpred, y = merge.pos.pred$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos.pred, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos.pred), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos.pred, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)


t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r2.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R2_pos_none.csv")

# Incident diabetes status
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos2 <- merge.pos.none[!is.na(merge.pos.none$DM_seguimiento),]
x2 <- merge.pos2[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

pls_result <- pls_stats(x = x2, y = merge.pos2$DM_seguimiento, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r3.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R3_pos_none.csv")

# incident diabetes without prediabetes at the end
merge.pos2 <- merge.pos.none[merge.pos.none$DM_seguimiento %in% c(0,2),
                             c(covars,"DM_seguimiento", "ID.sample", grep(pattern = "SOI", x = colnames(merge.pos.none), value = TRUE))]

x2 <- merge.pos2[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

merge.pos2$DM_seguimiento[merge.pos2$DM_seguimiento==2] <- 1

pls_result <- pls_stats(x = x2, 
                        y = merge.pos2$DM_seguimiento, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.pos2, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.pos2, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r4.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R4_pos_none.csv")

# Glycemia change
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")

# glycemia change
# normoglycemia at t=0
merge.pos2 <- merge.pos.none[!is.na(merge.pos.none$DM_seguimiento),]
merge.pos2$Glyc_change <- paste(merge.pos2$Prediabetes, merge.pos2$DM_seguimiento, sep = "->")
merge.pos2$Glyc_change <- as.factor(merge.pos2$Glyc_change)
table(merge.pos2$Glyc_change)

merge.pos20 <- merge.pos2[merge.pos2$Glyc_change %in% c("0->0", "0->1", "0->2"),]
merge.pos20$Glyc_change <- as.factor(as.character(merge.pos20$Glyc_change))
levels(merge.pos20$Glyc_change) <- c("0", "1", "2")
merge.pos20$Glyc_change <- as.numeric(as.character(merge.pos20$Glyc_change))

x2 <- merge.pos20[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos20), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.pos20$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos20, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos20), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)
t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos20, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r5.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R5_pos_none.csv")

# prediabetes at t=0
merge.pos21 <- merge.pos2[merge.pos2$Glyc_change %in% c("1->0", "1->1", "1->2"),]
merge.pos21$Glyc_change <- as.factor(as.character(merge.pos21$Glyc_change))
levels(merge.pos21$Glyc_change) <- c("0", "1", "2")
merge.pos21$Glyc_change <- as.numeric(as.character(merge.pos21$Glyc_change))

x2 <- merge.pos21[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos21), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.pos21$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.pos21, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos21), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.pos21, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r6.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R6_pos_none.csv")

# Glycemia at follow-up
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

stats_results <- statistical_analysis_lm(df = merge.pos.none, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos.none), value = TRUE), 
                                         covars = covars, outcome = "glucosa_FU",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

r7.pos.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R7_pos_none.csv")

################################################################################
############################## NEGATIVE IONIZATION #############################
################################################################################

# data loading - positive ionization
MS_neg <- read.csv("data/DM_FIS2018_Hilic_neg_results2023_filled.csv")
rownames(MS_neg) <- paste0("SOI", MS_neg$Qidx)
MS_neg_t <- as.data.frame(t(MS_neg[,-c(1:9)]))
MS_neg_t <- MS_neg_t[!(grepl(pattern = "Bi", rownames(MS_neg_t)) | grepl(pattern = "Bf", rownames(MS_neg_t))),]

# data loading - injection order
inj_order_neg <- read.csv("data/DidacMauricio_hilic_neg_injectionorder.csv")
rownames(inj_order_neg) <- inj_order_neg$SampleName
inj_order_neg <- inj_order_neg[,-1]

# list of SOIs
feat.cols <- grep(pattern = "SOI", colnames(MS_neg_t), value = TRUE)
feats.before80 <- length(feat.cols)

# 80% rule filter and imputation with minimum/2
proc_object <- filtering_imputing(df = MS_neg_t, features.names = feat.cols, 
                                  imp.method = "min/2")

# number of features after filtering
feats.after80 <- length(grep(pattern = "SOI", x = colnames(proc_object$Imputed), value = TRUE))

# log-transformation
proc_object$log_transf <- proc_object$Imputed

for (j in 1:ncol(proc_object$log_transf)){
  proc_object$log_transf[,j] <- log10(proc_object$log_transf[,j])
}

res_prep_viz <- data_prep_vis(proc_object = proc_object, inj_order = inj_order_neg)
pca.initial.neg <- ggpubr::ggarrange(res_prep_viz$pca_data, 
                                     res_prep_viz$pca_inj, 
                                     ncol = 2)

pca_neg_none <- res_prep_viz$pca_data + labs(title = "No correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))


# outliers removals
data_log <- proc_object$log_transf
outlier_list <- plyr::llply(1:ncol(data_log), function(j){
  idx.rows <- which(data_log[,j] %in% boxplot.stats(data_log[,j])$out)
  out.samps <- rownames(data_log)[idx.rows]
  return(out.samps)
}) %>% unlist()

freq_out_samps <- as.data.frame(table(outlier_list))
freq_out_samps$prop <- freq_out_samps$Freq/feats.before80
out_samps <- as.character(freq_out_samps$outlier_list[freq_out_samps$prop>0.2])

data_log <- data_log[!(rownames(data_log) %in% out_samps),]
proc_object$log_transf <- data_log

res_prep_viz <- data_prep_vis(proc_object = proc_object, inj_order = inj_order_neg)
pca.initial.neg.nooutliers <- ggpubr::ggarrange(res_prep_viz$pca_data, 
                                                res_prep_viz$pca_inj, 
                                                ncol = 2)

pca_neg_none <- res_prep_viz$pca_data + labs(title = "No correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))

proc_object <- res_prep_viz$proc_object
common_proc_object <- proc_object

###################################### CPCA ####################################

normalization_ns_pcas <- plyr::llply(seq_len(4), function(n){
  normalization_res <- data_normalization(proc_object = proc_object, ncomps = n)
  pca.data <- prcomp(normalization_res$norm_data[,-c(1,2)], center = TRUE, scale = TRUE)
  pcs <- pca.data$x[grepl("QC", rownames(pca.data$x)),c(1,2)]
  centroid_coords <- rearrr::centroid(pcs[,1], pcs[,2])
  dist.centroid <- unlist(plyr::llply(1:nrow(pcs), function(i) euclidean_dist(p1 = pcs[i,], p2 = centroid_coords)))
  
  df.dispersion <- data.frame(ncomps = n, dist.centroid = mean(dist.centroid))
  
  p <- normalization_res$pca_norm + 
    labs(title = paste0("Number of CPCs = ",n, " (dist=",round(mean(dist.centroid),2),")")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          legend.position = "bottom", legend.title = element_blank(), 
          legend.text = element_text(size = 12))
  list(df.dispersion = df.dispersion, pca = p)
})

df.dispersion.all <- ldply(normalization_ns_pcas, function(n) n$df.dispersion)
pcas_norms <- llply(normalization_ns_pcas, function(n) n$pca)

pcas_ncomps <- ggpubr::ggarrange(plotlist = pcas_norms, nrow = 2, 
                                 ncol = 2, common.legend = TRUE, legend = "bottom") +
  theme(plot.title = element_text(size = 10))

ggsave(plot=pcas_ncomps, filename = "figures/PCA_cpca_ns_neg.jpeg", 
       width = 10, height = 8, dpi = 500)

# looking at the PCAs, we select the first n cpcs where technical bias is 
# suposedly removed in the first two PCs (n = 4)
normalization_res <- data_normalization(proc_object = proc_object, ncomps = 4)
normalization_res$pca_norm

pca_neg_cpca <- normalization_res$pca_norm + labs(title = "CPCA correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))

proc_object$log_transf <- normalization_res$norm_data
proc_object$scaled <- proc_object$log_transf
proc_object$scaled[,-c(1,2)] <- scale(proc_object$scaled[,-c(1,2)], 
                                      center = TRUE, scale = TRUE)

feats.to.keep <- CV_filtering(proc_object = proc_object)
feats.to.keep <- names(feats.to.keep)[feats.to.keep]

feats.afterCV <- length(feats.to.keep)

feats.df <- data.frame(Feats0 = feats.before80, feats80 = feats.after80, featsCV = feats.afterCV)

proc_object$scaled <- proc_object$scaled[,c("sampleID", "class", feats.to.keep)]

scaled_dat <- proc_object$scaled

merge.neg.cpca <- merge(db,proc_object$scaled[,-c(1,2)], by = 0)
rownames(merge.neg.cpca) <- merge.neg.cpca$Row.names
merge.neg.cpca <- merge.neg.cpca[,-c(1,2)]

############################# Statistical analysis #############################

merge.neg.cpca <- merge.neg.cpca[!is.na(merge.neg.cpca$DM_seguimiento),]

write.csv(merge.neg.cpca, "processed_files/cpca_processed_neg.csv")

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.neg.cpca[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg.cpca), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

# incident diabetis using all subjects
pls_result <- pls_stats(x = x, y = merge.neg.cpca$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg.cpca, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg.cpca), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats) 
volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg.cpca, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15, psign = 0.05)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r1.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R1_neg_cpca.csv")

# incident diabetis with only prediabetics at t=0
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg.pred <- merge.neg.cpca[merge.neg.cpca$Prediabetes==1,]
xpred <- merge.neg.pred[,c(covars,grep(pattern = "SOI", x = colnames(merge.neg.pred), value = TRUE))]
xpred$Sex <- as.numeric(as.factor(xpred$Sex))
xpred$HT <- as.numeric(as.factor(xpred$HT))

pls_result <- pls_stats(x = xpred, y = merge.neg.pred$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg.pred, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg.pred), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg.pred, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r2.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R2_neg_cpca.csv")

# Incident diabetes status
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg2 <- merge.neg.cpca[!is.na(merge.neg.cpca$DM_seguimiento),]
x2 <- merge.neg2[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

pls_result <- pls_stats(x = x2, y = merge.neg2$DM_seguimiento, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r3.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R3_neg_cpca.csv")

# incident diabetes without prediabetes at the end
merge.neg2 <- merge.neg.cpca[merge.neg.cpca$DM_seguimiento %in% c(0,2),
                             c(covars,"DM_seguimiento", "ID.sample", grep(pattern = "SOI", x = colnames(merge.neg.cpca), value = TRUE))]

x2 <- merge.neg2[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

merge.neg2$DM_seguimiento[merge.neg2$DM_seguimiento==2] <- 1

pls_result <- pls_stats(x = x2, 
                        y = merge.neg2$DM_seguimiento, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg2, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg2, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r4.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R4_neg_cpca.csv")

# Glycemia change

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")

# glycemia change
# normoglycemia at t=0
merge.neg2 <- merge.neg.cpca[!is.na(merge.neg.cpca$DM_seguimiento),]
merge.neg2$Glyc_change <- paste(merge.neg2$Prediabetes, merge.neg2$DM_seguimiento, sep = "->")
merge.neg2$Glyc_change <- as.factor(merge.neg2$Glyc_change)
table(merge.neg2$Glyc_change)

merge.neg20 <- merge.neg2[merge.neg2$Glyc_change %in% c("0->0", "0->1", "0->2"),]
merge.neg20$Glyc_change <- as.factor(as.character(merge.neg20$Glyc_change))
levels(merge.neg20$Glyc_change) <- c("0", "1", "2")
merge.neg20$Glyc_change <- as.numeric(as.character(merge.neg20$Glyc_change))

x2 <- merge.neg20[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg20), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.neg20$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg20, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg20), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg20, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r5.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R5_neg_cpca.csv")

# prediabetes at t=0
merge.neg21 <- merge.neg2[merge.neg2$Glyc_change %in% c("1->0", "1->1", "1->2"),]
merge.neg21$Glyc_change <- as.factor(as.character(merge.neg21$Glyc_change))
levels(merge.neg21$Glyc_change) <- c("0", "1", "2")
merge.neg21$Glyc_change <- as.numeric(as.character(merge.neg21$Glyc_change))

x2 <- merge.neg21[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg21), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.neg21$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg21, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg21), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg21, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r6.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R6_neg_cpca.csv")

# Glycemia at follow-up
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

stats_results <- statistical_analysis_lm(df = merge.neg.cpca, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg.cpca), value = TRUE), 
                                         covars = covars, outcome = "glucosa_FU",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

r7.neg.cpca <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R7_neg_cpca.csv")

###################################### LOESS ###################################
proc_object <- common_proc_object

res <- proc_object$log_transf
res <- merge(inj_order_pos, res, by = 0)
rownames(res) <- res$Row.names
res <- res[,-1]
res2 <- res[,-c(1:4)]
res2 <- t(res2)
batch <- rep(1, times=ncol(res2))

sum(is.na(res2))

corrected_data <- pmp::QCRSC(df=res2, order=res$injectionOrder, 
                             batch=batch,
                             log = F, # log = F perque ja hem aplicat el log
                             classes=res$class,
                             spar_lim = c(-2,2),
                             minQC=10)

corrdata <- as.data.frame(t(corrected_data))
sum(is.na(corrdata))
corrdata <- data.frame(sampleID = rownames(corrdata), class = "sample", corrdata)
corrdata$class[grepl(pattern = "QC", x = corrdata$sampleID)] <- "QC"

proc_object$scaled <- corrdata

proc_object$scaled[,-c(1,2)] <- scale(proc_object$scaled[,-c(1,2)], 
                                      center = TRUE, scale = TRUE)

p.loess <- autoplot(prcomp(proc_object$scaled[,-c(1,2)], center = FALSE, scale. = FALSE), 
                    data = proc_object$scaled, colour = "class")

pca_neg_loess <- p.loess +
  theme_bw() + scale_color_brewer(palette = "Set1") + labs(title = "QC-RSC correction") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 13))

feats.to.keep <- CV_filtering(proc_object = proc_object)
feats.to.keep <- names(feats.to.keep)[feats.to.keep]

feats.afterCV <- length(feats.to.keep)

feats.df <- data.frame(Feats0 = feats.before80, feats80 = feats.after80, featsCV = feats.afterCV)

proc_object$scaled <- proc_object$scaled[,c("sampleID", "class", feats.to.keep)]

merge.neg.loess <- merge(db,proc_object$scaled[,-c(1,2)], by = 0)
rownames(merge.neg.loess) <- merge.neg.loess$Row.names
merge.neg.loess <- merge.neg.loess[,-c(1,2)]

######################## Statistical analysis ##################################

merge.neg.loess <- merge.neg.loess[!is.na(merge.neg.loess$DM_seguimiento),]

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.neg.loess[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg.loess), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

# incident diabetis using all subjects
pls_result <- pls_stats(x = x, y = merge.neg.loess$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

pls_result$plot

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg.loess, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg.loess), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats) 

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg.loess, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15, psign = 0.05)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r1.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R1_neg_loess.csv")

# incident diabetis with only prediabetics at t=0
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg.pred <- merge.neg.loess[merge.neg.loess$Prediabetes==1,]
xpred <- merge.neg.pred[,c(covars,grep(pattern = "SOI", x = colnames(merge.neg.pred), value = TRUE))]
xpred$Sex <- as.numeric(as.factor(xpred$Sex))
xpred$HT <- as.numeric(as.factor(xpred$HT))

pls_result <- pls_stats(x = xpred, y = merge.neg.pred$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg.pred, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg.pred), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg.pred, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)


t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r2.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R2_neg_loess.csv")

# Incident diabetes status
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg2 <- merge.neg.loess[!is.na(merge.neg.loess$DM_seguimiento),]
x2 <- merge.neg2[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

pls_result <- pls_stats(x = x2, y = merge.neg2$DM_seguimiento, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r3.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R3_neg_loess.csv")

# incident diabetes without prediabetes at the end
merge.neg2 <- merge.neg.loess[merge.neg.loess$DM_seguimiento %in% c(0,2),
                              c(covars,"DM_seguimiento", "ID.sample", grep(pattern = "SOI", x = colnames(merge.neg.loess), value = TRUE))]

x2 <- merge.neg2[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

merge.neg2$DM_seguimiento[merge.neg2$DM_seguimiento==2] <- 1

pls_result <- pls_stats(x = x2, 
                        y = merge.neg2$DM_seguimiento, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg2, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg2, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r4.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R4_neg_loess.csv")

# Glycemia change

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")

# glycemia change
# normoglycemia at t=0
merge.neg2 <- merge.neg.loess[!is.na(merge.neg.loess$DM_seguimiento),]
merge.neg2$Glyc_change <- paste(merge.neg2$Prediabetes, merge.neg2$DM_seguimiento, sep = "->")
merge.neg2$Glyc_change <- as.factor(merge.neg2$Glyc_change)
table(merge.neg2$Glyc_change)

merge.neg20 <- merge.neg2[merge.neg2$Glyc_change %in% c("0->0", "0->1", "0->2"),]
merge.neg20$Glyc_change <- as.factor(as.character(merge.neg20$Glyc_change))
levels(merge.neg20$Glyc_change) <- c("0", "1", "2")
merge.neg20$Glyc_change <- as.numeric(as.character(merge.neg20$Glyc_change))

x2 <- merge.neg20[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg20), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.neg20$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg20, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg20), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg20, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r5.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R5_neg_loess.csv")

# prediabetes at t=0
merge.neg21 <- merge.neg2[merge.neg2$Glyc_change %in% c("1->0", "1->1", "1->2"),]
merge.neg21$Glyc_change <- as.factor(as.character(merge.neg21$Glyc_change))
levels(merge.neg21$Glyc_change) <- c("0", "1", "2")
merge.neg21$Glyc_change <- as.numeric(as.character(merge.neg21$Glyc_change))

x2 <- merge.neg21[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg21), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.neg21$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg21, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg21), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg21, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r6.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R6_neg_loess.csv")

# Glycemia at follow-up
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

stats_results <- statistical_analysis_lm(df = merge.neg.loess, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg.loess), value = TRUE), 
                                         covars = covars, outcome = "glucosa_FU",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

r7.neg.loess <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R7_neg_loess.csv")

###################################### NONE ####################################

proc_object <- common_proc_object

proc_object$scaled <- proc_object$log_transf
proc_object$scaled[,-c(1,2)] <- scale(proc_object$scaled[,-c(1,2)], 
                                      center = TRUE, scale = TRUE)

feats.to.keep <- CV_filtering(proc_object = proc_object)
feats.to.keep <- names(feats.to.keep)[feats.to.keep]

feats.afterCV <- length(feats.to.keep)

feats.df <- data.frame(Feats0 = feats.before80, feats80 = feats.after80, featsCV = feats.afterCV)

proc_object$scaled <- proc_object$scaled[,c("sampleID", "class", feats.to.keep)]

merge.neg.none <- merge(db,proc_object$scaled[,-c(1,2)], by = 0)
rownames(merge.neg.none) <- merge.neg.none$Row.names
merge.neg.none <- merge.neg.none[,-c(1,2)]


######################## Statistical analysis ##################################

merge.neg.none <- merge.neg.none[!is.na(merge.neg.none$DM_seguimiento),]

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.neg.none[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg.none), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

# incident diabetis using all subjects
pls_result <- pls_stats(x = x, y = merge.neg.none$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

pls_result$plot

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg.none, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg.none), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats) 

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg.none, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 20, psign = 0.05)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r1.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R1_neg_none.csv")

# incident diabetis with only prediabetics at t=0
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg.pred <- merge.neg.none[merge.neg.none$Prediabetes==1,]
xpred <- merge.neg.pred[,c(covars,grep(pattern = "SOI", x = colnames(merge.neg.pred), value = TRUE))]
xpred$Sex <- as.numeric(as.factor(xpred$Sex))
xpred$HT <- as.numeric(as.factor(xpred$HT))

pls_result <- pls_stats(x = xpred, y = merge.neg.pred$DM_INCIDENTE, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg.pred, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg.pred), value = TRUE), 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg.pred, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_INCIDENTE",
                                      do.Par = TRUE, nCore = 15)


t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r2.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R2_neg_none.csv")

# Incident diabetes status
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg2 <- merge.neg.none[!is.na(merge.neg.none$DM_seguimiento),]
x2 <- merge.neg2[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

pls_result <- pls_stats(x = x2, y = merge.neg2$DM_seguimiento, ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "DM_seguimiento",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r3.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R3_neg_none.csv")

# incident diabetes without prediabetes at the end
merge.neg2 <- merge.neg.none[merge.neg.none$DM_seguimiento %in% c(0,2),
                             c(covars,"DM_seguimiento", "ID.sample", grep(pattern = "SOI", x = colnames(merge.neg.none), value = TRUE))]

x2 <- merge.neg2[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))

merge.neg2$DM_seguimiento[merge.neg2$DM_seguimiento==2] <- 1

pls_result <- pls_stats(x = x2, 
                        y = merge.neg2$DM_seguimiento, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis(df = merge.neg2, 
                                      feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis(df = merge.neg2, 
                                      feats.names = vip.feats, 
                                      covars = covars, outcome = "DM_seguimiento",
                                      do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r4.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R4_neg_none.csv")

# Glycemia change

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")

# glycemia change
# normoglycemia at t=0
merge.neg2 <- merge.neg.none[!is.na(merge.neg.none$DM_seguimiento),]
merge.neg2$Glyc_change <- paste(merge.neg2$Prediabetes, merge.neg2$DM_seguimiento, sep = "->")
merge.neg2$Glyc_change <- as.factor(merge.neg2$Glyc_change)
table(merge.neg2$Glyc_change)

merge.neg20 <- merge.neg2[merge.neg2$Glyc_change %in% c("0->0", "0->1", "0->2"),]
merge.neg20$Glyc_change <- as.factor(as.character(merge.neg20$Glyc_change))
levels(merge.neg20$Glyc_change) <- c("0", "1", "2")
merge.neg20$Glyc_change <- as.numeric(as.character(merge.neg20$Glyc_change))

x2 <- merge.neg20[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg20), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.neg20$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg20, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg20), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg20, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r5.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R5_neg_none.csv")

# prediabetes at t=0
merge.neg21 <- merge.neg2[merge.neg2$Glyc_change %in% c("1->0", "1->1", "1->2"),]
merge.neg21$Glyc_change <- as.factor(as.character(merge.neg21$Glyc_change))
levels(merge.neg21$Glyc_change) <- c("0", "1", "2")
merge.neg21$Glyc_change <- as.numeric(as.character(merge.neg21$Glyc_change))

x2 <- merge.neg21[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg21), value = TRUE))]
x2$Sex <- as.numeric(as.factor(x2$Sex))
x2$HT <- as.numeric(as.factor(x2$HT))
pls_result <- pls_stats(x = x2, 
                        y = merge.neg21$Glyc_change, 
                        ncomp = 10, 
                        plot.title = "All subjects", ncomp.perf = 3)

vip.feats <- pls_result$vip_df$Feats[pls_result$vip_df$comp>=1]
vip.feats <- grep(pattern = "SOI", x = vip.feats, value = TRUE)

stats_results <- statistical_analysis_lm(df = merge.neg21, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg21), value = TRUE), 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

stats_results <- statistical_analysis_lm(df = merge.neg21, 
                                         feats.names = vip.feats, 
                                         covars = covars, outcome = "Glyc_change",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

volcano_plot(stats_results = stats_results, pls_result = pls_result)

r6.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)

write.csv(stats_results$Pvals.df, "results/discovery/R6_neg_none.csv")

# Glycemia at follow-up
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

stats_results <- statistical_analysis_lm(df = merge.neg.none, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg.none), value = TRUE), 
                                         covars = covars, outcome = "glucosa_FU",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)


r7.neg.none <- sum(stats_results$Pvals.df$qvalue_trunc<0.05)
write.csv(stats_results$Pvals.df, "results/discovery/R7_neg_none.csv")


number.sign.feats <- data.frame(R1.pos.cpca = r1.pos.cpca, R2.pos.cpca = r2.pos.cpca,
                                R3.pos.cpca = r3.pos.cpca, R4.pos.cpca = r4.pos.cpca,
                                R5.pos.cpca = r5.pos.cpca, R6.pos.cpca = r6.pos.cpca,
                                R7.pos.cpca = r7.pos.cpca,
                                R1.pos.loess = r1.pos.loess, R2.pos.loess = r2.pos.loess,
                                R3.pos.loess = r3.pos.loess, R4.pos.loess = r4.pos.loess,
                                R5.pos.loess = r5.pos.loess, R6.pos.loess = r6.pos.loess,
                                R7.pos.loess = r7.pos.loess,
                                R1.pos.none = r1.pos.none, R2.pos.none = r2.pos.none,
                                R3.pos.none = r3.pos.none, R4.pos.none = r4.pos.none,
                                R5.pos.none = r5.pos.none, R6.pos.none = r6.pos.none,
                                R7.pos.none = r7.pos.none,
                                R1.neg.cpca = r1.neg.cpca, R2.neg.cpca = r2.neg.cpca,
                                R3.neg.cpca = r3.neg.cpca, R4.neg.cpca = r4.neg.cpca,
                                R5.neg.cpca = r5.neg.cpca, R6.neg.cpca = r6.neg.cpca,
                                R7.neg.cpca = r7.neg.cpca,
                                R1.neg.loess = r1.neg.loess, R2.neg.loess = r2.neg.loess,
                                R3.neg.loess = r3.neg.loess, R4.neg.loess = r4.neg.loess,
                                R5.neg.loess = r5.neg.loess, R6.neg.loess = r6.neg.loess,
                                R7.neg.loess = r7.neg.loess,
                                R1.neg.none = r1.neg.none, R2.neg.none = r2.neg.none,
                                R3.neg.none = r3.neg.none, R4.neg.none = r4.neg.none,
                                R5.neg.none = r5.neg.none, R6.neg.none = r6.neg.none,
                                R7.neg.none = r7.neg.none)

pcas_all <- ggpubr::ggarrange(pca_pos_none, pca_neg_none, 
                              pca_pos_cpca, pca_neg_cpca, 
                              pca_pos_loess, pca_neg_loess,
                              common.legend = TRUE, 
                              labels = c("A", "B", "", "", "", ""),
                              nrow = 3, ncol = 2, legend = "bottom")

ggsave(plot=pcas_all, filename = "figures/allPCAs.jpeg", 
       width = 10, height = 12, dpi = 500)


