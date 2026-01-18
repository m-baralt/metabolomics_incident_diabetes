library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)
library(magrittr)
source("Functions_validation.R")

files_results <- list.files("results/discovery/", full.names = TRUE)

files.pos <- files_results[grep("pos", files_results)]
files.neg <- files_results[grep("neg", files_results)]

### RESULTS WITH CPCA CORRECTION
files.pos.cpca <- grep(pattern = "cpca", x = files.pos, value = TRUE)

cpca_pos <- plyr::ldply(files.pos.cpca, .inform = TRUE, function(f){
  name.analysis <- gsub(pattern = ".csv", replacement = "", x = basename(f))
  df <- read.csv(f)
  analysis <- unlist(strsplit(name.analysis, split = "_"))[1]
  
  if (nrow(df)>0){
    if (sum(grepl(pattern = "OR", x = colnames(df)))==1){
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, OR = df$OR, CI2.5 = df$CI2.5, 
                 CI97.5 = df$CI97.5, Analysis = analysis, Ionization = "Positive")
    } else {
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, Beta = df$beta, 
                 Analysis = analysis, Ionization = "Positive")
    }
  }
})

cpca_pos$correction <- "CPCA"

files.neg.cpca <- grep(pattern = "cpca", x = files.neg, value = TRUE)

cpca_neg <- plyr::ldply(files.neg.cpca, .inform = TRUE, function(f){
  name.analysis <- gsub(pattern = ".csv", replacement = "", x = basename(f))
  df <- read.csv(f)
  analysis <- unlist(strsplit(name.analysis, split = "_"))[1]
  
  if (nrow(df)>0){
    if (sum(grepl(pattern = "OR", x = colnames(df)))==1){
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, OR = df$OR, CI2.5 = df$CI2.5, 
                 CI97.5 = df$CI97.5, Analysis = analysis, Ionization = "Negative")
    } else {
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, Beta = df$beta, 
                 Analysis = analysis, Ionization = "Negative")
    }
  }
  
})

cpca_neg$correction <- "CPCA"

## RESULTS WITHOUT CORRECTION
files.pos.none <- grep(pattern = "none", x = files.pos, value = TRUE)
none_pos <- plyr::ldply(files.pos.none, .inform = TRUE, function(f){
  name.analysis <- gsub(pattern = ".csv", replacement = "", x = basename(f))
  df <- read.csv(f)
  analysis <- unlist(strsplit(name.analysis, split = "_"))[1]
  
  if (nrow(df)>0){
    if (sum(grepl(pattern = "OR", x = colnames(df)))==1){
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, OR = df$OR, CI2.5 = df$CI2.5, 
                 CI97.5 = df$CI97.5, Analysis = analysis, Ionization = "Positive")
    } else {
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, Beta = df$beta, 
                 Analysis = analysis, Ionization = "Positive")
    }
  }
})

none_pos$correction <- "None"

files.neg.none <- grep(pattern = "none", x = files.neg, value = TRUE)

none_neg <- plyr::ldply(files.neg.none, .inform = TRUE, function(f){
  name.analysis <- gsub(pattern = ".csv", replacement = "", x = basename(f))
  df <- read.csv(f)
  analysis <- unlist(strsplit(name.analysis, split = "_"))[1]
  
  if (nrow(df)>0){
    if (sum(grepl(pattern = "OR", x = colnames(df)))==1){
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, OR = df$OR, CI2.5 = df$CI2.5, 
                 CI97.5 = df$CI97.5, Analysis = analysis, Ionization = "Negative")
    } else {
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, Beta = df$beta, 
                 Analysis = analysis, Ionization = "Negative")
    }
  }
})

none_neg$correction <- "None"

## RESULTS WITH LOESS CORRECTION
files.pos.loess <- grep(pattern = "loess", x = files.pos, value = TRUE)

loess_pos <- plyr::ldply(files.pos.loess, .inform = TRUE, function(f){
  name.analysis <- gsub(pattern = ".csv", replacement = "", x = basename(f))
  df <- read.csv(f)
  analysis <- unlist(strsplit(name.analysis, split = "_"))[1]
  
  if (nrow(df)>0){
    if (sum(grepl(pattern = "OR", x = colnames(df)))==1){
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, OR = df$OR, CI2.5 = df$CI2.5, 
                 CI97.5 = df$CI97.5, Analysis = analysis, Ionization = "Positive")
    } else {
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, Beta = df$beta, 
                 Analysis = analysis, Ionization = "Positive")
    }
  }
})

loess_pos$correction <- "QC-RSC"

files.neg.loess <- grep(pattern = "loess", x = files.neg, value = TRUE)
loess_neg <- plyr::ldply(files.neg.loess, .inform = TRUE, function(f){
  name.analysis <- gsub(pattern = ".csv", replacement = "", x = basename(f))
  df <- read.csv(f)
  analysis <- unlist(strsplit(name.analysis, split = "_"))[1]
  
  if (nrow(df)>0){
    if (sum(grepl(pattern = "OR", x = colnames(df)))==1){
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, OR = df$OR, CI2.5 = df$CI2.5, 
                 CI97.5 = df$CI97.5, Analysis = analysis, Ionization = "Negative")
    } else {
      data.frame(Feature = df$Feature, qvalue = df$qvalue_trunc, Beta = df$beta, 
                 Analysis = analysis, Ionization = "Negative")
    }
  }
  
})

loess_neg$correction <- "QC-RSC"
table_pos_all_corrections <- rbind(cpca_pos,loess_pos,none_pos)
table_neg_all_corrections <- rbind(cpca_neg, loess_neg, none_neg)

MS_pos <- read.csv("data/DM_FIS2018_Hilic_pos_results2023_filled.csv")
rownames(MS_pos) <- paste0("SOI", MS_pos$Qidx)

MS_neg <- read.csv("data/DM_FIS2018_Hilic_neg_results2023_filled.csv")
rownames(MS_neg) <- paste0("SOI", MS_neg$Qidx)

pos_info <- MS_pos[,c(1:9)]
pos_info <- data.frame(Feature = rownames(pos_info), pos_info)
merge_table_pos_all_corr <- merge(pos_info, table_pos_all_corrections, by = "Feature")

neg_info <- MS_neg[,c(1:9)]
neg_info <- data.frame(Feature = rownames(neg_info), neg_info)
merge_table_neg_all_corr <- merge(neg_info, table_neg_all_corrections, by = "Feature")

cmp_validation <- readxl::read_excel("data/Compostos_validació_DM.xlsx")
cmp_validation <- cmp_validation[1:6,1:5]
cmp_validation <- plyr::ldply(1:nrow(cmp_validation), function(i){
  dx <- cmp_validation[i,]
  qidxs <- unlist(strsplit(cmp_validation$Qidx[i], split = "[/]"))
  dx <- cbind(dx,rep(row.names(dx), each = length(qidxs)))
  dx$Qidx <- unlist(strsplit(cmp_validation$Qidx[i], split = "[/]"))
  dx[,1:5]
})

colnames(cmp_validation)[1] <- "compound"
cmp_validation_pos <- cmp_validation[cmp_validation$Polaritat %in% "Positiu",]
cmp_validation_neg <- cmp_validation[cmp_validation$Polaritat %in% "Negatiu",]

cmp_features_pos <- merge(cmp_validation_pos, 
                          merge_table_pos_all_corr[merge_table_pos_all_corr$Qidx %in% cmp_validation_pos$Qidx,], 
                          by = "Qidx")
cmp_features_neg <- merge(cmp_validation_neg,
                          merge_table_neg_all_corr[merge_table_neg_all_corr$Qidx %in% cmp_validation_neg$Qidx,],
                          by = "Qidx")

total_cmp_features <- rbind(cmp_features_pos, cmp_features_neg)
total_cmp_features <- total_cmp_features[!(total_cmp_features$Analysis %in% "R7"),]
total_cmp_features$compound <- paste(total_cmp_features$Qidx, total_cmp_features$compound, sep = " - ")

## R1, R2, R4
glm_cmp_features <- total_cmp_features[total_cmp_features$Analysis %in% c("R1", "R2", "R4"),]

# Complete all combinations of compound, correction, and Analysis, filling missing values with 0
glm_cmp_features <- glm_cmp_features %>%
  group_by(Analysis) %>%
  complete(compound, correction, fill = list(OR = 0))

glm_cmp_features$signif <- ""
glm_cmp_features$signif[glm_cmp_features$qvalue<0.05] <- "*"

glm_cmp_features <- glm_cmp_features[,c("compound", "correction", "OR", 
                                          "CI2.5", "CI97.5", "Analysis", "signif")]
glm_cmp_features$compound <- as.factor(glm_cmp_features$compound)

glm_cmp_features$correction <- as.factor(glm_cmp_features$correction)
levels(glm_cmp_features$correction) <- paste0("Discovery - ", glm_cmp_features$correction)

glm_validation <- read.csv("processed_files/validation_results_glm.csv")
glm_validation$correction <- "Validation"
glm_validation$signif <- ""
glm_validation$signif[glm_validation$pval<0.05] <- "*"
glm_validation$signif[glm_validation$pval<0.01] <- "**"
glm_validation$signif[glm_validation$pval<0.001] <- "***"
glm_validation$signif[glm_validation$pval<0.0001] <- "****"

glm_validation <- glm_validation[,c("Compound", "correction", "OR", 
                                    "CI2.5", "CI97.5", "analysis", "signif")]
colnames(glm_validation)[c(1,6)] <- c("compound", "Analysis")


glm_validation$compound <- as.factor(glm_validation$compound)

levels(glm_validation$compound) <- c("139 - Adenine", "449 - Citrulline", "799 - Ecgonine",
                                     "445 - Guanine", "430 - Phenyl sulfate", "1881 - Pregnenolone sulfate")

ecgonine_ids <- grep(pattern = "Ecgonine", x = glm_validation$compound)
guanine_ids <- grep(pattern = "Guanine", x = glm_validation$compound)
rep_df <- glm_validation[rep(c(ecgonine_ids, guanine_ids), each = 1), ]
rep_df$compound <- as.factor(as.character(rep_df$compound))
levels(rep_df$compound) <- c("447 - Guanine", "801 - Ecgonine")

glm_validation <- rbind(glm_validation, rep_df)

glm_cmp_features <- rbind(glm_cmp_features, glm_validation)

glm_cmp_features <- glm_cmp_features %>%
  group_by(Analysis) %>%
  complete(compound, correction, fill = list(OR = 0))

glm_cmp_features$signif[is.na(glm_cmp_features$signif)] <- ""

glm_plot <- ggplot(data=glm_cmp_features, aes(x=compound, y=OR, fill = correction)) + 
  geom_bar(stat="identity", width=0.6, alpha = 0.7, position=position_dodge(0.6)) + 
  geom_errorbar(aes(ymin = CI2.5, ymax = CI97.5), color = "gray10",
                width=0.2, position=position_dodge(0.6), alpha=0.8)+
  facet_grid(~Analysis)+
  scale_fill_manual(values = c("#E66100", "#5E3C99", "goldenrod1", "#1B9E77")) +
  theme_classic()+
  coord_flip() +
  theme(legend.position = c("bottom"), axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 13), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 13), 
        legend.text = element_text(size = 13), legend.title = element_blank(), axis.title.x = element_text(size = 14),
        plot.background = element_rect(fill = "transparent", color = NA))+
  geom_hline(yintercept = 1, linetype="dashed", color="firebrick")+
  geom_text(data = glm_cmp_features, aes(x=compound, y = CI97.5, label = signif), color = "gray10",
            vjust = 0.8, hjust = -0.5,size = 4.5, position=position_dodge(0.6)) +
  scale_y_continuous(limits = c(0, max(glm_cmp_features$CI97.5, na.rm = T)+0.1))


## R3, R5, R6
lm_cmp_features <- total_cmp_features[total_cmp_features$Analysis %in% c("R3", "R5", "R6"),]

# Complete all combinations of compound, correction, and Analysis, filling missing values with 0
lm_cmp_features <- lm_cmp_features %>%
  group_by(Analysis) %>%
  complete(compound, correction, fill = list(Beta = 0))

lm_cmp_features$signif <- ""
lm_cmp_features$signif[lm_cmp_features$qvalue<0.05] <- "*"

lm_cmp_features <- lm_cmp_features[,c("compound", "correction", "Beta", "Analysis", "signif")]
lm_cmp_features$compound <- as.factor(lm_cmp_features$compound)

lm_cmp_features$correction <- as.factor(lm_cmp_features$correction)
levels(lm_cmp_features$correction) <- paste0("Discovery - ", lm_cmp_features$correction)

lm_validation <- read.csv("processed_files/validation_results_lm.csv")
lm_validation$correction <- "Validation"
lm_validation$signif <- ""
lm_validation$signif[lm_validation$pval<0.05] <- "*"
lm_validation$signif[lm_validation$pval<0.01] <- "**"
lm_validation$signif[lm_validation$pval<0.001] <- "***"
lm_validation$signif[lm_validation$pval<0.0001] <- "****"

lm_validation <- lm_validation[,c("Compound", "correction", "beta", "analysis", "signif")]
colnames(lm_validation)[c(1,3,4)] <- c("compound", "Beta", "Analysis")

lm_validation$compound <- as.factor(lm_validation$compound)

levels(lm_validation$compound) <- c("139 - Adenine", "449 - Citrulline", "799 - Ecgonine",
                                     "445 - Guanine", "430 - Phenyl sulfate", "1881 - Pregnenolone sulfate")

ecgonine_ids <- grep(pattern = "Ecgonine", x = lm_validation$compound)
guanine_ids <- grep(pattern = "Guanine", x = lm_validation$compound)
rep_df <- lm_validation[rep(c(ecgonine_ids, guanine_ids), each = 1), ]
rep_df$compound <- as.factor(as.character(rep_df$compound))
levels(rep_df$compound) <- c("447 - Guanine", "801 - Ecgonine")

lm_validation <- rbind(lm_validation, rep_df)

lm_cmp_features <- rbind(lm_cmp_features, lm_validation)

lm_cmp_features <- lm_cmp_features %>%
  group_by(Analysis) %>%
  complete(compound, correction, fill = list(Beta = 0))

lm_cmp_features$signif[is.na(lm_cmp_features$signif)] <- ""

lm_plot <- ggplot(data=lm_cmp_features, aes(x=compound, y=Beta, fill = correction)) + 
  geom_bar(stat="identity", width=0.6, alpha = 0.7, position=position_dodge(0.6)) +
  facet_grid(~Analysis)+
  scale_fill_manual(values = c("#E66100", "#5E3C99", "goldenrod1", "#1B9E77"))+
  theme_classic()+
  coord_flip() +
  theme(legend.position = c("bottom"), axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 13), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 13), 
        legend.text = element_text(size = 13), legend.title = element_blank(), axis.title.x = element_text(size = 14),
        plot.background = element_rect(fill = "transparent", color = NA))+
  geom_hline(yintercept = 0, linetype="dashed", color="firebrick") +
  scale_y_continuous(limits = c(min(lm_cmp_features$Beta, na.rm = T)-0.15, max(lm_cmp_features$Beta, na.rm = T)+0.15)) +
  geom_text(data = lm_cmp_features, aes(x=compound, y = ifelse(Beta<0, Beta + (sign(Beta)*0.12), Beta + (sign(Beta)*0.1)), label = signif), color = "gray10",
            vjust = 0.8, size = 4.5, position=position_dodge(0.6)) + labs(y = "Regression coefficient")


lm_plot <- lm_plot + theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())

association_plot <- ggpubr::ggarrange(glm_plot, NULL, lm_plot, labels = c("A","B", ""),
                  ncol = 3, common.legend = TRUE, legend = "bottom", widths = c(0.54, 0.03, 0.43))

ggsave(association_plot, 
       filename = "figures/Figure2.jpg", 
       width = 16, height = 12, dpi = 900)


db <- read.csv("data/Updated_db.csv", row.names = 1)
matching_ids <- readxl::read_excel("data/match_ids.xlsx")
colnames(matching_ids)[1] <- "subName"
db <- merge(matching_ids[,-3], db, by = "CIP")

object_data <- readRDS("data/all_batches_alldata.RDS")
corrected.data <- object_data$CorrectedData

merge.df <- merge(db, corrected.data, by = "subName")

merge.df.validation <- merge.df[!(merge.df$FU_time<=0),]

merge.df.validation <- merge.df.validation[,c("event.DM2", "Sex", "Age", 
                        "BMI", "HT", "Smoking", "FU_time", "Prediabetes", object_data$metabs.used)]

merge.df.validation$event.DM2 <- as.factor(merge.df.validation$event.DM2)

covars <- c("Sex", "Age", "BMI", "HT", "Smoking", "FU_time", "Prediabetes")
f_regr <- as.formula(paste0("event.DM2~",paste(covars, collapse = "+")))
lm.model <- lm(MZ_395.1898_MZ_135.1208~Age+BMI+HT+Smoking+FU_time+Prediabetes, 
               data = merge.df.validation)
mm_merge <- merge(merge.df.validation, data.frame(lm.model$residuals), by = 0)

model_women <- glm(event.DM2~lm.model.residuals, data = mm_merge[mm_merge$Sex %in% "Women",], family = "binomial")
s_w <- summary(model_women)
pval_w <- s_w$coefficients["lm.model.residuals",4]

model_men <- glm(event.DM2~lm.model.residuals, data = mm_merge[mm_merge$Sex %in% "Men",], family = "binomial")
s_m <- summary(model_men)
pval_m <- s_m$coefficients["lm.model.residuals",4]

my_comparisons <- list( c("", "1"), c("1", "2"), c("0.5", "2") )

library(ggpubr)

p_validation <- ggplot(mm_merge, aes(x = Sex, y = lm.model.residuals))+ 
  geom_point(size=2, alpha = 0.5, aes(color = event.DM2),
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7))+
  geom_hpline(stat = "summary", fun = "mean", width = 0.25, aes(group = event.DM2), 
              position = position_dodge(width = 0.7))+
  theme_bw() +
  labs(color = "T2D event")+
  scale_colour_manual(values = c("#E66100", "#1B9E77")) +
  theme(axis.text.y = element_text(family = "Times New Roman", size = 17),
        axis.text.x = element_text(family = "Times New Roman", size = 17),
        plot.title = element_text(family = "Times New Roman", size = 14, hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.title = element_text(family = "Times New Roman", size = 19), 
        legend.text = element_text(family = "Times New Roman", size = 18)) +
  stat_compare_means(aes(group = event.DM2), label = "p.format", method = "t.test", size = 6, family = "Times New Roman", label.y = 3.5) +
  stat_compare_means(method = "t.test", label.y = 4, label.x = 1.35, size = 6, family = "Times New Roman") + ylim(-2.5, 4.2)


discovery_neg_cpca <- read.csv("processed_files/cpca_processed_neg.csv")

## Pregnenolone sulfate is soi 1881
discovery_neg_cpca$DM_INCIDENTE <- as.factor(discovery_neg_cpca$DM_INCIDENTE)
discovery_neg_cpca$Sex <- as.factor(discovery_neg_cpca$Sex)
levels(discovery_neg_cpca$Sex) <- c("Women", "Men")
discovery_neg_cpca$Sex <- as.factor(as.character(discovery_neg_cpca$Sex))

lm.model <- lm(SOI1881~Age+BMI+HT+HOMA_IR+Smoking+Family_history_DM+seguimiento+Glucose, 
               data = discovery_neg_cpca)
mm_merge <- merge(discovery_neg_cpca, data.frame(lm.model$residuals), by = 0)


p_discovery <- ggplot(mm_merge, aes(x = Sex, y = lm.model.residuals))+ 
  geom_point(size=2, alpha = 0.5, aes(color = DM_INCIDENTE),
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7))+
  geom_hpline(stat = "summary", fun = "mean", width = 0.25, aes(group = DM_INCIDENTE), 
              position = position_dodge(width = 0.7))+
  labs(color = "T2D event")+
  theme_bw() +
  scale_colour_manual(values = c("#E66100", "#1B9E77")) +
  theme(axis.text.y = element_text(family = "Times New Roman", size = 17),
        axis.text.x = element_text(family = "Times New Roman", size = 17),
        plot.title = element_text(family = "Times New Roman", size = 14, hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.title = element_text(family = "Times New Roman", size = 19), 
        legend.text = element_text(family = "Times New Roman", size = 18)) +
  stat_compare_means(aes(group = DM_INCIDENTE), label = "p.format", method = "t.test", size = 6, family = "Times New Roman", label.y = 2.5) +
  stat_compare_means(method = "t.test", label.y = 3, label.x = 1.35, size = 6, family = "Times New Roman") + ylim(-4.1, 3.2)



sex_pregS_plot <- ggarrange(p_discovery, p_validation, ncol = 2, nrow = 1, 
                            labels = c("A", "B"), common.legend = TRUE, legend = "bottom")

ggsave(sex_pregS_plot, filename = "figures/Figure3.jpg", 
       width = 12, height = 5.5, dpi = 900)


