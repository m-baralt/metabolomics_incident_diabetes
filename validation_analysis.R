library(ggfortify)
library(magrittr)
library(kableExtra)
library(ggplot2)
library(patchwork)
library(pROC)
library(caret)

path <- "/home/maria/Incident_type2_diabetes/"
setwd(path)
source("Functions_validation.R")

#####################
## Data processing ##
#####################

batch1_data <- read_prepare_data(path = "~/FIS2018_subestudi1/Data_072024/DM1_2024_export.rds")
batch2_data <- read_prepare_data(path = "~/FIS2018_subestudi1/Data_072024/DM2_2024_export.rds")
batch3_data <- read_prepare_data(path = "~/FIS2018_subestudi1/Data_072024/DM3_2024_export.rds")
batch4_data <- read_prepare_data(path = "~/FIS2018_subestudi1/Data_072024/DM4_2024_export.rds")
kableExtra::kable(batch1_data$compounds_info)
kableExtra::kable(batch1_data$metabs_names, col.names = "Features names")

## Batch 1
object_filter_impute_b1 <- filter_impute(metabs = batch1_data$metabs_names, 
                                         df = batch1_data$Long_data, 
                                         var.batch = "subBatch", 
                                         n.batch = length(unique(batch1_data$Long_data$subBatch)),
                                         perc.missings = 0.2)

kableExtra::kable(object_filter_impute_b1$NAs_per_batch)

pca_b1 <- pca_plot(df = object_filter_impute_b1$log_transf, 
                   metabs = object_filter_impute_b1$metabs.used, 
                   var.colour = "Class") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pca_b1 <- pca_b1 + theme_bw() + labs(title = "Batch 1") +
  theme(legend.position = "bottom", legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
        axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), legend.spacing.x = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
  scale_color_brewer(palette = "Set1") 

## Batch 2
object_filter_impute_b2 <- filter_impute(metabs = batch2_data$metabs_names, 
                                         df = batch2_data$Long_data, 
                                         var.batch = "subBatch", 
                                         n.batch = length(unique(batch2_data$Long_data$subBatch)),
                                         perc.missings = 0.2)

kableExtra::kable(object_filter_impute_b2$NAs_per_batch)

pca_b2 <- pca_plot(df = object_filter_impute_b2$log_transf, 
                   metabs = object_filter_impute_b2$metabs.used, 
                   var.colour = "Class") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pca_b2 <- pca_b2 + theme_bw() + labs(title = "Batch 2") +
  theme(legend.position = "bottom", legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
        axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), legend.spacing.x = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
  scale_color_brewer(palette = "Set1") 

## Batch 3
object_filter_impute_b3 <- filter_impute(metabs = batch3_data$metabs_names, 
                                         df = batch3_data$Long_data, 
                                         var.batch = "subBatch", 
                                         n.batch = length(unique(batch3_data$Long_data$subBatch)),
                                         perc.missings = 0.2)

kableExtra::kable(object_filter_impute_b3$NAs_per_batch)

pca_b3 <- pca_plot(df = object_filter_impute_b3$log_transf, 
                   metabs = object_filter_impute_b3$metabs.used, 
                   var.colour = "Class") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pca_b3 <- pca_b3 + theme_bw() + labs(title = "Batch 3") +
  theme(legend.position = "bottom", legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
        axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), legend.spacing.x = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
  scale_color_brewer(palette = "Set1") 

## Batch 4
object_filter_impute_b4 <- filter_impute(metabs = batch4_data$metabs_names, 
                                         df = batch4_data$Long_data, 
                                         var.batch = "subBatch", 
                                         n.batch = length(unique(batch4_data$Long_data$subBatch)),
                                         perc.missings = 0.2)

kableExtra::kable(object_filter_impute_b4$NAs_per_batch)

pca_b4 <- pca_plot(df = object_filter_impute_b4$log_transf, 
                   metabs = object_filter_impute_b4$metabs.used, 
                   var.colour = "Class") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pca_b4 <- pca_b4 + theme_bw() + labs(title = "Batch 4") +
  theme(legend.position = "bottom", legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
        axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), legend.spacing.x = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
  scale_color_brewer(palette = "Set1") 

combined_pcas <- ggpubr::ggarrange(pca_b1, pca_b2, pca_b3, pca_b4, ncol = 2, 
                                   nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(plot=combined_pcas, filename = "figures/PCA_batches.jpeg", 
       width = 10, height = 8, dpi = 500)

df1 <- as.data.frame(t(object_filter_impute_b1$filtered_imputed_data[,batch1_data$metabs_names]))
df2 <- as.data.frame(t(object_filter_impute_b2$filtered_imputed_data[,batch2_data$metabs_names]))
df3 <- as.data.frame(t(object_filter_impute_b3$filtered_imputed_data[,batch3_data$metabs_names]))
df4 <- as.data.frame(t(object_filter_impute_b4$filtered_imputed_data[,batch4_data$metabs_names]))

list.dfs <- list(data.frame(metabolites = rownames(df1), df1),
                 data.frame(metabolites = rownames(df2), df2),
                 data.frame(metabolites = rownames(df3), df3),
                 data.frame(metabolites = rownames(df4), df4))

list.injection <- list(batch1_data$injectionOrder, 
                       batch2_data$injectionOrder, 
                       batch3_data$injectionOrder,
                       batch4_data$injectionOrder)

all_batches <- joint_data(list.dfs = list.dfs, list.injection = list.injection)

distr_plots <- data_distribution(df = all_batches, 
                                 metabs = batch1_data$metabs_names)

ggpubr::ggarrange(plotlist = distr_plots[1:3], nrow = 3)
ggpubr::ggarrange(plotlist = distr_plots[4:6], nrow = 3)

metabolites <- batch1_data$metabs_names
all_batches[,metabolites] <- log(all_batches[,metabolites])

pca_all_before <- pca_plot(df = all_batches, 
                           metabs = metabolites, 
                           var.colour = "Batch") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pca_all_before <- pca_all_before + theme_bw() + labs(title = "Before ComBat correction") +
  theme(legend.position = "bottom", legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
        axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), legend.spacing.x = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
  scale_color_brewer(palette = "Set1") 

corrected.data <- combat_correction(df = all_batches, 
                                    metabs = batch1_data$metabs_names)

pca_all_after <- pca_plot(df = corrected.data, 
                          metabs = metabolites, 
                          var.colour = "Batch") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pca_all_after <- pca_all_after + theme_bw() + labs(title = "After ComBat correction") +
  theme(legend.position = "bottom", legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
        axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
        axis.text.y = element_text(size = 10), legend.spacing.x = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
  scale_color_brewer(palette = "Set1") 

combined_pcas_combat <- ggpubr::ggarrange(pca_all_before, pca_all_after, ncol = 2, 
                                          nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave(plot=combined_pcas_combat, filename = "figures/PCA_combat.jpeg", 
       width = 10, height = 5, dpi = 500)

# The distribution of each log-transformed metabolite is visualized once the data has been corrected. 
# The metabolite MZ_134.0472_MZ_135.1208 does not contain a mixture distribution (two peaks) anymore.
distr_plots <- data_distribution_corr(df = corrected.data, 
                                      metabs = batch1_data$metabs_names)
plots_metabs_hist <- ggpubr::ggarrange(plotlist = distr_plots, nrow = 3, ncol = 2)

ggsave(plot=plots_metabs_hist, filename = "figures/distrib_metabs.jpeg", 
       width = 10, height = 12, dpi = 500)

object_data <- list(AllBatches = all_batches, 
                    CorrectedData = corrected.data,
                    metabs.used = object_filter_impute_b1$metabs.used)
saveRDS(object_data, "processed_files/all_batches_alldata.RDS")

##########################
## Statistical Analysis ##
##########################

db <- read.csv("~/FIS2018_subestudi1/Updated_db.csv", row.names = 1)
matching_ids <- readxl::read_excel("~/FIS2018_subestudi1/match_ids.xlsx")
colnames(matching_ids)[1] <- "subName"
db <- merge(matching_ids[,-3], db, by = "CIP")

object_data <- readRDS("~/FIS2018_subestudi1/Data_072024/all_batches_alldata.RDS")
corrected.data <- object_data$CorrectedData
merge.df <- merge(db, corrected.data, by = "subName")
merge.df <- merge.df[!(merge.df$FU_time<=0),]
kable(table(merge.df$event.DM2), col.names = c("Event DM2", "Freq"))
kable(table(merge.df$event.prediab.dm2), col.names = c("Event prediabetes and DM2", "Freq"))

merge.df$transition <- paste(merge.df$Prediabetes, merge.df$event.prediab.dm2, sep = "->")
merge.df$transition <- as.factor(merge.df$transition)
kable(table(merge.df$transition), col.names = c("Glucemia change", "Freq"))

merge.df$Sex <- as.factor(merge.df$Sex)
merge.df$HT <- as.factor(merge.df$HT)
merge.df$Smoking <- as.factor(merge.df$Smoking)
merge.df$Prediabetes <- as.factor(merge.df$Prediabetes)

res1 <- compareGroups::compareGroups(event.DM2 ~ Sex+Age+BMI+HT+Smoking+
                                       FU_time+Prediabetes, 
                                     data = merge.df)
compTab <- compareGroups::createTable(res1, show.ratio = FALSE)

compareGroups::export2csv(compTab, "figures/comparegroups_validation.csv")

covars <- c("Sex", "Age", "BMI", "HT", "Smoking", "FU_time", "Prediabetes")
merge.df <- merge.df[,c("event.DM2", "event.prediab.dm2", "transition", "Sex", "Age", 
                        "BMI", "HT", "Smoking", "FU_time", "Prediabetes", object_data$metabs.used)]

merge.df[,c("Age", "BMI","FU_time", object_data$metabs.used)] <- 
  scale(merge.df[,c("Age", "BMI","FU_time",object_data$metabs.used)], center = TRUE, scale = TRUE)


# R1
stats_results_r1 <- statistical_analysis(df = merge.df,
                                         metabolites = object_data$metabs.used,
                                         covars = covars, outcome = "event.DM2",
                                         do.Par = TRUE, nCore = 15)

# R2 - only prediabetics
covars <- c("Sex", "Age", "BMI", "HT", "Smoking", "FU_time")
stats_results_r2 <- statistical_analysis(df = merge.df[merge.df$Prediabetes==1,],
                                         metabolites = object_data$metabs.used,
                                         covars = covars, outcome = "event.DM2",
                                         do.Par = TRUE, nCore = 15)

# R3
covars <- c("Sex", "Age", "BMI", "HT", "Smoking", "FU_time", "Prediabetes")
stats_results_r3 <- statistical_analysis_lm(df = merge.df, 
                                            metabolites = object_data$metabs.used,
                                            covars = covars, outcome = "event.prediab.dm2",
                                            do.Par = TRUE, nCore = 15)

# R4
stats_results_r4 <- statistical_analysis(df = merge.df[merge.df$event.prediab.dm2!=1,],
                                         metabolites = object_data$metabs.used,
                                         covars = covars, outcome = "event.DM2",
                                         do.Par = TRUE, nCore = 15)

# R5 glycemia transition with normoglycemia at t0
df_0 <- merge.df[merge.df$transition %in% c("0->0", "0->1","0->2"),]
df_0$transition <- as.factor(as.character(df_0$transition))
levels(df_0$transition) <- c("0", "1", "2")
df_0$transition <- as.numeric(as.character(df_0$transition))

covars <- c("Sex", "Age", "BMI", "HT", "Smoking", "FU_time")
stats_results_r5 <- statistical_analysis_lm(df = df_0, 
                                            metabolites = object_data$metabs.used,
                                            covars = covars, outcome = "transition",
                                            do.Par = TRUE, nCore = 15)

# R6 glycemia transition with prediabetes at t0
df_1 <- merge.df[merge.df$transition %in% c("1->0", "1->1","1->2"),]
df_1$transition <- as.factor(as.character(df_1$transition))
levels(df_1$transition) <- c("0", "1", "2")
df_1$transition <- as.numeric(as.character(df_1$transition))

covars <- c("Sex", "Age", "BMI", "HT", "Smoking", "FU_time")
stats_results_r6 <- statistical_analysis_lm(df = df_1, 
                                            metabolites = object_data$metabs.used,
                                            covars = covars, outcome = "transition",
                                            do.Par = TRUE, nCore = 15)

stats_results_r1 <- data.frame(analysis = "R1", stats_results_r1)
stats_results_r2 <- data.frame(analysis = "R2", stats_results_r2)
stats_results_r3 <- data.frame(analysis = "R3", stats_results_r3)
stats_results_r4 <- data.frame(analysis = "R4", stats_results_r4)
stats_results_r5 <- data.frame(analysis = "R5", stats_results_r5)
stats_results_r6 <- data.frame(analysis = "R6", stats_results_r6)

stats_results_glm <- rbind(stats_results_r1, stats_results_r2, stats_results_r4)
stats_results_lm <- rbind(stats_results_r3, stats_results_r5, stats_results_r6)

allinfo <- readRDS("~/FIS2018_subestudi1/Data_072024/DM1_2024_export.rds")
compounds_info <- allinfo[[2]] %>% as.data.frame()
compounds_info$Metabolite <- paste("MZ",compounds_info$`m/z`,"MZ_135.1208", sep = "_")
matching <- compounds_info[,c("Metabolite", "Compound")]
stats_results_glm <- merge(matching, stats_results_glm, by = "Metabolite")
colnames(matching)[1] <- "Feature"
stats_results_lm <- merge(matching, stats_results_lm, by = "Feature")

write.csv(stats_results_glm, "processed_files/validation_results_glm.csv")
write.csv(stats_results_lm, "processed_files/validation_results_lm.csv")

#######################
## predictive models ##
#######################

set.seed(2409)
dm2_subjects <- merge.df[merge.df$event.DM2 == 1,]
ngroups <- 5
idx_folds <- sample(cut(seq(1,nrow(dm2_subjects)), breaks=ngroups, labels=FALSE))
dm2_folds <- split(dm2_subjects, idx_folds)

ct_subjects <- merge.df[merge.df$event.DM2 == 0,]
ngroups <- 5
idx_folds <- sample(cut(seq(1,nrow(ct_subjects)), breaks=ngroups, labels=FALSE))
ct_folds <- split(ct_subjects, idx_folds)

complete_folds <- plyr::llply(1:5, function(i){
  fold_set <- rbind(dm2_folds[[i]], ct_folds[[i]])
  return(fold_set)
})

metrics_CV_alldata <- plyr::ldply(1:5, function(i){
  test_set <- complete_folds[[i]]
  train_set <- ldply(seq(1,5)[!(seq(1,5) %in% i)], function(j) complete_folds[[j]])
  
  class_weights <- ifelse(train_set$event.DM2 == 1, 
                          sum(train_set$event.DM2==0)/nrow(train_set), 
                          sum(train_set$event.DM2==1)/nrow(train_set))
  pred_model_basic <- glm(event.DM2~Sex+Age+BMI+HT+Smoking+Prediabetes, 
                          data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_basic, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  f <- as.formula(paste0("event.DM2~Sex+Age+BMI+HT+Smoking+Prediabetes+", 
                         paste(object_data$metabs.used, collapse = "+")))
  
  pred_model_metabs <- glm(f, data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_metabs, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic_metabs <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  f <- as.formula(paste0("event.DM2~Sex+Age+BMI+HT+Smoking+Prediabetes+MZ_152.0567_MZ_135.1208+MZ_395.1898_MZ_135.1208"))
  
  pred_model_metabs_sign <- glm(f, data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_metabs_sign, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic_metabs_sign <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  rbind(data.frame(Model = "Basic", t(cmat_basic$byClass)),
        data.frame(Model = "All_metabs", t(cmat_basic_metabs$byClass)),
        data.frame(Model = "Sign_metabs", t(cmat_basic_metabs_sign$byClass)))
  
})

summary_metrics_all <- plyr::ddply(metrics_CV_alldata, "Model", function(dx){
  dx_set <- dx[,c("Sensitivity","Specificity","F1","Balanced.Accuracy")] 
  allMeans <- colMeans(dx_set) %>% t() %>% as.data.frame()
  allSDs <- apply(dx_set, 2, function(j) sd(j)) %>% t() %>% as.data.frame()
  colnames(allSDs) <- paste0("sd_", colnames(allSDs))
  cbind(allMeans, allSDs)
})

set.seed(2409)
preds_data <- merge.df[merge.df$Prediabetes==1,]
dm2_subjects <- preds_data[preds_data$event.DM2 == 1,]
ngroups <- 5
idx_folds <- sample(cut(seq(1,nrow(dm2_subjects)), breaks=ngroups, labels=FALSE))
dm2_folds <- split(dm2_subjects, idx_folds)

ct_subjects <- preds_data[preds_data$event.DM2 == 0,]
ngroups <- 5
idx_folds <- sample(cut(seq(1,nrow(ct_subjects)), breaks=ngroups, labels=FALSE))
ct_folds <- split(ct_subjects, idx_folds)

complete_folds <- plyr::llply(1:5, function(i){
  fold_set <- rbind(dm2_folds[[i]], ct_folds[[i]])
  return(fold_set)
})

metrics_CV_only_pred <- plyr::ldply(1:5, function(i){
  test_set <- complete_folds[[i]]
  train_set <- ldply(seq(1,5)[!(seq(1,5) %in% i)], function(j) complete_folds[[j]])
  
  class_weights <- ifelse(train_set$event.DM2 == 1, 
                          sum(train_set$event.DM2==0)/nrow(train_set), 
                          sum(train_set$event.DM2==1)/nrow(train_set))
  pred_model_basic <- glm(event.DM2~Sex+Age+BMI+HT+Smoking, 
                          data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_basic, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  f <- as.formula(paste0("event.DM2~Sex+Age+BMI+HT+Smoking+", 
                         paste(object_data$metabs.used, collapse = "+")))
  
  pred_model_metabs <- glm(f, data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_metabs, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic_metabs <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  f <- as.formula(paste0("event.DM2~Sex+Age+BMI+HT+Smoking+MZ_152.0567_MZ_135.1208+MZ_395.1898_MZ_135.1208"))
  
  pred_model_metabs_sign <- glm(f, data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_metabs_sign, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic_metabs_sign <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  rbind(data.frame(Model = "Basic", t(cmat_basic$byClass)),
        data.frame(Model = "All_metabs", t(cmat_basic_metabs$byClass)),
        data.frame(Model = "Sign_metabs", t(cmat_basic_metabs_sign$byClass)))
  
})

summary_metrics_only_pred <- plyr::ddply(metrics_CV_only_pred, "Model", function(dx){
  dx_set <- dx[,c("Sensitivity","Specificity","F1","Balanced.Accuracy")] 
  allMeans <- colMeans(dx_set) %>% t() %>% as.data.frame()
  allSDs <- apply(dx_set, 2, function(j) sd(j)) %>% t() %>% as.data.frame()
  colnames(allSDs) <- paste0("sd_", colnames(allSDs))
  cbind(allMeans, allSDs)
})


set.seed(2409)
preds_data <- merge.df[!(merge.df$event.prediab.dm2==1),]
dm2_subjects <- preds_data[preds_data$event.DM2 == 1,]
ngroups <- 5
idx_folds <- sample(cut(seq(1,nrow(dm2_subjects)), breaks=ngroups, labels=FALSE))
dm2_folds <- split(dm2_subjects, idx_folds)

ct_subjects <- preds_data[preds_data$event.DM2 == 0,]
ngroups <- 5
idx_folds <- sample(cut(seq(1,nrow(ct_subjects)), breaks=ngroups, labels=FALSE))
ct_folds <- split(ct_subjects, idx_folds)

complete_folds <- plyr::llply(1:5, function(i){
  fold_set <- rbind(dm2_folds[[i]], ct_folds[[i]])
  return(fold_set)
})

metrics_CV_nopreds <- plyr::ldply(1:5, function(i){
  test_set <- complete_folds[[i]]
  train_set <- ldply(seq(1,5)[!(seq(1,5) %in% i)], function(j) complete_folds[[j]])
  class_weights <- ifelse(train_set$event.DM2 == 1, 
                          sum(train_set$event.DM2==0)/nrow(train_set), 
                          sum(train_set$event.DM2==1)/nrow(train_set))
  pred_model_basic <- glm(event.DM2~Sex+Age+BMI+HT+Smoking+Prediabetes, 
                          data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_basic, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  f <- as.formula(paste0("event.DM2~Sex+Age+BMI+HT+Smoking+Prediabetes+", 
                         paste(object_data$metabs.used, collapse = "+")))
  
  pred_model_metabs <- glm(f, data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_metabs, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic_metabs <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  f <- as.formula(paste0("event.DM2~Sex+Age+BMI+HT+Smoking+Prediabetes+MZ_152.0567_MZ_135.1208+MZ_395.1898_MZ_135.1208"))
  
  pred_model_metabs_sign <- glm(f, data = train_set, family = "binomial", weights = class_weights)
  
  probabilities <- predict(pred_model_metabs_sign, newdata = test_set, type = "response")
  
  roc_curve <- roc(test_set$event.DM2, probabilities)
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
  optimal_threshold <- optimal_threshold$threshold
  
  predictions <- as.factor(as.character(ifelse(probabilities > optimal_threshold, 1, 0)))
  
  cmat_basic_metabs_sign <- confusionMatrix(predictions, as.factor(test_set$event.DM2), positive = "1")
  
  rbind(data.frame(Model = "Basic", t(cmat_basic$byClass)),
        data.frame(Model = "All_metabs", t(cmat_basic_metabs$byClass)),
        data.frame(Model = "Sign_metabs", t(cmat_basic_metabs_sign$byClass)))
  
})

summary_metrics_nopreds <- plyr::ddply(metrics_CV_nopreds, "Model", function(dx){
  dx_set <- dx[,c("Sensitivity","Specificity","F1","Balanced.Accuracy")] 
  allMeans <- colMeans(dx_set) %>% t() %>% as.data.frame()
  allSDs <- apply(dx_set, 2, function(j) sd(j)) %>% t() %>% as.data.frame()
  colnames(allSDs) <- paste0("sd_", colnames(allSDs))
  cbind(allMeans, allSDs)
})

summary_all_models <- rbind(data.frame(Data = "All_data", summary_metrics_all), 
                            data.frame(Data = "Only prediabetics", summary_metrics_only_pred), 
                            data.frame(Data = "No prediabetics", summary_metrics_nopreds))
write.csv(summary_all_models, "processed_files/metrics_prediction.csv")

