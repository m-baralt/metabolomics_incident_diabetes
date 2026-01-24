library(mWISE)
library(ggplot2)

processed_data <- read.csv("processed_files/cpca_processed_pos.csv")
rownames(processed_data) <- processed_data$X
processed_features <- processed_data[,grepl(pattern = "SOI", x = colnames(processed_data))]
processed_features <-  as.data.frame(t(processed_features))

MS_pos <- read.csv("data/DM_FIS2018_Hilic_pos_results2023_filled.csv")
rownames(MS_pos) <- paste0("SOI", MS_pos$Qidx)

merged_dat <- merge(MS_pos[,c("mass", "rtmed")], processed_features, by = 0)
rownames(merged_dat) <- merged_dat$Row.names

Peak_list <- data.frame(Peak.Id = rownames(merged_dat), mz = merged_dat$mass, 
                        rt = merged_dat$rtmed, merged_dat[,4:ncol(merged_dat)])
Peak_list <- Peak_list[,!grepl(pattern = "QC", x = colnames(Peak_list))]

data("Info.Add")
positive_add <- Info.Add[Info.Add$polarity %in% "positive",]
selectedAdds <- positive_add$name[positive_add$Freq>0.001]

Cpdadd <- CpdaddPreparation(do.Par = TRUE, nClust = 20)

Annotated_List <- matchingStage(Peak.List = Peak_list, 
                                Cpd.Add = Cpdadd,
                                polarity = "positive", 
                                Add.List = selectedAdds, 
                                nClust = 20,
                                do.Par = TRUE)

Intensity.idx <- seq(4,ncol(Peak_list))

clustered <- featuresClustering(Peak.List = Peak_list, 
                                Intensity.idx = Intensity.idx, 
                                nClust = 20,
                                do.Par = TRUE)

Annotated.Tab <- merge(Annotated_List$Peak.Cpd,
                       clustered$Peak.List[,c("Peak.Id", "pcgroup")],
                       by = "Peak.Id")

MH.Tab <- clusterBased.filter(df = Annotated.Tab, 
                              polarity = "positive")

Input.diffusion <- diffusion.input(df = MH.Tab,
                                   input.type = "binary",
                                   Unique.Annotation = TRUE,
                                   do.Par = FALSE)

diff.Cpd <- set.diffusion(df = Input.diffusion,
                          scores = "raw",
                          graph.name = "fella",
                          do.Par = FALSE)

Diffusion.Results <- diff.Cpd$Diffusion.Results

MH.Tab <- recoveringPeaks(Annotated.Tab = Annotated.Tab, 
                          MH.Tab = MH.Tab)

Diff.Tab <- merge(x = MH.Tab, 
                  y = Diffusion.Results,
                  by = "Compound", 
                  all.x = TRUE)

Ranked.Tab <- finalResults(Diff.Tab = Diff.Tab, 
                           score = "raw", do.Par = FALSE)

Annotated.Tab2 <- modifiedTabs(df = Annotated.Tab, 
                               do.Par = FALSE)

MH.Tab2 <- modifiedTabs(df = MH.Tab, 
                        do.Par = FALSE)

Annotated.dataset <- list(Annotated.Tab = Annotated.Tab2,
                          Clustered.Tab = clustered,
                          MH.Tab = MH.Tab2,
                          Diff.Tab = Diff.Tab,
                          Ranked.Tab = Ranked.Tab)

write.csv(Annotated.dataset$Ranked.Tab, "processed_files/annotated_mwise_pos.csv")

### NEGATIVE ###
processed_data <- read.csv("processed_files/cpca_processed_neg.csv")
rownames(processed_data) <- processed_data$X
processed_features <- processed_data[,grepl(pattern = "SOI", x = colnames(processed_data))]
processed_features <-  as.data.frame(t(processed_features))

MS_neg <- read.csv("data/DM_FIS2018_Hilic_neg_results2023_filled.csv")
rownames(MS_neg) <- paste0("SOI", MS_neg$Qidx)

merged_dat <- merge(MS_neg[,c("mass", "rtmed")], processed_features, by = 0)
rownames(merged_dat) <- merged_dat$Row.names

Peak_list <- data.frame(Peak.Id = rownames(merged_dat), mz = merged_dat$mass, rt = merged_dat$rtmed, merged_dat[,4:ncol(merged_dat)])
Peak_list <- Peak_list[,!grepl(pattern = "QC", x = colnames(Peak_list))]

data("Info.Add")
negative_add <- Info.Add[Info.Add$polarity %in% "negative",]
selectedAdds <- negative_add$name[(negative_add$Freq>0.001) | (negative_add$name %in% "M+Cl")]

Annotated_List <- matchingStage(Peak.List = Peak_list, 
                                Cpd.Add = Cpdadd,
                                polarity = "negative", 
                                Add.List = selectedAdds, 
                                nClust = 20,
                                do.Par = TRUE)

Intensity.idx <- seq(4,ncol(Peak_list))

clustered <- featuresClustering(Peak.List = Peak_list, 
                                Intensity.idx = Intensity.idx, 
                                nClust = 20,
                                do.Par = TRUE)

Annotated.Tab <- merge(Annotated_List$Peak.Cpd,
                       clustered$Peak.List[,c("Peak.Id", "pcgroup")],
                       by = "Peak.Id")

MH.Tab <- clusterBased.filter(df = Annotated.Tab, 
                              polarity = "negative")

Input.diffusion <- diffusion.input(df = MH.Tab,
                                   input.type = "binary",
                                   Unique.Annotation = TRUE,
                                   do.Par = FALSE)

diff.Cpd <- set.diffusion(df = Input.diffusion,
                          scores = "raw",
                          graph.name = "fella",
                          do.Par = FALSE)

Diffusion.Results <- diff.Cpd$Diffusion.Results

MH.Tab <- recoveringPeaks(Annotated.Tab = Annotated.Tab, 
                          MH.Tab = MH.Tab)

Diff.Tab <- merge(x = MH.Tab, 
                  y = Diffusion.Results,
                  by = "Compound", 
                  all.x = TRUE)

Ranked.Tab <- finalResults(Diff.Tab = Diff.Tab, 
                           score = "raw", do.Par = FALSE)

Annotated.Tab2 <- modifiedTabs(df = Annotated.Tab, 
                               do.Par = FALSE)

MH.Tab2 <- modifiedTabs(df = MH.Tab, 
                        do.Par = FALSE)

Annotated.dataset <- list(Annotated.Tab = Annotated.Tab2,
                          Clustered.Tab = clustered,
                          MH.Tab = MH.Tab2,
                          Diff.Tab = Diff.Tab,
                          Ranked.Tab = Ranked.Tab)

write.csv(Annotated.dataset$Ranked.Tab, "processed_files/annotated_mwise_neg.csv")

## Results for enrichment ##

files_results <- list.files("results/discovery/", full.names = TRUE)

files.pos <- files_results[grep("pos", files_results)]
files.neg <- files_results[grep("neg", files_results)]

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

mwise_pos <- read.csv("processed_files/annotated_mwise_pos.csv", row.names = 1)
colnames(mwise_pos)[2] <- "Feature"
stats_cmps_pos <- merge(cpca_pos, mwise_pos, by = "Feature")

mwise_neg <- read.csv("processed_files/annotated_mwise_neg.csv", row.names = 1)
colnames(mwise_neg)[2] <- "Feature"
stats_cmps_neg <- merge(cpca_neg, mwise_neg, by = "Feature")

## R1 ##
stats_cmps_pos_r1 <- stats_cmps_pos[stats_cmps_pos$Analysis %in% "R1",]
stats_cmps_neg_r1 <- stats_cmps_neg[stats_cmps_neg$Analysis %in% "R1",]

stats_cmps_pos_r1 <- stats_cmps_pos_r1[stats_cmps_pos_r1$Ranking<=3,]
stats_cmps_neg_r1 <- stats_cmps_neg_r1[stats_cmps_neg_r1$Ranking<=3,]



r1_cmp <- c(stats_cmps_pos_r1$OR[stats_cmps_pos_r1$qvalue<0.05],
            stats_cmps_neg_r1$OR[stats_cmps_neg_r1$qvalue<0.05])

r1_cmp <- log(r1_cmp)

names(r1_cmp) <- c(stats_cmps_pos_r1$Compound[stats_cmps_pos_r1$qvalue<0.05],
                   stats_cmps_neg_r1$Compound[stats_cmps_neg_r1$qvalue<0.05])

## R2 ##
stats_cmps_pos_r2 <- stats_cmps_pos[stats_cmps_pos$Analysis %in% "R2",]
stats_cmps_neg_r2 <- stats_cmps_neg[stats_cmps_neg$Analysis %in% "R2",]

stats_cmps_pos_r2 <- stats_cmps_pos_r2[stats_cmps_pos_r2$Ranking<=3,]
stats_cmps_neg_r2 <- stats_cmps_neg_r2[stats_cmps_neg_r2$Ranking<=3,]

r2_cmp <- c(stats_cmps_pos_r2$OR[stats_cmps_pos_r2$qvalue<0.05],
            stats_cmps_neg_r2$OR[stats_cmps_neg_r2$qvalue<0.05])

r2_cmp <- log(r2_cmp)

names(r2_cmp) <- c(stats_cmps_pos_r2$Compound[stats_cmps_pos_r2$qvalue<0.05],
                   stats_cmps_neg_r2$Compound[stats_cmps_neg_r2$qvalue<0.05])

## R3 ##
stats_cmps_pos_r3 <- stats_cmps_pos[stats_cmps_pos$Analysis %in% "R3",]
stats_cmps_neg_r3 <- stats_cmps_neg[stats_cmps_neg$Analysis %in% "R3",]

stats_cmps_pos_r3 <- stats_cmps_pos_r3[stats_cmps_pos_r3$Ranking<=3,]
stats_cmps_neg_r3 <- stats_cmps_neg_r3[stats_cmps_neg_r3$Ranking<=3,]

r3_cmp <- c(stats_cmps_pos_r3$Beta[stats_cmps_pos_r3$qvalue<0.05],
            stats_cmps_neg_r3$Beta[stats_cmps_neg_r3$qvalue<0.05])

names(r3_cmp) <- c(stats_cmps_pos_r3$Compound[stats_cmps_pos_r3$qvalue<0.05],
                   stats_cmps_neg_r3$Compound[stats_cmps_neg_r3$qvalue<0.05])

## R4 ##
stats_cmps_pos_r4 <- stats_cmps_pos[stats_cmps_pos$Analysis %in% "R4",]
stats_cmps_neg_r4 <- stats_cmps_neg[stats_cmps_neg$Analysis %in% "R4",]

stats_cmps_pos_r4 <- stats_cmps_pos_r4[stats_cmps_pos_r4$Ranking<=3,]
stats_cmps_neg_r4 <- stats_cmps_neg_r4[stats_cmps_neg_r4$Ranking<=3,]

r4_cmp <- c(stats_cmps_pos_r4$OR[stats_cmps_pos_r4$qvalue<0.05],
            stats_cmps_neg_r4$OR[stats_cmps_neg_r4$qvalue<0.05])

r4_cmp <- log(r4_cmp)

names(r4_cmp) <- c(stats_cmps_pos_r4$Compound[stats_cmps_pos_r4$qvalue<0.05],
                   stats_cmps_neg_r4$Compound[stats_cmps_neg_r4$qvalue<0.05])

## R5 ##
stats_cmps_pos_r5 <- stats_cmps_pos[stats_cmps_pos$Analysis %in% "R5",]
stats_cmps_neg_r5 <- stats_cmps_neg[stats_cmps_neg$Analysis %in% "R5",]

stats_cmps_pos_r5 <- stats_cmps_pos_r5[stats_cmps_pos_r5$Ranking<=3,]
stats_cmps_neg_r5 <- stats_cmps_neg_r5[stats_cmps_neg_r5$Ranking<=3,]

r5_cmp <- c(stats_cmps_pos_r5$Beta[stats_cmps_pos_r5$qvalue<0.05],
            stats_cmps_neg_r5$Beta[stats_cmps_neg_r5$qvalue<0.05])

names(r5_cmp) <- c(stats_cmps_pos_r5$Compound[stats_cmps_pos_r5$qvalue<0.05],
                   stats_cmps_neg_r5$Compound[stats_cmps_neg_r5$qvalue<0.05])



### statistical analysis

# Function that applies a linear regression model for each LC-MS feature,
# corrects the p-values and outputs the resulting data frame and a summary table
statistical_analysis_lm <- function(df, feats.names, covars, outcome, 
                                    do.Par = TRUE, nCore = 10, psign = 0.05){
  doParallel::registerDoParallel(nCore)
  pvals_df <- plyr::ldply(feats.names, function(soi){
    f_regr <- as.formula(paste0(soi,"~",paste(covars, collapse = "+"),"+",outcome))
    lm.model <- lm(f_regr, data = df)
    s <- summary(lm.model)
    pval <- s$coefficients[outcome,4]
    data.frame(Feature = soi, pval = pval, 
               beta = s$coefficients[outcome,1],
               n = nobs(lm.model))
  },.parallel = do.Par)
  pvals_df$qvalue_trunc <- qvalue::qvalue_truncp(pvals_df$pval)$qvalues
  pvals_df$qvalue <- qvalue::qvalue(pvals_df$pval)$qvalues
  pvals_df$p_fdr <- p.adjust(pvals_df$pval, "fdr")
  summary.stats <- data.frame(Variable_interest = outcome,
                              Model_type = "linear regression",
                              sample_size = mean(pvals_df$n),
                              Number_features = length(feats.names),
                              Significant_features = sum(pvals_df$pval<psign),
                              Significant_features_fdr = sum(pvals_df$p_fdr<psign),
                              Significant_features_qvalue = sum(pvals_df$qvalue<psign),
                              Significant_features_qvalue_trunc = sum(pvals_df$qvalue_trunc<psign),
                              Min_qvalue = min(pvals_df$qvalue),
                              Sign_p = psign)
  p.hist <- ggplot2::ggplot(aes(x = pval), data = pvals_df) + geom_histogram()
  return(list(summary.stats = summary.stats,
              Pvals.df = pvals_df, p.hist = p.hist))
}


### Positive
merge.pos.cpca <- read.csv("processed_files/cpca_processed_pos.csv", row.names = 1)
covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.pos.cpca[,c(covars, grep(pattern = "SOI", x = colnames(merge.pos.cpca), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

### Negative
merge.neg.cpca <- read.csv("processed_files/cpca_processed_neg.csv", row.names = 1)

covars <- c("Sex","Age","BMI","HT","HOMA_IR","Smoking",
            "Family_history_DM","Glucose",
            "seguimiento")

x <- merge.neg.cpca[,c(covars, grep(pattern = "SOI", x = colnames(merge.neg.cpca), value = TRUE))]
x$Sex <- as.numeric(as.factor(x$Sex))
x$HT <- as.numeric(as.factor(x$HT))

### HOMA-IR
## positive
covars <- c("Sex","Age","BMI","HT","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.pos2 <- merge.pos.cpca[!is.na(merge.pos.cpca$HOMA_IR),]

stats_results <- statistical_analysis_lm(df = merge.pos2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.pos2), value = TRUE), 
                                         covars = covars, outcome = "HOMA_IR",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

pos_stats <- stats_results$Pvals.df

stats_cmps_pos <- merge(pos_stats, mwise_pos, by = "Feature")

covars <- c("Sex","Age","BMI","HT","Smoking",
            "Family_history_DM", "Glucose",
            "seguimiento")
merge.neg2 <- merge.neg.cpca[!is.na(merge.neg.cpca$HOMA_IR),]

stats_results <- statistical_analysis_lm(df = merge.neg2, 
                                         feats.names = grep(pattern = "SOI", x = colnames(merge.neg2), value = TRUE), 
                                         covars = covars, outcome = "HOMA_IR",
                                         do.Par = TRUE, nCore = 15)

t(stats_results$summary.stats)

neg_stats <- stats_results$Pvals.df

stats_cmps_neg <- merge(neg_stats, mwise_neg, by = "Feature")


stats_cmps_pos <- stats_cmps_pos[stats_cmps_pos$Ranking<=3,]
stats_cmps_neg <- stats_cmps_neg[stats_cmps_neg$Ranking<=3,]

cmp_ir <- c(stats_cmps_pos$beta[stats_cmps_pos$qvalue_trunc<0.05],
            stats_cmps_neg$beta[stats_cmps_neg$qvalue_trunc<0.05])

names(cmp_ir) <- c(stats_cmps_pos$Compound[stats_cmps_pos$qvalue_trunc<0.05],
                   stats_cmps_neg$Compound[stats_cmps_neg$qvalue_trunc<0.05])

compounds_analysis <- list(R1 = r1_cmp,
                           R3 = r3_cmp,
                           R4 = r4_cmp,
                           R5 = r5_cmp,
                           IR = cmp_ir)

saveRDS(compounds_analysis, "processed_files/list_sign_cmps.rds")

