# function to correct technical drift using common principal components analysis
cpca.driftRemoval <- function (df, nComps = 1) {
  if (!is.null(dim(df)) & !(ncol(as.matrix(df)) <= nComps)) {
    tmp <- df
    if (any(is.na(tmp))) 
      stop(paste("ERR: NA present in DataSet.", "Consider Imputation."))
    if (class(tmp) != "matrix") 
      tmp <- as.matrix(tmp)
    classes <- as.character(substr(rownames(tmp), 1, 2))
    classLabels <- unique(classes)
    tmp <- t(tmp)
    X <- t(tmp)[classes %in% classLabels, ]
    Y <- factor(classes[classes %in% classLabels])
    ngr <- nlevels(factor(Y))
    X <- scale(X, center = TRUE, scale = TRUE)
    rownames(X) <- NULL
    colnames(X) <- NULL
    C <- plyr::daply(data.frame(X, Y), "Y", function(x) cov(x[, 
                                                              -ncol(x)]))
    C <- aperm(C, c(2, 3, 1))
    cat(fill = TRUE)
    cat("Computing CPC model...")
    if (nComps <= 3) {
      cpcs <- cpca::cpc(C, k = 3)
    }
    else {
      cpcs <- cpca::cpc(C, k = nComps)
    }
    cat("Model done!", fill = TRUE)
    out <- cpcs$CPC
    cat(fill = TRUE)
    cat("Computing CPC explained variance...", fill = TRUE)
    rownames(out) <- colnames(X)
    colnames(out) <- paste("CPC", 1:ncol(out), sep = "")
    var.projected <- apply(out, 2, function(e) sum((X %*% 
                                                      e)^2))
    var.total <- sum(apply(X, 2, function(x) sum((x)^2)))
    var.cpc <- var.projected/var.total
    var.cpc <- round(var.cpc, 3)[c(1:nComps)]
    show(var.cpc)
    out <- cpcs$CPC[, 1:nComps]
    if (is.null(dim(out))) {
      out <- (as.matrix(out, ncol = 1))
    }
    aux <- list(out, var.cpc)
    names(aux) <- c("out", "variance")
    cat("Correcting drift...")
    dat <- scale(t(tmp), center = TRUE, scale = TRUE)
    diff.Mat <- ((dat %*% aux$out) %*% t(aux$out))
    Xc <- dat - ((dat %*% aux$out) %*% t(aux$out))
    XcRes <- Xc
    for (i in c(1:dim(dat)[2])) {
      XcRes[, i] <- (Xc[, i] * attributes(dat)$`scaled:scale`[i]) + 
        attributes(dat)$`scaled:center`[i]
    }
    Xc <- XcRes
    norm_tmp <- Xc
    cat("Correction done!", fill = TRUE)
    return(norm_df = norm_tmp)
  }
  else {
    return(norm_df = df)
  }
}

# function that applies 80% rule and imputes with 0s or minimum divided by 2
filtering_imputing <- function(df, features.names, imp.method = c("0s", "min/2")){
  nas.feats <- apply(df[,features.names], 2, function(c) sum(is.na(c)))
  feats.to.rm <- names(nas.feats[nas.feats > 0.2*nrow(df)])
  df_filt <- df[,!(colnames(df) %in% feats.to.rm)]
  new.names <- grep(pattern = "SOI", x = colnames(df_filt), value = TRUE)
  df_imp <- df_filt
  if (imp.method == "0s"){
    for (j in new.names){
      df_imp[is.na(df_imp[,j]),j] <- 0
    }
  } else {
    for (j in new.names){
      df_imp[is.na(df_imp[,j]),j] <- min(df_imp[,j], na.rm = TRUE)/2
    }
  }
  return(list(Original = df,
              Filtered = df_filt,
              Imputed = df_imp))
}

# Function that builds PCAs plots of the imputed data
data_prep_vis <- function(proc_object, inj_order){
  imp_data <- proc_object$log_transf
  imp_data <- data.frame(sampleID = rownames(imp_data), class = "sample", imp_data)
  imp_data$class[grepl(pattern = "QC", x = imp_data$sampleID)] <- "QC"
  imp_data$class[grepl(pattern = "Bf", x = imp_data$sampleID)] <- "Blank"
  imp_data$class[grepl(pattern = "Bi", x = imp_data$sampleID)] <- "Blank"
  p.pre <- autoplot(prcomp(imp_data[!(imp_data$class %in% "Blank"),-c(1,2)], center = TRUE, 
                           scale. = TRUE), 
                    data = imp_data[!(imp_data$class %in% "Blank"),], scale = FALSE,
                    colour = "class", size = 2)+
    theme(legend.position = "bottom", legend.title = element_text(size = 11), 
          legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
          axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
          axis.text.y = element_text(size = 1), legend.spacing.x = unit(0.4, 'cm')) + 
    theme_bw() + scale_color_brewer(palette = "Set1")
  proc_object$log_transf <- imp_data
  merged_inj <- merge(inj_order, imp_data, by = 0)
  rownames(merged_inj) <- merged_inj$Row.names
  merged_inj <- merged_inj[,-c(1,2)]
  p.inj <- autoplot(prcomp(merged_inj[!(merged_inj$class %in% "Blank"),-c(1:3)], center = TRUE, scale. = TRUE), 
                    data = merged_inj[!(merged_inj$class %in% "Blank"),], colour = "injectionOrder") +
    theme(legend.position = "bottom", legend.title = element_text(size = 11, vjust = 0.8), 
          legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
          axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
          axis.text.y = element_text(size = 1), legend.spacing.x = unit(0.6, 'cm')) + theme_bw()
  merged_inj$injectionOrder[merged_inj$class %in% "sample"] <- 600
  p.inj.qc <- autoplot(prcomp(merged_inj[!(merged_inj$class %in% "Blank"),-c(1:3)], center = TRUE, scale. = TRUE), 
                       data = merged_inj[!(merged_inj$class %in% "Blank"),], colour = "injectionOrder") +
    theme(legend.position = "bottom", legend.title = element_text(size = 11), 
          legend.text = element_text(size = 9), axis.title.x = element_text(size = 11), 
          axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 11), 
          axis.text.y = element_text(size = 1), legend.spacing.x = unit(0.4, 'cm')) + theme_bw()
  return(list(proc_object = proc_object, pca_data = p.pre, 
              pca_inj = p.inj, pca_inj_qc = p.inj.qc))
}

# function that applies CPCA correction to scaled data
data_normalization <- function(proc_object, ncomps = 3){
  imp_data <- proc_object$log_transf
  imp_data$class[imp_data$class %in% "sample"] <- "sa"
  imp_data$class[imp_data$class %in% "QC"] <- "qc"
  rownames(imp_data) <- paste(imp_data$class, imp_data$sampleID, sep = "_")
  data.proc <- cpca.driftRemoval(df = imp_data[,-c(1,2)], nComps = ncomps)
  data.proc <- data.frame(sampleID = rownames(data.proc),
                          class = rownames(data.proc), data.proc)
  data.proc$sampleID <- gsub(pattern = "sa_", replacement = "", 
                             x = data.proc$sampleID)
  data.proc$sampleID <- gsub(pattern = "qc_", replacement = "", 
                             x = data.proc$sampleID)
  data.proc$class[grepl("sa_", data.proc$class)] <- "sample"
  data.proc$class[grepl("qc_", data.proc$class)] <- "QC"
  rownames(data.proc) <- data.proc$sampleID
  pca.norm <- autoplot(prcomp(data.proc[,-c(1,2)], center = TRUE, scale. = TRUE), 
                       data = data.proc, colour = "class", scale = FALSE,
                       size = 2) +
    theme_bw() + scale_color_brewer(palette = "Set1")
  return(list(norm_data = data.proc, pca_norm = pca.norm))
}

# Function that returns the features to keep after CV-based filter
CV_filtering <- function(proc_object){
  data_imp <- proc_object$scaled
  cv.feats <- llply(3:ncol(data_imp), function(j){
    cv.sample <- raster::cv(data_imp[data_imp$class %in% "sample",j], na.rm = TRUE)
    cv.qc <- raster::cv(data_imp[data_imp$class %in% "QC",j], na.rm = TRUE)
    cv.sample>cv.qc
  }) %>% unlist
  names(cv.feats) <- colnames(data_imp)[-c(1,2)]
  return(cv.feats)
}

# Function that applies a logistic regression model for each LC-MS feature,
# corrects the p-values and outputs the resulting data frame and a summary table
statistical_analysis <- function(df, feats.names, covars, outcome, 
                                 do.Par = TRUE, nCore = 10, psign = 0.05){
  doParallel::registerDoParallel(nCore)
  pvals_df <- plyr::ldply(feats.names, function(soi){
    f_regr <- as.formula(paste0(outcome,"~",paste(covars, collapse = "+"),"+",soi))
    gl.model <- glm(f_regr, data = df, family = "binomial")
    s <- summary(gl.model)
    OR <- exp(cbind(OR = coef(gl.model), confint(gl.model)))[soi,]
    #OR <- exp(s$coefficients[j,1])
    pval <- s$coefficients[soi,4]
    data.frame(Feature = soi, pval = pval, 
               OR = OR[1], CI2.5 = OR[2], 
               CI97.5 = OR[3], n = nobs(gl.model))
  },.parallel = do.Par)
  pvals_df$qvalue_trunc <- qvalue_truncp(pvals_df$pval)$qvalues
  pvals_df$qvalue <- qvalue(pvals_df$pval)$qvalues
  pvals_df$p_fdr <- p.adjust(pvals_df$pval, "fdr")
  summary.stats <- data.frame(Variable_interest = outcome,
                              Model_type = "logistic regression",
                              sample_size = mean(pvals_df$n),
                              Number_features = length(feats.names),
                              Significant_features = sum(pvals_df$pval<psign),
                              Significant_features_fdr = sum(pvals_df$p_fdr<psign),
                              Significant_features_qvalue = sum(pvals_df$qvalue<psign),
                              Significant_features_qvalue_trunc = sum(pvals_df$qvalue_trunc<psign),
                              Min_qvalue = min(pvals_df$qvalue),
                              Sign_p = psign)
  p.hist <- ggplot(aes(x = pval), data = pvals_df) + geom_histogram()
  return(list(summary.stats = summary.stats,
              Pvals.df = pvals_df, p.hist = p.hist))
}

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
  pvals_df$qvalue_trunc <- qvalue_truncp(pvals_df$pval)$qvalues
  pvals_df$qvalue <- qvalue(pvals_df$pval)$qvalues
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
  p.hist <- ggplot(aes(x = pval), data = pvals_df) + geom_histogram()
  return(list(summary.stats = summary.stats,
              Pvals.df = pvals_df, p.hist = p.hist))
}

# function that applies a plsda model and outputs the pls plot and a data frame
# with a VIP for each feature
pls_stats <- function(x, y, ncomp = 10, plot.title, ncomp.perf){
  pls_model <- mixOmics::plsda(X = x, Y = y, ncomp = ncomp)
  pall <- mixOmics::plotIndiv(pls_model, ind.names = FALSE,  legend = TRUE)
  df.all <- pall$df
  pall.gg <- ggplot(df.all, aes(x = x, y = y, colour = group)) + 
    geom_point(size = 0.8) + 
    stat_ellipse(geom = "polygon", aes(fill = group), alpha = 0.25, show.legend = F) +
    labs(colour = "Incident DM", x = pall$graph$labels$x, 
         y = pall$graph$labels$y, title = plot.title) +
    scale_color_manual(values = c("cornflowerblue", "coral"))+
    scale_fill_manual(values = c("cornflowerblue", "coral")) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  vipAll <- as.data.frame(mixOmics::vip(pls_model))
  vipAll <- data.frame(Feats = rownames(vipAll), vipAll)
  vip.comp <- vipAll[,c(1,ncomp.perf+1)]
  colnames(vip.comp)[2] <- "comp"
  vip.comp <- arrange(vip.comp, -comp)
  return(list(plot = pall.gg,
              vip_df = vip.comp))
}

# Function that returns a plot to decide how many components should be considered
# in the plsda model
perf_ncomp <- function(x, y, ncomp = 10){
  pls_model <- mixOmics::plsda(X = x, Y = y, ncomp = ncomp)
  perf.plsda <- mixOmics::perf(pls_model, validation = "Mfold", folds = 5,
                               progressBar = FALSE, auc = TRUE, nrepeat = 10)
  return(plot(perf.plsda))
}

volcano_plot <- function(stats_results, pls_result){
  pvals <- stats_results$Pvals.df
  vips <- pls_result$vip_df
  colnames(vips)[1] <- "Feature"
  merge.df <- merge(pvals, vips, by = "Feature")
  merge.df$discriminant <- "Non-discriminant"
  merge.df$discriminant[merge.df$qvalue<0.05] <- "q<0.05"
  merge.df$discriminant[merge.df$comp>1] <- "VIP>1"
  merge.df$discriminant[merge.df$comp>1 & merge.df$qvalue<0.05] <- "VIP>1 and q<0.05"
  merge.df$discriminant <- as.factor(merge.df$discriminant)
  ggplot(data = merge.df, aes(x = -log10(qvalue), y = comp, colour = discriminant)) + 
    geom_point(alpha = 0.5) +
    labs(y = "VIP") +
    theme_bw() + 
    scale_colour_manual(values = c("skyblue", 
                                   "goldenrod1","mediumpurple1"))
}