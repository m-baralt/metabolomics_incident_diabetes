# function to read and prepare the data
read_prepare_data <- function(path){
  dm <- readRDS(path)
  injectionOrder <- dm[[1]]
  injectionOrder <- injectionOrder[,c("SampleName", "injectionOrder", 
                                      "Batch", "Class", "subBatch", 
                                      "subName")]
  injectionOrder$Batch <- unlist(strsplit(unique(injectionOrder$Batch), 
                                          split = "_"))[3]
  injectionOrder$Batch <- substring(text = injectionOrder$Batch, 
                                    first = nchar(unique(injectionOrder$Batch)), 
                                    last = nchar(unique(injectionOrder$Batch)))
  
  injectionOrder$subBatch <- paste(injectionOrder$Batch, 
                                   injectionOrder$subBatch, sep = "_")
  subjects.samps <- injectionOrder$SampleName[injectionOrder$Class=="Sample"]
  compounds.df <- dm[[2]]
  norm.df <- dm[[5]]
  #norm.df <- norm.df[,colnames(norm.df) %in% paste0("X",subjects.samps)]
  metabs.names <- rownames(norm.df)
  norm.t.df <- as.data.frame(t(norm.df))
  norm.t.df <- data.frame(SampleName=rownames(norm.t.df), norm.t.df)
  norm.t.df$SampleName <- gsub(pattern = "X", replacement = "", 
                               norm.t.df$SampleName)
  norm.t.df <- merge(injectionOrder, norm.t.df, by = "SampleName")
  norm.t.df <- norm.t.df[!(norm.t.df$Class %in% "Blank"),]
  norm.t.df <- norm.t.df[!(norm.t.df$Class %in% "Cleaning"),]
  norm.t.df$SampleName[norm.t.df$Class %in% "QC"] <- 
    paste0("B",unique(norm.t.df$Batch),norm.t.df$SampleName[norm.t.df$Class %in% "QC"])
  rownames(norm.t.df) <- paste0("X",norm.t.df$SampleName)
  norm.df <- as.data.frame(t(norm.t.df[,-c(1:6)]))
  norm.df <- data.frame(metabolites = rownames(norm.df), norm.df)
  result.object <- list(
    Normalized_data = norm.df,
    Long_data = norm.t.df,
    injectionOrder = injectionOrder,
    compounds_info = compounds.df,
    metabs_names = metabs.names,
    samples = norm.t.df$subName
  )
  return(result.object)
}

# function to filter and impute the data
filter_impute <- function(metabs, df, var.batch, n.batch, perc.missings=0.2){
  raw.df <- df
  nas_batch <- plyr::llply(unique(df[,var.batch]), function(b){
    batch.df <- df[df[,var.batch] %in% b, metabs]
    dd <- data.frame(apply(batch.df, 2, function(c) sum(is.na(c))/nrow(batch.df)))
    colnames(dd) <- paste("NAs", b, sep = "_")
    dd <- data.frame(metabolites = rownames(dd), dd)
    return(dd)
  })
  nas_per_batch <- Reduce(function(x, y)
    merge(x, y, by = "metabolites"), nas_batch)
  metabs.filter <- plyr::ldply(nas_per_batch$metabolites, function(m){
    data.frame(metabolites = m, 
               NAs_total = sum(nas_per_batch[nas_per_batch$metabolites %in% m,-1]<perc.missings))
  })
  # Metabolites with less than 20% of missing values in all batches are treated together
  metabs.keep <- metabs.filter$metabolites[metabs.filter$NAs_total==n.batch]
  df <- df[,colnames(df) %in% c("subName","SampleName","injectionOrder",
                                "Batch","Class","subBatch",metabs.keep)]
  # impute minimum/2
  imputed_data <- df
  for (m in metabs.keep){
    imputed_data[is.na(imputed_data[,m]),m] <- min(imputed_data[,m], na.rm = T)/2
  }
  # logarithm
  log_data <- imputed_data
  log_data[,metabs.keep] <- log(log_data[,metabs.keep])
  
  result.object <- list(
    NAs_per_batch = nas_per_batch,
    raw_data = raw.df,
    filtered_data = df,
    filtered_imputed_data = imputed_data,
    log_transf = log_data,
    metabs.used = metabs.keep
  )
  return(result.object)
}

# function to plot pcas
pca_plot <- function(df, metabs, var.colour, text.title = NULL){
  pca <- prcomp(df[,metabs], center = TRUE, scale. = TRUE)
  
  if (is.null(text.title)){
    autoplot(pca, data = df, colour = var.colour, scale = FALSE, size = 2)
  } else {
    autoplot(pca, data = df, colour = var.colour, scale = FALSE, size = 2)+
      labs(title = text.title) + theme(plot.title = element_text(hjust = 0.5))
  }
}

# function to plot the distribution of the data
data_distribution <- function(df, metabs){
  df.log <- df
  df.log[,metabs] <- log(df.log[,metabs])
  df.sqrt <- df
  df.sqrt[,metabs] <- sqrt(df.sqrt[,metabs])
  hist.plots <- plyr::llply(metabs, function(m){
    p <- ggplot(df, aes_string(x=m)) + 
      geom_histogram(color="black", fill="white", bins = 50) + 
      theme(axis.title.y = element_blank())
    
    p.log <- ggplot(df.log, aes_string(x=m)) + 
      geom_histogram(color="black", fill="white", bins = 50) +
      theme(axis.title.y = element_blank()) +
      labs(x = paste0("log(",m,")"))
    
    p.sqrt <- ggplot(df.sqrt, aes_string(x=m)) +
      geom_histogram(color="black", fill="white", bins = 50) +
      theme(axis.title.y = element_blank()) +
      labs(x = paste0("sqrt(",m,")"))
    
    ggpubr::ggarrange(p,p.log, ncol=2)
    
  })
  return(hist.plots)
}

# function to plot the distribution of the corrected data
data_distribution_corr <- function(df, metabs){
  hist.plots <- plyr::llply(metabs, function(m){
    p <- ggplot(df, aes_string(x=m)) + 
      geom_histogram(color="black", fill="white", bins = 50) + 
      theme(axis.title.y = element_blank())
  })
  return(hist.plots)
}

# function for batch correction using combat
combat_correction <- function(df, metabs){
  corrected.data <- sva::ComBat(dat = t(df[,metabs]), 
                                batch = as.numeric(as.factor(df$subBatch)))
  corrected.data <- as.data.frame(t(corrected.data))
  corrected.data <- data.frame(SampleName = rownames(corrected.data), corrected.data)
  corrected.data$SampleName <- gsub(pattern = "X", replacement = "", x = corrected.data$SampleName)
  corrected.data <- merge(df[,!(colnames(df) %in% metabs)], corrected.data, by = "SampleName")
  rownames(corrected.data) <- corrected.data$SampleName
  return(corrected.data)
}

# function to joint all batches
joint_data <- function(list.dfs, list.injection){
  merged.data <- Reduce(function(x, y)
    merge(x, y, by = "metabolites"), list.dfs)
  injectionOrder <- do.call("rbind", list.injection)
  rownames(merged.data) <- merged.data$metabolites
  merged.t.data <- as.data.frame(t(merged.data[,-1]))
  merged.t.data <- data.frame(subName=rownames(merged.t.data), merged.t.data)
  merged.t.data$subName <- gsub(pattern = "X", replacement = "", 
                                merged.t.data$subName)
  injectionOrder <- injectionOrder[injectionOrder$subName %in% merged.t.data$subName,]
  merged.t.data <- merge(injectionOrder, merged.t.data, by = "subName")
  rownames(merged.t.data) <- paste0("X",merged.t.data$SampleName)
  return(merged.t.data)
}

# functions for statistical analysis
statistical_analysis <- function(df, metabolites, covars, outcome, 
                                 do.Par = TRUE, nCore = 10, psign = 0.05){
  doParallel::registerDoParallel(nCore)
  pvals_df <- plyr::ldply(metabolites, function(soi){
    f_regr <- as.formula(paste0(outcome,"~",paste(covars, collapse = "+"),"+",soi))
    gl.model <- glm(f_regr, data = df, family = "binomial")
    s <- summary(gl.model)
    OR <- exp(cbind(OR = coef(gl.model), confint(gl.model)))[soi,]
    #OR <- exp(s$coefficients[j,1])
    pval <- s$coefficients[soi,4]
    data.frame(Metabolite = soi, pval = pval, 
               OR = OR[1], CI2.5 = OR[2], 
               CI97.5 = OR[3], n = nobs(gl.model))
  },.parallel = do.Par)
  
  return(pvals_df)
}

statistical_analysis_lm <- function(df, metabolites, covars, outcome, 
                                    do.Par = TRUE, nCore = 10, psign = 0.05){
  doParallel::registerDoParallel(nCore)
  pvals_df <- plyr::ldply(metabolites, function(soi){
    f_regr <- as.formula(paste0(soi,"~",paste(covars, collapse = "+"),"+",outcome))
    lm.model <- lm(f_regr, data = df)
    s <- summary(lm.model)
    pval <- s$coefficients[outcome,4]
    data.frame(Feature = soi, pval = pval, 
               beta = s$coefficients[outcome,1],
               n = nobs(lm.model))
  },.parallel = do.Par)
  
  return(pvals_df)
}

geom_hpline <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHpline,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

GeomHpline <- ggproto("GeomHpline", GeomSegment,
                      required_aes = c("x", "y"),
                      non_missing_aes = c("size", "colour", "linetype", "width"),
                      default_aes = aes(
                        width = 0.5, colour = "black", size = 0.5, linetype = 1,
                        alpha = NA
                      ),
                      
                      draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                                            lineend = "butt", linejoin = "round", na.rm = FALSE) {
                        data <- mutate(data, x = x - width/2, xend = x + width, yend = y)
                        ggproto_parent(GeomSegment, self)$draw_panel(
                          data, panel_params, coord, arrow = arrow, arrow.fill = arrow.fill,
                          lineend = lineend, linejoin = linejoin, na.rm = na.rm
                        )
                      }
)
