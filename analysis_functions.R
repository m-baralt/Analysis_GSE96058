
pheno_data_preparation <- function(pheno_data){
  # detection of column names without variability
  cols_to_keep <- plyr::llply(1:ncol(pheno_data), .inform = T, function(j){
    col_var <- pheno_data[,j]
    if (is.character(col_var) | is.factor(col_var)){
      col_var <- as.factor(pheno_data[,j])
      if (length(levels(col_var))>1) return(colnames(pheno_data)[j])
    } else {
      if (var(col_var, na.rm = TRUE)>0) return(colnames(pheno_data)[j])
    }
  }) %>% unlist()
  
  # characteristics and prediction columns are removed
  cols_to_keep <- cols_to_keep[!grepl(pattern = "characteristics", 
                                      x = cols_to_keep)]
  cols_to_keep <- cols_to_keep[!grepl(pattern = "prediction", 
                                      x = cols_to_keep)]
  cols_to_keep <- cols_to_keep[!(cols_to_keep %in% c("geo_accession", "last_update_date", 
                                                     "instrument_model", "relation", 
                                                     "relation.1", "scan-b external id:ch1",
                                                     "platform_id"))]
  
  # Variables are converted to their proper classes
  pheno_data <- pheno_data[,cols_to_keep]
  rownames(pheno_data) <- pheno_data$title
  pheno_data <- pheno_data[,-1]
  
  colnames_numeric <- c("overall survival days:ch1", 
                        "age at diagnosis:ch1",
                        "tumor size:ch1")
  
  for (col in colnames_numeric){
    pheno_data[,col] <- as.numeric(as.character(pheno_data[,col]))
  }
  
  colnames_factor <- c("instrument model:ch1", "lymph node group:ch1",
                       "lymph node status:ch1", "nhg:ch1",
                       "pam50 subtype:ch1")
  
  for (col in colnames_factor){
    pheno_data[,col][pheno_data[,col]=="NA"] <- NA
    pheno_data[,col] <- as.factor(as.character(pheno_data[,col]))
  }
  
  colnames_rest <- colnames(pheno_data)[!(colnames(pheno_data) %in% c(colnames_numeric, colnames_factor))]
  for (col in colnames_rest){
    pheno_data[,col] <- as.factor(as.character(as.numeric(pheno_data[,col])))
  }
  colnames(pheno_data) <- make.names(colnames(pheno_data))
  return(pheno_data)
  
}


# expr data is a gene expression matrix with samples as rows and genes as columns
variability_based_filtering <- function(expr_data, percentile = 0.1){
  gene_vars <- apply(expr_data, 2, var)
  threshold <- quantile(gene_vars, probs = percentile)
  filtered_expr <- expr_data[, gene_vars > threshold]
  return(filtered_expr)
}

pca_comput <- function(expr_data, pheno_data, center = TRUE, scale = TRUE, k = 10){
  pca_result <- RSpectra::svds(scale(as.matrix(expr_data), 
                                     center = center, scale = scale), 
                               k = k)
  
  pca_u <- pca_result$u
  rownames(pca_u) <- rownames(expr_data)
  colnames(pca_u) <- paste0("PC", 1:k)
  pca_u <- as.data.frame(pca_u)
  
  singular_values <- pca_result$d
  explained_variance <- singular_values^2 / sum(singular_values^2)
  pcs_var <- data.frame(PC = paste0("PC",1:10), var = explained_variance[1:10])
  pcs_var$PC <- factor(pcs_var$PC, levels = unique(pcs_var$PC))
  
  merge_pca <- merge(pheno_data, pca_u, by = 0)
  rownames(merge_pca) <- merge_pca$Row.names
  merge_pca <- merge_pca[,-1]
  merge_pca <- data.frame(Replicates = as.factor(grepl(pattern = "repl", x = rownames(merge_pca))*1), merge_pca)
  
  return(list(PCA_pheno = merge_pca,
              Explained_variance = pcs_var))
}

pca_plot <- function(pca_results, color_variable, variable_nameplot, variable_type){
  
  PCA_pheno <- pca_results$PCA_pheno
  Explained_variance <- pca_results$Explained_variance
  
  # scale color
  custom_colors <- c(
    "#FFCA00",
    "#FF4500",  
    "#7B1FA2",  
    "#0288D1",  
    "#7A0000"  
  )
  
  if (variable_type=="numeric") {
    color_scale <- scale_color_gradientn(colors = c("yellow", "orange", "red", "darkred"))
  } else {
    color_scale <- scale_color_manual(values = custom_colors)
  }
  
  p <- ggplot(aes(x = PC1, y = PC2, color = !! sym(color_variable)), data = PCA_pheno) +
    geom_point(alpha = 0.7, size = 1.5) +
    theme_bw() +
    labs(x = paste0("PC1 (", round(Explained_variance$var[1]*100, 2), "%)"),
         y = paste0("PC2 (", round(Explained_variance$var[2]*100, 2), "%)"),
         color = variable_nameplot) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.25, "cm"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 14, 
                                  margin = margin(r = 18, b = 10))
    ) +
    color_scale
  return(p)
}

# stattistical analysis
linear_model_analysis <- function(expr_data, pheno_data, covars, 
                                  var_interest, corr_method = "fdr",
                                  do.Par = TRUE){
  merged_data <- merge(pheno_data, expr_data, by = 0)
  rownames(merged_data) <- merged_data$Row.names
  merged_data <- merged_data[,-1]
  
  pvals_df <- plyr::ldply(colnames(expr_data), function(gene){
    formula_str <- as.formula(paste0("`", gene, "`", "~",
                                     paste(c(covars, var_interest), 
                                           collapse = "+")))
    model_lm <- lm(formula_str, data = merged_data)
    s <- summary(model_lm)
    stats_gene <- s$coefficients[grepl(pattern = var_interest, x = rownames(s$coefficients)),c(1,4)]
    data.frame(Gene = gene, t(stats_gene))
  },.parallel = do.Par)
  colnames(pvals_df)[3] <- "pvalue"
  pvals_df$p_adj <- p.adjust(pvals_df$pvalue, corr_method)
  return(pvals_df)
}

# survival analysis
library(survival)
library(survminer)
is_valid_name <- function(x) {
  make.names(x) == x
}
cox_model_analysis <- function(expr_data, pheno_data, covars, 
                               corr_method = "fdr",
                               do.Par = TRUE){
  merged_data <- merge(pheno_data, expr_data, by = 0)
  rownames(merged_data) <- merged_data$Row.names
  merged_data <- merged_data[,-1]
  merged_data$overall.survival.event.ch1 <- as.numeric(as.character(merged_data$overall.survival.event.ch1))
  
  pvals_cox <- plyr::ldply(colnames(expr_data), .inform = TRUE, function(gene){
    right_part <- paste(covars, collapse = "+")
    formula_part <- as.formula(paste0("Surv(overall.survival.days.ch1, overall.survival.event.ch1) ~","`",
                                      gene, "` + ", right_part))
    res.cox <- coxph(formula_part, data = merged_data)
    s <- summary(res.cox)
    if (!((sum(is.infinite(s$coefficients))>0) | (sum(is.infinite(s$conf.int))>0))){
      PH <- cox.zph(res.cox)
      PH.assum <- t(PH$table[,3])
      colnames(PH.assum)[1] <- "gene"
      
      #gene_name <- rownames(s$coefficients)[grepl(pattern = gene, rownames(s$coefficients))]
      gene_name <- ifelse(is_valid_name(gene), gene, paste0("`",gene,"`"))
      nSample <- s$n
      pval <- s$coefficients[gene_name,5]
      HR <- s$coefficients[gene_name,2]
      SE <- s$coefficients[gene_name,3]
      data.frame(Feature = gene, pval = pval, HR = HR, SE = SE, nSample = nSample, PH.assum)
    }
  },.parallel = TRUE)
  
  pvals_cox$p_adj <- p.adjust(pvals_cox$pval, corr_method)
  return(pvals_cox)
}

# pathway enrichment analysis
library(gprofiler2)
enrichment_analysis <- function(gene_list){
  gostres <- gost(
    query = gene_list,
    organism = "hsapiens",     
    sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"),  
    significant = TRUE, 
    correction_method = "gSCS" # built-in multiple testing correction
  )
  p <- gostplot(gostres, capped = TRUE, interactive = TRUE)
  return(list(enrichment_result = gostres$result, goplot = p))
}

# pathview - does not work with the vpn
# library(pathview)
# data(gse16873.d)
# pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
#                    species = "hsa", out.suffix = "gse16873")


# biclustering
# library(biclust)
# filtered_expr_data <- variability_based_filtering(expr_data = expr_data_t, percentile = 0.95)
# expr_data_biclust <- na.omit(filtered_expr_data)
# expr_data_biclust <- t(expr_data_biclust)
# res <- biclust(expr_data_biclust, method = BCBimax(), minr = 5, minc = 5)
# res@Number
# genes_in_i <- which(res@RowxNumber[, 70])
# gene_names <- rownames(expr_data_biclust)[genes_in_i]
# 
# View(res@RowxNumber)
# View(res@RowxNumber)

# random survival forest

#expr_data <- variability_based_filtering(expr_data = expr_data, percentile = 0.5)

library(ranger)
random_survival_forest <- function(expr_data, pheno_data, covars, path, 
                                   num_threads = 10, seed = 5){
  set.seed(seed)
  merged_data <- merge(pheno_data, expr_data, by = 0)
  rownames(merged_data) <- merged_data$Row.names
  merged_data <- merged_data[,-1]
  survival_merged_expr <- merged_data[,c(covars, colnames(expr_data))]
  survival_merged_expr <- as.data.frame(na.omit(survival_merged_expr))
  
  train_data <- survival_merged_expr[sample(1:nrow(survival_merged_expr), 
                                            nrow(survival_merged_expr)*0.8),]
  test_data <- survival_merged_expr[!(rownames(survival_merged_expr) %in% rownames(train_data)),]
  
  x_train <- train_data[, -which(names(train_data) %in% c("overall.survival.days.ch1", "overall.survival.event.ch1"))]
  y_train <- Surv(train_data$overall.survival.days.ch1, train_data$overall.survival.event.ch1)
  
  model_surv <- ranger(
    dependent.variable.name = NULL,
    y = y_train,
    x = x_train,
    num.trees = 500,
    importance = "impurity",
    splitrule = "logrank",
    num.threads = num_threads
  )
  x_test <- test_data[, -which(names(test_data) %in% 
                                 c("overall.survival.days.ch1", 
                                   "overall.survival.event.ch1"))]
  y_test <- Surv(test_data$overall.survival.days.ch1, 
                 test_data$overall.survival.event.ch1)
  
  pred <- predict(model_surv, data = x_test)
  risk_score <- pred$chf[, ncol(pred$chf)]
  c_index <- survcomp::concordance.index(x = risk_score, 
                               surv.time = test_data$overall.survival.days.ch1,
                               surv.event = test_data$overall.survival.event.ch1,
                               method = "noether")
  importance_vals <- model_surv$variable.importance
  importance_vars <- sort(importance_vals, decreasing = TRUE)
  
  surv_formula <- as.formula(paste0("Surv(overall.survival.days.ch1, overall.survival.event.ch1) ~ ", 
                             paste(covars[-which(covars %in% c("overall.survival.days.ch1", "overall.survival.event.ch1"))], 
                                   collapse = "+")))
  model_surv_baseline <- ranger(
    formula = surv_formula,
    data = train_data,
    num.trees = 500,
    importance = "impurity",
    splitrule = "logrank",
    num.threads = num_threads
  )
  pred <- predict(model_surv_baseline, data = x_test)
  risk_score <- pred$chf[, ncol(pred$chf)]
  c_index_baseline <- survcomp::concordance.index(x = risk_score, 
                                         surv.time = test_data$overall.survival.days.ch1,
                                         surv.event = test_data$overall.survival.event.ch1,
                                         method = "noether")
  
  total_result <- list(model = model_surv,
                       concordance_index = c_index, 
                       VIP = importance_vars, 
                       concordance_index_baseline = c_index_baseline)
  
  saveRDS(total_result, file = path)
  return(list(concordance_index = c_index, VIP = importance_vars, 
              concordance_index_baseline = c_index_baseline))
}









