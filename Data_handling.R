library(GEOquery)
library(magrittr)

working_path <- "/home/maria/BERN/"
setwd(working_path)
source("analysis_functions.R")

# GEO download and preparation
setwd("data/")
gds <- getGEO("GSE96058")
data1 <- gds$`GSE96058-GPL11154_series_matrix.txt.gz`
data2 <- gds$`GSE96058-GPL18573_series_matrix.txt.gz`

# phenotypic data
pheno_data <- data1@phenoData@data
nrow(pheno_data)
pheno_data2 <- data2@phenoData@data
nrow(pheno_data2)
# binding of platforms data
pheno_data_total <- rbind(pheno_data, pheno_data2)
# pheno data preparation
pheno_data_total <- pheno_data_preparation(pheno_data = pheno_data_total)
# save cleaned and prepare data
write.csv(pheno_data_total, "/home/maria/BERN/data/pheno_data.csv")

# Gene expression data - saving to RDS to speed up the process - only need to run it once

# system.time({
# expr_data <- read.csv(gzfile("/home/maria/BERN/data/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"), row.names = 1)
# })
# user  system elapsed 
# 83.374   0.053  83.395 
#saveRDS(expr_data, "/home/maria/BERN/data/GSE96058_expr_data.rds")

system.time({
  expr_data <- readRDS("/home/maria/BERN/data/GSE96058_expr_data.rds")
})
# user  system elapsed 
# 2.617   0.058   2.675

#pheno_data <- read.csv("/home/maria/BERN/data/pheno_data.csv", row.names = 1)

expr_data_t <- as.data.frame(t(expr_data))
# filter based on gene variability
filtered_expr_data <- variability_based_filtering(expr_data = expr_data_t, percentile = 0.2)

# PCA
pca_results <- pca_comput(expr_data = filtered_expr_data, pheno_data = pheno_data_total,
           center = TRUE, scale = TRUE, k = 10)

# QC PCA plots
pca_plot(pca_results = pca_results, color_variable = "Replicates",
         variable_nameplot = "Replicates", variable_type = "factor")

# PCA plots

## REMOVE REPLICATES

pca_plot(pca_results = pca_results, color_variable = "age.at.diagnosis.ch1",
         variable_nameplot = "Age at diagnosis", variable_type = "numeric")

pca_plot(pca_results = pca_results, color_variable = "chemo.treated.ch1",
         variable_nameplot = "Chemotherapy", variable_type = "factor")


# analysis
doParallel::registerDoParallel(20)
pvals_lymph_node <- linear_model_analysis(expr_data = filtered_expr_data, 
                                          pheno_data = pheno_data_total, 
                                          covars = c("chemo.treated.ch1","endocrine.treated.ch1", 
                                                     "age.at.diagnosis.ch1", "tumor.size.ch1"),
                                          var_interest = "lymph.node.status.ch1", 
                                          corr_method = "fdr", 
                                          do.Par = TRUE)

pvals_survival_cox <- cox_model_analysis(expr_data = filtered_expr_data,
                                         pheno_data = pheno_data_total,
                                         covars = c("chemo.treated.ch1","endocrine.treated.ch1", 
                                                    "age.at.diagnosis.ch1", "tumor.size.ch1",
                                                    "er.status.ch1", "her2.status.ch1",
                                                    "ki67.status.ch1", "nhg.ch1", "pgr.status.ch1",
                                                    "lymph.node.group.ch1"), 
                                         corr_method = "fdr", do.Par = TRUE)

# enrichment
enrichment_results <- enrichment_analysis(gene_list = pvals_lymph_node$Gene[pvals_lymph_node$p_adj<0.05])
enrichment_results$goplot


# random survival forest
expr_data_t <- expr_data_t[!grepl(pattern = "repl", x = rownames(expr_data_t)),]
pheno_data_total <- pheno_data_total[!grepl(pattern = "repl", x = rownames(pheno_data_total)),]

pheno_data_total$overall.survival.event.ch1 <- as.numeric(as.character(pheno_data_total$overall.survival.event.ch1))

iter_rsf_results <- plyr::llply(c(0.1,0.25,0.50,0.75,0.95, 0.99), function(q){
  filtered_expr_data <- variability_based_filtering(expr_data = expr_data_t, percentile = q)
  name <- unlist(strsplit(as.character(q), split = "[.]"))[2]
  rsf_results <- random_survival_forest(expr_data = filtered_expr_data, 
                                        pheno_data = pheno_data_total, 
                                        covars = c("age.at.diagnosis.ch1", "chemo.treated.ch1", 
                                                   "endocrine.treated.ch1","lymph.node.group.ch1", 
                                                   "tumor.size.ch1", "overall.survival.days.ch1", 
                                                   "overall.survival.event.ch1"), 
                                        path = paste0("/home/maria/BERN/ranger_survival_model_",
                                                      name,".rds"),
                                        num_threads = 15, seed = 2607)
  rsf_results
})


iter_rsf_results[[1]]$concordance_index$c.index
iter_rsf_results[[2]]$concordance_index$c.index
iter_rsf_results[[3]]$concordance_index$c.index
iter_rsf_results[[4]]$concordance_index$c.index

iter_rsf_results[[1]]$concordance_index_baseline$c.index
iter_rsf_results[[2]]$concordance_index_baseline$c.index
iter_rsf_results[[3]]$concordance_index_baseline$c.index
iter_rsf_results[[4]]$concordance_index_baseline$c.index

filtered_expr_data <- variability_based_filtering(expr_data = expr_data_t, percentile = 0.1)
pvals_survival_cox <- cox_model_analysis(expr_data = filtered_expr_data,
                                         pheno_data = pheno_data,
                                         covars = c("age.at.diagnosis.ch1", "chemo.treated.ch1", 
                                                    "endocrine.treated.ch1","lymph.node.group.ch1", 
                                                    "tumor.size.ch1"), 
                                         corr_method = "fdr", do.Par = TRUE)
pvals_survival_cox <- pvals_survival_cox[!is.na(pvals_survival_cox$p_adj),]

write.csv(pvals_survival_cox, "/home/maria/BERN/cox_results.csv")

features_sign <- pvals_survival_cox$Feature[(pvals_survival_cox$p_adj<0.05) & (pvals_survival_cox$GLOBAL>0.05)]
filtered_expr_data <- expr_data_t[,features_sign]

rsf_results <- random_survival_forest(expr_data = filtered_expr_data, pheno_data = pheno_data, 
                       covars = c("age.at.diagnosis.ch1", "chemo.treated.ch1", 
                                  "endocrine.treated.ch1","lymph.node.group.ch1", 
                                  "tumor.size.ch1", "overall.survival.days.ch1", 
                                  "overall.survival.event.ch1"), 
                       path = "/home/maria/BERN/ranger_survival_model_cox.rds",
                       num_threads = 15, seed = 2607)

rsf_results$concordance_index$c.index
rsf_results$VIP[1:100]
rsf_results$concordance_index_baseline$c.index



