###############################################################
# Project : Differential Gene Expression Analysis in Lung Cancer
###############################################################

###############################################################
# STEP 1 : Load Required Packages
###############################################################

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)

###############################################################
# STEP 2 : Download GEO Dataset
###############################################################
###############################################################
# STEP 2 : Download GEO Dataset
###############################################################

if (!exists("gse")) {
  gse <- getGEO(
    "GSE19804",
    GSEMatrix = TRUE,
    getGPL = FALSE
  )
}

if (!exists("lung_data")) {
  lung_data <- gse[[1]]
}

###############################################################
# STEP 3 : Extract ExpressionSet
###############################################################

lung_data <- gse[[1]]

###############################################################
# STEP 4 : View Dataset Information
###############################################################

lung_data

sampleNames(lung_data)

pData(lung_data)
###############################################################
# STEP 5 : Extract Sample Information
###############################################################

sample_info <- pData(lung_data)

head(sample_info$title)

tail(sample_info$title)
###############################################################
# STEP 6 : Extract Gene Expression Matrix
###############################################################

expr_matrix <- exprs(lung_data)
###############################################################
# STEP 7 : Identify Cancer and Normal Samples
###############################################################

sample_titles <- sample_info$title

head(sample_titles)

tail(sample_titles)
###############################################################
# STEP 8 : Create Sample Groups
###############################################################

group <- ifelse(
  grepl("Cancer", sample_titles),
  "Cancer",
  "Normal"
)

group <- factor(group)

group
###############################################################
# STEP 9 : Count Cancer and Normal Samples
###############################################################

table(group)
###############################################################
# STEP 10 : Create Design Matrix
###############################################################

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

design

head(design)
###############################################################
# STEP 11 : Create Contrast Matrix
###############################################################

contrast.matrix <- makeContrasts(
  Cancer_vs_Normal = Cancer - Normal,
  levels = design
)

contrast.matrix
###############################################################
# STEP 12 : Fit Linear Model
###############################################################

fit <- lmFit(expr_matrix, design)

fit
###############################################################
# STEP 13 : Apply Contrast
###############################################################

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2
###############################################################
# STEP 14 : Apply Empirical Bayes
###############################################################

fit2 <- eBayes(fit2)

fit2
###############################################################
# STEP 15 : Extract Differentially Expressed Genes
###############################################################

results <- topTable(
  fit2,
  adjust = "fdr",
  number = Inf
)

head(results)

dim(results)
###############################################################
# STEP 16 : Save Differential Expression Results
###############################################################

write.csv(
  results,
  file = "results/Differential_Expression_Results.csv",
  row.names = TRUE
)
###############################################################
# STEP 17 : Filter Significant Genes
###############################################################

significant_genes <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
]

head(significant_genes)

dim(significant_genes)
###############################################################
# STEP 18 : Save Significant Genes
###############################################################

write.csv(
  significant_genes,
  file = "results/Significant_Genes.csv",
  row.names = TRUE
)
###############################################################
# STEP 19 : Count Upregulated and Downregulated Genes
###############################################################

upregulated <- significant_genes[
  significant_genes$logFC > 1,
]

downregulated <- significant_genes[
  significant_genes$logFC < -1,
]

cat("Upregulated genes :", nrow(upregulated), "\n")
cat("Downregulated genes :", nrow(downregulated), "\n")
###############################################################
# STEP 20 : Create Volcano Plot
###############################################################

results$Significance <- "Not Significant"

results$Significance[
  results$adj.P.Val < 0.05 &
    results$logFC > 1
] <- "Upregulated"

results$Significance[
  results$adj.P.Val < 0.05 &
    results$logFC < -1
] <- "Downregulated"

volcano_plot <- ggplot(
  results,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Significance
  )
) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

volcano_plot
###############################################################
# STEP 21 : Heatmap of Top 50 Differentially Expressed Genes
###############################################################

top50 <- rownames(results)[1:50]

heatmap_data <- expr_matrix[top50, ]

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)
###############################################################
# STEP 22 : Analysis Summary
###############################################################

cat("=====================================\n")
cat(" Lung Cancer Differential Expression\n")
cat("=====================================\n\n")

cat("Total genes analyzed :", nrow(results), "\n")
cat("Significant genes    :", nrow(significant_genes), "\n")
cat("Upregulated genes    :", nrow(upregulated), "\n")
cat("Downregulated genes  :", nrow(downregulated), "\n")
###############################################################
# STEP 23 : Save Analysis Summary
###############################################################

summary_text <- c(
  "Lung Cancer Differential Expression Analysis",
  "",
  paste("Total genes analyzed :", nrow(results)),
  paste("Significant genes :", nrow(significant_genes)),
  paste("Upregulated genes :", nrow(upregulated)),
  paste("Downregulated genes :", nrow(downregulated))
)

writeLines(
  summary_text,
  "results/Analysis_Summary.txt"
)

cat("Analysis summary saved successfully!\n")