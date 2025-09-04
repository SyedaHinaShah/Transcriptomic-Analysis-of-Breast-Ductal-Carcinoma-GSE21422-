# Transcriptomic-Analysis-of-Breast-Ductal-Carcinoma-GSE21422-
This repository provides an open-source transcriptomic analysis workflow of the GEO dataset GSE21422, profiling ductal carcinoma in situ (DCIS) and invasive ductal carcinoma (IDC). Using R and Bioconductor, the project covers data preprocessing, normalization, differential gene expression (limma), and visualization (volcano plot)
str(GSE21422)
# Install if not installed
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("preprocessCore")
}
library(preprocessCore)

# Your dataset: GSE21422_series_matrix (tibble with ID_REF + 20 columns)

# Extract expression matrix (remove ID_REF)
expr_matrix <- as.matrix(GSE21422[ , -1])  # numeric matrix of expression values

# Check if log2 transform is needed (usually microarray data are log2 scale)
# Basic heuristic: if max value > 100, probably not log-transformed yet
max_val <- max(expr_matrix, na.rm=TRUE)
if (max_val > 100) {
  expr_matrix <- log2(expr_matrix + 1)  # add 1 to avoid log(0)
  message("Data log2 transformed")
} else {
  message("Data looks already log2 transformed")
}

# Quantile Normalization
expr_norm <- normalize.quantiles(expr_matrix)

# Add back row names and column names
rownames(expr_norm) <- GSE21422$ID_REF
colnames(expr_norm) <- colnames(GSE21422)[-1]

# Convert to data frame for convenience
expr_norm_df <- as.data.frame(expr_norm)

# Optional: add ID_REF as column again
expr_norm_df$ID_REF <- rownames(expr_norm_df)

# View normalized data summary
summary(expr_norm_df[ , 1:5])  # view summary of first 5 samples

# Now expr_norm_df contains normalized expression values for downstream analysis
###boxplot
# Extract expression matrix before normalization
expr_raw <- as.matrix(GSE21422[ , -1])
rownames(expr_raw) <- GSE21422$ID_REF

# Log2 transform if needed (same logic as before)
max_val <- max(expr_raw, na.rm=TRUE)
if (max_val > 100) {
  expr_raw <- log2(expr_raw + 1)
}

# Boxplot before normalization
boxplot(expr_raw,
        main = "Boxplot Before Normalization",
        las=2,
        col = "red",
        outline=FALSE)

# Quantile normalization
expr_norm <- normalize.quantiles(expr_raw)
rownames(expr_norm) <- rownames(expr_raw)
colnames(expr_norm) <- colnames(expr_raw)

# Boxplot after normalization
boxplot(expr_norm,
        main = "Boxplot After Quantile Normalization",
        las=2,
        col = "green",
        outline=FALSE)
###updatedbefore working 
# -------------------------------
# Differential Expression Script
# -------------------------------

# 1. Load necessary packages
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("preprocessCore")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}
library(preprocessCore)
library(limma)
library(ggplot2)

# 2. Load your dataset (already loaded in your case)
# If not, load it using:
# library(readxl)
# GSE21422_series_matrix <- read_excel("GSE21422.xlsx")

# 3. Extract expression matrix
expr_raw <- as.matrix(GSE21422[ , -1])  # Remove ID_REF
rownames(expr_raw) <- GSE21422$ID_REF

# 4. Log2 transform (if needed)
if (max(expr_raw, na.rm = TRUE) > 100) {
  expr_raw <- log2(expr_raw + 1)
  message("Log2 transformation applied")
} else {
  message("Data already log-transformed")
}

# 5. Plot boxplot before normalization
boxplot(expr_raw,
        main = "Boxplot Before Normalization",
        las = 2,
        col = "lightblue",
        outline = FALSE)

# 6. Quantile normalization
expr_norm <- normalize.quantiles(expr_raw)
rownames(expr_norm) <- rownames(expr_raw)
colnames(expr_norm) <- colnames(expr_raw)

# 7. Plot boxplot after normalization
boxplot(expr_norm,
        main = "Boxplot After Quantile Normalization",
        las = 2,
        col = "lightgreen",
        outline = FALSE)

# 
#####DGE
# 1. Load libraries
library(preprocessCore)
library(limma)
library(ggplot2)
library(pheatmap)

# 2. Prepare expression matrix
expr_raw <- as.matrix(GSE21422_series_matrix[ , -1])  # remove first column (IDs)
rownames(expr_raw) <- GSE21422_series_matrix$ID_REF

# 3. Log2 transform if needed
if (max(expr_raw, na.rm = TRUE) > 100) {
  expr_raw <- log2(expr_raw + 1)
  message("Log2 transformation applied")
}

# 4. Quantile normalization
expr_norm <- normalize.quantiles(expr_raw)
rownames(expr_norm) <- rownames(expr_raw)
colnames(expr_norm) <- colnames(expr_raw)

# 5. Define group labels (adjust if needed)
group <- factor(c(rep("Control", 9), rep("Treatment", 10)))

# 6. Differential expression using limma
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_norm, design)
contrast.matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 7. Get DEGs
deg_results <- topTable(fit2, number = Inf, adjust.method = "BH")
write.csv(deg_results, "DEG_results_GSE21422.csv")

# 8. Volcano Plot
deg_results$Significant <- with(deg_results, ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: GSE21422", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")

# 9. Heatmap of top 50 DEGs
top_genes <- rownames(deg_results)[1:50]
heatmap_data <- expr_norm[top_genes, ]
heatmap_data <- t(scale(t(heatmap_data)))  # z-score scaling

pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = data.frame(Group = group),
         show_rownames = FALSE,
         main = "Top 50 DEGs Heatmap")

######
deg_results_clean <- deg_results[!is.na(deg_results$logFC) & !is.na(deg_results$adj.P.Val), ]
deg_results_clean$Significant <- with(deg_results_clean, ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))

ggplot(deg_results_clean, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: GSE21422", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")
###heatmap
annotation_col = data.frame(Group = group)
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr_norm)  # set rownames to sample names

pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 DEGs Heatmap")
##
# Example: save boxplot before normalization
png("boxplot_before_norm.png", width = 1200, height = 800, res = 150)
boxplot(expr_raw,
        main = "Boxplot Before Normalization",
        las = 2,
        col = "lightblue",
        outline = FALSE)
dev.off()

# Example: save boxplot after normalization
png("boxplot_after_norm.png", width = 1200, height = 800, res = 150)
boxplot(expr_norm,
        main = "Boxplot After Quantile Normalization",
        las = 2,
        col = "lightgreen",
        outline = FALSE)
dev.off()

# Save Volcano Plot
png("volcano_plot_GSE21422.png", width = 1200, height = 1000, res = 150)
ggplot(deg_results_clean, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: GSE21422", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")
dev.off()

# Save Heatmap
png("heatmap_top50_GSE21422.png", width = 1200, height = 1000, res = 150)
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 DEGs Heatmap")
dev.off()
