# Transcriptomic-Analysis-of-Breast-Ductal-Carcinoma-GSE21422-
This repository provides an open-source transcriptomic analysis workflow of the GEO dataset GSE21422, profiling ductal carcinoma in situ (DCIS) and invasive ductal carcinoma (IDC). Using R and Bioconductor, the project covers data preprocessing, normalization, differential gene expression (limma), and visualization (volcano plot)

# Install Packages
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

 library(readxl)
#GSE21422_series_matrix <- read_excel("GSE21422.xlsx")

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
##### DGE
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
  labs(title = "Volcano Plot: GSE21422", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")
  theme_minimal() +

  
<img width="1200" height="1000" alt="volcano_plot_GSE21422" src="https://github.com/user-attachments/assets/d4539880-1794-4fcc-961e-313e6adc038d" />
<img width="1200" height="800" alt="boxplot_after_norm" src="https://github.com/user-attachments/assets/78f45ff4-fc48-4cde-a713-02c1cb3c605c" />
<img width="1200" height="800" alt="boxplot_before_norm" src="https://github.com/user-attachments/assets/f445284c-6902-4292-8494-acb6dc8fd4bd" />

  labs(title = "Volcano Plot: GSE21422", x = "Log 2 Fold Change", y = "-log10 Adjusted P-Value")
