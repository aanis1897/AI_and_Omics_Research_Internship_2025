# install / load packages
# -------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# bioconductor packages needed (install if missing)
bioc_pkgs <- c("limma", "AnnotationDbi", "hugene10sttranscriptcluster.db", "GEOquery")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

# CRAN packages
cran_pkgs <- c("dplyr", "tibble", "ggplot2", "pheatmap")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# load libraries
library(GEOquery)
library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)    # annotation for Affymetrix Human Gene 1.0 ST
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# 2. load normalized expression data #
# use  pre-processed data object from previous steps

getwd()
list.files()

save(processed_data, pheno6244, file = "processed_GSE16558.RData")

load("processed_GSE16558.RData")

dim(processed_data)
head(rownames(processed_data))

# 3. Probe ID → Gene Symbol Mapping #
probe_ids <- rownames(processed_data)

# reinstall & reload

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi", force = TRUE)
library(AnnotationDbi)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hugene10sttranscriptcluster.db", force = TRUE)

library(hugene10sttranscriptcluster.db)

library(AnnotationDbi)

probe_ids <- rownames(processed_data)

gene_symbols <- mapIds(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# add gene symbols to expression data
# reinstall & reload
library(dplyr)
install.packages("tidyverse")
library(tidyverse)

processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  filter(!is.na(SYMBOL)) %>%
  relocate(SYMBOL, .after = PROBEID)

# 4. handle duplicate gene mappings #
# multiple probes may map to the same gene
duplicate_summary <- processed_data_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

# collapse duplicates using average expression
expr_matrix <- processed_data_df %>%
  select(-PROBEID, -SYMBOL)
averaged_data <- limma::avereps(as.matrix(expr_matrix),
                                ID = processed_data_df$SYMBOL)

# 5. define experimental groups #
# based on metadata: bone marrow, healthy (control) vs Bone marrow, MM (disease)
groups <- factor(pheno6244$source_name_ch1,
                 levels = c("Bone marrow, healthy", "Bone marrow, MM"),
                 labels = c("control", "myeloma"))
table(groups)

# 6. design matrix and contrast #
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)

fit_1 <- lmFit(averaged_data, design)
contrast_matrix <- makeContrasts(myeloma_vs_control = myeloma - control,
                                 levels = design)
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)
fit_2 <- eBayes(fit_contrast)

#### 7. extract DEGs ####
deg_results <- topTable(fit_2,
                        coef = "myeloma_vs_control",
                        number = Inf,
                        adjust.method = "BH")


# add regulation status
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
))

# split subsets
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

# save results
dir.create("Results", showWarnings = FALSE)
write.csv(deg_results, "results/DEG_all_results.csv")
write.csv(upregulated, "results/DEG_upregulated.csv")
write.csv(downregulated, "results/DEG_downregulated.csv")

cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# 8. volcano plot #
dir.create("plots", showWarnings = FALSE)

png("plots/volcano_plot_GSE16558.png", width = 2000, height = 1500, res = 300)
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Myeloma vs Healthy Plasma Cells",
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "Regulation")
dev.off()

# 9. heatmap of top 25 DEGs #

install.packages("pheatmap")
library(pheatmap)

top_genes <- head(rownames(deg_results[order(deg_results$adj.P.Val), ]), 25)
heatmap_data <- averaged_data[top_genes, ]

png("plots/heatmap_top25_GSE16558.png", width = 2000, height = 1500, res = 300)
pheatmap(heatmap_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 6,
         main = "Top 25 Differentially Expressed Genes – GSE16558")
dev.off()

