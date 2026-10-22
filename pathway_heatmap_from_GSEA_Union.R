# ------ Pathway heatmap from GSEA result-----

# Install if needed
# install.packages(c("tidyverse", "readxl", "pheatmap", "openxlsx"))
# BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi"))

library(tidyverse)
library(readxl)
library(openxlsx)
library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)

# ----1. Set working directory--------
setwd("./RNAseq/Yanlin")

# -------2. Load count and metadata--------
count_data <- read_excel("gene_count.xlsx")
meta_data  <- read_xlsx("RNASeq_sPD_sample_info.xlsx")

meta_data$individual_ID <- as.factor(meta_data$individual_ID)
meta_data$treatment <- ifelse(grepl("DMSO", meta_data[[5]]), "DMSO", "TF")
meta_data$group <- factor(meta_data$donor, levels = c("control", "patient"))
meta_data$treatment <- factor(meta_data$treatment, levels = c("DMSO", "TF"))

# nested donor variable
meta_data$ind.n <- NA
patient_idx <- meta_data$group == "patient"
meta_data$ind.n[patient_idx] <- rep(1:4, each = 2)
control_idx <- meta_data$group == "control"
meta_data$ind.n[control_idx] <- rep(1:3, each = 2)
meta_data$ind.n <- factor(meta_data$ind.n)

meta_data$label <- factor(paste(meta_data$group, meta_data$treatment, sep = "_"))

# -------3. Build count matrix---------
gene_ids <- count_data$gene_id
count_matrix <- count_data[, 2:15] %>% as.data.frame()
rownames(count_matrix) <- count_data[[1]]

count_matrix <- count_matrix[, match(meta_data$`sample name`, colnames(count_matrix))]
stopifnot(identical(colnames(count_matrix), meta_data$`sample name`))

# ----4. Create DESeq2 object and vst----------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta_data,
  design = ~ ind.n + label
)

rownames(dds) <- gene_ids
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
batch <- colData(dds)$ind.n
design <- model.matrix(~ colData(vsd)$label)

assay(vsd) <- limma::removeBatchEffect(
  assay(vsd),
  batch = batch,
  design = design
)

expr_mat <- assay(vsd)
p <- plotPCA(vsd, intgroup = "treatment")

# -------5. Load GSEA workbook------------
gsea_file <- "GSEA_results_KEGG.xlsx"

# See sheet names
sheet_names <- getSheetNames(gsea_file)
print(sheet_names)

# -------6. USER INPUT---------

# Two comparisons:
sheet_DMSO <- "pat_TF_vs_DMSO"
sheet_TF   <- "ctl_TF_vs_DMSO"

selected_pathway <- "Spliceosome"

output_prefix <- "Heatmap_splice"

# ------7. Read selected GSEA sheet--------------
gsea_DMSO <- read.xlsx(gsea_file, sheet = sheet_DMSO)
gsea_TF   <- read.xlsx(gsea_file, sheet = sheet_TF)

# Check that pathway exists in both sheets
pathway_DMSO <- gsea_DMSO %>%
  filter(Description == selected_pathway)

pathway_TF <- gsea_TF %>%
  filter(Description == selected_pathway)

if (nrow(pathway_DMSO) == 0) {
  stop(paste("Pathway not found in sheet:", sheet_DMSO))
}

if (nrow(pathway_TF) == 0) {
  stop(paste("Pathway not found in sheet:", sheet_TF))
}


# -------8. Extract core enrichment genes-----------
# clusterProfiler gseKEGG stores genes in core_enrichment as ENTREZ IDs separated by "/"
core_DMSO <- unlist(
  strsplit(pathway_DMSO$core_enrichment[1], "/")
)

core_TF <- unlist(
  strsplit(pathway_TF$core_enrichment[1], "/")
)

core_DMSO <- unique(core_DMSO)
core_TF   <- unique(core_TF)

# Union of the two core-enrichment gene sets
core_genes_entrez <- union(core_DMSO, core_TF)

cat("DMSO core genes:", length(core_DMSO), "\n")
cat("TF core genes:", length(core_TF), "\n")
cat("Union core genes:", length(core_genes_entrez), "\n")

# Genes shared between the two comparisons
shared_genes <- intersect(core_DMSO, core_TF)

cat("Shared core genes:", length(shared_genes), "\n")

# -----9. Map ENTREZ -> ENSEMBL and SYMBOL--------
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = core_genes_entrez,
  columns = c("ENTREZID", "ENSEMBL", "SYMBOL"),
  keytype = "ENTREZID"
)

gene_map <- gene_map %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(ENTREZID, ENSEMBL, .keep_all = TRUE)

gene_map <- gene_map %>%
  filter(!grepl("^LOC[0-9]+", SYMBOL))


# Remove version if present in expression matrix rownames? 
# In your pipeline rownames seem like raw Ensembl IDs already,
# but we make a safe cleaned version:
expr_rownames_clean <- sub("\\..*", "", rownames(expr_mat))

# Match genes in expression matrix
gene_map$ENSEMBL_clean <- sub("\\..*", "", gene_map$ENSEMBL)

matched_idx <- match(gene_map$ENSEMBL_clean, expr_rownames_clean)
gene_map$expr_row <- matched_idx

gene_map <- gene_map %>%
  filter(!is.na(expr_row))

cat("Number of pathway genes found in expression matrix:", nrow(gene_map), "\n")

if (nrow(gene_map) == 0) {
  stop("No pathway genes matched the expression matrix.")
}

# --------10. Subset expression matrix-------------
heatmap_mat <- expr_mat[gene_map$expr_row, , drop = FALSE]

# Set rownames to gene symbols where possible
row_labels <- gene_map$SYMBOL
row_labels[is.na(row_labels) | row_labels == ""] <- gene_map$ENSEMBL[is.na(row_labels) | row_labels == ""]
rownames(heatmap_mat) <- make.unique(row_labels)

# ------11. Optional sample ordering-----
sample_order_names <- meta_data %>%
  arrange(treatment, group) %>%
  pull(`sample name`)

# sample_order_names <- meta_data %>%
#   arrange(group, treatment) %>%
#   pull(`sample name`)


sample_order_names <- sample_order_names[sample_order_names %in% colnames(heatmap_mat)]
heatmap_mat <- heatmap_mat[, sample_order_names, drop = FALSE]

meta_data_ordered <- meta_data[match(colnames(heatmap_mat), meta_data$`sample name`), ]
stopifnot(identical(colnames(heatmap_mat), meta_data_ordered$`sample name`))


# -------12. Scale by gene (row z-score)-------
# heatmap_scaled <- t(scale(t(heatmap_mat)))

#do if want to switch row and column                    
heatmap_scaled <- t(scale(t(heatmap_mat)))
heatmap_scaled <- t(heatmap_scaled)

# Remove rows with all NA after scaling
heatmap_scaled <- heatmap_scaled[complete.cases(heatmap_scaled), , drop = FALSE]

# -------13. Sample annotation---------
annotation_col <- data.frame(
  Treatment = meta_data_ordered$treatment,
  Group = meta_data_ordered$group
)
rownames(annotation_col) <- meta_data_ordered$`sample name`

ann_colors <- list(
  Group = c(control = "#4C78A8", patient = "#E45756"),
  Treatment = c(DMSO = "#72B7B2", TF = "#F58518"),
  Individual = setNames(RColorBrewer::brewer.pal(length(unique(meta_data_ordered$ind.n)), "Set2"),
                        levels(meta_data_ordered$ind.n))
)

# -------14. Draw and save heatmap------------
png(
  filename = paste0(output_prefix, ".png"),
  width = 2200,
  height = 3200,
  res = 300
)

pheatmap(
  heatmap_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = paste0(selected_pathway, "\n")
)

dev.off()

cat("Heatmap saved as:", paste0(output_prefix, ".png"), "\n")

#Or
png(
  filename = paste0(output_prefix, ".png"),
  width = 5000,
  height = 1500,
  res = 300
)

pheatmap(
  heatmap_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 315,
  main = paste(selected_pathway, "", sep = "\n")
)

dev.off()

cat("Heatmap saved as:", paste0(output_prefix, ".png"), "\n")

