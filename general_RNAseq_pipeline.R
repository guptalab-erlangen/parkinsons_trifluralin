
# SET PATH
setwd("./RNAseq/Yanlin")

#BiocManager::install("DESeq2")
# Load libraries
library(tidyverse)
library(edgeR)
library(DESeq2)
library(readxl)
library(limma)
library(biomaRt)
library(ggrepel)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(openxlsx)
library(enrichplot)


# Read count data
count_data <- read_excel ("gene_count.xlsx")
meta_data <- read_xlsx("RNASeq_sPD_sample_info.xlsx")
#meta_data$sample_id <- gsub(pattern = '_Isogenic| |_|_Mutant',replacement = '',\paste(meta_data$`cell_line`, meta_data$condition,sep='_'))

meta_data$individual_ID <- as.factor(meta_data$individual_ID)


# Separate gene IDs from count data
gene_ids <- count_data$gene_id  # Store gene IDs (first column)
count_matrix <- count_data[, 2:15]
count_matrix <- as.data.frame(count_matrix)
rownames(count_matrix) <- count_data[[1]] 

# Check that the dimensions now match
dim(count_matrix) 
dim(meta_data) 

count_matrix <- count_matrix[,match(meta_data$'sample name',colnames(count_matrix))]
identical(colnames(count_matrix),meta_data$'sample name')

meta_data$treatment <- ifelse(grepl("DMSO", meta_data[[5]]), "DMSO", "TF")

meta_data$group <- factor(meta_data$donor, levels = c("control", "patient"))
meta_data$treatment <- factor(meta_data$treatment, levels = c("DMSO", "TF"))

#Deal with nested data
meta_data$ind.n <- NA
patient_idx <- meta_data$group == "patient"
meta_data$ind.n[patient_idx] <- rep(1:4, each = 2)
control_idx <- meta_data$group == "control"
meta_data$ind.n[control_idx] <- rep(1:3, each = 2)
meta_data$ind.n <- factor(meta_data$ind.n)

#Create the combined label
meta_data$label <- paste(meta_data$group, meta_data$treatment, sep = "_")
meta_data$label <- factor(meta_data$label)

write.csv(meta_data,
          file = "sample_metadata.csv",
          row.names = FALSE) 


# Create DESeq2 object, can only infer "patient-biased treatment genes"
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta_data,
                              design = ~ ind.n + label)

rownames(dds) <- gene_ids

hist(res$pvalue, breaks=50, col="grey")

# PCA plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="group")
plotPCA(vsd, intgroup="label")
pcaData <-plotPCA(vsd, intgroup="treatment")
plotPCA(vsd, intgroup="individual_ID")


# Extract the batch information from colData
batch <- colData(dds)$ind.n
design <- model.matrix(~ colData(vsd)$label)

# Remove the batch effect using limma's removeBatchEffect function
assay(vsd) <- removeBatchEffect(assay(vsd), batch = batch, design = design)

# Plot PCA on batch-corrected data
plotPCA(vsd, intgroup = "treatment")
plotPCA(vsd, intgroup="group")
plotPCA(vsd, intgroup="individual_ID")
plotPCA(vsd, intgroup="ind.n")
plotPCA(vsd, intgroup="label")

# Filter out low counts and run DESeq2
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)

# Get DESeq2 results
res <- results(dds)

res_df <- as.data.frame(res)

resultsNames(dds)


#-----1. Patient vs Control in DMSO ------
#1. Patient vs Control in DMSO
res_pat_vs_ctl_DMSO <- results(dds, name = "label_patient_DMSO_vs_control_DMSO")
res_pat_vs_ctl_DMSO_df <- as.data.frame(res_pat_vs_ctl_DMSO)
sum(res_pat_vs_ctl_DMSO_df$padj<0.05, na.rm = TRUE)

res_pat_vs_ctl_DMSO_df$significant <- ifelse(res_pat_vs_ctl_DMSO_df$padj < 0.05, "significant", "not significant")

res_pat_vs_ctl_DMSO_df$regulation <- ifelse(res_pat_vs_ctl_DMSO_df$padj < 0.05 & res_pat_vs_ctl_DMSO_df$log2FoldChange > 1, "upregulated",
                                         ifelse(res_pat_vs_ctl_DMSO_df$padj < 0.05 & res_pat_vs_ctl_DMSO_df$log2FoldChange < -1, "downregulated", "not significant"))

res_pat_vs_ctl_DMSO_df$ENSEMBL_ID <- rownames(res_pat_vs_ctl_DMSO_df)

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes1 <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
               filters = 'ensembl_gene_id',
               values = res_pat_vs_ctl_DMSO_df$ENSEMBL_ID,
               mart = mart)

res_pat_vs_ctl_DMSO_df <- merge(res_pat_vs_ctl_DMSO_df, genes1, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(res_pat_vs_ctl_DMSO_df)[which(colnames(res_pat_vs_ctl_DMSO_df) == "external_gene_name")] <- "Gene_ID"

write.csv(res_pat_vs_ctl_DMSO_df,
          file = "DEGs_DMSO_pat_vs_ctrl.csv",
          row.names = TRUE)

## Creating Volcano Plot
res_pat_vs_ctl_DMSO_df$color <- ifelse(res_pat_vs_ctl_DMSO_df$padj < 0.05 & res_pat_vs_ctl_DMSO_df$log2FoldChange > 0.585, "Upregulated",
                                    ifelse(res_pat_vs_ctl_DMSO_df$padj < 0.05 & res_pat_vs_ctl_DMSO_df$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

print(head(subset(res_pat_vs_ctl_DMSO_df, significant == "significant")))

# Rebuild the volcano plot
label_df <- res_pat_vs_ctl_DMSO_df %>%
  filter(color %in% c("Upregulated", "Downregulated"),
         !is.na(Gene_ID), Gene_ID != "") %>%
  group_by(Gene_ID) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# Base volcano plot
volcano_plot <- ggplot(res_pat_vs_ctl_DMSO_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color), alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) + 
  labs(title = "DMSO: patient VS control",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") + 
  theme_bw() +  # white background
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_df,
    aes(label = Gene_ID),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  )

# Save the plot
ggsave(
  filename = "Volcano_Patient_vs_Control_DMSO.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300
)
# Connect to Ensembl and retrieve the correct dataset
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl IDs to ENTREZ IDs
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = res_pat_vs_ctl_DMSO_df$ENSEMBL_ID,
                         mart = mart)

# Merge conversion results with your DEGs data
res_pat_vs_ctl_DMSO_df <- merge(res_pat_vs_ctl_DMSO_df, gene_conversion, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)

# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # For human genes



#-----2. Patient vs Control in TF------
#2. Patient vs Control in TF
res_pat_vs_ctl_TF <- results(dds, contrast = list( "label_patient_TF_vs_control_DMSO", "label_control_TF_vs_control_DMSO"))
res_pat_vs_ctl_TF_df <- as.data.frame(res_pat_vs_ctl_TF)
sum(res_pat_vs_ctl_TF_df$padj<0.05, na.rm = TRUE)

res_pat_vs_ctl_TF_df$significant <- ifelse(res_pat_vs_ctl_TF_df$padj < 0.05, "significant", "not significant")

res_pat_vs_ctl_TF_df$regulation <- ifelse(res_pat_vs_ctl_TF_df$padj < 0.05 & res_pat_vs_ctl_TF_df$log2FoldChange > 1, "upregulated",
                                            ifelse(res_pat_vs_ctl_TF_df$padj < 0.05 & res_pat_vs_ctl_TF_df$log2FoldChange < -1, "downregulated", "not significant"))

res_pat_vs_ctl_TF_df$ENSEMBL_ID <- rownames(res_pat_vs_ctl_TF_df)

genes2 <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                filters = 'ensembl_gene_id',
                values = res_pat_vs_ctl_TF_df$ENSEMBL_ID,
                mart = mart)

res_pat_vs_ctl_TF_df <- merge(res_pat_vs_ctl_TF_df, genes2, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(res_pat_vs_ctl_TF_df)[which(colnames(res_pat_vs_ctl_TF_df) == "external_gene_name")] <- "Gene_ID"

write.csv(res_pat_vs_ctl_TF_df,
          file = "DEGs_TF_pat_vs_ctrl.csv",
          row.names = TRUE)

## Creating Volcano Plot
res_pat_vs_ctl_TF_df$color <- ifelse(res_pat_vs_ctl_TF_df$padj < 0.05 & res_pat_vs_ctl_TF_df$log2FoldChange > 0.585, "Upregulated",
                                       ifelse(res_pat_vs_ctl_TF_df$padj < 0.05 & res_pat_vs_ctl_TF_df$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

print(head(subset(res_pat_vs_ctl_TF_df, significant == "significant")))

# Base volcano plot
label_df <- res_pat_vs_ctl_TF_df %>%
  filter(color %in% c("Upregulated", "Downregulated"),
         !is.na(Gene_ID), Gene_ID != "") %>%
  group_by(Gene_ID) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup()

volcano_plot <- ggplot(res_pat_vs_ctl_TF_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color), alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) + 
  labs(title = "TF: patient VS control",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") + 
  theme_bw() +  # white background
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_df,
    aes(label = Gene_ID),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  )

# Save the plot
ggsave(
  filename = "Volcano_Patient_vs_Control_TF.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Connect to Ensembl and retrieve the correct dataset
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl IDs to ENTREZ IDs
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = res_pat_vs_ctl_TF_df$ENSEMBL_ID,
                         mart = mart)

# Merge conversion results with your DEGs data
res_pat_vs_ctl_TF_df <- merge(res_pat_vs_ctl_TF_df, gene_conversion, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)


#-----3. Patient TF vs DMSO------
#3. Patient TF vs DMSO
res_patient_TF_vs_DMSO <- results(dds, contrast = list( "label_patient_TF_vs_control_DMSO","label_patient_DMSO_vs_control_DMSO"))
res_patient_TF_vs_DMSO_df <- as.data.frame(res_patient_TF_vs_DMSO)
sum(res_patient_TF_vs_DMSO_df$padj<0.05, na.rm = TRUE)

res_patient_TF_vs_DMSO_df$significant <- ifelse(res_patient_TF_vs_DMSO_df$padj < 0.05, "significant", "not significant")

res_patient_TF_vs_DMSO_df$regulation <- ifelse(res_patient_TF_vs_DMSO_df$padj < 0.05 & res_patient_TF_vs_DMSO_df$log2FoldChange > 1, "upregulated",
                                          ifelse(res_patient_TF_vs_DMSO_df$padj < 0.05 & res_patient_TF_vs_DMSO_df$log2FoldChange < -1, "downregulated", "not significant"))

res_patient_TF_vs_DMSO_df$ENSEMBL_ID <- rownames(res_patient_TF_vs_DMSO_df)

genes3 <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                filters = 'ensembl_gene_id',
                values = res_patient_TF_vs_DMSO_df$ENSEMBL_ID,
                mart = mart)

res_patient_TF_vs_DMSO_df <- merge(res_patient_TF_vs_DMSO_df, genes3, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(res_patient_TF_vs_DMSO_df)[which(colnames(res_patient_TF_vs_DMSO_df) == "external_gene_name")] <- "Gene_ID"

write.csv(res_patient_TF_vs_DMSO_df,
          file = "DEGs_pat_TF_vs_DMSO.csv",
          row.names = TRUE)

## Creating Volcano Plot
res_patient_TF_vs_DMSO_df$color <- ifelse(res_patient_TF_vs_DMSO_df$padj < 0.05 & res_patient_TF_vs_DMSO_df$log2FoldChange > 0.585, "Upregulated",
                                       ifelse(res_patient_TF_vs_DMSO_df$padj < 0.05 & res_patient_TF_vs_DMSO_df$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

print(head(subset(res_patient_TF_vs_DMSO_df, significant == "significant")))

# Base volcano plot
label_df <- res_patient_TF_vs_DMSO_df %>%
  filter(color %in% c("Upregulated", "Downregulated"),
         !is.na(Gene_ID), Gene_ID != "") %>%
  group_by(Gene_ID) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup()

plot_df <- res_patient_TF_vs_DMSO_df %>%
  filter(
    !is.na(padj),
    !is.na(log2FoldChange)
  )

volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color), alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) + 
  labs(title = "Patient: TF vs DMSO",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") + 
  theme_bw() +  # white background
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_df,
    aes(label = Gene_ID),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  )

# Save the plot
ggsave(
  filename = "Volcano_Patient_TF_vs_DMSO.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300
)


# Convert Ensembl IDs to ENTREZ IDs
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = res_patient_TF_vs_DMSO_df$ENSEMBL_ID,
                         mart = mart)

# Merge conversion results with your DEGs data
res_patient_TF_vs_DMSO_df <- merge(res_patient_TF_vs_DMSO_df, gene_conversion, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)


#-----4. Control TF vs DMSO------
#4. Control TF vs DMSO
res_control_TF_vs_DMSO <- results(dds, name = "label_control_TF_vs_control_DMSO")
res_control_TF_vs_DMSO_df <- as.data.frame(res_control_TF_vs_DMSO)
sum(res_control_TF_vs_DMSO_df$padj<0.05, na.rm = TRUE)

res_control_TF_vs_DMSO_df$significant <- ifelse(res_control_TF_vs_DMSO_df$padj < 0.05, "significant", "not significant")

res_control_TF_vs_DMSO_df$regulation <- ifelse(res_control_TF_vs_DMSO_df$padj < 0.05 & res_control_TF_vs_DMSO_df$log2FoldChange > 1, "upregulated",
                                               ifelse(res_control_TF_vs_DMSO_df$padj < 0.05 & res_control_TF_vs_DMSO_df$log2FoldChange < -1, "downregulated", "not significant"))

res_control_TF_vs_DMSO_df$ENSEMBL_ID <- rownames(res_control_TF_vs_DMSO_df)

genes4 <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                filters = 'ensembl_gene_id',
                values = res_control_TF_vs_DMSO_df$ENSEMBL_ID,
                mart = mart)

res_control_TF_vs_DMSO_df <- merge(res_control_TF_vs_DMSO_df, genes4, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(res_control_TF_vs_DMSO_df)[which(colnames(res_control_TF_vs_DMSO_df) == "external_gene_name")] <- "Gene_ID"

write.csv(res_control_TF_vs_DMSO_df,
          file = "DEGs_ctrl_TF_vs_DMSO.csv",
          row.names = TRUE)

## Creating Volcano Plot
res_control_TF_vs_DMSO_df$color <- ifelse(res_control_TF_vs_DMSO_df$padj < 0.05 & res_control_TF_vs_DMSO_df$log2FoldChange > 0.585, "Upregulated",
                                          ifelse(res_control_TF_vs_DMSO_df$padj < 0.05 & res_control_TF_vs_DMSO_df$log2FoldChange < -0.585, "Downregulated", "Not Significant"))

print(head(subset(res_control_TF_vs_DMSO_df, significant == "significant")))

# Base volcano plot
label_df <- res_control_TF_vs_DMSO_df %>%
  filter(color %in% c("Upregulated", "Downregulated"),
         !is.na(Gene_ID), Gene_ID != "") %>%
  group_by(Gene_ID) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup()

plot_df <- res_control_TF_vs_DMSO_df %>%
  filter(
    !is.na(padj),
    !is.na(log2FoldChange)
  )

volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color), alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) + 
  labs(title = "Control: DMSO VS TF",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") + 
  theme_bw() +  # white background
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_df,
    aes(label = Gene_ID),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  )

# Save the plot
ggsave(
  filename = "Volcano_Control_TF_vs_DMSO.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Convert Ensembl IDs to ENTREZ IDs
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = res_control_TF_vs_DMSO_df$ENSEMBL_ID,
                         mart = mart)

# Merge conversion results with your DEGs data
res_control_TF_vs_DMSO_df <- merge(res_control_TF_vs_DMSO_df, gene_conversion, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)


#-----5. patient-biased treatment genes with ~ ind.n + label ------
# significant in patients, not in controls
#TF vs DMSO is clearly present in patients, but absent or weak in controls
res_diff_of_diff <- results(
  dds,
  contrast = list(
    "label_patient_TF_vs_control_DMSO",
    c("label_control_TF_vs_control_DMSO",
      "label_patient_DMSO_vs_control_DMSO")
  )
)

res_diff_df <- as.data.frame(res_diff_of_diff)
sum(res_diff_df$padj<0.05)
res_diff_df$ENSEMBL_ID <- rownames(res_diff_df)
genes_diff <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                    filters = 'ensembl_gene_id',
                    values = res_diff_df$ENSEMBL_ID,
                    mart = mart)
res_diff_df <- merge(res_diff_df, genes_diff, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(res_diff_df)[colnames(res_diff_df) == "external_gene_name"] <- "Gene_ID"

write.csv(res_diff_df,
          file = "DEGs_interaction.csv",
          row.names = TRUE)

res_diff_df$color <- "Not Significant"
res_diff_df$color[res_diff_df$padj < 0.05 & res_diff_df$log2FoldChange > 0] <- "Patient-Specific Increase"
res_diff_df$color[res_diff_df$padj < 0.05 & res_diff_df$log2FoldChange < 0] <- "Control-Specific Increase"

volcano_diff <- ggplot(res_diff_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color), alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("Patient-Specific Increase" = "darkorange", 
                                "Control-Specific Increase" = "purple", 
                                "Not Significant" = "grey80")) + 
  labs(title = "Interaction Effect: (Pat TF-DMSO) vs (Ctl TF-DMSO)",
       subtitle = "Genes responding differently to treatment between groups",
       x = "Interaction Log2 Fold Change",
       y = "-Log10 Adjusted P-value") + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # Label top 15 most significant interaction genes
  geom_text_repel(
    data = subset(res_diff_df, padj < 0.05) %>% slice_min(padj, n = 15),
    aes(label = Gene_ID),
    size = 3.5,
    max.overlaps = 10
  )

ggsave("Volcano_Interaction_Effect.png", plot = volcano_diff, width = 10, height = 7)

#---- Perform GSEA enrichment analysis ------
prepare_gsea_input <- function(res_df) {
  
  # Remove NA pvalues
  res_df <- res_df[!is.na(res_df$pvalue), ]
  
  # Remove ENSEMBL version numbers if present
  res_df$ENSEMBL_ID <- sub("\\..*", "", res_df$ENSEMBL_ID)
  
  # Map to ENTREZ
  gene_map <- bitr(
    res_df$ENSEMBL_ID,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Merge mapping back
  res_mapped <- merge(res_df, gene_map,
                      by.x = "ENSEMBL_ID",
                      by.y = "ENSEMBL",
                      all = FALSE)
  
  # Rank by statistic (BEST choice for GSEA)
  res_mapped <- res_mapped %>%
    arrange(desc(stat))
  
  # Remove duplicated ENTREZ IDs (keep highest stat)
  res_mapped <- res_mapped[!duplicated(res_mapped$ENTREZID), ]
  
  # Create named numeric vector
  gene_list <- res_mapped$stat
  names(gene_list) <- res_mapped$ENTREZID
  
  # Sort decreasing
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

gsea_results <- list()

gsea_results$pat_vs_ctl_DMSO <- gseKEGG(
  geneList = prepare_gsea_input(res_pat_vs_ctl_DMSO_df),
  organism = "hsa"
)

gsea_results$pat_vs_ctl_TF <- gseKEGG(
  geneList = prepare_gsea_input(res_pat_vs_ctl_TF_df),
  organism = "hsa"
)

gsea_results$pat_TF_vs_DMSO <- gseKEGG(
  geneList = prepare_gsea_input(res_patient_TF_vs_DMSO_df),
  organism = "hsa"
)

gsea_results$ctl_TF_vs_DMSO <- gseKEGG(
  geneList = prepare_gsea_input(res_control_TF_vs_DMSO_df),
  organism = "hsa"
)

gsea_results$interaction_like <- gseKEGG(
  geneList = prepare_gsea_input(res_diff_df),
  organism = "hsa"
)

wb <- createWorkbook()
addWorksheet(wb, "pat_vs_ctl_DMSO")
writeData(wb, "pat_vs_ctl_DMSO", as.data.frame(gsea_results$pat_vs_ctl_DMSO))

addWorksheet(wb, "pat_vs_ctl_TF")
writeData(wb, "pat_vs_ctl_TF", as.data.frame(gsea_results$pat_vs_ctl_TF))

addWorksheet(wb, "pat_TF_vs_DMSO")
writeData(wb, "pat_TF_vs_DMSO", as.data.frame(gsea_results$pat_TF_vs_DMSO))

addWorksheet(wb, "ctl_TF_vs_DMSO")
writeData(wb, "ctl_TF_vs_DMSO", as.data.frame(gsea_results$ctl_TF_vs_DMSO))

addWorksheet(wb, "interaction_like")
writeData(wb, "interaction_like", as.data.frame(gsea_results$interaction_like))

saveWorkbook(wb, "GSEA_results_KEGG.xlsx", overwrite = TRUE)

gseaplot2(
  gsea_results$pat_TF_vs_DMSO,
  geneSetID = 1,
  title = "Top pathway: Patient TF vs DMSO"
)


gsea_list <- gsea_results

for (name in names(gsea_list)) {
  
  gsea_obj <- gsea_list[[name]]
  
  if (is.null(gsea_obj)) {
    message(name, " is NULL")
    next
  }
  
  if (nrow(gsea_obj@result) == 0) {
    message(name, " has no results")
    next
  }

  gsea_obj <- gsea_obj
  gsea_obj@result <- gsea_obj@result[gsea_obj@result$p.adjust < 0.05, ]
  
  # skip if filtering removes everything
  if (nrow(gsea_obj@result) == 0) {
    message(name, " has no significant pathways after filtering")
    next
  }
  
  p_dot <- dotplot(gsea_obj, showCategory = 15) +
    ggtitle(paste0(name, " - Dotplot"))
  
  print(p_dot)
  
  ggsave(
    filename = paste0(name, "_dotplot.pdf"),
    plot = p_dot,
    width = 6,
    height = 5
  )
  
  p_cnet <- cnetplot(
    gsea_obj,
    showCategory = 5,
    foldChange = gsea_obj@geneList
  ) +
    ggtitle(paste0(name, " - Cnetplot"))
  
  print(p_cnet)
  
  ggsave(
    filename = paste0(name, "_cnetplot.pdf"),
    plot = p_cnet,
    width = 7,
    height = 6
  )
  
  gsea_obj2 <- pairwise_termsim(gsea_obj)
  
  p_emap <- emapplot(gsea_obj2, showCategory = 15) +
    ggtitle(paste0(name, " - Emapplot"))
  
  print(p_emap)
  
  ggsave(
    filename = paste0(name, "_emapplot.pdf"),
    plot = p_emap,
    width = 7,
    height = 6
  )
}




