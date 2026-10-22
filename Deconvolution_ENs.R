# Initialization
rm(list = ls())
graphics.off()
Sys.setenv(lang = "EN")

#Deconvolution Analysis based on MuSiC Pipeline (REF: https://xuranw.github.io/MuSiC/articles/MuSiC.html)
Sys.setenv(lang = "EN")
library(MuSiC)
library(Biobase)
library(SingleCellExperiment)
library(biomaRt)
library(readxl)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(ggpubr)


setwd(".Yanlin/Deconvolution")
bulk_eset <- readRDS("bulkRNAseq.rds")
bulk_eset

# Check First Values of gene_idm Column 
head(exprs(bulk_eset)[, 1])
# Extract gene_id correctly
gene_ids <- as.character(exprs(bulk_eset)[, 1])  # Make sure it is set as 'character'

# Remove gene_id column out of the Expression Matrix
bulk.mtx <- exprs(bulk_eset)

# Set gene_ids as rownames
rownames(bulk.mtx) <- gene_ids

# Check rownames once again
head(rownames(bulk.mtx))  # It should now show gene identifiers and not numbers

# Remove lines with NAs or empty cells
valid_gene_ids <- !is.na(gene_ids) & gene_ids != ""
bulk.mtx <- bulk.mtx[valid_gene_ids, ]
gene_ids <- gene_ids[valid_gene_ids]

bulk.mtx <- as.matrix(exprs(bulk_eset)[, -1])
mode(bulk.mtx) <- "numeric"
rownames(bulk.mtx) <- gene_ids


# Also update gene_ids to match the filtered matrix
gene_ids <- rownames(bulk.mtx)

# Set gene_ids as rownames
rownames(bulk.mtx) <- gene_ids

# Check rownames once again
head(rownames(bulk.mtx))

EMTAB.sce <- readRDS('sce_ref.rds')
EMTAB.sce
EMTAB.sce$annotation <- as.factor(EMTAB.sce$annotation)
levels(EMTAB.sce$annotation)

# Produce the first step information
ENSsub.basis = music_basis(EMTAB.sce, clusters = 'annotation', samples = 'orig.ident', 
                           select.ct = c("PEMN_1", "PEMN_2", "PEMN_3", "PEMN_4", "PIMN_1", "PIMN_2",  
                                         "PIMN_3", "PIMN_4", "PIMN_5", "PIN_1", "PIN_2", "PSVN_1", "PSVN_2", 'PSN'))
# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(ENSsub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(ENSsub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

clusters.type = list(C1 = c('PEMN_1', 'PEMN_2', 'PEMN_3','PEMN_4'), C2 = c('PIMN_1', 'PIMN_2', 'PIMN_3', 'PIMN_4', 'PIMN_5'),
                     C3 = c('PIN_1', 'PIN_2'),  C4 = c('PSVN_1', 'PSVN_2'), C5 = 'PSN')

cl.type = as.character(EMTAB.sce$annotation)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
EMTAB.sce$clusterType = factor(cl.type, levels = c(names(clusters.type)))

# 14 selected cell types
s.type = unlist(clusters.type)
s.type
#C11      C12      C13      C14      C21      C22      C23      C24      C25      C31      C32      C41      C42       C5 
#"PEMN_1" "PEMN_2" "PEMN_3" "PEMN_4" "PIMN_1" "PIMN_2" "PIMN_3" "PIMN_4" "PIMN_5"  "PIN_1"  "PIN_2" "PSVN_1" "PSVN_2"    "PSN"  

marker_table <- read_excel("ENSMarkers.xlsx")
group.markers <- split(marker_table$gene, marker_table$ident)
str(group.markers)
ENSmarkers <- group.markers

#this was the one that worked

clusters.type = list(
  C1 = c('PEMN_1', 'PEMN_2', 'PEMN_3', 'PEMN_4'),
  C2 = c('PIMN_1', 'PIMN_2', 'PIMN_3', 'PIMN_4', 'PIMN_5'),
  C3 = c('PIN_1', 'PIN_2'),
  C4 = c('PSVN_1', 'PSVN_2'),
  C5 = c('PSN')
)
group.markers <- split(marker_table$gene, marker_table$ident)

ENSmarkers <- lapply(clusters.type, function(cluster_cells) {
  # Combine marker genes from all subtypes within this cluster
  unique(unlist(group.markers[cluster_cells]))
})

str(ENSmarkers)
# Should show List of 5, names: C1, C2, C3, C4, C5

#If the code gives out an error saying it is not numeric: run these two lines and try again
bulk.mtx <- apply(bulk.mtx, 2, as.numeric)
rownames(bulk.mtx) <- gene_ids

# ConnecbiomaRt# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"   # ou "asia", "www"
)

# Retrieve mapping
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(bulk.mtx),
  mart = ensembl
)

# Remove unmapped or empty gene symbols
mapping <- mapping[mapping$hgnc_symbol != "", ]
mapping <- mapping[!duplicated(mapping$ensembl_gene_id), ]

# Subset bulk.mtx to those with mapped symbols
bulk.mtx <- bulk.mtx[rownames(bulk.mtx) %in% mapping$ensembl_gene_id, ]

# Replace Ensembl IDs with HGNC symbols
rownames(bulk.mtx) <- mapping$hgnc_symbol[match(rownames(bulk.mtx), mapping$ensembl_gene_id)]
length(intersect(rownames(bulk.mtx), rownames(EMTAB.sce)))

rownames(counts) <- sub("\\..*$", "", counts$gene_id)
rownames(EMTAB.sce) <- sub("\\..*$", "", rownames(EMTAB.sce))
ENSmarkers <- lapply(ENSmarkers, function(x) sub("\\..*$", "", x))

length(intersect(rownames(bulk.mtx), rownames(EMTAB.sce)))

head(rownames(bulk.mtx))
head(rownames(EMTAB.sce))


Est.entericneurons.bulk = music_prop.cluster(
  bulk.mtx = bulk.mtx, 
  sc.sce = EMTAB.sce, 
  clusters = 'annotation', 
  group.markers = ENSmarkers,
  groups = 'clusterType', 
  samples = 'orig.ident', 
  clusters.type = clusters.type
)

#Save the Data

write.csv(Est.entericneurons.bulk, 'ENSEstimation_sPD.csv')

names(Est.entericneurons.bulk)


################### Extract proportions To Perform a Stacked Plot
prop_df <- as.data.frame(Est.entericneurons.bulk$Est.prop.weighted.cluster)
prop_df$Sample <- rownames(prop_df)

# Convert to long format for ggplot
prop_long <- melt(prop_df, id.vars = "Sample", variable.name = "Cell_Type", value.name = "Proportion")

# Plot
ggplot(prop_long, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Estimated Cell Type Proportions per Sample", x = "Sample", y = "Proportion") +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),  # Increases legend text size and makes it bold
    legend.title = element_text(size = 16, face = "bold")  # Increases legend title size and makes it bold
  )

# First: Ensure Condition column is correctly assigned
# Adjust based on your actual sample names
prop_long$Condition <- ifelse(grepl("Ctl", prop_long$Sample), "Ctl", "Pat")
prop_long$Treatment <- ifelse(grepl("DMSO", prop_long$Sample), "DMSO", "TF")

# Group variable for convenience
prop_long$Group <- paste(prop_long$Condition, prop_long$Treatment, sep = "_")


compare_groups <- function(data, condition1, condition2,
                           treatment_filter = NULL,
                           group_col = "Condition",
                           paired = FALSE,
                           test_method = "wilcox.test") {
  
  df <- data
  
  if (!is.null(treatment_filter)) {
    df <- df[df$Treatment == treatment_filter, ]
  }
  
  df <- df[df[[group_col]] %in% c(condition1, condition2), ]
  df[[group_col]] <- factor(df[[group_col]], levels = c(condition1, condition2))
  
  stat_results <- df %>%
    group_by(Cell_Type) %>%
    summarise(
      p_value = tryCatch({
        if (test_method == "t.test") {
          t.test(
            Proportion ~ .data[[group_col]],
            data = pick(everything()),
            paired = FALSE
          )$p.value
        } else {
          wilcox.test(
            Proportion ~ .data[[group_col]],
            data = pick(everything()),
            paired = FALSE
          )$p.value
        }
      }, error = function(e) NA),
      .groups = "drop"
    )
  
  df %>%
    left_join(stat_results, by = "Cell_Type")
}

ctl_vs_pat_DMSO <- compare_groups(
  prop_long,
  "Ctl", "Pat",
  treatment_filter = "DMSO",
  paired = FALSE
)

ctl_vs_pat_TF <- compare_groups(
  prop_long,
  "Ctl", "Pat",
  treatment_filter = "TF",
  paired = FALSE
)

pat_TF_vs_DMSO <- compare_groups(
  prop_long[prop_long$Condition == "Pat", ],
  "TF", "DMSO",
  group_col = "Treatment",
  paired = FALSE   # ✅ corrigido
)

ctl_TF_vs_DMSO <- compare_groups(
  prop_long[prop_long$Condition == "Ctl", ],
  "TF", "DMSO",
  group_col = "Treatment",
  paired = FALSE   # ✅ corrigido
)

ggplot(ctl_vs_pat_DMSO, aes(x = Cell_Type, y = Proportion, fill = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Control vs sPD in DMSO", x = "Cell Type", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")) +   # <-- this adds back axis lines
  stat_compare_means(
    method = "wilcox.test",   # teste não paramétrico
    label = "p.format"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Ctl" = "steelblue", "Pat" = "firebrick"))

ggplot(ctl_vs_pat_TF , aes(x = Cell_Type, y = Proportion, fill = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Control vs sPD in TF", x = "Cell Type", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")) +   # <-- this adds back axis lines
  stat_compare_means(
    method = "wilcox.test",   # teste não paramétrico
    label = "p.format"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Ctl" = "steelblue", "Pat" = "firebrick"))

ggplot(ctl_TF_vs_DMSO, aes(x = Cell_Type, y = Proportion, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = c("DMSO" = "steelblue", "TF" = "firebrick")) +
  labs(title = "Control: TF vs DMSO", x = "Cell Type", y = "Proportion") +
  theme_classic() +
  stat_compare_means(
    method = "wilcox.test",   # teste não paramétrico
    label = "p.format"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(pat_TF_vs_DMSO, aes(x = Cell_Type, y = Proportion, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = c("DMSO" = "steelblue", "TF" = "firebrick")) +
  labs(title = "sPD: TF vs DMSO", x = "Cell Type", y = "Proportion") +
  theme_classic() +
  stat_compare_means(
    method = "wilcox.test",   # teste não paramétrico
    label = "p.format"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




