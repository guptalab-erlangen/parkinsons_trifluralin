#Deconvolution Analysis based on MuSiC Pipeline (REF: https://xuranw.github.io/MuSiC/articles/MuSiC.html)
Sys.setenv(lang = "EN")


library(MuSiC)
library(Biobase)
library(SingleCellExperiment)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)


setwd("./Yanlin/Deconvolution")
bulk_eset <- readRDS("bulkRNAseq.rds")
bulk_eset

# Check First Values of gene_idm Column 
head(exprs(bulk_eset)[, 1])

# Remove gene_id column out of the Expression Matrix
bulk.mtx <- exprs(bulk_eset)

# Check rownames once again
head(rownames(bulk.mtx))  # It should now show gene identifiers and not numbers

valid_gene_ids <- !is.na(rownames(bulk.mtx)) & rownames(bulk.mtx) != ""

bulk.mtx <- bulk.mtx[valid_gene_ids, , drop = FALSE]
head(rownames(bulk.mtx))
dim(bulk.mtx)

EMTAB.sce <- readRDS('sce_ref_glia.rds')
EMTAB.sce
EMTAB.sce$annotation <- as.factor(EMTAB.sce$annotation)
levels(EMTAB.sce$annotation)

# Produce the first step information
ENSsub.basis = music_basis(EMTAB.sce, clusters = 'annotation', samples = 'orig.ident', 
                           select.ct = c('Glia_1', 'Glia_2', 'Glia_3', 'Glia_4', 'Glia_5'))
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

clusters.type = list(C1 = 'Glia_5', C2 = 'Glia_4',
                     C3 = 'Glia_3', C4 = c('Glia_1', 'Glia_2'))

cl.type = as.character(EMTAB.sce$annotation)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
EMTAB.sce$clusterType = factor(cl.type, levels = c(names(clusters.type)))

# 14 selected cell types
s.mouse = unlist(clusters.type)
s.mouse
#C1       C2       C3      C41      C42 
#"Glia_5" "Glia_4" "Glia_3" "Glia_1" "Glia_2"

marker_table <- readxl::read_xlsx("ENSGlialMarkers.xlsx")
group.markers <- split(marker_table$gene, marker_table$ident)
str(group.markers)
ENSGlialmarkers <- group.markers

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(bulk.mtx),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)
head(mapping)

# Remove unmapped or empty gene symbols
mapping <- mapping[
  !is.na(mapping$SYMBOL) & mapping$SYMBOL != "",
]
mapping <- mapping[
  !duplicated(mapping$ENSEMBL),
]

# Subset bulk.mtx to those with mapped symbols
bulk.mtx <- bulk.mtx[
  rownames(bulk.mtx) %in% mapping$ENSEMBL,
  ,
  drop = FALSE
]
# Replace Ensembl IDs with HGNC symbols
rownames(bulk.mtx) <- mapping$SYMBOL[
  match(rownames(bulk.mtx), mapping$ENSEMBL)
]
bulk.mtx <- bulk.mtx[
  !duplicated(rownames(bulk.mtx)),
  ,
  drop = FALSE
]


common_genes <- intersect(
  rownames(bulk.mtx),
  rownames(EMTAB.sce)
)


names(ENSGlialmarkers) = c('C4', 'C4', 'C3', 'C2', 'C1')
Est.Glial.bulk = music_prop.cluster(bulk.mtx = bulk.mtx, 
                                    sc.sce = EMTAB.sce, clusters = 'annotation', group.markers = ENSGlialmarkers,
                                    groups = 'clusterType', samples = 'orig.ident', clusters.type = clusters.type)


names(Est.Glial.bulk)
write.csv(Est.Glial.bulk, 'GliaEstimation.csv')

# Extract proportions
prop_df <- as.data.frame(Est.Glial.bulk$Est.prop.weighted.cluster)
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
prop_long$Condition <- ifelse(grepl("^A", prop_long$Sample), "Adult",
                              ifelse(grepl("^Neo", prop_long$Sample), "Neonate", NA))

# Remove rows with NA in Condition
prop_long <- prop_long[!is.na(prop_long$Condition), ]

# Now define the comparison function
compare_groups <- function(data, condition1, condition2, group_col = "Condition",
                           paired = FALSE, test_method = "t.test") {
  
  df <- data[data[[group_col]] %in% c(condition1, condition2), ]  # <- this is the key line
  df[[group_col]] <- factor(df[[group_col]], levels = c(condition1, condition2))
  
  stat_results <- df %>%
    group_by(Cell_Type) %>%
    summarise(
      p_value = tryCatch({
        if (test_method == "t.test") {
          t.test(Proportion ~ get(group_col), data = ., paired = paired)$p.value
        } else {
          wilcox.test(Proportion ~ get(group_col), data = ., paired = paired)$p.value
        }
      }, error = function(e) NA)
    )
  
  df <- merge(df, stat_results, by = "Cell_Type")
  return(df)
}

# ✅ Now run the comparison
adult_vs_neonate <- compare_groups(prop_long, "Adult", "Neonate", paired = FALSE)

ggplot(adult_vs_neonate, aes(x = Cell_Type, y = Proportion, fill = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PD patients vs Healthy controls - enteric glia", x = "Cell Type", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")) +   # <-- this adds back axis lines
  stat_compare_means(aes(group = Condition), method = "t.test", label = "p.signif") +
  scale_fill_manual(values = c("Adult" = "steelblue", "Neonate" = "firebrick"))



prop_df <- as.data.frame(Est.Glial.bulk$Est.prop.weighted.cluster)
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

