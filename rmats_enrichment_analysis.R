# clusterProfiler enrichment with shared study-wide universe
# all significant AS genes by comparison
# Universe = union of all tested AS genes across comparisons

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(tibble)
})

# -----1. Paths------
base_dir <- "./rmats"
results_dir <- file.path(base_dir, "results")

out_dir <- file.path(results_dir, "clusterProfiler_enrichment")
out_as_dir <- file.path(out_dir, "all_significant_AS_genes")
out_overlap_dir <- file.path(out_dir, "DEG_AS_overlap")
out_universe_dir <- file.path(out_dir, "universes")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_as_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_overlap_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_universe_dir, showWarnings = FALSE, recursive = TRUE)

# -----Helper functions-------
clean_genes <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    .[!is.na(.)] %>%
    .[. != ""] %>%
    unique()
}

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  return(NA)
}

convert_symbols_to_entrez <- function(gene_symbols) {
  gene_symbols <- clean_genes(gene_symbols)
  
  conv <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = org.Hs.eg.db
  )
  
  conv %>% distinct(SYMBOL, .keep_all = TRUE)
}

save_enrich_result <- function(enrich_obj, out_folder, prefix, comparison_name, name_suffix) {
  if (!is.null(enrich_obj) && nrow(as.data.frame(enrich_obj)) > 0) {
    df <- as.data.frame(enrich_obj)
    
    write.csv(
      df,
      file.path(out_folder, paste0(prefix, "_", comparison_name, "_", name_suffix, ".csv")),
      row.names = FALSE
    )
    
    pdf(
      file.path(out_folder, paste0(prefix, "_", comparison_name, "_", name_suffix, "_dotplot.pdf")),
      width = 10, height = 8
    )
    print(
      dotplot(enrich_obj, showCategory = 10) +
        ggtitle(paste(prefix, comparison_name, name_suffix, sep = " | ")) +
        theme(axis.text.y = element_text(size = 10))
    )
    dev.off()
    
    pdf(
      file.path(out_folder, paste0(prefix, "_", comparison_name, "_", name_suffix, "_barplot.pdf")),
      width = 10, height = 10
    )
    print(
      barplot(enrich_obj, showCategory = 15) +
        ggtitle(paste(prefix, comparison_name, name_suffix, sep = " | ")) +
        theme(axis.text.y = element_text(size = 13))
    )
    dev.off()
  }
}

read_gene_symbols_from_file <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- make.names(colnames(df))
  
  gene_symbol_col <- find_col(df, c("GeneSymbol", "geneSymbol", "gene.symbol", "Symbol", "symbol", "GeneSymbols"))
  gene_id_col <- find_col(df, c("Gene", "gene", "GeneID", "geneID"))
  
  if (!is.na(gene_symbol_col)) {
    genes <- df[[gene_symbol_col]]
  } else if (!is.na(gene_id_col)) {
    genes <- df[[gene_id_col]]
  } else {
    stop("Could not find gene symbol column in file: ", file)
  }
  
  clean_genes(genes)
}

run_cp_enrichment <- function(genes, universe_genes, comparison_name, out_folder, prefix) {
  genes <- clean_genes(genes)
  universe_genes <- clean_genes(universe_genes)
  
  if (length(genes) < 2) {
    message("Skipping ", comparison_name, " (", prefix, "): fewer than 2 foreground genes.")
    return(NULL)
  }
  
  if (length(universe_genes) < 2) {
    message("Skipping ", comparison_name, " (", prefix, "): fewer than 2 universe genes.")
    return(NULL)
  }
  
  conv_fg <- tryCatch({
    convert_symbols_to_entrez(genes)
  }, error = function(e) {
    message("Foreground ID conversion failed for ", comparison_name, " (", prefix, "): ", e$message)
    return(NULL)
  })
  
  conv_bg <- tryCatch({
    convert_symbols_to_entrez(universe_genes)
  }, error = function(e) {
    message("Universe ID conversion failed for ", comparison_name, " (", prefix, "): ", e$message)
    return(NULL)
  })
  
  if (is.null(conv_fg) || nrow(conv_fg) == 0) {
    message("No Entrez IDs mapped for foreground: ", comparison_name, " (", prefix, ").")
    return(NULL)
  }
  
  if (is.null(conv_bg) || nrow(conv_bg) == 0) {
    message("No Entrez IDs mapped for universe: ", comparison_name, " (", prefix, ").")
    return(NULL)
  }
  
  entrez_fg <- unique(conv_fg$ENTREZID)
  entrez_bg <- unique(conv_bg$ENTREZID)
  
  entrez_fg <- intersect(entrez_fg, entrez_bg)
  
  if (length(entrez_fg) < 2) {
    message("Skipping ", comparison_name, " (", prefix, "): fewer than 2 foreground Entrez IDs after intersecting with universe.")
    return(NULL)
  }
  
  if (length(entrez_bg) < 2) {
    message("Skipping ", comparison_name, " (", prefix, "): fewer than 2 universe Entrez IDs.")
    return(NULL)
  }
  
  write.csv(
    conv_fg,
    file.path(out_folder, paste0(prefix, "_", comparison_name, "_foreground_geneID_mapping.csv")),
    row.names = FALSE
  )
  
  write.csv(
    conv_bg,
    file.path(out_folder, paste0(prefix, "_", comparison_name, "_universe_geneID_mapping.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_fg),
    file.path(out_folder, paste0(prefix, "_", comparison_name, "_foreground_entrez_ids_used.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_bg),
    file.path(out_folder, paste0(prefix, "_", comparison_name, "_universe_entrez_ids_used.csv")),
    row.names = FALSE
  )
  
  message(
    "Running enrichment for ", comparison_name, " (", prefix, ") with ",
    length(entrez_fg), " foreground genes and ",
    length(entrez_bg), " universe genes."
  )
  
  ego_bp <- tryCatch({
    enrichGO(
      gene = entrez_fg,
      universe = entrez_bg,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }, error = function(e) {
    message("GO BP failed for ", comparison_name, " (", prefix, "): ", e$message)
    NULL
  })
  
  ego_cc <- tryCatch({
    enrichGO(
      gene = entrez_fg,
      universe = entrez_bg,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }, error = function(e) {
    message("GO CC failed for ", comparison_name, " (", prefix, "): ", e$message)
    NULL
  })
  
  ego_mf <- tryCatch({
    enrichGO(
      gene = entrez_fg,
      universe = entrez_bg,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }, error = function(e) {
    message("GO MF failed for ", comparison_name, " (", prefix, "): ", e$message)
    NULL
  })
  
  ekegg <- tryCatch({
    enrichKEGG(
      gene = entrez_fg,
      universe = entrez_bg,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    message("KEGG failed for ", comparison_name, " (", prefix, "): ", e$message)
    NULL
  })
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  save_enrich_result(ego_bp, out_folder, prefix, comparison_name, "GO_BP")
  save_enrich_result(ego_cc, out_folder, prefix, comparison_name, "GO_CC")
  save_enrich_result(ego_mf, out_folder, prefix, comparison_name, "GO_MF")
  save_enrich_result(ekegg, out_folder, prefix, comparison_name, "KEGG")
  
  return(list(
    foreground_mapping = conv_fg,
    universe_mapping = conv_bg,
    GO_BP = ego_bp,
    GO_CC = ego_cc,
    GO_MF = ego_mf,
    KEGG = ekegg
  ))
}

# -------3. Detect comparison folders---------
comparison_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
comparison_dirs <- comparison_dirs[basename(comparison_dirs) != "clusterProfiler_enrichment"]

if (length(comparison_dirs) == 0) {
  stop("No comparison subfolders found in results_splicing_deg.")
}

comparison_names <- basename(comparison_dirs)
print(comparison_names)

# --------4. Build per-comparison tested gene lists----------
universe_list <- list()
as_sig_list <- list()
overlap_list <- list()

for (comp_dir in comparison_dirs) {
  comp <- basename(comp_dir)
  
  tested_file <- list.files(
    comp_dir,
    pattern = "_rMATS_tested_AS_events\\.csv$",
    full.names = TRUE
  )
  
  sig_file <- list.files(
    comp_dir,
    pattern = "_rMATS_significant_AS_events\\.csv$",
    full.names = TRUE
  )
  
  overlap_file <- list.files(
    comp_dir,
    pattern = "_rMATS_DEG_overlap_gene_level\\.csv$",
    full.names = TRUE
  )
  
  if (length(tested_file) == 0) {
    message("No tested AS file found for ", comp)
    next
  }
  
  if (length(sig_file) == 0) {
    message("No significant AS file found for ", comp)
    next
  }
  
  universe_genes <- read_gene_symbols_from_file(tested_file[1])
  sig_genes <- read_gene_symbols_from_file(sig_file[1])
  
  universe_list[[comp]] <- universe_genes
  as_sig_list[[comp]] <- sig_genes
  
  write.csv(
    data.frame(GeneSymbol = universe_genes),
    file.path(out_universe_dir, paste0(comp, "_tested_AS_genes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(GeneSymbol = sig_genes),
    file.path(out_as_dir, paste0(comp, "_foreground_significant_AS_genes.csv")),
    row.names = FALSE
  )
  
  if (length(overlap_file) > 0) {
    overlap_df <- read.csv(overlap_file[1], stringsAsFactors = FALSE, check.names = FALSE)
    colnames(overlap_df) <- make.names(colnames(overlap_df))
    
    overlap_symbol_col <- find_col(overlap_df, c("GeneSymbol", "geneSymbol", "gene.symbol", "Symbol", "symbol"))
    
    if (!is.na(overlap_symbol_col)) {
      overlap_genes <- clean_genes(overlap_df[[overlap_symbol_col]])
      overlap_list[[comp]] <- overlap_genes
      
      write.csv(
        data.frame(GeneSymbol = overlap_genes),
        file.path(out_overlap_dir, paste0(comp, "_foreground_overlap_genes.csv")),
        row.names = FALSE
      )
    } else {
      message("No GeneSymbol column found in overlap file for ", comp)
    }
  }
}

# ----5. Build shared global universe---------
global_universe <- unique(unlist(universe_list))
global_universe <- clean_genes(global_universe)

if (length(global_universe) < 2) {
  stop("Global universe has fewer than 2 genes. Cannot run enrichment.")
}

write.csv(
  data.frame(GeneSymbol = global_universe),
  file.path(out_universe_dir, "global_universe_all_tested_AS_genes_union.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    universe_type = "global_union_of_all_tested_AS_genes",
    gene_count = length(global_universe)
  ),
  file.path(out_universe_dir, "global_universe_summary.csv"),
  row.names = FALSE
)

cat("Global universe size:", length(global_universe), "\n")

# ------6. Summary tables----------
universe_counts <- tibble(
  comparison = names(universe_list),
  tested_AS_gene_count = sapply(universe_list, length)
)

write.csv(
  universe_counts,
  file.path(out_universe_dir, "tested_AS_gene_counts_by_comparison.csv"),
  row.names = FALSE
)

as_counts <- tibble(
  comparison = names(as_sig_list),
  foreground_gene_count = sapply(as_sig_list, length),
  foreground_in_global_universe = sapply(as_sig_list, function(x) sum(clean_genes(x) %in% global_universe)),
  foreground_not_in_global_universe = sapply(as_sig_list, function(x) sum(!clean_genes(x) %in% global_universe))
)

write.csv(
  as_counts,
  file.path(out_as_dir, "all_significant_AS_gene_counts_by_comparison.csv"),
  row.names = FALSE
)

if (length(overlap_list) > 0) {
  overlap_counts <- tibble(
    comparison = names(overlap_list),
    foreground_gene_count = sapply(overlap_list, length),
    foreground_in_global_universe = sapply(overlap_list, function(x) sum(clean_genes(x) %in% global_universe)),
    foreground_not_in_global_universe = sapply(overlap_list, function(x) sum(!clean_genes(x) %in% global_universe))
  )
  
  write.csv(
    overlap_counts,
    file.path(out_overlap_dir, "overlap_gene_counts_by_comparison.csv"),
    row.names = FALSE
  )
}

# -------7. Run enrichment for all significant AS genes-----------
for (comp in names(as_sig_list)) {
  run_cp_enrichment(
    genes = as_sig_list[[comp]],
    universe_genes = global_universe,
    comparison_name = comp,
    out_folder = out_as_dir,
    prefix = "allAS"
  )
}

