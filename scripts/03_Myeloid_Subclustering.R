# =========================================================================
# 03_Myeloid_Subclustering.R
# =========================================================================
# Purpose: Extract myeloid cells from Stage 2 annotated object, re-process
#          (normalize, PCA, Harmony, UMAP), subcluster, identify markers,
#          annotate subclusters, and explore Control vs APAP changes.
# Input:  data/processed/02_Global_Annotation/seurat_annotated.rds
# Output: data/processed/03_Myeloid_Subclustering/myeloid_subclustered.rds
#         results/figures/03_Myeloid_Subclustering/*.pdf
#         results/tables/03_Myeloid_Subclustering/*.csv
#         notes/03_Myeloid_Subclustering_note.md
# =========================================================================
#
# MANUAL INTERVENTION POINTS (search for "MANUAL_REVIEW"):
#   1. Cluster resolution choice (Section C)
#   2. Myeloid subcluster annotation labels (Section F)
#
# The script will run end-to-end with a TEMPORARY default resolution,
# produce preliminary annotations, and save all outputs.  You MUST
# inspect the figures/tables and update the manual_review sections
# before trusting the annotation or proceeding to Stage 4.
# =========================================================================

# -------------------------------------------------------------------------
# 0. Setup
# -------------------------------------------------------------------------

# --- Check conda environment ---
current_env <- Sys.getenv("CONDA_DEFAULT_ENV")
if (current_env != "singlecell") {
  stop("Current conda environment is '", current_env,
       "'. Please activate 'singlecell' environment:\n",
       "  conda activate singlecell")
}
message("conda environment: ", current_env)

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(future)

plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)

# --- project root from script path ---
script_args <- commandArgs(trailingOnly = FALSE)
file_arg    <- grep("^--file=", script_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg)
if (length(script_path) == 0) {
  PROJ_ROOT <- getwd()
} else {
  PROJ_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
}
setwd(PROJ_ROOT)
message("Project root: ", PROJ_ROOT)

# --- directories ---
PROC_DIR   <- file.path(PROJ_ROOT, "data", "processed", "03_Myeloid_Subclustering")
FIG_DIR    <- file.path(PROJ_ROOT, "results", "figures",   "03_Myeloid_Subclustering")
TAB_DIR    <- file.path(PROJ_ROOT, "results", "tables",    "03_Myeloid_Subclustering")
NOTE_DIR   <- file.path(PROJ_ROOT, "notes")
STAGE2_RDS <- file.path(PROJ_ROOT, "data", "processed", "02_Global_Annotation",
                        "seurat_annotated.rds")

dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(NOTE_DIR, recursive = TRUE, showWarnings = FALSE)

NPCS <- 25

# -------------------------------------------------------------------------
# A. Extract myeloid subset from Stage-2 object
# -------------------------------------------------------------------------

if (!file.exists(STAGE2_RDS)) {
  stop("Stage-2 RDS not found: ", STAGE2_RDS)
}
message("Reading: ", STAGE2_RDS)
merged <- readRDS(STAGE2_RDS)
message("Loaded object: ", ncol(merged), " cells x ", nrow(merged), " genes")

# --- validate metadata ---
req_meta <- c("celltype_manual", "sample_id", "group")
missing_meta <- setdiff(req_meta, colnames(merged@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ",
       paste(missing_meta, collapse = ", "))
}

# --- myeloid cell types to keep ---
myeloid_keep <- c(
  "Kupffer cell",
  "Monocyte-derived macrophage / DC-like",
  "Inflammatory neutrophil",
  "Activated neutrophil",
  "Neutrophil"
)

# Check each type exists
avail_types <- unique(merged$celltype_manual)
missing_types <- setdiff(myeloid_keep, avail_types)
if (length(missing_types) > 0) {
  stop("The following myeloid types not found in celltype_manual:\n",
       paste(missing_types, collapse = ", "),
       "\nAvailable types:\n",
       paste(sort(avail_types), collapse = ", "))
}

# Select myeloid cells
myeloid_barcodes <- colnames(merged)[merged$celltype_manual %in% myeloid_keep]
message("Myeloid barcodes identified: ", length(myeloid_barcodes))

if (length(myeloid_barcodes) < 500) {
  stop("Too few myeloid cells for subclustering: ", length(myeloid_barcodes),
       " (minimum 500 required)")
}

myeloid <- subset(merged, cells = myeloid_barcodes)
rm(merged)
gc()
message("Myeloid subset: ", ncol(myeloid), " cells x ", nrow(myeloid), " genes")

# --- validate per-sample coverage ---
sample_check <- myeloid@meta.data %>%
  count(sample_id, group, celltype_manual, name = "n_cells")

# Check each sample has myeloid cells
samples_with_myeloid <- unique(sample_check$sample_id)
expected_samples <- c("con1", "con2", "con3", "AP-300", "AP-300-2", "AP-300-3")
missing_samples <- setdiff(expected_samples, samples_with_myeloid)
if (length(missing_samples) > 0) {
  warning("The following samples have NO myeloid cells: ",
          paste(missing_samples, collapse = ", "))
}

# Samples with very few myeloid cells
low_samples <- sample_check %>%
  group_by(sample_id) %>%
  summarise(total = sum(n_cells), .groups = "drop") %>%
  filter(total < 50)
if (nrow(low_samples) > 0) {
  warning("Samples with <50 myeloid cells: ",
          paste(low_samples$sample_id, "(", low_samples$total, "cells)", collapse = ", "))
}

# Output input cell counts
write_csv(sample_check, file.path(TAB_DIR, "myeloid_input_cell_counts.csv"))
message("Saved: myeloid_input_cell_counts.csv")

message("Myeloid input summary:")
print(sample_check %>% group_by(celltype_manual) %>%
        summarise(n = sum(n_cells), .groups = "drop"))


# -------------------------------------------------------------------------
# B. Re-process myeloid subset
# -------------------------------------------------------------------------

message("=== Re-processing myeloid subset ===")

DefaultAssay(myeloid) <- "RNA"

# 1. NormalizeData
myeloid <- NormalizeData(myeloid, normalization.method = "LogNormalize",
                         scale.factor = 10000)
message("NormalizeData done")

# 2. FindVariableFeatures
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst",
                                nfeatures = 2000)
message("FindVariableFeatures done")

# 3. ScaleData
myeloid <- ScaleData(myeloid, features = rownames(myeloid))
gc()
message("ScaleData done")

# 4. RunPCA (compute 30 PCs, use NPCS=25 for downstream)
myeloid <- RunPCA(myeloid, features = VariableFeatures(myeloid), npcs = 30,
                  verbose = FALSE)
message("RunPCA done")

# 5. Harmony batch/sample integration with fallback
harmony_success <- FALSE
harmony_note    <- ""

tryCatch({
  myeloid <- RunHarmony(
    myeloid,
    group.by.vars     = "sample_id",
    reduction.use     = "pca",
    dims.use          = 1:NPCS,
    theta             = 2,
    plot_convergence  = FALSE
  )
  harmony_success <- TRUE
  message("RunHarmony done")
}, error = function(e) {
  harmony_note <<- paste("Harmony failed:", e$message,
                         "- Falling back to PCA-based UMAP")
  message(harmony_note)
})

# 6. RunUMAP
if (harmony_success) {
  myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:NPCS,
                     verbose = FALSE)
} else {
  myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:NPCS,
                     verbose = FALSE)
}
message("RunUMAP done")

# 7. FindNeighbors
if (harmony_success) {
  myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:NPCS,
                           verbose = FALSE)
} else {
  myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:NPCS,
                           verbose = FALSE)
}
message("FindNeighbors done")

# 8. FindClusters at multiple resolutions
resolutions <- c(0.2, 0.4, 0.6, 0.8)

for (res in resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  myeloid <- FindClusters(myeloid, resolution = res, verbose = FALSE)
  # Rename to consistent naming
  if (col_name %in% colnames(myeloid@meta.data)) {
    # Already named correctly by Seurat
  }
  message("FindClusters done: resolution = ", res,
          " (", length(unique(myeloid@meta.data[[col_name]])), " clusters)")
}
gc()

# --- Resolution UMAPs ---
for (res in resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  n_clust  <- length(unique(myeloid@meta.data[[col_name]]))

  p <- DimPlot(myeloid, reduction = "umap", group.by = col_name,
               label = TRUE, repel = TRUE, pt.size = 0.5, shuffle = TRUE) +
    ggtitle(paste0("Myeloid UMAP: resolution ", res, " (", n_clust, " clusters)")) +
    theme(aspect.ratio = 1)

  pdf(file.path(FIG_DIR, paste0("myeloid_umap_res_", res, ".pdf")),
      width = 9, height = 8)
  print(p)
  dev.off()
  message("Saved: myeloid_umap_res_", res, ".pdf (", n_clust, " clusters)")
}

# --- Resolution summary table ---
res_summary <- data.frame(
  resolution  = resolutions,
  n_clusters  = sapply(resolutions, function(r) {
    length(unique(myeloid@meta.data[[paste0("RNA_snn_res.", r)]]))
  }),
  stringsAsFactors = FALSE
)
write_csv(res_summary, file.path(TAB_DIR, "myeloid_resolution_summary.csv"))
message("Saved: myeloid_resolution_summary.csv")
print(res_summary)

# =========================================================================
# MANUAL_REVIEW #1 — Cluster resolution choice
# =========================================================================
# TEMPORARY DEFAULT: resolution = 0.4
# 
# This is a PRELIMINARY choice.  You MUST inspect the four UMAP plots
# (myeloid_umap_res_0.2/0.4/0.6/0.8.pdf), the resolution summary table,
# and the DotPlot / FeaturePlot outputs in Section E before finalizing.
#
# Guidelines:
#   - 0.2: Coarse, may merge biologically distinct myeloid states
#   - 0.4: Moderate granularity; typically 6-10 clusters for ~13K myeloid cells
#   - 0.6: Higher granularity; may reveal transitional states
#   - 0.8: Fine granularity; risk of over-fragmentation
#
# Choose the resolution that best separates known myeloid subtypes
# (Kupffer vs monocyte-derived vs neutrophil states) without fragmenting
# homogeneous populations into spurious clusters.
# =========================================================================

DEFAULT_RES <- 0.4
myeloid$myeloid_clusters <- myeloid@meta.data[[paste0("RNA_snn_res.", DEFAULT_RES)]]
Idents(myeloid) <- myeloid$myeloid_clusters
message("Using temporary resolution: ", DEFAULT_RES,
        " (", length(unique(myeloid$myeloid_clusters)), " clusters)")
message("MANUAL_REVIEW #1: Verify this resolution before trusting annotations.")

# -------------------------------------------------------------------------
# C. (consolidated into Section B — resolution choice above)
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# D. Myeloid subcluster marker identification
# -------------------------------------------------------------------------

message("=== Myeloid subcluster marker identification ===")

# Seurat v5: join sample-level layers before DE testing
myeloid <- JoinLayers(myeloid)

# --- Marker identification: try presto first, fall back to Seurat ---
presto_available <- requireNamespace("presto", quietly = TRUE)

if (presto_available) {
  message("Using presto::wilcoxauc() for fast marker detection ...")
  library(presto)

  # Input: log-normalised data matrix + cluster assignments
  expr_mat  <- GetAssayData(myeloid, assay = "RNA", layer = "data")
  clust_vec <- as.character(myeloid$myeloid_clusters)

  presto_raw <- wilcoxauc(expr_mat, clust_vec)
  gc()

  # Save raw presto output (full table with AUC)
  write_csv(presto_raw, file.path(TAB_DIR, "myeloid_cluster_markers_presto.csv"))
  message("Saved: myeloid_cluster_markers_presto.csv (", nrow(presto_raw), " rows)")

  # Format to Seurat-like marker table (positive LFC, padj < 0.05, logFC > 0.25)
  all_markers <- presto_raw %>%
    filter(padj < 0.05, logFC > 0.25) %>%
    rename(
      cluster    = group,
      gene       = feature,
      avg_log2FC = logFC,
      p_val_adj  = padj,
      p_val      = pval,
      pct.1      = pct_in,
      pct.2      = pct_out
    ) %>%
    select(cluster, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, auc)

  message("presto markers after filtering: ", nrow(all_markers), " rows")

} else {
  message("presto not available, falling back to Seurat::FindAllMarkers() ...")
  message("To speed up marker detection, consider installing presto:")
  message("  R -e \"install.packages('presto', repos='https://cloud.r-project.org')\"")

  all_markers <- FindAllMarkers(
    myeloid,
    only.pos       = TRUE,
    min.pct        = 0.25,
    logfc.threshold = 0.25,
    test.use       = "wilcox"
  )
  gc()
}

# Save filtered marker table
write_csv(all_markers, file.path(TAB_DIR, "myeloid_cluster_markers.csv"))
message("Saved: myeloid_cluster_markers.csv (", nrow(all_markers), " rows)")

if (nrow(all_markers) == 0) {
  stop("Marker detection returned 0 rows. Check JoinLayers and assay setup.")
}

# Top markers per cluster (top 20 by avg_log2FC)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  ungroup()
write_csv(top_markers, file.path(TAB_DIR, "myeloid_top_markers_per_cluster.csv"))
message("Saved: myeloid_top_markers_per_cluster.csv")


# -------------------------------------------------------------------------
# E. Annotation marker visualization
# -------------------------------------------------------------------------

message("=== Annotation marker visualization ===")

# Canonical marker sets for myeloid subcluster annotation
myeloid_marker_sets <- list(
  "Resident_Kupffer"        = c("Clec4f", "Vsig4", "Cd5l", "Timd4", "Marco",
                                "Folr2", "Cd163", "C1qa", "C1qb", "C1qc"),
  "Monocyte_derived_Mph"    = c("Lyz2", "Ccr2", "Ly6c2", "Fcgr1", "Ms4a6c",
                                "Cx3cr1", "Plac8", "Chil3"),
  "Inflammatory_Mph"        = c("Il1b", "Tnf", "Nfkbia", "Ccl2", "Ccl3",
                                "Ccl4", "Cxcl2", "Ptgs2", "Nos2"),
  "Injury_response_Mph"     = c("Lgals3", "Spp1", "Gpnmb", "Trem2", "Apoe",
                                "Lpl", "Ctsb", "Ctsd"),
  "Neutrophil"              = c("S100a8", "S100a9", "Ly6g", "Cxcr2", "Retnlg",
                                "Camp", "Ltf", "Ngp", "Mpo", "Csf3r"),
  "Activated_neutrophil"    = c("Il1b", "Trem1", "Cxcl2", "Ccl3", "Ccl4",
                                "Slc16a3", "Hcar2", "Cd14"),
  "Interferon_response"     = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1",
                                "Oas1a", "Irf7", "Stat1"),
  "Antigen_presentation"    = c("H2-Aa", "H2-Ab1", "Cd74", "H2-Eb1", "Flt3",
                                "Itgax", "Ccr7"),
  "Cycling"                 = c("Mki67", "Top2a", "Stmn1", "Pclaf", "Ccna2"),
  "Contamination_warning"   = c("Siglech", "Ccr9", "Bst2",          # pDC-like
                                "Xcr1", "Clec9a",                    # cDC1
                                "Pecam1", "Kdr", "Cdh5", "Lyve1", "Stab2", # endothelial
                                "Alb", "Apoa1", "Ttr")              # hepatocyte ambient
)

# --- DotPlot of canonical markers ---
all_canonical <- unique(unlist(myeloid_marker_sets))
genes_avail <- intersect(all_canonical, rownames(myeloid))
message("Canonical markers available for DotPlot: ", length(genes_avail), "/",
        length(all_canonical))

p_dot <- DotPlot(myeloid, features = genes_avail,
                 group.by = "myeloid_clusters", assay = "RNA") +
  RotatedAxis() +
  ggtitle("Myeloid canonical markers by subcluster") +
  theme(axis.text.x = element_text(size = 6))

pdf(file.path(FIG_DIR, "myeloid_dotplot_canonical_markers.pdf"),
    width = max(18, length(genes_avail) * 0.3), height = 7)
print(p_dot)
dev.off()
message("Saved: myeloid_dotplot_canonical_markers.pdf")

# --- FeaturePlot of key markers in batches of 6 ---
key_markers_myeloid <- c(
  # Resident Kupffer
  "Clec4f", "Vsig4", "Cd5l", "Timd4", "Cd163",
  # Monocyte-derived
  "Lyz2", "Ccr2", "Fcgr1", "Ms4a6c", "Cx3cr1", "Plac8",
  # Inflammatory macrophage
  "Il1b", "Tnf", "Ccl2", "Ptgs2", "Nos2",
  # Injury-response
  "Lgals3", "Spp1", "Gpnmb", "Trem2", "Apoe", "Ctsd",
  # Neutrophil
  "S100a8", "S100a9", "Ly6g", "Cxcr2", "Camp", "Ltf", "Mpo", "Csf3r",
  # Activated neutrophil
  "Trem1", "Hcar2", "Cd14",
  # IFN response
  "Ifit1", "Ifit3", "Isg15", "Mx1", "Irf7", "Stat1",
  # Antigen presentation
  "H2-Aa", "H2-Ab1", "Cd74",
  # Cycling
  "Mki67", "Top2a", "Stmn1",
  # Contamination check
  "Siglech", "Pecam1", "Alb"
)
key_markers_myeloid <- intersect(key_markers_myeloid, rownames(myeloid))

batch_size <- 6
batches <- split(key_markers_myeloid,
                 ceiling(seq_along(key_markers_myeloid) / batch_size))

for (bname in names(batches)) {
  genes_batch <- batches[[bname]]
  p <- FeaturePlot(myeloid, features = genes_batch, reduction = "umap",
                   pt.size = 0.5, ncol = min(3, length(genes_batch)),
                   order = TRUE, combine = TRUE) +
    plot_annotation(title = paste0("Myeloid key markers (batch ", bname, ")"))
  pdf(file.path(FIG_DIR, paste0("myeloid_featureplot_key_markers_batch", bname, ".pdf")),
      width = 14, height = 4 * ceiling(length(genes_batch) / 3))
  print(p)
  dev.off()
}
message("Saved: myeloid_featureplot_key_markers_batch*.pdf (",
        length(batches), " files)")

# --- UMAP by sample ---
p_samp <- DimPlot(myeloid, reduction = "umap", group.by = "sample_id",
                  pt.size = 0.5, shuffle = TRUE) +
  ggtitle("Myeloid UMAP by sample") + theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "myeloid_umap_by_sample.pdf"), width = 10, height = 8)
print(p_samp)
dev.off()
message("Saved: myeloid_umap_by_sample.pdf")

# --- UMAP by group ---
p_grp <- DimPlot(myeloid, reduction = "umap", group.by = "group",
                 pt.size = 0.5, shuffle = TRUE) +
  ggtitle("Myeloid UMAP by group (Control vs APAP)") + theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "myeloid_umap_by_group.pdf"), width = 9, height = 8)
print(p_grp)
dev.off()
message("Saved: myeloid_umap_by_group.pdf")

# --- UMAP by original celltype ---
p_orig <- DimPlot(myeloid, reduction = "umap", group.by = "celltype_manual",
                  pt.size = 0.5, shuffle = TRUE, label = TRUE, repel = TRUE) +
  ggtitle("Myeloid UMAP by original cell type") + theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "myeloid_umap_by_original_celltype.pdf"),
    width = 10, height = 8)
print(p_orig)
dev.off()
message("Saved: myeloid_umap_by_original_celltype.pdf")

# --- UMAP by subcluster ---
p_sub <- DimPlot(myeloid, reduction = "umap", group.by = "myeloid_clusters",
                 label = TRUE, repel = TRUE, pt.size = 0.5, shuffle = TRUE) +
  ggtitle("Myeloid UMAP by subcluster") + theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "myeloid_umap_by_subcluster.pdf"), width = 9, height = 8)
print(p_sub)
dev.off()
message("Saved: myeloid_umap_by_subcluster.pdf")


# -------------------------------------------------------------------------
# F. Preliminary myeloid subcluster annotation
# -------------------------------------------------------------------------

message("=== Preliminary subcluster annotation ===")

# --- Score each cluster against each marker set ---
# Use AverageExpression to get mean expression per cluster
avg_expr <- AverageExpression(myeloid, assays = "RNA", group.by = "myeloid_clusters",
                              layer = "data", verbose = FALSE)$RNA

# AverageExpression prefixes cluster names with "g" (e.g., "g0", "g1")
# Strip prefix for consistent lookup against myeloid_clusters metadata
clusters_raw <- colnames(avg_expr)
clusters     <- gsub("^g", "", clusters_raw)
score_matrix <- data.frame(cluster = clusters, stringsAsFactors = FALSE)

for (set_name in names(myeloid_marker_sets)) {
  genes_in_set <- intersect(myeloid_marker_sets[[set_name]], rownames(avg_expr))
  if (length(genes_in_set) == 0) {
    score_matrix[[set_name]] <- NA_real_
    next
  }
  # Mean expression of all genes in the set per cluster
  if (length(genes_in_set) == 1) {
    set_scores <- avg_expr[genes_in_set, ]
  } else {
    set_scores <- colMeans(as.matrix(avg_expr[genes_in_set, ]))
  }
  score_matrix[[set_name]] <- round(set_scores, 4)
}

# Determine best-scoring set per cluster (excluding Contamination_warning and Cycling)
annotation_sets <- setdiff(names(myeloid_marker_sets),
                           c("Contamination_warning", "Cycling"))

best_set <- apply(score_matrix[, annotation_sets, drop = FALSE], 1, function(x) {
  valid <- !is.na(x)
  if (!any(valid)) return("Unknown myeloid state")
  names(which.max(x[valid]))
})

score_matrix$preliminary_annotation <- best_set

# Add contamination warning scores
if ("Contamination_warning" %in% names(myeloid_marker_sets)) {
  cw_genes <- intersect(myeloid_marker_sets[["Contamination_warning"]],
                        rownames(avg_expr))
  if (length(cw_genes) > 0) {
    if (length(cw_genes) == 1) {
      cw_score <- avg_expr[cw_genes, ]
    } else {
      cw_score <- colMeans(as.matrix(avg_expr[cw_genes, ]))
    }
    score_matrix$contamination_score <- round(cw_score, 4)
  }
}

write_csv(score_matrix,
          file.path(TAB_DIR, "myeloid_preliminary_subcluster_annotation.csv"))
message("Saved: myeloid_preliminary_subcluster_annotation.csv")
print(score_matrix)

# =========================================================================
# MANUAL_REVIEW #2 — Myeloid subcluster annotation (COMPLETED)
# =========================================================================
# After reviewing DotPlot, FeaturePlots, and top markers per cluster,
# the following manual working annotation is applied.
# Key manual corrections relative to auto-scoring:
#   - Cluster 1: auto-scored as Monocyte_derived_Mph; manually corrected to
#     Resident Kupffer cell based on Clec4f/Vsig4/Timd4 co-expression
#   - Cluster 9: auto-scored as Neutrophil; manually corrected to
#     Possible doublet / hepatocyte ambient-contaminated neutrophil
#     (high Alb/Apoa1 ambient signal + S100a8/S100a9)
#   - Cluster 13: auto-scored as Neutrophil; manually corrected to
#     Cycling neutrophil (Mki67/Top2a/Stmn1 positive)
# These labels are WORKING ANNOTATIONS — not final publication-grade labels.
# =========================================================================

# --- Save auto-scored preliminary annotation for traceability ---
myeloid$preliminary_annotation <- unname(best_set[as.character(myeloid$myeloid_clusters)])
message("Saved auto-scored preliminary_annotation to metadata for reference")

# --- Apply manual working annotation ---
manual_myeloid_annotation <- c(
  "0"  = "Inflammatory neutrophil",
  "1"  = "Resident Kupffer cell",
  "2"  = "Antigen-presenting macrophage / Mo-Mac",
  "3"  = "DC-like antigen-presenting myeloid",
  "4"  = "Inflammatory neutrophil",
  "5"  = "Activated neutrophil",
  "6"  = "Mature neutrophil",
  "7"  = "Monocyte-derived macrophage",
  "8"  = "Activated inflammatory neutrophil",
  "9"  = "Possible doublet / hepatocyte ambient-contaminated neutrophil",
  "10" = "IFN-responsive inflammatory neutrophil",
  "11" = "Neutrophil with erythroid/ambient signal",
  "12" = "Antigen-presenting myeloid with lymphoid-like signal",
  "13" = "Cycling neutrophil",
  "14" = "CCR7+ DC-like antigen-presenting myeloid",
  "15" = "Rare injury-associated macrophage-like cells"
)

# Map cluster IDs to manual annotation
myeloid$myeloid_subtype <- unname(
  manual_myeloid_annotation[as.character(myeloid$myeloid_clusters)]
)

# Flag unannotated clusters
unannotated <- setdiff(unique(as.character(myeloid$myeloid_clusters)),
                       names(manual_myeloid_annotation))
if (length(unannotated) > 0) {
  warning("The following clusters have NO manual annotation: ",
          paste(unannotated, collapse = ", "))
  myeloid$myeloid_subtype[myeloid$myeloid_clusters %in% unannotated] <- "Unannotated"
}

message("Applied manual working annotation. Subtypes: ",
        paste(sort(unique(myeloid$myeloid_subtype)), collapse = "\n  "))

# --- Add manual annotation to score_matrix for cross-reference ---
score_matrix$manual_annotation <- manual_myeloid_annotation[score_matrix$cluster]
write_csv(score_matrix,
          file.path(TAB_DIR, "myeloid_preliminary_subcluster_annotation.csv"))
message("Updated: myeloid_preliminary_subcluster_annotation.csv (with manual_annotation column)")


# -------------------------------------------------------------------------
# G. Group comparison (Control vs APAP)
# -------------------------------------------------------------------------

message("=== Group comparison ===")

# --- Counts per sample per subcluster ---
count_df <- myeloid@meta.data %>%
  count(sample_id, group, myeloid_subtype, name = "n_cells")

write_csv(count_df, file.path(TAB_DIR, "myeloid_subcluster_counts_by_sample.csv"))
message("Saved: myeloid_subcluster_counts_by_sample.csv")

# --- Proportions per sample ---
prop_df <- count_df %>%
  group_by(sample_id) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(prop_df, file.path(TAB_DIR, "myeloid_subcluster_proportions_by_sample.csv"))
message("Saved: myeloid_subcluster_proportions_by_sample.csv")

# --- Stacked barplot ---
p_stack <- ggplot(prop_df, aes(x = sample_id, y = proportion,
                                fill = myeloid_subtype)) +
  geom_bar(stat = "identity", position = "fill", colour = "grey20",
           linewidth = 0.1) +
  facet_wrap(~ group, scales = "free_x") +
  labs(x = "Sample", y = "Proportion", fill = "Myeloid subtype",
       title = "Myeloid subcluster proportions by sample") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

pdf(file.path(FIG_DIR, "myeloid_subcluster_proportion_stacked.pdf"),
    width = 12, height = 7)
print(p_stack)
dev.off()
message("Saved: myeloid_subcluster_proportion_stacked.pdf")

# --- Group-level barplot ---
prop_group <- prop_df %>%
  group_by(group, myeloid_subtype) %>%
  summarise(
    mean_prop = mean(proportion),
    sd_prop   = sd(proportion),
    n_samples = n(),
    .groups   = "drop"
  )

p_group <- ggplot(prop_group, aes(x = myeloid_subtype, y = mean_prop,
                                   fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, colour = "grey20", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Myeloid subtype", y = "Mean proportion", fill = "Group",
       title = "Mean myeloid subcluster proportions: Control vs APAP") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(FIG_DIR, "myeloid_subcluster_proportion_grouped.pdf"),
    width = 12, height = 7)
print(p_group)
dev.off()
message("Saved: myeloid_subcluster_proportion_grouped.pdf")

# --- Exploratory Wilcoxon test per subcluster ---
# WARNING: n = 3 per group — p-values are EXPLORATORY only, NOT confirmatory.
subtypes <- unique(prop_df$myeloid_subtype)
wilcox_results <- data.frame(
  myeloid_subtype = character(0),
  W_stat          = numeric(0),
  p_value         = numeric(0),
  mean_ctrl       = numeric(0),
  mean_apap       = numeric(0),
  direction       = character(0),
  stringsAsFactors = FALSE
)

for (st in subtypes) {
  ctrl_vals <- prop_df$proportion[prop_df$myeloid_subtype == st &
                                  prop_df$group == "Control"]
  apap_vals <- prop_df$proportion[prop_df$myeloid_subtype == st &
                                  prop_df$group == "APAP"]

  if (length(ctrl_vals) < 2 || length(apap_vals) < 2) next
  if (all(ctrl_vals == 0) || all(apap_vals == 0)) next

  wt <- wilcox.test(apap_vals, ctrl_vals, exact = FALSE)

  wilcox_results <- rbind(wilcox_results, data.frame(
    myeloid_subtype = st,
    W_stat          = wt$statistic,
    p_value         = round(wt$p.value, 4),
    mean_ctrl       = round(mean(ctrl_vals), 4),
    mean_apap       = round(mean(apap_vals), 4),
    direction       = ifelse(mean(apap_vals) > mean(ctrl_vals),
                             "Up in APAP", "Down in APAP"),
    stringsAsFactors = FALSE
  ))
}

write_csv(wilcox_results,
          file.path(TAB_DIR, "myeloid_subcluster_proportion_wilcox.csv"))
message("Saved: myeloid_subcluster_proportion_wilcox.csv (EXPLORATORY — n=3 per group)")
print(wilcox_results)

# --- Cluster-level counts/proportions (by myeloid_clusters, not merged subtype) ---
count_cluster_df <- myeloid@meta.data %>%
  count(sample_id, group, myeloid_clusters, name = "n_cells")

write_csv(count_cluster_df,
          file.path(TAB_DIR, "myeloid_cluster_counts_by_sample.csv"))
message("Saved: myeloid_cluster_counts_by_sample.csv")

prop_cluster_df <- count_cluster_df %>%
  group_by(sample_id) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(prop_cluster_df,
          file.path(TAB_DIR, "myeloid_cluster_proportions_by_sample.csv"))
message("Saved: myeloid_cluster_proportions_by_sample.csv")

# --- Wilcoxon by myeloid_clusters (exploratory, n=3 per group) ---
clusters_unique <- sort(unique(prop_cluster_df$myeloid_clusters))
wilcox_cluster_results <- data.frame(
  myeloid_cluster = character(0),
  W_stat          = numeric(0),
  p_value         = numeric(0),
  mean_ctrl       = numeric(0),
  mean_apap       = numeric(0),
  direction       = character(0),
  stringsAsFactors = FALSE
)

for (cl in clusters_unique) {
  ctrl_vals <- prop_cluster_df$proportion[prop_cluster_df$myeloid_clusters == cl &
                                           prop_cluster_df$group == "Control"]
  apap_vals <- prop_cluster_df$proportion[prop_cluster_df$myeloid_clusters == cl &
                                           prop_cluster_df$group == "APAP"]

  if (length(ctrl_vals) < 2 || length(apap_vals) < 2) next
  if (all(ctrl_vals == 0) && all(apap_vals == 0)) next

  wt <- wilcox.test(apap_vals, ctrl_vals, exact = FALSE)

  wilcox_cluster_results <- rbind(wilcox_cluster_results, data.frame(
    myeloid_cluster = cl,
    W_stat          = wt$statistic,
    p_value         = round(wt$p.value, 4),
    mean_ctrl       = round(mean(ctrl_vals), 4),
    mean_apap       = round(mean(apap_vals), 4),
    direction       = ifelse(mean(apap_vals) > mean(ctrl_vals),
                             "Up in APAP", "Down in APAP"),
    stringsAsFactors = FALSE
  ))
}

write_csv(wilcox_cluster_results,
          file.path(TAB_DIR, "myeloid_cluster_proportion_wilcox.csv"))
message("Saved: myeloid_cluster_proportion_wilcox.csv (EXPLORATORY — n=3 per group)")

# --- QC table: myeloid_clusters × sample_id × group ---
qc_sample <- myeloid@meta.data %>%
  count(myeloid_clusters, sample_id, group, name = "n_cells") %>%
  arrange(myeloid_clusters, group, sample_id)

write_csv(qc_sample,
          file.path(TAB_DIR, "myeloid_cluster_by_sample_qc.csv"))
message("Saved: myeloid_cluster_by_sample_qc.csv")

# --- QC table: myeloid_clusters × original celltype ---
qc_celltype <- myeloid@meta.data %>%
  count(myeloid_clusters, celltype_manual, name = "n_cells") %>%
  arrange(myeloid_clusters, desc(n_cells))

write_csv(qc_celltype,
          file.path(TAB_DIR, "myeloid_cluster_by_original_celltype_qc.csv"))
message("Saved: myeloid_cluster_by_original_celltype_qc.csv")


# -------------------------------------------------------------------------
# H. Module scores for myeloid functional states
# -------------------------------------------------------------------------

message("=== Module scores ===")

# Gene sets for module scoring
module_sets <- list(
  "Kupffer_resident_score"      = c("Clec4f", "Vsig4", "Cd5l", "Timd4",
                                     "Marco", "Folr2", "Cd163", "C1qa",
                                     "C1qb", "C1qc"),
  "Monocyte_derived_Mph_score"  = c("Lyz2", "Ccr2", "Ly6c2", "Fcgr1",
                                     "Ms4a6c", "Cx3cr1", "Plac8", "Chil3"),
  "Inflammatory_Mph_score"      = c("Il1b", "Tnf", "Nfkbia", "Ccl2", "Ccl3",
                                     "Ccl4", "Cxcl2", "Ptgs2", "Nos2"),
  "Injury_response_Mph_score"   = c("Lgals3", "Spp1", "Gpnmb", "Trem2",
                                     "Apoe", "Lpl", "Ctsb", "Ctsd"),
  "Neutrophil_score"            = c("S100a8", "S100a9", "Ly6g", "Cxcr2",
                                     "Retnlg", "Camp", "Ltf", "Ngp", "Mpo",
                                     "Csf3r"),
  "Activated_neutrophil_score"  = c("Il1b", "Trem1", "Cxcl2", "Ccl3", "Ccl4",
                                     "Slc16a3", "Hcar2", "Cd14"),
  "Interferon_response_score"   = c("Ifit1", "Ifit2", "Ifit3", "Isg15",
                                     "Mx1", "Oas1a", "Irf7", "Stat1"),
  "Antigen_presentation_score"  = c("H2-Aa", "H2-Ab1", "Cd74", "H2-Eb1",
                                     "Flt3", "Itgax", "Ccr7")
)

# AddModuleScore one at a time for clean column names
for (score_name in names(module_sets)) {
  genes_in_set <- intersect(module_sets[[score_name]], rownames(myeloid))
  n_genes <- length(genes_in_set)

  if (n_genes < 3) {
    warning("Skipping ", score_name, ": only ", n_genes,
            " genes detected (minimum 3 required)")
    next
  }

  message("Computing AddModuleScore for ", score_name, " (", n_genes, " genes)")

  myeloid <- AddModuleScore(myeloid, features = list(genes_in_set),
                            name = score_name, assay = "RNA", ctrl = 100)
  gc()
}

# Identify score columns (AddModuleScore appends "1" to name)
score_cols <- grep("_score1$", colnames(myeloid@meta.data), value = TRUE)
message("Module score columns: ", paste(score_cols, collapse = ", "))

if (length(score_cols) == 0) {
  warning("No module score columns found. Skipping score visualizations.")
} else {
  # --- Violin plots by group (combined) ---
  plot_list <- list()
  for (sc in score_cols) {
    label <- gsub("_score1$", "", sc)
    p <- VlnPlot(myeloid, features = sc, group.by = "group",
                 pt.size = 0.05, assay = "RNA") +
      ggtitle(label) + NoLegend() + theme(aspect.ratio = 1)
    plot_list[[sc]] <- p
  }

  ncol_vln <- min(3, length(plot_list))
  nrow_vln <- ceiling(length(plot_list) / ncol_vln)

  p_vln <- wrap_plots(plot_list, ncol = ncol_vln) +
    plot_annotation(title = "Myeloid module scores by group")

  pdf(file.path(FIG_DIR, "myeloid_module_scores_violin_by_group.pdf"),
      width = 5 * ncol_vln, height = 4 * nrow_vln)
  print(p_vln)
  dev.off()
  message("Saved: myeloid_module_scores_violin_by_group.pdf")

  # --- FeaturePlot of module scores on UMAP ---
  ncol_feat <- min(3, length(score_cols))
  nrow_feat <- ceiling(length(score_cols) / ncol_feat)

  p_feat <- FeaturePlot(myeloid, features = score_cols, reduction = "umap",
                        pt.size = 0.5, ncol = ncol_feat,
                        order = TRUE, combine = TRUE) +
    plot_annotation(title = "Module scores on myeloid UMAP")

  pdf(file.path(FIG_DIR, "myeloid_module_scores_featureplot.pdf"),
      width = 5 * ncol_feat, height = 4 * nrow_feat)
  print(p_feat)
  dev.off()
  message("Saved: myeloid_module_scores_featureplot.pdf")

  # --- Module scores per cell ---
  score_df <- myeloid@meta.data %>%
    select(sample_id, group, myeloid_clusters, myeloid_subtype, all_of(score_cols))
  write_csv(score_df, file.path(TAB_DIR, "myeloid_module_scores_per_cell.csv"))
  message("Saved: myeloid_module_scores_per_cell.csv (", nrow(score_df), " cells)")

  # --- Module scores aggregated by sample_id / group / myeloid_clusters ---
  score_agg <- score_df %>%
    group_by(sample_id, group, myeloid_clusters, myeloid_subtype) %>%
    summarise(
      n_cells = n(),
      across(all_of(score_cols), mean, .names = "{.col}_mean"),
      .groups = "drop"
    )
  write_csv(score_agg, file.path(TAB_DIR, "myeloid_module_scores_aggregated.csv"))
  message("Saved: myeloid_module_scores_aggregated.csv (", nrow(score_agg), " rows)")
}


# -------------------------------------------------------------------------
# I. Save final object and sessionInfo
# -------------------------------------------------------------------------

message("=== Saving outputs ===")

saveRDS(myeloid, file.path(PROC_DIR, "myeloid_subclustered.rds"))
message("Saved: myeloid_subclustered.rds (", ncol(myeloid), " cells)")

writeLines(capture.output(sessionInfo()),
           file.path(TAB_DIR, "sessionInfo_03_Myeloid_Subclustering.txt"))
message("Saved: sessionInfo_03_Myeloid_Subclustering.txt")

# -------------------------------------------------------------------------
# J. Write stage-3 notes (Chinese)
# -------------------------------------------------------------------------

message("=== Writing stage-3 notes ===")

n_clust <- length(unique(myeloid$myeloid_clusters))
n_subtypes <- length(unique(myeloid$myeloid_subtype))
total_myeloid <- ncol(myeloid)

# Count cells per original celltype
orig_counts <- myeloid@meta.data %>%
  count(celltype_manual, name = "n") %>%
  arrange(desc(n))

orig_lines <- paste0("- ", orig_counts$celltype_manual, "：",
                     orig_counts$n, " cells")
orig_summary <- paste(orig_lines, collapse = "\n")

# Resolution summary lines
res_lines <- apply(res_summary, 1, function(r) {
  paste0("- ", r["resolution"], "：", r["n_clusters"], " clusters")
})
res_summary_text <- paste(res_lines, collapse = "\n")

# Current subtypes
subtype_list <- paste(unique(myeloid$myeloid_subtype), collapse = "、")

note_lines <- c(
  "# 第三阶段 Myeloid Subclustering 分析笔记",
  "",
  paste0("**日期**：", Sys.Date()),
  "**脚本**：`scripts/03_Myeloid_Subclustering.R`",
  "**版本**：人工 working annotation（16 clusters，基于 res=0.4）",
  "",
  "---",
  "",
  "## 1. 本阶段为什么选择髓系细胞",
  "",
  "髓系细胞是 APAP 急性肝损伤中最早响应且最核心的效应细胞群体。",
  "文献表明，APAP 过量后 Kupffer cell 来源的促炎因子驱动初始损伤，",
  "随后 monocyte-derived macrophage 浸润并参与炎症消退和组织修复。",
  "中性粒细胞在早期炎症中也发挥重要作用。",
  "",
  "Stage 2 全局注释识别了 5 个髓系相关 cluster：Kupffer cell、",
  "Monocyte-derived macrophage/DC-like、Inflammatory neutrophil、",
  "Neutrophil、Activated neutrophil，合计 13,556 cells，",
  "数量充足且 marker 表达清晰，具备 subclustering 条件。",
  "",
  "本阶段目标：在髓系内部发现更精细的亚群结构，",
  "评估 Control vs APAP 的组成和功能状态变化趋势。",
  "",
  "---",
  "",
  "## 2. 输入对象、纳入和排除的细胞类型",
  "",
  "**输入**：`data/processed/02_Global_Annotation/seurat_annotated.rds`",
  paste0("- ", ncol(myeloid), " myeloid cells x ", nrow(myeloid), " genes"),
  "- 6 samples（3 Control + 3 APAP）",
  "",
  "**纳入的细胞类型**：",
  orig_summary,
  "",
  "**排除的细胞类型**：",
  "- pDC-like（虽为髓系来源但属于 pDC lineage）",
  "- cDC1（属于 DC lineage）",
  "- Mast / basophil-like（细胞数少，身份需确认）",
  "- Endothelial（肝脏结构细胞）",
  "- Low-quality / mixed（质量控制排除）",
  "- B cell / possible doublet",
  "- Hepatocyte contamination",
  "- Cholangiocyte / epithelial",
  "",
  "---",
  "",
  "## 3. 重新聚类使用的参数",
  "",
  "- NormalizeData：LogNormalize，scale factor = 10,000",
  "- FindVariableFeatures：vst，nfeatures = 2,000",
  "- RunPCA：npcs = 30",
  paste0("- Harmony：group.by.vars = \"sample_id\"，dims = 1:", NPCS,
         if (harmony_success) "（成功）" else paste0("（失败，原因：", harmony_note, "）")),
  paste0("- RunUMAP：reduction = ",
         if (harmony_success) "\"harmony\"" else "\"pca\"",
         "，dims = 1:", NPCS),
  paste0("- FindNeighbors：reduction = ",
         if (harmony_success) "\"harmony\"" else "\"pca\"",
         "，dims = 1:", NPCS),
  "- FindClusters：resolution = 0.2, 0.4, 0.6, 0.8",
  "- presto::wilcoxauc() 用于 marker 检测（padj < 0.05, logFC > 0.25）",
  "",
  "---",
  "",
  "## 4. Resolution 选择依据",
  "",
  "测试了 4 个 resolution：",
  res_summary_text,
  "",
  paste0("最终选择 **resolution = ", DEFAULT_RES, "**（", n_clust, " clusters）。"),
  "",
  "resolution = 0.4 在髓系内部提供了适中的粒度，",
  "能够区分 Kupffer cell、monocyte-derived macrophage、",
  "neutrophil 的不同激活状态和功能亚群，",
  "同时避免 resolution = 0.6/0.8 的过度碎片化。",
  "",
  "---",
  "",
  "## 5. 每个亚群的 marker 和人工命名（working annotation）",
  "",
  "**重要声明：以下标签为 working annotation，不是最终发表级别的细胞注释。**",
  "注释基于 DotPlot、FeaturePlot、top marker 和文献中已知的髓系 marker 综合判断。",
  "",
  "当前 16 个 cluster 的人工注释：",
  paste0("- Cluster ", names(manual_myeloid_annotation), "：", manual_myeloid_annotation, collapse = "\n"),
  "",
  "**重点人工修正（与自动 marker 评分不同的 cluster）**：",
  "- **Cluster 1**：自动评分判为 Monocyte_derived_Mph，人工复核发现高表达 Clec4f/Vsig4/Timd4，",
  "  修正为 Resident Kupffer cell。Kupffer 标记物在此 cluster 的表达模式与 cluster 7 (真正的 Mo-Mac) 明显不同。",
  "- **Cluster 9**：自动评分判为 Neutrophil，人工复核发现同时携带高 Alb/Apoa1 环境信号 + S100a8/S100a9，",
  "  修正为 Possible doublet / hepatocyte ambient-contaminated neutrophil。",
  "  这些细胞可能并非真实的生物学亚群，需在后续分析中审慎对待。",
  "- **Cluster 13**：自动评分判为 Neutrophil，人工复核发现 Mki67/Top2a/Stmn1 阳性，",
  "  修正为 Cycling neutrophil。代表了处于增殖状态的中性粒细胞前体或活化群体。",
  "",
  "详细 marker 和评分见：",
  "- `results/tables/03_Myeloid_Subclustering/myeloid_top_markers_per_cluster.csv`",
  "- `results/tables/03_Myeloid_Subclustering/myeloid_preliminary_subcluster_annotation.csv`",
  "- `results/tables/03_Myeloid_Subclustering/myeloid_cluster_markers_presto.csv`（完整 presto 输出含 AUC）",
  "- `results/figures/03_Myeloid_Subclustering/myeloid_dotplot_canonical_markers.pdf`",
  "- `results/figures/03_Myeloid_Subclustering/myeloid_featureplot_key_markers_batch*.pdf`",
  "",
  "---",
  "",
  "## 6. 哪些亚群最不确定，需要人工复核",
  "",
  "以下 cluster 的身份不确定性较高，建议后续实验验证或谨慎解释：",
  "",
  "- **Cluster 9**（Possible doublet / hepatocyte ambient-contaminated neutrophil）：",
  "  高 contamination_score，可能为实验 artifact 或双细胞。不建议用于生物学结论。",
  "- **Cluster 11**（Neutrophil with erythroid/ambient signal）：",
  "  可能携带红细胞或环境 RNA 污染信号，需确认其 marker 表达是否独立于 contamination。",
  "- **Cluster 12**（Antigen-presenting myeloid with lymphoid-like signal）：",
  "  同时表达髓系和淋系 marker，需排除 doublet 可能性。",
  "- **Cluster 15**（Rare injury-associated macrophage-like cells）：",
  "  细胞数极少的稀有群体，可能为真实稀有亚群或 outlier。",
  "",
  "一般注意事项：",
  "- 如果一个 cluster 同时高表达多个 lineage marker，可能为 doublet 或过渡状态",
  "- 如果 contamination_score 较高，需检查是否为实验 artifact",
  "- 如果某个 cluster 细胞数很少（< 50 cells），可能为 outlier cluster",
  "- 如果某 cluster 仅在单个 sample 中出现，可能为 sample-specific batch effect",
  "",
  "---",
  "",
  "## 7. Control vs APAP 的比例变化，谨慎描述",
  "",
  "**重要声明**：n = 3 Control vs 3 APAP，Wilcoxon 秩和检验 p 值仅作探索性参考，",
  "**不能作为强统计结论**。需要更大样本量或独立的验证队列。",
  "",
  "分别输出了按 myeloid_subtype 和按 myeloid_clusters 的 Wilcoxon 结果：",
  "- `myeloid_subcluster_proportion_wilcox.csv`",
  "- `myeloid_cluster_proportion_wilcox.csv`",
  "- `myeloid_cluster_counts_by_sample.csv`",
  "- `myeloid_cluster_proportions_by_sample.csv`",
  "",
  "主要趋势（待人工根据实际输出填写）：",
  "- 待填写...",
  "",
  "---",
  "",
  "## 8. Module score 的趋势，谨慎描述",
  "",
  paste0("计算了 ", length(score_cols), " 个基因集的 AddModuleScore。"),
  "Module score 反映的是基因集在单个细胞中的平均表达水平，",
  "不代表通路激活状态，不能作为机制证明。",
  "",
  "输出文件：",
  "- `myeloid_module_scores_per_cell.csv`：每个细胞的 module score（per-cell level）",
  "- `myeloid_module_scores_aggregated.csv`：按 sample_id/group/myeloid_clusters 聚合的均值表",
  "",
  "主要观察（待人工根据实际输出填写）：",
  "- 待填写...",
  "",
  "---",
  "",
  "## 9. 是否建议第四阶段进入 pseudotime",
  "",
  "**当前决策：不自动启动 Monocle3 pseudotime。**",
  "",
  "原因：髓系 subcluster 之间主要是功能状态差异（Kupffer vs Mo-Mac vs neutrophil），",
  "而非同一 lineage 的发育连续体。中性粒细胞各亚群（0/4/5/6/8/10/13）之间可能存在",
  "一定的激活梯度，但并非严格的 pseudotime 轨迹。",
  "",
  "替代方案：",
  "- 使用 module score（Inflammatory/Activated/IFN-response）在 cluster 间比较功能状态趋势",
  "- 使用 top marker 的表达模式描述激活/分化梯度",
  "",
  "**仅在以下情况下启动 Monocle3**：",
  "- UMAP 上中性粒细胞各亚群呈现明显的连续弧线结构",
  "- 或 Kupffer→Mo-Mac 之间存在清晰的过渡路径",
  "- 选择 root 需要明确的生物学依据（如选择 resident Kupffer 作为起始点）",
  "",
  "如果第四阶段确认进入 pseudotime，建议仅对中性粒细胞 lineage 或 Kupffer/Mo-Mac lineage",
  "分别独立运行，而非将所有髓系细胞混合跑 pseudotime。",
  "",
  "---",
  "",
  "## 10. 声明",
  "",
  "本阶段注释基于自动 marker 富集评分作为参考，经人工复核后应用 working annotation。",
  "标签不代表最终发表级别的细胞注释，部分 cluster 的身份需要在独立数据集中验证。",
  "module score 和 Wilcoxon 结果仅为探索性描述，不是机制证明。",
  "Cluster 1/9/13 为重点人工修正，与自动评分结果不同，已在第 5 节详细说明了修正依据。",
  "",
  "**此版本为 working annotation，不视为最终发表级别的细胞注释。**"
)

writeLines(note_lines, file.path(NOTE_DIR, "03_Myeloid_Subclustering_note.md"))
message("Saved: notes/03_Myeloid_Subclustering_note.md")

# -------------------------------------------------------------------------
# Completion
# -------------------------------------------------------------------------

message("")
message("=== 03_Myeloid_Subclustering.R complete ===")
message("")
message("Output summary:")
message("  Object:   data/processed/03_Myeloid_Subclustering/myeloid_subclustered.rds")
message("  Figures:  results/figures/03_Myeloid_Subclustering/ (", 
        length(list.files(FIG_DIR, pattern = "\\.pdf$")), " PDFs)")
message("  Tables:   results/tables/03_Myeloid_Subclustering/ (",
        length(list.files(TAB_DIR, pattern = "\\.csv$")), " CSVs)")
message("  Notes:    notes/03_Myeloid_Subclustering_note.md")
message("")
message("Next steps:")
message("  1. Review myeloid_umap_res_0.2/0.4/0.6/0.8.pdf and choose resolution")
message("  2. Review marker expression in DotPlot and FeaturePlot outputs")
message("  3. Fill in manual_myeloid_annotation in Section F")
message("  4. Re-run proportion comparison after annotation is stable")
message("  5. Review module score trends before deciding on Stage 4 (pseudotime)")
message("  6. If trajectory structure is unclear, prefer module scores / marker trends over Monocle3")

