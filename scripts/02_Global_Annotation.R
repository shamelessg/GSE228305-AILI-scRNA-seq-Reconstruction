# =========================================================================
# 02_Global_Annotation.R
# =========================================================================
# Purpose: Global clustering, marker identification, coarse cell-type
#          annotation, and Control-vs-APAP cell proportion comparison.
# Input:  data/processed/01_QC_and_Integration/seurat_integrated.rds
# Output: data/processed/02_Global_Annotation/seurat_annotated.rds
#         results/figures/02_Global_Annotation/*.pdf
#         results/tables/02_Global_Annotation/*.csv
# =========================================================================
#
# MANUAL INTERVENTION POINTS (search for "MANUAL_REVIEW"):
#   1. Cluster resolution choice (Section B)
#   2. Cell-type annotation labels (Section D)
#   3. Stage-3 myeloid vs. other lineage decision (Section G)
#
# The script will run end-to-end with a TEMPORARY default resolution,
# produce preliminary annotations, and save all outputs.  You MUST
# inspect the figures/tables and update the manual_review sections
# before trusting the annotation or proceeding to Stage 3.
# =========================================================================

# -------------------------------------------------------------------------
# 0. Setup
# -------------------------------------------------------------------------

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

PROJ_ROOT <- getwd()
setwd(PROJ_ROOT)

# --- directories ---
PROC_DIR   <- file.path(PROJ_ROOT, "data", "processed", "02_Global_Annotation")
FIG_DIR    <- file.path(PROJ_ROOT, "results", "figures",   "02_Global_Annotation")
TAB_DIR    <- file.path(PROJ_ROOT, "results", "tables",    "02_Global_Annotation")
STAGE1_RDS <- file.path(PROJ_ROOT, "data", "processed", "01_QC_and_Integration",
                        "seurat_integrated.rds")

dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# A. Read & validate stage-1 object
# -------------------------------------------------------------------------

if (!file.exists(STAGE1_RDS)) {
  stop("Stage-1 RDS not found: ", STAGE1_RDS)
}
message("Reading: ", STAGE1_RDS)
merged <- readRDS(STAGE1_RDS)
message("Loaded object: ", ncol(merged), " cells x ", nrow(merged), " genes")

# --- validate metadata ---
req_meta <- c("sample_id", "orig.ident", "percent.mt")
missing_meta <- setdiff(req_meta, colnames(merged@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ",
       paste(missing_meta, collapse = ", "))
}

# Defensive: unify group/condition column to "group"
has_group     <- "group"     %in% colnames(merged@meta.data)
has_condition <- "condition" %in% colnames(merged@meta.data)

if (!has_group && !has_condition) {
  stop("No 'group' or 'condition' column found in metadata.")
}
if (!has_group && has_condition) {
  merged$group <- merged$condition
  message("Renamed 'condition' -> 'group'")
}
if (has_group && has_condition) {
  # If both exist but differ, prefer the one with Control/APAP values
  grp_vals <- unique(merged$group)
  if (!any(c("Control", "APAP") %in% grp_vals)) {
    merged$group <- merged$condition
    message("Switched to 'condition' column (contains Control/APAP labels)")
  }
}

# --- validate group values ---
grp_levels <- unique(merged$group)
message("Group levels found: ", paste(grp_levels, collapse = ", "))
if (!all(c("Control", "APAP") %in% grp_levels)) {
  stop("Expected 'Control' and 'APAP' in group column. Found: ",
       paste(grp_levels, collapse = ", "))
}

# --- validate sample count ---
n_samples <- length(unique(merged$sample_id))
message("Number of samples: ", n_samples)
if (n_samples != 6) {
  warning("Expected 6 samples (3 Control + 3 APAP). Found: ", n_samples,
          ". Check sample_id naming.")
}

n_ctrl <- sum(grepl("^con", unique(merged$sample_id)))
n_apap <- sum(grepl("^AP",  unique(merged$sample_id)))
message("Detected: ", n_ctrl, " Control + ", n_apap, " APAP samples")
if (n_ctrl != 3 || n_apap != 3) {
  warning("Expected 3 Control + 3 APAP. Check sample_id naming patterns.")
}

# --- validate reductions ---
req_reduc <- c("pca", "harmony", "umap")
missing_reduc <- setdiff(req_reduc, names(merged@reductions))
if (length(missing_reduc) > 0) {
  stop("Missing required reductions: ",
       paste(missing_reduc, collapse = ", "))
}

# --- validate assays ---
if (!"RNA" %in% names(merged@assays)) {
  stop("Missing 'RNA' assay in Seurat object")
}

message("=== Stage-1 object validation passed ===")

# -------------------------------------------------------------------------
# B. Global clustering (multiple resolutions)
# -------------------------------------------------------------------------

NPCS <- 30  # same as stage 1

# Use Harmony reduction (already present)
message("Running FindNeighbors on harmony reduction, dims 1:", NPCS)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:NPCS)
gc()

# Cluster at multiple resolutions for evaluation
resolutions <- c(0.2, 0.4, 0.6, 0.8)
message("Clustering at resolutions: ", paste(resolutions, collapse = ", "))
merged <- FindClusters(merged, resolution = resolutions)
gc()

# --- UMAP and cluster count for each resolution ---
for (res in resolutions) {
  res_col <- paste0("RNA_snn_res.", res)
  if (!res_col %in% colnames(merged@meta.data)) {
    warning("Column ", res_col, " not found; skipping resolution ", res)
    next
  }

  n_clusters <- length(unique(merged@meta.data[[res_col]]))
  message("Resolution ", res, ": ", n_clusters, " clusters")

  # UMAP coloured by cluster
  p <- DimPlot(merged, reduction = "umap", group.by = res_col,
               label = TRUE, repel = TRUE, pt.size = 0.2, shuffle = TRUE) +
    ggtitle(paste0("Resolution = ", res, " (", n_clusters, " clusters)")) +
    theme(aspect.ratio = 1) +
    NoLegend()

  pdf(file.path(FIG_DIR, paste0("umap_res_", res, ".pdf")), width = 9, height = 8)
  print(p)
  dev.off()
  message("Saved: umap_res_", res, ".pdf")
}

# --- Resolution summary table ---
res_summary <- data.frame(
  resolution = resolutions,
  n_clusters = sapply(resolutions, function(r) {
    length(unique(merged@meta.data[[paste0("RNA_snn_res.", r)]]))
  })
)
write_csv(res_summary, file.path(TAB_DIR, "resolution_summary.csv"))
message("Saved: resolution_summary.csv")
print(res_summary)

# =========================================================================
# MANUAL_REVIEW #1 — Cluster resolution
# =========================================================================
# Resolution 0.4 is a MANUAL preliminary choice for coarse global annotation.
# It is NOT an algorithmic optimum — it balances interpretable cluster numbers
# (27 clusters) against over-fragmentation.
# Resolution summary (informational):
#   0.2 = 21 clusters
#   0.4 = 27 clusters  <-- current choice
#   0.6 = 34 clusters
#   0.8 = 39 clusters
# If subsequent DotPlot / FeaturePlot / marker review reveals over-merging
# (e.g., distinct cell types collapsed into one cluster), re-evaluate at
# resolution 0.6 and update DEFAULT_RES accordingly.
# =========================================================================

DEFAULT_RES <- 0.4

# Set active clusters to the chosen default resolution
default_col <- paste0("RNA_snn_res.", DEFAULT_RES)
merged$seurat_clusters <- merged@meta.data[[default_col]]
Idents(merged) <- "seurat_clusters"
message("Active clusters set to resolution = ", DEFAULT_RES,
        " (", length(unique(merged$seurat_clusters)), " clusters)")

# UMAP by cluster (default resolution) — main reference plot
p_cluster <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters",
                     label = TRUE, repel = TRUE, pt.size = 0.2, shuffle = TRUE) +
  ggtitle(paste0("UMAP by cluster (resolution = ", DEFAULT_RES, ")")) +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_by_cluster.pdf"), width = 10, height = 8)
print(p_cluster)
dev.off()
message("Saved: umap_by_cluster.pdf")

# -------------------------------------------------------------------------
# C. Global marker identification
# -------------------------------------------------------------------------

# Seurat v5: join sample-level layers before DE testing
merged <- JoinLayers(merged)

# --- Marker identification: try presto first, fall back to Seurat ---
# presto::wilcoxauc() is orders of magnitude faster than Seurat::FindAllMarkers()
# for large datasets.  We try it first; if unavailable, we fall back to the
# standard Seurat Wilcoxon test.

presto_available <- requireNamespace("presto", quietly = TRUE)

if (presto_available) {
  message("Using presto::wilcoxauc() for fast marker detection ...")
  library(presto)

  # Input: log-normalised data matrix + cluster assignments
  expr_mat  <- GetAssayData(merged, assay = "RNA", layer = "data")
  clust_vec <- as.character(merged$seurat_clusters)

  presto_raw <- wilcoxauc(expr_mat, clust_vec)
  gc()

  # Save raw presto output (full table with AUC)
  write_csv(presto_raw, file.path(TAB_DIR, "global_cluster_markers_presto.csv"))
  message("Saved: global_cluster_markers_presto.csv (", nrow(presto_raw), " rows)")

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

  all_markers <- FindAllMarkers(
    merged,
    only.pos       = TRUE,
    min.pct        = 0.25,
    logfc.threshold = 0.25,
    test.use       = "wilcox"
  )
  gc()
}

# Full table (from whichever method was used)
write_csv(all_markers, file.path(TAB_DIR, "global_cluster_markers.csv"))
message("Saved: global_cluster_markers.csv (", nrow(all_markers), " rows)")

if (nrow(all_markers) == 0) {
  stop("Marker detection returned 0 rows. Check JoinLayers and assay setup.")
}

# Top markers per cluster (top 20 by avg_log2FC)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  ungroup()
write_csv(top_markers, file.path(TAB_DIR, "top_markers_per_cluster.csv"))
message("Saved: top_markers_per_cluster.csv")

# Top 5 per cluster for DotPlot
top5 <- top_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  ungroup()
top5_genes <- unique(top5$gene)

p_dot_top <- DotPlot(merged, features = top5_genes, group.by = "seurat_clusters",
                     assay = "RNA") +
  RotatedAxis() +
  ggtitle("Top 5 markers per cluster") +
  theme(axis.text.x = element_text(size = 6))

pdf(file.path(FIG_DIR, "dotplot_top5_per_cluster.pdf"),
    width = max(14, length(top5_genes) * 0.35), height = 7)
print(p_dot_top)
dev.off()
message("Saved: dotplot_top5_per_cluster.pdf")

# -------------------------------------------------------------------------
# D. Coarse cell-type annotation
# -------------------------------------------------------------------------

# Canonical marker sets for mouse liver non-parenchymal cells
canonical_markers <- list(
  "T cell"                = c("Cd3d", "Cd3e", "Trac", "Cd4", "Cd8a", "Lef1", "Tcf7", "Il7r"),
  "B cell"                = c("Ms4a1", "Cd79a", "Cd79b", "Bank1", "Pax5", "Fcer2a", "Ighd"),
  "NK cell"               = c("Ncr1", "Nkg7", "Klrb1c", "Gzma", "Gzmb", "Prf1", "Xcl1", "Klra9"),
  "Neutrophil"            = c("S100a8", "S100a9", "Mpo", "Cxcr2", "Csf3r", "Ltf", "Camp", "Ngp", "Cd177", "Ly6g"),
  "Endothelial"           = c("Pecam1", "Kdr", "Cdh5", "Lyve1", "Stab2", "Robo4", "Mmrn2", "Ptprb", "Thbd", "Nos3"),
  "Kupffer cell"          = c("Cd5l", "Clec4f", "Folr2", "Fcna", "Cd163", "Vsig4", "C6", "Adgre1", "C1qa"),
  "Monocyte-derived Mph/DC" = c("Cx3cr1", "Lyz2", "Fcgr1", "Clec4b1", "Ms4a6c", "Cd209a", "Clec4a3"),
  "Cholangiocyte / Epithelial" = c("Epcam", "Krt8", "Krt18", "Krt19", "Sox9", "Dmbt1"),
  "Hepatocyte"            = c("Alb", "Apoa1", "Ttr", "Hnf4a"),
  "Fibroblast / Stellate" = c("Col1a1", "Col1a2", "Dcn", "Lum", "Dpep1", "Dpt", "Pdgfra", "Tcf21", "Mmp2"),
  "Mesothelial"           = c("Msln", "Muc16", "Upk3b", "Rspo1"),
  "Plasma cell"           = c("Jchain", "Mzb1", "Xbp1", "Tnfrsf17", "Igha"),
  "pDC"                   = c("Siglech", "Ccr9", "Flt3", "Cd300c"),
  "cDC1"                  = c("Xcr1", "Clec9a", "Flt3", "Itgax"),
  "Cycling"               = c("Mki67", "Top2a", "Stmn1", "Pclaf", "Ccna2", "Hist1h3c"),
  "Mast / Basophil"       = c("Cpa3", "Gata2", "Kit", "Prtn3"),
  "VSMC / Pericyte"       = c("Lmod1", "Cnn1", "Myh11", "Tagln")
)

# --- Preliminary annotation by marker enrichment ---
# For each cluster, compute the fraction of cells expressing each canonical
# marker, then score each cell-type signature by mean expression of detected
# markers.  Assign the cell type with the highest score.
# THIS IS PRELIMINARY — must be reviewed manually.

# Get average expression per cluster
DefaultAssay(merged) <- "RNA"
avg_expr <- AverageExpression(
  merged,
  assays       = "RNA",
  group.by     = "seurat_clusters",
  slot         = "data"
)[["RNA"]]  # genes x clusters matrix
avg_expr <- as.matrix(avg_expr)

# Strip "g" prefix that Seurat v5 AverageExpression adds to numeric cluster names
colnames(avg_expr) <- gsub("^g", "", colnames(avg_expr))

# Score each cluster for each cell type (mean expression of marker genes present in data)
cluster_ids <- colnames(avg_expr)
type_scores <- matrix(NA, nrow = length(cluster_ids),
                      ncol = length(canonical_markers))
rownames(type_scores) <- cluster_ids
colnames(type_scores) <- names(canonical_markers)

for (ct in names(canonical_markers)) {
  genes_present <- intersect(canonical_markers[[ct]], rownames(avg_expr))
  if (length(genes_present) == 0) {
    type_scores[, ct] <- 0
    next
  }
  # Mean expression of the marker set
  sub <- avg_expr[genes_present, , drop = FALSE]
  type_scores[, ct] <- colMeans(sub)
}

# Assign each cluster to the cell type with the highest score
prelim_assignment <- apply(type_scores, 1, function(x) names(which.max(x)))
prelim_df <- data.frame(
  cluster              = cluster_ids,
  preliminary_label    = prelim_assignment,
  stringsAsFactors     = FALSE
)
# Round scores for readability
scores_rounded <- as.data.frame(round(type_scores, 3))
scores_rounded$cluster <- rownames(scores_rounded)
prelim_df <- merge(prelim_df, scores_rounded, by = "cluster")

write_csv(prelim_df, file.path(TAB_DIR, "preliminary_cluster_annotation.csv"))
message("Saved: preliminary_cluster_annotation.csv")
print(prelim_df)

# =========================================================================
# MANUAL_REVIEW #2 — Cell-type annotation
# =========================================================================
# The manual annotation below is based on manual review of:
#   1. global_cluster_markers.csv (or global_cluster_markers_presto.csv)
#   2. dotplot_top5_per_cluster.pdf
#   3. dotplot_canonical_markers.pdf
#   4. featureplot_key_markers*.pdf
#   5. umap_by_cluster.pdf
#
# Each label is assigned by inspecting cluster-enriched genes against known
# mouse liver NPC markers.  Key marker evidence is noted inline.
# Labels marked "uncertain" require further validation (DotPlot / FeaturePlot
# re-check), and some may be revised after subclustering in Stage 3.
#
# NOTE: This is a manual first-pass annotation — NOT a computational
#       assignment.  Do NOT replace with automatic scoring output.
# =========================================================================

# --- Manual cluster annotation (resolution = 0.4, 27 clusters) ---
# Key marker rationale per cluster:
#
#  0  — T cell:
#       Lef1, Tcf7, Il7r, Trac, Cd28  => naive / conventional T cell
#  1  — B cell:
#       Pax5, Ms4a1, Cd79a, Cd79b, Fcer2a, Bank1
#  2  — Inflammatory neutrophil:
#       Cxcr2, Il1f9, Acod1, Ptgs2, Cxcl2, S100a8/S100a9 signal
#  3  — NK cell:
#       Ncr1, Klra9, Klrb1b, Gzma, Prf1, Gzmb, Xcl1
#  4  — Cholangiocyte / epithelial:
#       Epcam, Krt8, Krt18, Krt19, Sox9 / Dmbt1 (not hepatocyte)
#  5  — Endothelial:
#       Robo4, Lyve1, Mmrn2, Ptprb, Kdr, Cdh5, Stab2, Pecam1
#  6  — Hepatocyte contamination:
#       Alb, Apoa1, Ttr, Cyp genes, Hnf4a
#  7  — Low-quality / mixed:
#       mt-Rnr2, mt-Rnr1, Gm42418, Malat1, Cst3; NOT a clean biological type
#  8  — Monocyte-derived macrophage / DC-like:
#       Clec4b1, Ms4a6c, Cd209a, Clec4a genes, Cx3cr1, Lyz2, Fcgr1
#  9  — Mesothelial:
#       Msln, Muc16, Upk3b, Rspo1
# 10  — Cytotoxic T / NK-like:
#       Gzmk, Ifit1, Ifit3, Mx1, Cd8a, Cd8b1, Ccl5 (IFN-responsive)
# 11  — pDC-like:
#       Siglech, Ccr9, Flt3, Cd300c
# 12  — Cycling cell:
#       Mki67, Top2a, Stmn1, Hist1h3c, Pclaf, Ccna2
# 13  — B cell:
#       Pax5, Ms4a1, Cd79a, Cd79b, Fcer2a, Bank1 (similar to cluster 1)
# 14  — Kupffer cell:
#       Cd5l, Clec4f, Folr2, Fcna, Cd163, Vsig4, C6
# 15  — cDC1:
#       Xcr1, Clec9a, Flt3, Itgax (NOT monocyte)
# 16  — Plasma cell:
#       Jchain, Tnfrsf17, Igha, Ighg genes, Mzb1/Xbp1
# 17  — Activated neutrophil:
#       Il1r2, Trem1, Csf3r, Cxcr2, Il1b
# 18  — Neutrophil:
#       Ltf, Camp, Ngp, Cd177, Ly6g
# 19  — VSMC / pericyte:
#       Lmod1, Cnn1, Myh11, Tagln
# 20  — B cell / possible doublet:
#       Ms4a1, Pax5, Fcer2a, Ighd, Bank1, Cd79a + Alb/Apoa1/Ttr background
# 21  — Fibroblast / stellate:
#       Dpep1, Dpt, Gdf10, Mmp2, Tcf21, Pdgfra, Lum
# 22  — Mast / basophil-like:
#       Cpa3, Gata2, Kit, Prtn3 (NOT monocyte)
# 23  — Cholangiocyte / epithelial:
#       Epcam, Krt8, Krt18, Krt19, Sox9 (similar to cluster 4)
# 24  — Hepatocyte contamination:
#       Alb, Apoa1, Ttr, Cyp genes, Hnf4a (similar to cluster 6)
# 25  — Endothelial:
#       Thbd, Nos3, Lyve1, Stab2, Cdh5, Kdr, Pecam1 (NOT neutrophil)
# 26  — pDC-like:
#       Siglech, Ccr9, Flt3, Cd300c (similar to cluster 11)
#
# Clusters needing further review:
#   7  — Low-quality / mixed; check mt% and nFeature
#  10  — Cytotoxic T / NK-like; check Cd8 vs Ncr1 co-expression
#  20  — B cell / possible doublet; consider SoupX / CellBender ambient RNA removal
#  22  — Mast / basophil-like; confirm with Cpa3/Kit co-expression

manual_cluster_annotation <- c(
  "0"  = "T cell",
  "1"  = "B cell",
  "2"  = "Inflammatory neutrophil",
  "3"  = "NK cell",
  "4"  = "Cholangiocyte / epithelial",
  "5"  = "Endothelial",
  "6"  = "Hepatocyte contamination",
  "7"  = "Low-quality / mixed",
  "8"  = "Monocyte-derived macrophage / DC-like",
  "9"  = "Mesothelial",
  "10" = "Cytotoxic T / NK-like",
  "11" = "pDC-like",
  "12" = "Cycling cell",
  "13" = "B cell",
  "14" = "Kupffer cell",
  "15" = "cDC1",
  "16" = "Plasma cell",
  "17" = "Activated neutrophil",
  "18" = "Neutrophil",
  "19" = "VSMC / pericyte",
  "20" = "B cell / possible doublet",
  "21" = "Fibroblast / stellate",
  "22" = "Mast / basophil-like",
  "23" = "Cholangiocyte / epithelial",
  "24" = "Hepatocyte contamination",
  "25" = "Endothelial",
  "26" = "pDC-like"
)

# Apply annotation to Seurat object
# unname() required for Seurat v5 which validates cell-name overlap on named vectors
merged$celltype_manual <- unname(manual_cluster_annotation[as.character(merged$seurat_clusters)])
message("Applied manual annotation. Cell types: ",
        paste(unique(merged$celltype_manual), collapse = ", "))

# -------------------------------------------------------------------------
# E. Annotation figures
# -------------------------------------------------------------------------

# --- UMAP by sample ---
p_samp <- DimPlot(merged, reduction = "umap", group.by = "sample_id",
                  pt.size = 0.2, shuffle = TRUE) +
  ggtitle("UMAP by sample") + theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_by_sample.pdf"), width = 10, height = 8)
print(p_samp)
dev.off()
message("Saved: umap_by_sample.pdf")

# --- UMAP by group ---
p_grp <- DimPlot(merged, reduction = "umap", group.by = "group",
                 pt.size = 0.2, shuffle = TRUE) +
  ggtitle("UMAP by group (Control vs APAP)") + theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_by_group.pdf"), width = 9, height = 8)
print(p_grp)
dev.off()
message("Saved: umap_by_group.pdf")

# --- UMAP by annotated cell type ---
p_ct <- DimPlot(merged, reduction = "umap", group.by = "celltype_manual",
                label = TRUE, repel = TRUE, pt.size = 0.2, shuffle = TRUE) +
  ggtitle("UMAP by annotated cell type") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_by_celltype.pdf"), width = 12, height = 8)
print(p_ct)
dev.off()
message("Saved: umap_by_celltype.pdf")

# --- DotPlot of canonical markers ---
all_canonical <- unique(unlist(canonical_markers))
genes_avail <- intersect(all_canonical, rownames(merged))

p_dot_canon <- DotPlot(merged, features = genes_avail,
                       group.by = "seurat_clusters", assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical markers by cluster") +
  theme(axis.text.x = element_text(size = 7))

pdf(file.path(FIG_DIR, "dotplot_canonical_markers.pdf"),
    width = max(16, length(genes_avail) * 0.35), height = 7)
print(p_dot_canon)
dev.off()
message("Saved: dotplot_canonical_markers.pdf")

# --- FeaturePlot of key markers ---
# Select a representative subset that are likely to be present
key_markers_feat <- c(
  # T cell
  "Cd3d", "Lef1", "Tcf7", "Cd8a",
  # B cell
  "Ms4a1", "Pax5", "Cd79a", "Bank1",
  # NK
  "Ncr1", "Gzma", "Gzmb", "Xcl1",
  # Neutrophil
  "S100a8", "Cxcr2", "Ltf", "Camp",
  # Endothelial
  "Pecam1", "Kdr", "Lyve1", "Robo4",
  # Kupffer cell
  "Cd5l", "Clec4f", "Vsig4",
  # Monocyte-derived Mph/DC
  "Cx3cr1", "Fcgr1", "Clec4b1",
  # Cholangiocyte
  "Epcam", "Krt18", "Sox9",
  # Hepatocyte
  "Alb", "Hnf4a",
  # Fibroblast / Stellate
  "Dpep1", "Dpt", "Pdgfra",
  # Mesothelial
  "Msln", "Upk3b",
  # pDC
  "Siglech", "Ccr9",
  # cDC1
  "Xcr1", "Clec9a",
  # Plasma cell
  "Jchain", "Tnfrsf17",
  # Mast cell
  "Cpa3", "Kit",
  # VSMC / pericyte
  "Cnn1", "Myh11",
  # Cycling
  "Mki67", "Top2a"
)
key_markers_feat <- intersect(key_markers_feat, rownames(merged))

# Plot in batches of 6 to keep files manageable
batch_size <- 6
batches <- split(key_markers_feat, ceiling(seq_along(key_markers_feat) / batch_size))

for (bname in names(batches)) {
  genes_batch <- batches[[bname]]
  p <- FeaturePlot(merged, features = genes_batch, reduction = "umap",
                   pt.size = 0.2, ncol = min(3, length(genes_batch)),
                   order = TRUE, combine = TRUE) +
    plot_annotation(title = paste0("Key markers (batch ", bname, ")"))
  pdf(file.path(FIG_DIR, paste0("featureplot_key_markers_batch", bname, ".pdf")),
      width = 14, height = 4 * ceiling(length(genes_batch) / 3))
  print(p)
  dev.off()
}
message("Saved: featureplot_key_markers_batch*.pdf (", length(batches), " files)")

# -------------------------------------------------------------------------
# F. Cell proportion comparison (Control vs APAP)
# -------------------------------------------------------------------------

# --- Counts per sample per cell type ---
count_df <- merged@meta.data %>%
  count(sample_id, group, celltype_manual, name = "n_cells")

write_csv(count_df, file.path(TAB_DIR, "celltype_counts_by_sample.csv"))
message("Saved: celltype_counts_by_sample.csv")

# --- Proportions per sample ---
prop_df <- count_df %>%
  group_by(sample_id) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(prop_df, file.path(TAB_DIR, "celltype_proportions_by_sample.csv"))
message("Saved: celltype_proportions_by_sample.csv")

# --- Stacked barplot ---
p_stack <- ggplot(prop_df, aes(x = sample_id, y = proportion, fill = celltype_manual)) +
  geom_bar(stat = "identity", position = "fill", colour = "grey20", linewidth = 0.1) +
  facet_wrap(~ group, scales = "free_x") +
  labs(x = "Sample", y = "Proportion", fill = "Cell type",
       title = "Cell-type proportions by sample") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

pdf(file.path(FIG_DIR, "celltype_proportion_stacked.pdf"), width = 12, height = 7)
print(p_stack)
dev.off()
message("Saved: celltype_proportion_stacked.pdf")

# --- Group-level barplot ---
prop_group <- prop_df %>%
  group_by(group, celltype_manual) %>%
  summarise(
    mean_prop = mean(proportion),
    sd_prop   = sd(proportion),
    n_samples = n(),
    .groups   = "drop"
  )

p_group <- ggplot(prop_group, aes(x = celltype_manual, y = mean_prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7,
           colour = "grey20", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Cell type", y = "Mean proportion", fill = "Group",
       title = "Mean cell-type proportions: Control vs APAP") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(FIG_DIR, "celltype_proportion_grouped.pdf"), width = 12, height = 7)
print(p_group)
dev.off()
message("Saved: celltype_proportion_grouped.pdf")

# --- Exploratory Wilcoxon test per cell type ---
# WARNING: n = 3 per group — p-values are EXPLORATORY only, NOT confirmatory.
celltypes <- unique(prop_df$celltype_manual)
wilcox_results <- data.frame(
  celltype   = character(0),
  W_stat     = numeric(0),
  p_value    = numeric(0),
  mean_ctrl  = numeric(0),
  mean_apap  = numeric(0),
  direction  = character(0),
  stringsAsFactors = FALSE
)

for (ct in celltypes) {
  ctrl_vals <- prop_df$proportion[prop_df$celltype_manual == ct &
                                  prop_df$group == "Control"]
  apap_vals <- prop_df$proportion[prop_df$celltype_manual == ct &
                                  prop_df$group == "APAP"]

  if (length(ctrl_vals) < 2 || length(apap_vals) < 2) next

  wt <- wilcox.test(apap_vals, ctrl_vals, exact = FALSE)

  wilcox_results <- rbind(wilcox_results, data.frame(
    celltype   = ct,
    W_stat     = wt$statistic,
    p_value    = round(wt$p.value, 4),
    mean_ctrl  = round(mean(ctrl_vals), 4),
    mean_apap  = round(mean(apap_vals), 4),
    direction  = ifelse(mean(apap_vals) > mean(ctrl_vals), "Up in APAP", "Down in APAP"),
    stringsAsFactors = FALSE
  ))
}

write_csv(wilcox_results, file.path(TAB_DIR, "celltype_proportion_wilcox.csv"))
message("Saved: celltype_proportion_wilcox.csv (EXPLORATORY — n=3 per group)")
print(wilcox_results)

# -------------------------------------------------------------------------
# G. Myeloid assessment for stage-3 decision
# -------------------------------------------------------------------------

# Summary table: per cell type — cell count, group means, direction, myeloid relevance
summary_ct <- count_df %>%
  group_by(celltype_manual) %>%
  summarise(total_cells = sum(n_cells), .groups = "drop")

summary_prop <- prop_df %>%
  group_by(celltype_manual, group) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_prop, values_fill = 0)

myeloid_assessment <- summary_ct %>%
  left_join(summary_prop, by = "celltype_manual") %>%
  left_join(wilcox_results[, c("celltype", "direction", "p_value")],
            by = c("celltype_manual" = "celltype"))

# Flag myeloid-relevant cell types
myeloid_types <- c("Kupffer cell", "Monocyte-derived macrophage / DC-like",
                   "Inflammatory neutrophil", "Neutrophil", "Activated neutrophil")
myeloid_assessment$myeloid_relevant <- myeloid_assessment$celltype_manual %in% myeloid_types
myeloid_assessment$direction[is.na(myeloid_assessment$direction)] <- "N/A"
myeloid_assessment$p_value[is.na(myeloid_assessment$p_value)] <- NA_real_

write_csv(myeloid_assessment,
          file.path(TAB_DIR, "myeloid_assessment_for_stage3.csv"))
message("Saved: myeloid_assessment_for_stage3.csv")
print(myeloid_assessment)

# =========================================================================
# MANUAL_REVIEW #3 — Stage-3 direction decision
# =========================================================================
# Review myeloid_assessment_for_stage3.csv and the cell-type proportion plots.
# Key questions:
#   - Do macrophage/monocyte/neutrophil clusters have sufficient cells (>500)?
#   - Are myeloid marker genes clearly expressed and cluster-specific?
#   - Is there a visible APAP-induced shift in myeloid proportions?
#   - Are non-myeloid populations (endothelial, T/NK, B) showing stronger changes?
#
# Default recommendation: proceed to 03_Myeloid_Subclustering IF:
#   - Total myeloid cells > 1000
#   - Marker expression is clean
#   - There is a visible APAP effect in proportion or marker expression
#
# Otherwise, review which lineage shows the strongest APAP response and
# redirect stage 3 accordingly.
# =========================================================================

# -------------------------------------------------------------------------
# H. Save final object, sessionInfo
# -------------------------------------------------------------------------

saveRDS(merged, file.path(PROC_DIR, "seurat_annotated.rds"))
message("Saved: seurat_annotated.rds (", ncol(merged), " cells)")

writeLines(capture.output(sessionInfo()),
           file.path(TAB_DIR, "sessionInfo_02_Global_Annotation.txt"))
message("Saved: sessionInfo_02_Global_Annotation.txt")

message("=== 02_Global_Annotation.R complete ===")
message("Next steps:")
message("  1. Review manual annotation against DotPlot / FeaturePlot outputs")
message("  2. Re-check uncertain clusters (7, 10, 20, 22)")
message("  3. Review myeloid_assessment_for_stage3.csv for Stage-3 subclustering decisions")
message("  4. Proceed to 03_Myeloid_Subclustering.R when annotation is stable")
