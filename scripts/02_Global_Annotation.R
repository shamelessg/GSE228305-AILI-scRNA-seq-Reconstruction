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
# TEMPORARY DEFAULT: resolution = 0.4
# Review the UMAP plots for each resolution (umap_res_0.2.pdf, etc.)
# and the resolution_summary.csv table.  Choose the resolution that gives
# biologically interpretable clusters without over-fragmentation.
# Then update DEFAULT_RES below and re-run from this section.
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

message("Running FindAllMarkers (this may take several minutes for ~88k cells)...")

# Seurat v5: join sample-level layers before DE testing
merged <- JoinLayers(merged)

all_markers <- FindAllMarkers(
  merged,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25,
  test.use       = "wilcox"
)
gc()

# Full table
write_csv(all_markers, file.path(TAB_DIR, "global_cluster_markers.csv"))
message("Saved: global_cluster_markers.csv (", nrow(all_markers), " rows)")

if (nrow(all_markers) == 0) {
  stop("FindAllMarkers returned 0 rows. Check JoinLayers and assay setup.")
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
  "Macrophage / Kupffer" = c("Adgre1", "C1qa", "C1qb", "C1qc", "Lyz2",
                              "Cd68", "Marco", "Vsig4"),
  "Monocyte / Inflammatory Mph" = c("Ly6c2", "S100a8", "S100a9", "Ccr2",
                                     "Lyz2", "Itgam"),
  "Neutrophil"           = c("S100a8", "S100a9", "Retnlg", "Mpo", "Cxcr2"),
  "Endothelial"          = c("Pecam1", "Kdr", "Cdh5", "Vwf", "Lyve1", "Stab2"),
  "Stellate / Fibroblast" = c("Col1a1", "Col1a2", "Dcn", "Lum", "Acta2"),
  "T cell"               = c("Cd3d", "Cd3e", "Trac", "Cd4", "Cd8a"),
  "NK cell"              = c("Nkg7", "Klrb1c", "Gzmb", "Prf1"),
  "B cell"               = c("Ms4a1", "Cd79a", "Cd79b", "Bank1"),
  "Plasma cell"          = c("Jchain", "Mzb1", "Xbp1"),
  "Cycling"              = c("Mki67", "Top2a", "Stmn1"),
  "Hepatocyte"           = c("Alb", "Apoa1", "Ttr")
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
# The preliminary annotation above is COMPUTATIONAL ONLY.
# You MUST review it against:
#   1. global_cluster_markers.csv — full DEG table per cluster
#   2. dotplot_top5_per_cluster.pdf — expression pattern of top markers
#   3. dotplot_canonical_markers.pdf — expression of canonical markers
#   4. featureplot_key_markers.pdf — spatial distribution of markers on UMAP
#   5. umap_by_cluster.pdf — cluster positions on UMAP
#
# After review, fill in the manual_cluster_annotation vector below.
# Use "Unknown" for clusters you cannot confidently assign.
# =========================================================================

# --- MANUAL ANNOTATION TABLE (EDIT AFTER REVIEW) ---
# For now, the preliminary labels are used as a starting point.
# Replace the right-hand side with your final labels after reviewing outputs.
manual_cluster_annotation <- setNames(
  prelim_df$preliminary_label,
  prelim_df$cluster
)
# If you want to start from scratch, uncomment and edit:
# manual_cluster_annotation <- c(
#   "0"  = "Macrophage / Kupffer",
#   "1"  = "Endothelial",
#   "2"  = "T cell",
#   ...
# )

# Apply annotation to Seurat object
merged$celltype_manual <- manual_cluster_annotation[as.character(merged$seurat_clusters)]
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
  ggtitle("UMAP by annotated cell type (preliminary)") +
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
  # Macrophage
  "Adgre1", "C1qa", "Cd68",
  # Monocyte / inflammatory
  "Ly6c2", "Ccr2",
  # Neutrophil
  "S100a8", "Mpo",
  # Endothelial
  "Pecam1", "Kdr",
  # Stellate
  "Col1a1", "Dcn",
  # T cell
  "Cd3d", "Cd4",
  # NK
  "Nkg7", "Gzmb",
  # B cell
  "Ms4a1", "Cd79a",
  # Plasma
  "Jchain", "Mzb1",
  # Cycling
  "Mki67", "Top2a",
  # Hepatocyte
  "Alb", "Ttr"
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
myeloid_types <- c("Macrophage / Kupffer", "Monocyte / Inflammatory Mph",
                   "Neutrophil")
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
message("  1. Review resolution UMAPs and choose final resolution")
message("  2. Review marker tables and fill in manual_cluster_annotation")
message("  3. Review myeloid_assessment_for_stage3.csv to decide stage-3 direction")
