# =========================================================================
# 05_Cell_Communication.R
# =========================================================================
# Purpose: Targeted exploratory cell-cell communication analysis.
#          Uses CellChatDB.mouse ligand-receptor pairs as reference, but
#          computes interaction scores manually (avg expr L × avg expr R)
#          to avoid CellChat full-pipeline fragility with small cell counts.
#          This is NOT a full CellChat discovery analysis.
#          n=3 Control vs 3 APAP — results are weak supplementary evidence.
#
# Targeted axes (Stage 4 informed):
#   1. Neutrophil / Activated neutrophil -> Endothelial
#   2. Kupffer / macrophage -> neutrophil or endothelial
#   3. Antigen-presenting myeloid -> T/NK
#   4. Fibroblast / stellate -> endothelial/myeloid
#
# Strategy:
#   - Subset seurat_global to focus celltypes only
#   - Per group x celltype, max 500 cells (memory & disk control)
#   - Compute mean log-norm expression per gene per celltype per group
#   - Score targeted LR pairs as mean(L_src) * mean(R_tgt) for each group
#   - Compare APAP vs Control by delta score
#   - No full CellChat objects saved
#
# Input:  data/processed/04_APAP_Control_Functional_Comparison/seurat_functional_comparison.rds
# Output: data/processed/05_Cell_Communication/cellchat_targeted_summary.rds
#         results/tables/05_Cell_Communication/*.csv
#         results/figures/05_Cell_Communication/*.pdf
# =========================================================================

# -------------------------------------------------------------------------
# 0. Setup
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
})

PROJ_ROOT <- getwd()
setwd(PROJ_ROOT)

# --- directories ---
PROC_DIR <- file.path(PROJ_ROOT, "data", "processed", "05_Cell_Communication")
FIG_DIR  <- file.path(PROJ_ROOT, "results", "figures",   "05_Cell_Communication")
TAB_DIR  <- file.path(PROJ_ROOT, "results", "tables",    "05_Cell_Communication")

dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Human-editable parameters
# -------------------------------------------------------------------------

TARGET_CELLTYPES <- c(
  "Inflammatory neutrophil",
  "Activated neutrophil",
  "Neutrophil",
  "Kupffer cell",
  "Monocyte-derived macrophage / DC-like",
  "Endothelial",
  "T cell",
  "NK cell",
  "Fibroblast / stellate"
)

MAX_CELLS_PER_GROUP_CELLTYPE <- 500L

# Targeted source -> target axes from Stage 4 interpretation
TARGET_SOURCE_TARGET_PAIRS <- list(
  c("Inflammatory neutrophil", "Endothelial"),
  c("Activated neutrophil",    "Endothelial"),
  c("Neutrophil",              "Endothelial"),
  c("Kupffer cell",            "Neutrophil"),
  c("Kupffer cell",            "Inflammatory neutrophil"),
  c("Kupffer cell",            "Endothelial"),
  c("Monocyte-derived macrophage / DC-like", "T cell"),
  c("Monocyte-derived macrophage / DC-like", "NK cell"),
  c("Fibroblast / stellate",   "Endothelial"),
  c("Fibroblast / stellate",   "Kupffer cell"),
  c("Fibroblast / stellate",   "Monocyte-derived macrophage / DC-like")
)

# Only report LR pairs where at least one group has score > this
MIN_INTERACTION_SCORE <- 0.01

# -------------------------------------------------------------------------
# 2. Load and subset data
# -------------------------------------------------------------------------

message("\n=== Loading Stage 4 data ===")
x <- readRDS(file.path(PROJ_ROOT, "data", "processed",
                       "04_APAP_Control_Functional_Comparison",
                       "seurat_functional_comparison.rds"))
seu <- x$seurat_global
rm(x); gc()

message("Full seurat_global: ", ncol(seu), " cells x ", nrow(seu), " genes")

# Keep only target celltypes
cells_keep <- colnames(seu)[seu$celltype_manual %in% TARGET_CELLTYPES]
seu <- subset(seu, cells = cells_keep)
message("After celltype filter: ", ncol(seu), " cells")

# Keep only Control and APAP
seu <- subset(seu, subset = group %in% c("Control", "APAP"))
message("After group filter: ", ncol(seu), " cells")

# -------------------------------------------------------------------------
# 3. Downsample per group x celltype
# -------------------------------------------------------------------------

set.seed(20250515)
meta <- seu@meta.data
meta$cell_barcode <- rownames(meta)

cells_keep <- meta %>%
  group_by(group, celltype_manual) %>%
  slice_sample(n = MAX_CELLS_PER_GROUP_CELLTYPE) %>%
  pull(cell_barcode)

seu <- subset(seu, cells = cells_keep)
message("After downsampling: ", ncol(seu), " cells")

ds_summary <- meta %>%
  filter(cell_barcode %in% cells_keep) %>%
  count(group, celltype_manual) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0L)
message("Downsampled composition:\n",
        paste(capture.output(print(as.data.frame(ds_summary))), collapse = "\n"))

# -------------------------------------------------------------------------
# 4. Load CellChatDB.mouse to get ligand-receptor pairs
# -------------------------------------------------------------------------

message("\n=== Loading CellChatDB.mouse ===")
data("CellChatDB.mouse", package = "CellChat")
ccdb <- CellChatDB.mouse
message("interaction entries: ", nrow(ccdb$interaction))

# Build a clean LR table with unique ligand-receptor pairs per pathway
lr_db <- ccdb$interaction %>%
  select(ligand, receptor, pathway_name, annotation) %>%
  distinct()

message("Unique LR pairs in DB: ", nrow(lr_db))

# Also get gene info for symbol validation
db_symbols <- unique(c(lr_db$ligand, lr_db$receptor))

# -------------------------------------------------------------------------
# 5. Compute mean expression per gene per celltype per group
# -------------------------------------------------------------------------

message("\n=== Computing mean expression per celltype x group ===")

seu <- JoinLayers(seu)

# Get log-normalized data
expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
meta_df  <- seu@meta.data

# Only keep genes that are in the CellChatDB
genes_use <- intersect(rownames(expr_mat), db_symbols)
message("Genes overlapping with CellChatDB: ", length(genes_use),
        " / ", length(db_symbols))
expr_mat <- expr_mat[genes_use, , drop = FALSE]

# Compute mean expression per celltype x group
celltype_groups <- interaction(meta_df$celltype_manual, meta_df$group, sep = "||")
mean_expr_list <- list()

for (cg in levels(celltype_groups)) {
  cells_in_cg <- rownames(meta_df)[celltype_groups == cg]
  if (length(cells_in_cg) == 0) next
  parts <- strsplit(cg, "\\|\\|")[[1]]
  ct  <- parts[1]
  grp <- parts[2]

  mean_vec <- Matrix::rowMeans(expr_mat[, cells_in_cg, drop = FALSE])
  mean_expr_list[[cg]] <- data.frame(
    gene      = names(mean_vec),
    celltype  = ct,
    group     = grp,
    mean_expr = as.numeric(mean_vec),
    n_cells   = length(cells_in_cg),
    stringsAsFactors = FALSE
  )
}

mean_expr_df <- bind_rows(mean_expr_list)
message("Mean expression table: ", nrow(mean_expr_df), " rows")

# -------------------------------------------------------------------------
# 6. Score targeted LR interactions
# -------------------------------------------------------------------------

message("\n=== Scoring targeted LR interactions ===")

score_lr_for_group <- function(mean_df, lr_db, grp_name) {
  # Build lookup: celltype x gene -> mean_expr
  grp_mean <- mean_df[mean_df$group == grp_name, ]

  results <- list()
  for (pair in TARGET_SOURCE_TARGET_PAIRS) {
    src <- pair[1]; tgt <- pair[2]

    # Get relevant LR pairs for this source-target
    # We score ALL known LR pairs, not just a subset
    for (i in seq_len(nrow(lr_db))) {
      L <- lr_db$ligand[i]
      R <- lr_db$receptor[i]

      expr_L <- grp_mean$mean_expr[grp_mean$gene == L & grp_mean$celltype == src]
      expr_R <- grp_mean$mean_expr[grp_mean$gene == R & grp_mean$celltype == tgt]

      if (length(expr_L) == 0 || length(expr_R) == 0) next
      if (expr_L == 0 && expr_R == 0) next

      score <- expr_L * expr_R
      if (is.na(score)) next

      results[[length(results) + 1]] <- data.frame(
        source       = src,
        target       = tgt,
        ligand       = L,
        receptor     = R,
        pathway_name = lr_db$pathway_name[i],
        annotation   = lr_db$annotation[i],
        expr_ligand  = expr_L,
        expr_receptor = expr_R,
        interaction_score = score,
        group        = grp_name,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- bind_rows(results)
  if (nrow(out) == 0) return(out)
  out <- out[order(out$interaction_score, decreasing = TRUE), ]
  out
}

message("Scoring Control...")
lr_control <- score_lr_for_group(mean_expr_df, lr_db, "Control")
message("  ", nrow(lr_control), " LR pairs scored")

message("Scoring APAP...")
lr_apap <- score_lr_for_group(mean_expr_df, lr_db, "APAP")
message("  ", nrow(lr_apap), " LR pairs scored")

# -------------------------------------------------------------------------
# 7. APAP vs Control comparison
# -------------------------------------------------------------------------

message("\n=== Comparing APAP vs Control ===")

# Merge on L+R+source+target
lr_compare <- merge(
  lr_control %>% rename(score_control = interaction_score,
                        exprL_control = expr_ligand,
                        exprR_control = expr_receptor) %>%
    select(source, target, ligand, receptor, pathway_name, annotation,
           score_control, exprL_control, exprR_control),

  lr_apap %>% rename(score_apap = interaction_score,
                     exprL_apap = expr_ligand,
                     exprR_apap = expr_receptor) %>%
    select(source, target, ligand, receptor, score_apap, exprL_apap, exprR_apap),

  by = c("source", "target", "ligand", "receptor"),
  all = TRUE
)

# Fill missing
lr_compare$score_control[is.na(lr_compare$score_control)] <- 0
lr_compare$score_apap[is.na(lr_compare$score_apap)]       <- 0
lr_compare$exprL_control[is.na(lr_compare$exprL_control)] <- 0
lr_compare$exprR_control[is.na(lr_compare$exprR_control)] <- 0
lr_compare$exprL_apap[is.na(lr_compare$exprL_apap)]       <- 0
lr_compare$exprR_apap[is.na(lr_compare$exprR_apap)]       <- 0

lr_compare$delta_score <- lr_compare$score_apap - lr_compare$score_control
lr_compare$max_score   <- pmax(lr_compare$score_control, lr_compare$score_apap)
lr_compare$axis <- paste(lr_compare$source, "->", lr_compare$target)

# Filter: at least one group has score > threshold
lr_compare <- lr_compare[lr_compare$max_score > MIN_INTERACTION_SCORE, ]
lr_compare <- lr_compare[order(abs(lr_compare$delta_score), decreasing = TRUE), ]

message("LR pairs passing score > ", MIN_INTERACTION_SCORE, ": ", nrow(lr_compare))
message("Top 10 by |delta|:")
top10 <- head(lr_compare[, c("axis", "ligand", "receptor", "score_control",
                              "score_apap", "delta_score")], 10)
print(top10, row.names = FALSE)

# -------------------------------------------------------------------------
# 8. Pathway-level summary across targeted axes
# -------------------------------------------------------------------------

pathway_summary <- lr_compare %>%
  group_by(axis, pathway_name) %>%
  summarise(
    n_lr_pairs        = n(),
    mean_score_control = mean(score_control, na.rm = TRUE),
    mean_score_apap    = mean(score_apap, na.rm = TRUE),
    mean_delta         = mean(delta_score, na.rm = TRUE),
    max_delta          = max(delta_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(mean_delta)))

# -------------------------------------------------------------------------
# 9. Decision summary
# -------------------------------------------------------------------------

message("\n=== Building decision summary ===")

# Count meaningful signals
strong_delta <- lr_compare[abs(lr_compare$delta_score) > 0.1, ]
n_strong <- nrow(strong_delta)

decision_notes <- c(
  "Stage 5 is a targeted exploratory analysis, not a full CellChat discovery.",
  "Method: manual LR scoring = mean(logNormExpr[L] in source) * mean(logNormExpr[R] in target)",
  "LR database: CellChatDB.mouse (Secreted Signaling + Cell-Cell Contact)",
  paste("Number of targeted source->target axes:", length(TARGET_SOURCE_TARGET_PAIRS)),
  paste("Max cells per group x celltype:", MAX_CELLS_PER_GROUP_CELLTYPE),
  paste("Sample size: 3 Control vs 3 APAP"),
  "All results are weak supplementary evidence only.",
  "No full CellChat objects are saved.",
  paste("Genes overlapping with CellChatDB:", length(genes_use)),
  paste("LR pairs scored (Control):", nrow(lr_control)),
  paste("LR pairs scored (APAP):", nrow(lr_apap)),
  paste("LR pairs with |delta| > 0.1:", n_strong)
)

if (n_strong > 0) {
  decision_notes <- c(decision_notes,
    paste(n_strong, "LR pairs show noticeable APAP vs Control difference (|delta| > 0.1)."),
    "Due to n=3 and manual scoring, these are hypotheses, not findings.",
    "Top altered axes:",
    paste("  ", head(strong_delta$axis, 5), "|",
          head(strong_delta$ligand, 5), "-", head(strong_delta$receptor, 5),
          " delta:", round(head(strong_delta$delta_score, 5), 3))
  )
} else {
  decision_notes <- c(decision_notes,
    "NO LR pairs show |delta| > 0.1 between APAP and Control in targeted axes.",
    "No signal detected — results stop at supplementary-table level."
  )
}

decision_summary <- data.frame(
  parameter = c(
    "analysis_type", "scoring_method", "lr_database",
    "targeted_axes_n", "max_cells_per_group_celltype",
    "n_control", "n_apap", "genes_overlap_db",
    "n_lr_pairs_control", "n_lr_pairs_apap",
    "n_lr_pairs_compared", "n_strong_delta",
    "files_written"
  ),
  value = c(
    "targeted exploratory (manual LR scoring)",
    "mean(L_expr_in_src) * mean(R_expr_in_tgt)",
    "CellChatDB.mouse",
    as.character(length(TARGET_SOURCE_TARGET_PAIRS)),
    as.character(MAX_CELLS_PER_GROUP_CELLTYPE),
    "3", "3",
    as.character(length(genes_use)),
    as.character(nrow(lr_control)),
    as.character(nrow(lr_apap)),
    as.character(nrow(lr_compare)),
    as.character(n_strong),
    paste(
      "targeted_lr_pairs_control.csv",
      "targeted_lr_pairs_apap.csv",
      "targeted_lr_pairs_apap_vs_control.csv",
      "targeted_pathway_summary.csv",
      "stage5_decision_summary.csv",
      "targeted_bubble_control.pdf",
      "targeted_bubble_apap.pdf",
      "targeted_interaction_delta_summary.pdf",
      sep = "; "
    )
  ),
  stringsAsFactors = FALSE
)

decision_summary <- rbind(
  decision_summary,
  data.frame(parameter = paste0("note_", seq_along(decision_notes)),
             value = decision_notes,
             stringsAsFactors = FALSE)
)

# -------------------------------------------------------------------------
# 10. Figures
# -------------------------------------------------------------------------

message("\n=== Generating figures ===")

group_cols <- c("Control" = "#377eb8", "APAP" = "#e41a1c")

# 10a. Bubble plot: targeted axes, top LR pairs per group
make_bubble <- function(lr_df, grp_title, score_col = "interaction_score",
                        max_n = 30) {
  if (nrow(lr_df) == 0) return(NULL)

  df <- lr_df[order(lr_df[[score_col]], decreasing = TRUE), ]
  if (nrow(df) > max_n) df <- df[1:max_n, ]
  df$axis <- paste(df$source, "->", df$target)
  df$lr_label <- paste0(df$ligand, " - ", df$receptor)

  ggplot(df, aes(x = lr_label, y = axis, size = .data[[score_col]])) +
    geom_point(color = if(grp_title == "Control") "#377eb8" else "#e41a1c",
               alpha = 0.7) +
    scale_size_continuous(name = "score", range = c(1, 6)) +
    labs(title = paste("Targeted LR Scores —", grp_title),
         x = "Ligand - Receptor", y = "Source -> Target") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}

fig_control <- make_bubble(lr_control, "Control")
fig_apap    <- make_bubble(lr_apap, "APAP")

pdf(file.path(FIG_DIR, "targeted_bubble_control.pdf"), width = 12, height = 7)
if (!is.null(fig_control)) print(fig_control) else
  { plot.new(); text(0.5, 0.5, "Control: no targeted LR pairs pass score threshold") }
dev.off()

pdf(file.path(FIG_DIR, "targeted_bubble_apap.pdf"), width = 12, height = 7)
if (!is.null(fig_apap)) print(fig_apap) else
  { plot.new(); text(0.5, 0.5, "APAP: no targeted LR pairs pass score threshold") }
dev.off()

# 10b. Delta heatmap
pdf(file.path(FIG_DIR, "targeted_interaction_delta_summary.pdf"), width = 12, height = 8)
if (nrow(lr_compare) > 0) {
  # Take top 50 by |delta|
  heat_df <- lr_compare[order(abs(lr_compare$delta_score), decreasing = TRUE), ]
  if (nrow(heat_df) > 50) heat_df <- heat_df[1:50, ]

  mat <- matrix(heat_df$delta_score, ncol = 1)
  rownames(mat) <- paste0(heat_df$ligand, "|", heat_df$receptor, " (",
                          heat_df$source, "->", heat_df$target, ")")
  colnames(mat) <- "APAP vs Control"
  mat[mat > 0.5]  <- 0.5
  mat[mat < -0.5] <- -0.5

  pheatmap(mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("#377eb8", "white", "#e41a1c"))(100),
           main = "Targeted LR Score Delta (APAP - Control)",
           fontsize_row = 7, fontsize = 9,
           display_numbers = FALSE)
} else {
  plot.new()
  text(0.5, 0.5, "No LR pairs pass score threshold in either group.")
}
dev.off()

# -------------------------------------------------------------------------
# 11. Write output tables
# -------------------------------------------------------------------------

message("\n=== Writing output tables ===")

write.csv(lr_control,
          file.path(TAB_DIR, "targeted_lr_pairs_control.csv"),
          row.names = FALSE)
write.csv(lr_apap,
          file.path(TAB_DIR, "targeted_lr_pairs_apap.csv"),
          row.names = FALSE)
write.csv(lr_compare,
          file.path(TAB_DIR, "targeted_lr_pairs_apap_vs_control.csv"),
          row.names = FALSE)
write.csv(pathway_summary,
          file.path(TAB_DIR, "targeted_pathway_summary.csv"),
          row.names = FALSE)
write.csv(decision_summary,
          file.path(TAB_DIR, "stage5_decision_summary.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 12. Save lightweight summary RDS (tables + params only, NO full objects)
# -------------------------------------------------------------------------

message("\n=== Saving summary RDS ===")

cellchat_targeted_summary <- list(
  targeted_lr_control  = lr_control,
  targeted_lr_apap     = lr_apap,
  lr_compare           = lr_compare,
  pathway_summary      = pathway_summary,
  decision_summary     = decision_summary,

  params = list(
    target_celltypes              = TARGET_CELLTYPES,
    max_cells_per_group_celltype  = MAX_CELLS_PER_GROUP_CELLTYPE,
    target_source_target_pairs    = TARGET_SOURCE_TARGET_PAIRS,
    min_interaction_score         = MIN_INTERACTION_SCORE,
    genes_overlap_db              = length(genes_use),
    n_cells_after_downsample      = ncol(seu),
    lr_database                   = "CellChatDB.mouse",
    scoring_method                = "mean(L_expr) * mean(R_expr)",
    n_control_samples             = 3L,
    n_apap_samples                = 3L
  ),

  session_info = sessionInfo()
)

saveRDS(cellchat_targeted_summary,
        file.path(PROC_DIR, "cellchat_targeted_summary.rds"))

message("\n=== Stage 5 complete ===")
message("Output directory: ", PROC_DIR)
message("Figure directory:  ", FIG_DIR)
message("Table directory:   ", TAB_DIR)
