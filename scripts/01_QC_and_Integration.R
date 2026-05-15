# =========================================================================
# 01_QC_and_Integration.R
# =========================================================================
# Purpose: Read 6 mouse liver scRNA-seq samples (3 Control, 3 APAP),
#          perform QC filtering, normalize, and integrate with Harmony.
# Input:  data/external/GSM*/  (10X format: barcodes.tsv, genes.tsv, matrix.mtx)
# Output: data/processed/01_QC_and_Integration/seurat_integrated.rds
#         results/figures/01_QC_and_Integration/*.pdf
#         results/tables/01_QC_and_Integration/sample_qc_summary.csv
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
library(Matrix)
library(future)

# Use sequential plan to avoid macOS fork memory duplication on 16GB RAM
plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)

PROJ_ROOT <- getwd()
setwd(PROJ_ROOT)

DATA_DIR   <- file.path(PROJ_ROOT, "data", "external")
PROC_DIR   <- file.path(PROJ_ROOT, "data", "processed", "01_QC_and_Integration")
FIG_DIR    <- file.path(PROJ_ROOT, "results", "figures", "01_QC_and_Integration")
TAB_DIR    <- file.path(PROJ_ROOT, "results", "tables", "01_QC_and_Integration")

dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Sample metadata (explicit, not inferred from folder names)
# -------------------------------------------------------------------------

sample_meta <- data.frame(
  group     = c("Control", "Control", "Control", "APAP", "APAP", "APAP"),
  sample_id = c("con1", "con2", "con3", "AP-300", "AP-300-2", "AP-300-3"),
  input_dir = c(
    "GSM7118537_con1_matrix",
    "GSM7118537_con2_matrix",
    "GSM7118537_con3_matrix",
    "GSM7118538_AP-300_matrix",
    "GSM7118538_AP-300-2_matrix",
    "GSM7118538_AP-300-3_matrix"
  ),
  stringsAsFactors = FALSE
)
sample_meta$condition <- sample_meta$group  # Control vs APAP
rownames(sample_meta) <- sample_meta$sample_id

# -------------------------------------------------------------------------
# 2. Read each sample and build individual Seurat objects
# -------------------------------------------------------------------------

seurat_list <- list()

for (i in seq_len(nrow(sample_meta))) {
  sid    <- sample_meta$sample_id[i]
  grp    <- sample_meta$group[i]
  cond   <- sample_meta$condition[i]
  in_dir <- file.path(DATA_DIR, sample_meta$input_dir[i])

  if (!dir.exists(in_dir)) {
    stop("Input directory does not exist: ", in_dir)
  }

  message("Reading sample: ", sid, " from ", in_dir)
  counts <- Read10X(data.dir = in_dir)

  obj <- CreateSeuratObject(
    counts  = counts,
    project = sid,
    assay   = "RNA",
    min.cells = 3,
    min.features = 200
  )

  obj$sample_id  <- sid
  obj$group      <- grp
  obj$condition  <- cond
  obj$orig.ident <- sid

  seurat_list[[sid]] <- obj
  message("  -> ", ncol(obj), " cells after per-sample filtering")
}

# -------------------------------------------------------------------------
# 3. Merge all samples
# -------------------------------------------------------------------------

merged <- merge(
  seurat_list[[1]],
  y        = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

message("Merged object: ", ncol(merged), " cells x ", nrow(merged), " genes")

# -------------------------------------------------------------------------
# 4. Calculate percent.mt (mouse: ^mt-)
# -------------------------------------------------------------------------

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")

# -------------------------------------------------------------------------
# 5. QC metrics — before filtering
#    Inspect these plots to choose thresholds; do not blindly apply defaults.
# -------------------------------------------------------------------------

# Violin of key metrics
p1 <- VlnPlot(merged, features = c("nFeature_RNA"), group.by = "sample_id", pt.size = 0.1) +
  NoLegend() + ggtitle("Genes per cell (before filter)")

p2 <- VlnPlot(merged, features = c("nCount_RNA"), group.by = "sample_id", pt.size = 0.1) +
  NoLegend() + ggtitle("UMIs per cell (before filter)")

p3 <- VlnPlot(merged, features = c("percent.mt"), group.by = "sample_id", pt.size = 0.1) +
  NoLegend() + ggtitle("% MT (before filter)")

qc_violin_before <- wrap_plots(p1, p2, p3, ncol = 1)

pdf(file.path(FIG_DIR, "qc_violin_before_filter.pdf"), width = 10, height = 14)
print(qc_violin_before)
dev.off()
message("Saved: qc_violin_before_filter.pdf")

# Feature scatter for nCount vs nFeature vs percent.mt (by sample)
p4 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "sample_id", pt.size = 0.3) +
  ggtitle("nCount vs nFeature (before filter)")

p5 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt",
                     group.by = "sample_id", pt.size = 0.3) +
  ggtitle("nCount vs %MT (before filter)")

pdf(file.path(FIG_DIR, "qc_scatter_before_filter.pdf"), width = 12, height = 6)
print(p4 + p5)
dev.off()

# -------------------------------------------------------------------------
# 6. QC thresholds
#    Adjust manually after inspecting the pre-filter violin/scatter plots.
#    - nFeature_RNA: low = potential empty/dying cells; high = potential doublets
#    - nCount_RNA:   low = poor quality; high = potential doublets
#    - percent.mt:   high = dying/damaged cells (mouse liver: usually <20-30%)
#
#    Initial exploratory thresholds (conservative, keep broad):
# -------------------------------------------------------------------------

QC_MIN_FEATURES <- 300
QC_MAX_FEATURES <- 6000
QC_MAX_COUNT    <- 40000
QC_MAX_PCT_MT   <- 20

# Record pre-filter cell counts
pre_qc_n <- table(merged$sample_id)

# Apply filters
merged <- subset(merged,
  subset = nFeature_RNA > QC_MIN_FEATURES &
           nFeature_RNA < QC_MAX_FEATURES &
           nCount_RNA   < QC_MAX_COUNT    &
           percent.mt   < QC_MAX_PCT_MT
)

message("After QC filter: ", ncol(merged), " cells (removed ",
        sum(pre_qc_n) - ncol(merged), " cells)")

# -------------------------------------------------------------------------
# 7. QC violin — after filtering
# -------------------------------------------------------------------------

p6 <- VlnPlot(merged, features = c("nFeature_RNA"), group.by = "sample_id", pt.size = 0.1) +
  NoLegend() + ggtitle("Genes per cell (after filter)")

p7 <- VlnPlot(merged, features = c("nCount_RNA"), group.by = "sample_id", pt.size = 0.1) +
  NoLegend() + ggtitle("UMIs per cell (after filter)")

p8 <- VlnPlot(merged, features = c("percent.mt"), group.by = "sample_id", pt.size = 0.1) +
  NoLegend() + ggtitle("% MT (after filter)")

qc_violin_after <- wrap_plots(p6, p7, p8, ncol = 1)

pdf(file.path(FIG_DIR, "qc_violin_after_filter.pdf"), width = 10, height = 14)
print(qc_violin_after)
dev.off()
message("Saved: qc_violin_after_filter.pdf")

# -------------------------------------------------------------------------
# 8. QC summary table
# -------------------------------------------------------------------------

post_qc_n <- table(merged$sample_id)

qc_summary <- data.frame(
  sample_id      = names(pre_qc_n),
  group          = sample_meta[names(pre_qc_n), "group"],
  cells_before_QC = as.integer(pre_qc_n),
  cells_after_QC  = as.integer(post_qc_n[names(pre_qc_n)]),
  cells_removed   = as.integer(pre_qc_n - post_qc_n[names(pre_qc_n)]),
  pct_removed     = round((1 - as.numeric(post_qc_n[names(pre_qc_n)]) /
                                as.numeric(pre_qc_n)) * 100, 2),
  row.names       = NULL
)

write_csv(qc_summary, file.path(TAB_DIR, "sample_qc_summary.csv"))
message("Saved: sample_qc_summary.csv")
print(qc_summary)

# -------------------------------------------------------------------------
# 9. Normalization, variable features, scaling, PCA
# -------------------------------------------------------------------------

merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 10000)

merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = 2000)

# Scale only variable features to reduce memory usage; regress out percent.mt.
# For the initial run we regress only percent.mt to reduce technical noise
# without over-correcting biology.
merged <- ScaleData(merged, features = VariableFeatures(merged),
                    vars.to.regress = "percent.mt")
gc()

merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = 50)
gc()

# Elbow plot to guide dimensionality choice
pdf(file.path(FIG_DIR, "pca_elbow.pdf"), width = 8, height = 5)
ElbowPlot(merged, ndims = 50) + ggtitle("PCA Elbow Plot")
dev.off()

# -------------------------------------------------------------------------
# 10. UMAP before Harmony
#     Use 30 PCs as a reasonable default; inspect elbow plot to adjust.
# -------------------------------------------------------------------------

NPCS <- 30

merged <- RunUMAP(merged, dims = 1:NPCS, reduction = "pca")
gc()

# By sample, shuffled to avoid overplotting bias
p9 <- DimPlot(merged, reduction = "umap", group.by = "sample_id",
              pt.size = 0.3, shuffle = TRUE) +
  ggtitle("UMAP before Harmony (by sample)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_before_harmony_by_sample.pdf"), width = 10, height = 8)
print(p9)
dev.off()
message("Saved: umap_before_harmony_by_sample.pdf")

# Facet highlight: one panel per sample, that sample's cells in color, rest in grey
p9_facet <- DimPlot(merged, reduction = "umap", group.by = "sample_id",
                    split.by = "sample_id", pt.size = 0.1, ncol = 3, shuffle = TRUE) +
  ggtitle("UMAP before Harmony (per-sample highlight)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_before_harmony_by_sample_facet.pdf"), width = 18, height = 10)
print(p9_facet)
dev.off()

# By group, shuffled
p10 <- DimPlot(merged, reduction = "umap", group.by = "group",
               pt.size = 0.3, shuffle = TRUE) +
  ggtitle("UMAP before Harmony (by group)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_before_harmony_by_group.pdf"), width = 9, height = 8)
print(p10)
dev.off()

# Facet by group
p10_facet <- DimPlot(merged, reduction = "umap", group.by = "group",
                     split.by = "group", pt.size = 0.1, ncol = 2, shuffle = TRUE) +
  ggtitle("UMAP before Harmony (per-group highlight)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_before_harmony_by_group_facet.pdf"), width = 14, height = 7)
print(p10_facet)
dev.off()

# -------------------------------------------------------------------------
# 11. Harmony integration
#     Harmony should reduce sample-level technical variation while preserving
#     genuine biological differences between Control and APAP.
#     After running, inspect UMAP: samples should mix within the same condition,
#     but Control vs APAP should retain some separation if biology is real.
# -------------------------------------------------------------------------

message("Starting Harmony integration (this may be memory-intensive)...")
merged <- RunHarmony(
  merged,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use      = 1:NPCS,
  theta         = 2,         # Higher theta = stronger correction; 2 is default
  plot_convergence = FALSE
)
gc()
message("Harmony complete.")

# -------------------------------------------------------------------------
# 12. UMAP after Harmony
# -------------------------------------------------------------------------

merged <- RunUMAP(merged, reduction = "harmony", dims = 1:NPCS)
gc()

p11 <- DimPlot(merged, reduction = "umap", group.by = "sample_id",
               pt.size = 0.3, shuffle = TRUE) +
  ggtitle("UMAP after Harmony (by sample)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_after_harmony_by_sample.pdf"), width = 10, height = 8)
print(p11)
dev.off()
message("Saved: umap_after_harmony_by_sample.pdf")

# Facet highlight: one panel per sample
p11_facet <- DimPlot(merged, reduction = "umap", group.by = "sample_id",
                     split.by = "sample_id", pt.size = 0.1, ncol = 3, shuffle = TRUE) +
  ggtitle("UMAP after Harmony (per-sample highlight)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_after_harmony_by_sample_facet.pdf"), width = 18, height = 10)
print(p11_facet)
dev.off()

# By group, shuffled
p12 <- DimPlot(merged, reduction = "umap", group.by = "group",
               pt.size = 0.3, shuffle = TRUE) +
  ggtitle("UMAP after Harmony (by group)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_after_harmony_by_group.pdf"), width = 9, height = 8)
print(p12)
dev.off()

# Facet by group
p12_facet <- DimPlot(merged, reduction = "umap", group.by = "group",
                     split.by = "group", pt.size = 0.1, ncol = 2, shuffle = TRUE) +
  ggtitle("UMAP after Harmony (per-group highlight)") +
  theme(aspect.ratio = 1)

pdf(file.path(FIG_DIR, "umap_after_harmony_by_group_facet.pdf"), width = 14, height = 7)
print(p12_facet)
dev.off()

# -------------------------------------------------------------------------
# 13. Save integrated object
# -------------------------------------------------------------------------

saveRDS(merged, file.path(PROC_DIR, "seurat_integrated.rds"))
message("Saved: seurat_integrated.rds")
message("=== 01_QC_and_Integration.R complete ===\n")
