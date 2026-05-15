# =========================================================================
# 04_APAP_Control_Functional_Comparison.R
# =========================================================================
# Purpose: Integrate Stage 2 (global) and Stage 3 (myeloid subcluster) results
#          to answer: which non-parenchymal or myeloid subsets show composition
#          and functional state changes in APAP-injured liver?  Are these
#          changes more consistent with damage amplification, inflammatory
#          recruitment, antigen presentation, loss of homeostasis, or
#          repair/clearance responses?
#
# Four modules:
#   1. Unified composition delta summary (global + myeloid)
#   2. Functional gene-set module scores (global + myeloid)
#   3. Sample-level pseudobulk differential expression (conservative)
#   4. Interpretation summary table for README
#
# Input:  data/processed/02_Global_Annotation/seurat_annotated.rds
#         data/processed/03_Myeloid_Subclustering/myeloid_subclustered.rds
#         results/tables/02_Global_Annotation/celltype_proportions_by_sample.csv
#         results/tables/02_Global_Annotation/celltype_proportion_wilcox.csv
#         results/tables/03_Myeloid_Subclustering/myeloid_subcluster_proportions_by_sample.csv
#         results/tables/03_Myeloid_Subclustering/myeloid_subcluster_proportion_wilcox.csv
#
# Output: data/processed/04_APAP_Control_Functional_Comparison/seurat_functional_comparison.rds
#         results/figures/04_APAP_Control_Functional_Comparison/*.pdf
#         results/tables/04_APAP_Control_Functional_Comparison/*.csv
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
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(edgeR)

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
PROC_DIR <- file.path(PROJ_ROOT, "data", "processed", "04_APAP_Control_Functional_Comparison")
FIG_DIR  <- file.path(PROJ_ROOT, "results", "figures",   "04_APAP_Control_Functional_Comparison")
TAB_DIR  <- file.path(PROJ_ROOT, "results", "tables",    "04_APAP_Control_Functional_Comparison")

dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

# --- colour palettes ---
group_cols <- c("Control" = "#377eb8", "APAP" = "#e41a1c")

# -------------------------------------------------------------------------
# 0a. Input file validation
# -------------------------------------------------------------------------

# RDS inputs (required)
STAGE2_RDS  <- file.path(PROJ_ROOT, "data", "processed", "02_Global_Annotation", "seurat_annotated.rds")
STAGE3_RDS  <- file.path(PROJ_ROOT, "data", "processed", "03_Myeloid_Subclustering", "myeloid_subclustered.rds")

# Table inputs (from Stage 2 & 3)
GLOB_PROP_CSV  <- file.path(PROJ_ROOT, "results", "tables", "02_Global_Annotation", "celltype_proportions_by_sample.csv")
GLOB_WILCOX_CSV <- file.path(PROJ_ROOT, "results", "tables", "02_Global_Annotation", "celltype_proportion_wilcox.csv")
MYE_PROP_CSV   <- file.path(PROJ_ROOT, "results", "tables", "03_Myeloid_Subclustering", "myeloid_subcluster_proportions_by_sample.csv")
MYE_WILCOX_CSV <- file.path(PROJ_ROOT, "results", "tables", "03_Myeloid_Subclustering", "myeloid_subcluster_proportion_wilcox.csv")
MYE_SCORE_CSV  <- file.path(PROJ_ROOT, "results", "tables", "03_Myeloid_Subclustering", "myeloid_module_scores_aggregated.csv")

required_files <- c(
  "Stage2 RDS"                   = STAGE2_RDS,
  "Stage3 RDS"                   = STAGE3_RDS,
  "Global proportions"           = GLOB_PROP_CSV,
  "Global Wilcoxon"              = GLOB_WILCOX_CSV,
  "Myeloid proportions"          = MYE_PROP_CSV,
  "Myeloid Wilcoxon"             = MYE_WILCOX_CSV,
  "Myeloid module scores"        = MYE_SCORE_CSV
)

missing_files <- c()
for (nm in names(required_files)) {
  fpath <- required_files[nm]
  if (!file.exists(fpath)) {
    missing_files <- c(missing_files, sprintf("  [%s] %s", nm, fpath))
  }
}
if (length(missing_files) > 0) {
  stop("The following required input files are missing:\n",
       paste(missing_files, collapse = "\n"),
       "\nPlease run earlier stages first.")
}
message("All ", length(required_files), " required input files found.")

# --- Load RDS objects ---
message("Loading Stage 2 global object ...")
merged <- readRDS(STAGE2_RDS)
message("  ", ncol(merged), " cells x ", nrow(merged), " genes")

message("Loading Stage 3 myeloid object ...")
myeloid <- readRDS(STAGE3_RDS)
message("  ", ncol(myeloid), " cells x ", nrow(myeloid), " genes")

# Validate metadata columns in global object
global_req_meta <- c("celltype_manual", "sample_id", "group")
global_missing <- setdiff(global_req_meta, colnames(merged@meta.data))
if (length(global_missing) > 0) {
  stop("Missing metadata in global object: ", paste(global_missing, collapse = ", "))
}

# Validate metadata columns in myeloid object
mye_req_meta <- c("myeloid_subtype", "sample_id", "group")
mye_missing <- setdiff(mye_req_meta, colnames(myeloid@meta.data))
if (length(mye_missing) > 0) {
  stop("Missing metadata in myeloid object: ", paste(mye_missing, collapse = ", "))
}

message("Input validation complete.\n")

# .........................................................................
# Parameters (centralised)
# .........................................................................

# Module score: minimum genes per set
MIN_GENES_PER_SET <- 3

# Pseudobulk: minimum cells per sample per group to include a cell type
MIN_CELLS_PER_SAMPLE_PSEUDOBULK <- 10

# =========================================================================
# MODULE 1 — Unified composition delta summary
# =========================================================================
# Combine Stage 2 (global celltype) and Stage 3 (myeloid subtype) proportion
# data into a single interpretable table with effect-size emphasis.
# n = 3 per group → p-values are exploratory only.

message("===== Module 1: Unified composition delta summary =====")

# --- Load existing tables ---
glob_prop  <- read_csv(GLOB_PROP_CSV,  show_col_types = FALSE)
glob_wilcox <- read_csv(GLOB_WILCOX_CSV, show_col_types = FALSE)
mye_prop   <- read_csv(MYE_PROP_CSV,   show_col_types = FALSE)
mye_wilcox <- read_csv(MYE_WILCOX_CSV, show_col_types = FALSE)

# --- Global celltype summary ---
glob_summary <- glob_prop %>%
  group_by(group, celltype_manual) %>%
  summarise(mean_prop = mean(proportion), sd_prop = sd(proportion), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = c(mean_prop, sd_prop), values_fill = 0)

# Build global delta table
glob_delta <- glob_summary %>%
  mutate(
    level             = "global_celltype",
    mean_control      = mean_prop_Control,
    mean_apap         = mean_prop_APAP,
    delta_apap_minus_control = mean_prop_APAP - mean_prop_Control,
    fold_change_apap_vs_control = ifelse(mean_prop_Control > 0,
                                         mean_prop_APAP / mean_prop_Control,
                                         NA_real_),
    direction         = ifelse(delta_apap_minus_control > 0, "Up in APAP", "Down in APAP")
  ) %>%
  left_join(glob_wilcox %>% select(celltype, p_value),
            by = c("celltype_manual" = "celltype")) %>%
  rename(cell_group = celltype_manual, p_value_exploratory = p_value) %>%
  select(level, cell_group, mean_control, mean_apap,
         delta_apap_minus_control, fold_change_apap_vs_control,
         direction, p_value_exploratory)

# --- Myeloid subtype summary ---
mye_summary <- mye_prop %>%
  group_by(group, myeloid_subtype) %>%
  summarise(mean_prop = mean(proportion), sd_prop = sd(proportion), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = c(mean_prop, sd_prop), values_fill = 0)

mye_delta <- mye_summary %>%
  mutate(
    level             = "myeloid_subtype",
    mean_control      = mean_prop_Control,
    mean_apap         = mean_prop_APAP,
    delta_apap_minus_control = mean_prop_APAP - mean_prop_Control,
    fold_change_apap_vs_control = ifelse(mean_prop_Control > 0,
                                         mean_prop_APAP / mean_prop_Control,
                                         NA_real_),
    direction         = ifelse(delta_apap_minus_control > 0, "Up in APAP", "Down in APAP")
  ) %>%
  left_join(mye_wilcox %>% select(myeloid_subtype, p_value),
            by = "myeloid_subtype") %>%
  rename(cell_group = myeloid_subtype, p_value_exploratory = p_value) %>%
  select(level, cell_group, mean_control, mean_apap,
         delta_apap_minus_control, fold_change_apap_vs_control,
         direction, p_value_exploratory)

# --- Combine and add interpretation notes ---
composition_delta <- bind_rows(glob_delta, mye_delta) %>%
  mutate(
    interpretation_note = case_when(
      # Global trends
      cell_group == "Inflammatory neutrophil" & level == "global_celltype"
        ~ "Markedly expanded in APAP; suggests neutrophil recruitment/activation in injured liver.",
      cell_group == "Activated neutrophil" & level == "global_celltype"
        ~ "Increased in APAP; consistent with neutrophil activation in tissue damage.",
      cell_group == "Kupffer cell" & level == "global_celltype"
        ~ "Decreased proportion in APAP; likely reflects dilution by infiltrating cells and/or Kupffer loss.",
      cell_group == "Monocyte-derived macrophage / DC-like" & level == "global_celltype"
        ~ "Proportion relatively stable; may reflect balanced monocyte recruitment and differentiation.",
      cell_group == "Endothelial" & level == "global_celltype"
        ~ "Possible decrease; may reflect endothelial activation or loss in injured liver.",
      cell_group == "NK cell" & level == "global_celltype"
        ~ "Modest decrease; NK dynamics in APAP injury are context-dependent.",
      cell_group == "T cell" & level == "global_celltype"
        ~ "Relatively stable proportion; adaptive immune changes may be subtle at this time point.",
      cell_group == "B cell" & level == "global_celltype"
        ~ "Unchanged; B cells not expected primary responders in acute APAP.",
      cell_group == "Fibroblast / stellate" & level == "global_celltype"
        ~ "Possible increase; may indicate early stellate activation or ECM remodelling.",
      # Myeloid subtype trends
      cell_group == "Resident Kupffer cell" & level == "myeloid_subtype"
        ~ "Resident Kupffer proportion decreased in APAP; dilution by infiltrating myeloid or activation-induced marker loss.",
      cell_group == "Inflammatory neutrophil" & level == "myeloid_subtype"
        ~ "Expanded in APAP; major neutrophil-driven inflammatory response.",
      cell_group == "IFN-responsive inflammatory neutrophil" & level == "myeloid_subtype"
        ~ "Expanded in APAP; IFN-responsive neutrophil subset potentially involved in inflammatory signalling.",
      cell_group == "Antigen-presenting macrophage / Mo-Mac" & level == "myeloid_subtype"
        ~ "Decreased proportion; may be outcompeted by inflammatory neutrophil expansion rather than functional loss.",
      cell_group == "DC-like antigen-presenting myeloid" & level == "myeloid_subtype"
        ~ "Decreased; antigen-presenting myeloid compartments relatively contracted in acute injury.",
      cell_group == "Monocyte-derived macrophage" & level == "myeloid_subtype"
        ~ "Stable or slightly decreased; Mo-Mac may differentiate into other states post-injury.",
      cell_group == "Activated neutrophil" & level == "myeloid_subtype"
        ~ "Modestly increased; neutrophil activation states elevated in APAP.",
      cell_group == "Mature neutrophil" & level == "myeloid_subtype"
        ~ "Modestly increased; mature neutrophil pool expanded.",
      cell_group == "CCR7+ DC-like antigen-presenting myeloid" & level == "myeloid_subtype"
        ~ "Decreased; CCR7+ DC-like cells may migrate out or be diluted.",
      TRUE ~ ""
    )
  ) %>%
  arrange(level, desc(abs(delta_apap_minus_control)))

# Save
write_csv(composition_delta, file.path(TAB_DIR, "composition_delta_summary.csv"))
message("Saved: composition_delta_summary.csv (", nrow(composition_delta), " rows)")

# --- Composition delta summary plot ---
comp_plot_df <- composition_delta %>%
  filter(!is.na(fold_change_apap_vs_control), is.finite(fold_change_apap_vs_control)) %>%
  mutate(
    neg_log10_p = -log10(p_value_exploratory + 1e-300),
    label = ifelse(abs(delta_apap_minus_control) > 0.015 | neg_log10_p > 1,
                   cell_group, "")
  )

p_comp <- ggplot(comp_plot_df, aes(x = delta_apap_minus_control,
                                    y = neg_log10_p, color = level)) +
  geom_point(aes(size = abs(fold_change_apap_vs_control - 1)), alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 25, box.padding = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_color_manual(values = c("global_celltype" = "#1b9e77", "myeloid_subtype" = "#d95f02")) +
  labs(x = expression(Delta ~ "(APAP - Control proportion)"),
       y = expression(-log[10](p[exploratory])),
       color = "Level", size = "|FC - 1|",
       title = "Composition change: APAP vs Control",
       subtitle = "Global cell types (green) + myeloid subtypes (orange); n = 3 vs 3, p-values exploratory only") +
  theme_minimal(base_size = 12)

pdf(file.path(FIG_DIR, "composition_delta_summary.pdf"), width = 12, height = 8)
print(p_comp)
dev.off()
message("Saved: composition_delta_summary.pdf")

# =========================================================================
# MODULE 2 — Functional gene-set module scores
# =========================================================================
# Build small, interpretable mouse gene sets. Compute AddModuleScore on
# BOTH the global object (for key NPC types) and the myeloid object
# (to complement Stage 3 existing scores with additional functional axes).
# Aggregate by sample before comparing Control vs APAP.

message("\n===== Module 2: Functional module scores =====")

# --- Define mouse gene sets ---
# Each set is curated to be small, interpretable, and relevant to APAP injury biology.
functional_gene_sets <- list(
  "Neutrophil_inflammation" = c(
    "S100a8", "S100a9", "Mpo", "Ngp", "Camp", "Ltf", "Cxcr2", "Csf3r",
    "Retnlg", "Ly6g", "Trem1", "Il1b", "Cxcl2", "Ccl3", "Ccl4"
  ),
  "Chemokine_recruitment" = c(
    "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cxcl1", "Cxcl2",
    "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10", "Cxcl16", "Ccr2", "Ccr5"
  ),
  "IL1_TNF_inflammatory" = c(
    "Il1a", "Il1b", "Il1rn", "Tnf", "Tnfrsf1a", "Tnfrsf1b", "Nfkbia",
    "Nfkb1", "Rela", "Ptgs2", "Nos2", "Myd88", "Tlr4", "Tlr2", "Il6"
  ),
  "IFN_response" = c(
    "Ifit1", "Ifit2", "Ifit3", "Ifitm3", "Isg15", "Mx1", "Mx2",
    "Oas1a", "Oasl2", "Irf7", "Irf9", "Stat1", "Stat2", "Ifi44"
  ),
  "Antigen_presentation" = c(
    "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "H2-D1", "H2-K1",
    "B2m", "Ciita", "Flt3", "Itgax", "Ccr7", "Cd80", "Cd86"
  ),
  "Kupffer_resident_identity" = c(
    "Clec4f", "Vsig4", "Cd5l", "Timd4", "Marco", "Folr2",
    "Cd163", "C1qa", "C1qb", "C1qc", "C6", "Fcna", "Adgre1"
  ),
  "Phagocytosis_efferocytosis" = c(
    "Mertk", "Axl", "Tyro3", "Gas6", "Pros1", "Mfge8", "Timd4",
    "Stab2", "Itgav", "Itgb5", "Cd36", "Lrp1", "Elmo1", "Dock1"
  ),
  "Tissue_repair_ECM" = c(
    "Tgfb1", "Tgfb2", "Tgfbr1", "Col1a1", "Col1a2", "Col3a1",
    "Fn1", "Lox", "Timp1", "Timp2", "Mmp2", "Mmp9", "Mmp12",
    "Spp1", "Lgals3", "Vim", "Acta2", "Pdgfra", "Pdgfrb"
  ),
  "Endothelial_activation" = c(
    "Pecam1", "Cdh5", "Kdr", "Vwf", "Sele", "Selp", "Vcam1",
    "Icam1", "Icam2", "Nos3", "Thbd", "Angpt2", "Tek", "Robo4"
  ),
  "ROS_degranulation" = c(
    "Ncf1", "Ncf2", "Ncf4", "Cyba", "Cybb", "Nox1", "Sod1", "Sod2",
    "Cat", "Gpx1", "Prdx1", "Mpo", "Elane", "Prtn3", "Ctsg"
  ),
  "T_NK_activation" = c(
    "Cd3d", "Cd3e", "Cd3g", "Nkg7", "Gzma", "Gzmb", "Gzmk",
    "Prf1", "Ifng", "Tnfa", "Fasl", "Klrb1c", "Ncr1", "Ccl5"
  ),
  "Anti_inflammatory_repair" = c(
    "Il10", "Tgfb1", "Il1rn", "Arg1", "Chil3", "Retnla",
    "Mrc1", "Cd163", "Vsig4", "C1qa", "C1qb", "C1qc",
    "Lgals3", "Axl", "Mertk"
  )
)

message("Defined ", length(functional_gene_sets), " functional gene sets.")

# --- Validate gene sets against BOTH objects ---
gene_set_report <- data.frame(
  gene_set = character(0),
  n_defined = integer(0),
  n_in_global = integer(0),
  n_in_myeloid = integer(0),
  usable = logical(0),
  note = character(0),
  stringsAsFactors = FALSE
)

global_genes <- rownames(merged)
myeloid_genes <- rownames(myeloid)

for (gs_name in names(functional_gene_sets)) {
  defined <- functional_gene_sets[[gs_name]]
  in_global <- intersect(defined, global_genes)
  in_myeloid <- intersect(defined, myeloid_genes)
  n_global <- length(in_global)
  n_mye    <- length(in_myeloid)
  
  usable <- (n_global >= MIN_GENES_PER_SET) || (n_mye >= MIN_GENES_PER_SET)
  note <- ""
  if (n_global < MIN_GENES_PER_SET && n_mye < MIN_GENES_PER_SET) {
    note <- sprintf("SKIPPED: only %d genes in global, %d in myeloid (need >= %d)",
                    n_global, n_mye, MIN_GENES_PER_SET)
  } else if (n_global < MIN_GENES_PER_SET) {
    note <- sprintf("Global: %d genes (will skip global scoring)", n_global)
  }
  
  gene_set_report <- rbind(gene_set_report, data.frame(
    gene_set = gs_name, n_defined = length(defined),
    n_in_global = n_global, n_in_myeloid = n_mye,
    usable = usable, note = note,
    stringsAsFactors = FALSE
  ))
}

write_csv(gene_set_report, file.path(TAB_DIR, "gene_set_validation.csv"))
message("Saved: gene_set_validation.csv")
print(gene_set_report)

usable_sets <- gene_set_report$gene_set[gene_set_report$usable]
message("Usable gene sets: ", length(usable_sets), "/", length(functional_gene_sets))

# .........................................................................
# 2a. Module scores on GLOBAL object (for key NPC types)
# .........................................................................

message("\n--- Computing module scores on global object ---")

# Define the global-celltype gene sets to compute (skip if too few genes)
global_sets_to_score <- list()
global_skipped <- c()

for (gs_name in usable_sets) {
  in_genes <- intersect(functional_gene_sets[[gs_name]], global_genes)
  if (length(in_genes) >= MIN_GENES_PER_SET) {
    global_sets_to_score[[gs_name]] <- in_genes
  } else {
    global_skipped <- c(global_skipped, gs_name)
  }
}

if (length(global_skipped) > 0) {
  message("Skipping on global object (too few genes): ", paste(global_skipped, collapse = ", "))
}

# Compute AddModuleScore on global object
for (gs_name in names(global_sets_to_score)) {
  genes_in_set <- global_sets_to_score[[gs_name]]
  n_genes <- length(genes_in_set)
  message(sprintf("  %s: %d genes", gs_name, n_genes))
  
  merged <- AddModuleScore(merged, features = list(genes_in_set),
                           name = gs_name, assay = "RNA", ctrl = 100)
  gc()
}

# Identify module score columns (AddModuleScore appends "1")
global_score_cols <- grep("_score1$|1$", colnames(merged@meta.data), value = TRUE)
# Filter to only our functional gene sets
global_score_cols <- global_score_cols[grepl(
  paste0("^(", paste(names(functional_gene_sets), collapse = "|"), ")"), global_score_cols)]
message("Global module score columns: ", paste(global_score_cols, collapse = ", "))

# .........................................................................
# 2b. Module scores on MYELOID object (complement Stage 3)
# .........................................................................

message("\n--- Computing module scores on myeloid object ---")

# Stage 3 already has: Kupffer_resident, Monocyte_derived_Mph, Inflammatory_Mph,
#   Injury_response_Mph, Neutrophil, Activated_neutrophil, Interferon_response,
#   Antigen_presentation
# We add new functional axes that were NOT covered by Stage 3.

existing_mye_scores <- c(
  "Kupffer_resident_score1", "Monocyte_derived_Mph_score1",
  "Inflammatory_Mph_score1", "Injury_response_Mph_score1",
  "Neutrophil_score1", "Activated_neutrophil_score1",
  "Interferon_response_score1", "Antigen_presentation_score1"
)

myeloid_sets_to_score <- list()
myeloid_skipped <- c()

for (gs_name in usable_sets) {
  # Skip if already covered by Stage 3 equivalent
  in_genes <- intersect(functional_gene_sets[[gs_name]], myeloid_genes)
  if (length(in_genes) >= MIN_GENES_PER_SET) {
    myeloid_sets_to_score[[gs_name]] <- in_genes
  } else {
    myeloid_skipped <- c(myeloid_skipped, gs_name)
  }
}

if (length(myeloid_skipped) > 0) {
  message("Skipping on myeloid object (too few genes): ", paste(myeloid_skipped, collapse = ", "))
}

for (gs_name in names(myeloid_sets_to_score)) {
  genes_in_set <- myeloid_sets_to_score[[gs_name]]
  n_genes <- length(genes_in_set)
  message(sprintf("  %s: %d genes", gs_name, n_genes))
  
  myeloid <- AddModuleScore(myeloid, features = list(genes_in_set),
                            name = gs_name, assay = "RNA", ctrl = 100)
  gc()
}

# Identify new myeloid score columns
myeloid_new_score_cols <- grep("_score1$|1$", colnames(myeloid@meta.data), value = TRUE)
myeloid_new_score_cols <- myeloid_new_score_cols[grepl(
  paste0("^(", paste(names(functional_gene_sets), collapse = "|"), ")"), myeloid_new_score_cols)]
message("New myeloid module score columns: ", paste(myeloid_new_score_cols, collapse = ", "))

# --- Per-cell module scores output ---
# Global object: extract scores for cells of interest
key_celltypes_for_scoring <- c(
  "Inflammatory neutrophil", "Kupffer cell", "Neutrophil", "Activated neutrophil",
  "Monocyte-derived macrophage / DC-like", "Endothelial", "T cell", "B cell", "NK cell",
  "Fibroblast / stellate", "Cytotoxic T / NK-like"
)

global_cells_of_interest <- colnames(merged)[merged$celltype_manual %in% key_celltypes_for_scoring]

if (length(global_score_cols) > 0) {
  global_score_per_cell <- merged@meta.data[global_cells_of_interest, ] %>%
    select(sample_id, group, celltype_manual, all_of(global_score_cols))
  
  write_csv(global_score_per_cell,
            file.path(TAB_DIR, "module_scores_per_cell.csv"))
  message("Saved: module_scores_per_cell.csv (", nrow(global_score_per_cell), " cells)")
}

# Myeloid per-cell scores — combine Stage 3 existing + new
myeloid_all_score_cols <- c(existing_mye_scores[existing_mye_scores %in% colnames(myeloid@meta.data)],
                            myeloid_new_score_cols)
myeloid_score_per_cell <- myeloid@meta.data %>%
  select(sample_id, group, myeloid_subtype, all_of(myeloid_all_score_cols))

write_csv(myeloid_score_per_cell,
          file.path(TAB_DIR, "myeloid_module_scores_per_cell.csv"))
message("Saved: myeloid_module_scores_per_cell.csv (", nrow(myeloid_score_per_cell), " cells)")

# --- Per-sample aggregated module scores (global) ---
if (length(global_score_cols) > 0) {
  global_score_by_sample <- global_score_per_cell %>%
    group_by(sample_id, group, celltype_manual) %>%
    summarise(
      n_cells = n(),
      across(all_of(global_score_cols), mean, .names = "{.col}_mean"),
      .groups = "drop"
    )
  
  write_csv(global_score_by_sample,
            file.path(TAB_DIR, "module_scores_by_sample.csv"))
  message("Saved: module_scores_by_sample.csv (", nrow(global_score_by_sample), " rows)")
}

# --- Per-sample aggregated module scores (myeloid) ---
myeloid_score_by_sample <- myeloid_score_per_cell %>%
  group_by(sample_id, group, myeloid_subtype) %>%
  summarise(
    n_cells = n(),
    across(all_of(myeloid_all_score_cols), mean, .names = "{.col}_mean"),
    .groups = "drop"
  )

write_csv(myeloid_score_by_sample,
          file.path(TAB_DIR, "myeloid_module_scores_by_sample.csv"))
message("Saved: myeloid_module_scores_by_sample.csv (", nrow(myeloid_score_by_sample), " rows)")

message("Module 2 score computation complete.\n")
# =========================================================================
# MODULE 2 (continued) — Figures
# =========================================================================

message("\n--- Module 2 figures ---")

# .........................................................................
# Figure: myeloid_function_scores_by_group.pdf
# Compare module scores in myeloid subtypes between Control vs APAP,
# aggregated at sample level first.
# .........................................................................

if (length(myeloid_all_score_cols) > 0 && nrow(myeloid_score_by_sample) > 0) {
  
  # Use per-sample aggregation for group comparison
  myeloid_score_plot_data <- myeloid_score_by_sample %>%
    pivot_longer(cols = ends_with("_mean"), names_to = "module", values_to = "score") %>%
    mutate(module = gsub("_mean$", "", module)) %>%
    mutate(module = gsub("1$", "", module))

  # Select key myeloid subtypes for display (exclude low-confidence clusters)
  myeloid_subtypes_to_show <- c(
    "Inflammatory neutrophil", "Activated neutrophil", "Mature neutrophil",
    "Activated inflammatory neutrophil", "IFN-responsive inflammatory neutrophil",
    "Resident Kupffer cell", "Monocyte-derived macrophage",
    "Antigen-presenting macrophage / Mo-Mac",
    "DC-like antigen-presenting myeloid"
  )
  
  myeloid_score_plot_data <- myeloid_score_plot_data %>%
    filter(myeloid_subtype %in% myeloid_subtypes_to_show)
  
  p_mye_scores <- ggplot(myeloid_score_plot_data,
                          aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
    facet_grid(module ~ myeloid_subtype, scales = "free_y") +
    scale_fill_manual(values = group_cols) +
    labs(x = "", y = "Module score (per-sample mean)", fill = "Group",
         title = "Myeloid functional module scores: Control vs APAP",
         subtitle = "Each point = one sample; scores aggregated per sample before comparison") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(size = 7),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom")
  
  n_modules <- length(unique(myeloid_score_plot_data$module))
  n_subtypes <- length(unique(myeloid_score_plot_data$myeloid_subtype))
  
  pdf(file.path(FIG_DIR, "myeloid_function_scores_by_group.pdf"),
      width = max(14, n_subtypes * 1.8), height = max(8, n_modules * 1.2))
  print(p_mye_scores)
  dev.off()
  message("Saved: myeloid_function_scores_by_group.pdf")
}

# .........................................................................
# Figure: key_celltype_function_scores_by_group.pdf
# Compare module scores in key global cell types (Control vs APAP)
# .........................................................................

if (length(global_score_cols) > 0 && exists("global_score_by_sample")) {
  
  global_score_plot_data <- global_score_by_sample %>%
    pivot_longer(cols = ends_with("_mean"), names_to = "module", values_to = "score") %>%
    mutate(module = gsub("_mean$", "", module)) %>%
    mutate(module = gsub("1$", "", module))

  # Select a focused set of cell types and modules for clean display
  display_celltypes <- c(
    "Inflammatory neutrophil", "Kupffer cell", "Neutrophil",
    "Monocyte-derived macrophage / DC-like", "Endothelial",
    "T cell", "NK cell", "Fibroblast / stellate"
  )
  
  display_modules <- c(
    "Neutrophil_inflammation", "Chemokine_recruitment",
    "IL1_TNF_inflammatory", "IFN_response", "Antigen_presentation",
    "Kupffer_resident_identity", "Phagocytosis_efferocytosis",
    "Tissue_repair_ECM", "Endothelial_activation"
  )
  
  global_score_plot_data <- global_score_plot_data %>%
    filter(celltype_manual %in% display_celltypes,
           module %in% display_modules)
  
  if (nrow(global_score_plot_data) > 0) {
    p_glob_scores <- ggplot(global_score_plot_data,
                             aes(x = group, y = score, fill = group)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
      facet_grid(module ~ celltype_manual, scales = "free_y") +
      scale_fill_manual(values = group_cols) +
      labs(x = "", y = "Module score (per-sample mean)", fill = "Group",
           title = "Key cell-type functional module scores: Control vs APAP",
           subtitle = "Each point = one sample; scores aggregated per sample before comparison") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y = element_text(size = 7),
            strip.text.x = element_text(size = 8),
            legend.position = "bottom")
    
    n_mod  <- length(unique(global_score_plot_data$module))
    n_ct   <- length(unique(global_score_plot_data$celltype_manual))
    
    pdf(file.path(FIG_DIR, "key_celltype_function_scores_by_group.pdf"),
        width = max(14, n_ct * 1.8), height = max(8, n_mod * 1.2))
    print(p_glob_scores)
    dev.off()
    message("Saved: key_celltype_function_scores_by_group.pdf")
  }
}

# .........................................................................
# Figure: function_heatmap_selected_celltypes.pdf
# Heatmap of mean module scores (APAP - Control) for selected cell types
# .........................................................................

if (exists("global_score_by_sample") && nrow(global_score_by_sample) > 0) {
  
  heatmap_celltypes <- c(
    "Inflammatory neutrophil", "Kupffer cell", "Neutrophil", "Activated neutrophil",
    "Monocyte-derived macrophage / DC-like", "Endothelial",
    "T cell", "NK cell", "B cell", "Fibroblast / stellate"
  )
  
  heatmap_modules <- intersect(display_modules, gsub("1$", "", global_score_cols))
  
  # Compute delta (APAP mean - Control mean) per celltype per module
  heatmap_data <- global_score_by_sample %>%
    filter(celltype_manual %in% heatmap_celltypes) %>%
    pivot_longer(cols = ends_with("_mean"), names_to = "module", values_to = "score") %>%
    mutate(module = gsub("_mean$", "", module)) %>%
    mutate(module = gsub("1$", "", module)) %>%
    filter(module %in% heatmap_modules) %>%
    group_by(celltype_manual, module, group) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_score, values_fill = 0) %>%
    mutate(delta = APAP - Control)
  
  if (nrow(heatmap_data) > 0) {
    # Reshape to matrix
    heatmap_mat <- heatmap_data %>%
      select(celltype_manual, module, delta) %>%
      pivot_wider(names_from = celltype_manual, values_from = delta, values_fill = 0) %>%
      tibble::column_to_rownames("module") %>%
      as.matrix()
    
    # Compute per-sample delta for MSD error bars (optional: just show mean delta)
    # Cap extreme values for colour scale
    max_abs <- max(abs(heatmap_mat), na.rm = TRUE)
    breaks <- seq(-max_abs, max_abs, length.out = 101)
    
    pdf(file.path(FIG_DIR, "function_heatmap_selected_celltypes.pdf"),
        width = 12, height = 8)
    pheatmap(heatmap_mat,
             color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
             breaks = breaks,
             cluster_rows = TRUE, cluster_cols = TRUE,
             display_numbers = TRUE,
             number_format = "%.3f",
             fontsize_number = 8,
             main = "Module score delta (APAP - Control): selected cell types",
             angle_col = 45,
             border_color = NA)
    dev.off()
    message("Saved: function_heatmap_selected_celltypes.pdf")
  }
}

message("Module 2 figures complete.\n")

# =========================================================================
# MODULE 3 — Sample-level pseudobulk DE with edgeR
# =========================================================================
# For each cell type / subtype:
#   1. Build gene × sample count matrix from raw counts
#   2. edgeR pipeline: DGEList → filterByExpr → calcNormFactors (within-group)
#      → estimateDisp → glmQLFit → glmQLFTest
#   3. Output: logFC, logCPM, PValue, FDR
#
# Key principle: each sample = one biological replicate.
# Normalisation factors computed within each cell group (TMM), not globally.
# =========================================================================

message("===== Module 3: Pseudobulk DE with edgeR =====")

# --- Helper: build gene x sample pseudobulk count matrix for one cell group ---
build_pb_matrix <- function(seurat_obj, cell_group, group_by_col,
                            sample_col = "sample_id", group_col = "group") {
  meta <- seurat_obj@meta.data
  samples <- unique(meta[[sample_col]])

  count_list <- list()
  sample_groups <- c()
  n_cells_vec <- c()

  for (s in samples) {
    cells <- colnames(seurat_obj)[meta[[group_by_col]] == cell_group & meta[[sample_col]] == s]
    if (length(cells) < MIN_CELLS_PER_SAMPLE_PSEUDOBULK) next

    if (length(cells) == 1) {
      pb <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
    } else {
      pb <- Matrix::rowSums(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE])
    }
    count_list[[s]] <- as.numeric(pb)
    sample_groups <- c(sample_groups, meta[[group_col]][meta[[sample_col]] == s][1])
    n_cells_vec <- c(n_cells_vec, length(cells))
  }

  if (length(count_list) == 0) return(NULL)

  count_mat <- do.call(cbind, count_list)
  rownames(count_mat) <- rownames(seurat_obj)
  colnames(count_mat) <- names(count_list)

  list(
    counts       = count_mat,
    sample_id    = names(count_list),
    group        = sample_groups,
    n_cells      = n_cells_vec
  )
}

# --- edgeR pseudobulk DE for one cell group ---
run_edger_pseudobulk <- function(pb_obj, cell_group_name) {
  counts <- pb_obj$counts
  group  <- factor(pb_obj$group, levels = c("Control", "APAP"))

  # Remove genes with zero counts across all samples
  keep_genes <- rowSums(counts) > 0
  counts <- counts[keep_genes, , drop = FALSE]

  if (nrow(counts) < 10) {
    return(list(result = NULL, skip_reason = sprintf("Too few expressed genes: %d", nrow(counts))))
  }

  # Check sample counts per group
  n_ctrl <- sum(group == "Control")
  n_apap <- sum(group == "APAP")
  if (n_ctrl < 2 || n_apap < 2) {
    return(list(result = NULL,
                skip_reason = sprintf("Insufficient samples: %d Control, %d APAP", n_ctrl, n_apap)))
  }

  # edgeR pipeline
  dge <- DGEList(counts = counts, group = group)

  # filterByExpr within this cell group
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  if (nrow(dge) < 10) {
    return(list(result = NULL, skip_reason = sprintf("Too few genes after filterByExpr: %d", nrow(dge))))
  }

  # TMM normalisation within this cell group
  dge <- calcNormFactors(dge, method = "TMM")

  # Design matrix
  design <- model.matrix(~ group)

  # Estimate dispersion and fit QL model
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)  # coef=2 tests APAP vs Control

  # Extract results
  tt <- topTags(qlf, n = nrow(dge), sort.by = "none")$table
  tt$gene <- rownames(tt)
  tt$cell_group <- cell_group_name
  tt <- tt[, c("gene", "logFC", "logCPM", "PValue", "FDR", "cell_group")]
  list(result = tt, skip_reason = NULL)
}

# --- Skip groups (low-confidence / contamination) ---
skip_groups_global <- c(
  "Low-quality / mixed", "Hepatocyte contamination",
  "B cell / possible doublet", "Mast / basophil-like",
  "pDC-like", "cDC1", "Mesothelial", "Cycling cell", "Cholangiocyte / epithelial",
  "VSMC / pericyte", "Plasma cell"
)

skip_groups_myeloid <- c(
  "Possible doublet / hepatocyte ambient-contaminated neutrophil",
  "Neutrophil with erythroid/ambient signal",
  "Antigen-presenting myeloid with lymphoid-like signal",
  "Rare injury-associated macrophage-like cells"
)

skip_log <- data.frame(
  cell_group = character(0),
  level      = character(0),
  reason     = character(0),
  stringsAsFactors = FALSE
)

# =========================================================================
# 3a. Pseudobulk DE: Global cell types
# =========================================================================

message("\n--- edgeR pseudobulk: global cell types ---")

global_cell_groups <- setdiff(unique(merged$celltype_manual), skip_groups_global)
global_de_results <- NULL

for (cg in global_cell_groups) {
  pb <- build_pb_matrix(merged, cg, "celltype_manual")

  if (is.null(pb)) {
    skip_log <- rbind(skip_log, data.frame(
      cell_group = cg, level = "global_celltype",
      reason = "No sample with sufficient cells",
      stringsAsFactors = FALSE
    ))
    next
  }

  res <- run_edger_pseudobulk(pb, cg)

  if (!is.null(res$skip_reason)) {
    skip_log <- rbind(skip_log, data.frame(
      cell_group = cg, level = "global_celltype",
      reason = res$skip_reason,
      stringsAsFactors = FALSE
    ))
  } else {
    res$result$level <- "global_celltype"
    global_de_results <- if (is.null(global_de_results)) res$result else rbind(global_de_results, res$result)
  }
}

n_global_tested <- if (is.null(global_de_results)) 0L else length(unique(global_de_results$cell_group))
n_global_rows   <- if (is.null(global_de_results)) 0L else nrow(global_de_results)
message("Global pseudobulk: ", n_global_tested,
        " cell groups tested, ",
        nrow(skip_log %>% filter(level == "global_celltype")), " skipped.")

if (!is.null(global_de_results)) {
  write_csv(global_de_results,
            file.path(TAB_DIR, "pseudobulk_de_by_celltype.csv"))
  message("Saved: pseudobulk_de_by_celltype.csv (", n_global_rows, " rows)")
} else {
  message("No global pseudobulk results to save.")
}

# =========================================================================
# 3b. Pseudobulk DE: Myeloid subtypes
# =========================================================================

message("\n--- edgeR pseudobulk: myeloid subtypes ---")

myeloid <- JoinLayers(myeloid)
myeloid_cell_groups <- setdiff(unique(myeloid$myeloid_subtype), skip_groups_myeloid)
myeloid_de_results <- NULL

for (cg in myeloid_cell_groups) {
  pb <- build_pb_matrix(myeloid, cg, "myeloid_subtype")

  if (is.null(pb)) {
    skip_log <- rbind(skip_log, data.frame(
      cell_group = cg, level = "myeloid_subtype",
      reason = "No sample with sufficient cells",
      stringsAsFactors = FALSE
    ))
    next
  }

  res <- run_edger_pseudobulk(pb, cg)

  if (!is.null(res$skip_reason)) {
    skip_log <- rbind(skip_log, data.frame(
      cell_group = cg, level = "myeloid_subtype",
      reason = res$skip_reason,
      stringsAsFactors = FALSE
    ))
  } else {
    res$result$level <- "myeloid_subtype"
    myeloid_de_results <- if (is.null(myeloid_de_results)) res$result else rbind(myeloid_de_results, res$result)
  }
}

n_myeloid_tested <- if (is.null(myeloid_de_results)) 0L else length(unique(myeloid_de_results$cell_group))
n_myeloid_rows   <- if (is.null(myeloid_de_results)) 0L else nrow(myeloid_de_results)
message("Myeloid pseudobulk: ", n_myeloid_tested,
        " subtypes tested, ",
        nrow(skip_log %>% filter(level == "myeloid_subtype")), " skipped.")

if (!is.null(myeloid_de_results)) {
  write_csv(myeloid_de_results,
            file.path(TAB_DIR, "pseudobulk_de_by_myeloid_subtype.csv"))
  message("Saved: pseudobulk_de_by_myeloid_subtype.csv (", n_myeloid_rows, " rows)")
} else {
  message("No myeloid pseudobulk results to save.")
}

# Save skip log
write_csv(skip_log, file.path(TAB_DIR, "pseudobulk_skip_log.csv"))
message("Saved: pseudobulk_skip_log.csv (", nrow(skip_log), " rows)")
if (nrow(skip_log) > 0) print(skip_log)

# =========================================================================
# 3c. Volcano plots (FDR-based significance)
# =========================================================================

message("\n--- Pseudobulk volcano plots ---")

FDR_THRESHOLD <- 0.05

make_volcano <- function(de_df, cell_group_name, fdr_thresh = FDR_THRESHOLD, n_label = 20) {
  sub <- de_df %>% filter(cell_group == cell_group_name)
  if (nrow(sub) < 10) return(NULL)

  sub <- sub %>%
    mutate(
      neg_log10_fdr = -log10(FDR + 1e-300),
      is_sig        = !is.na(FDR) & FDR < fdr_thresh,
      abs_logFC     = abs(logFC)
    )

  top_genes <- sub %>%
    filter(is_sig) %>%
    slice_max(n = n_label, order_by = abs_logFC) %>%
    pull(gene)

  sub$label <- ifelse(sub$gene %in% top_genes, sub$gene, "")

  n_sig <- sum(sub$is_sig)

  p <- ggplot(sub, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(aes(color = is_sig), size = 0.8, alpha = 0.6) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                    box.padding = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dotted", colour = "grey50") +
    scale_color_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "grey60")) +
    labs(x = "log2(APAP / Control)", y = expression(-log[10](FDR)),
         title = paste("Pseudobulk DE:", cell_group_name),
         subtitle = sprintf("edgeR QLF; FDR < %s = %d genes",
                            fdr_thresh, n_sig)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  p
}

volcano_celltypes <- intersect(
  c("Inflammatory neutrophil", "Kupffer cell", "Monocyte-derived macrophage / DC-like",
    "Endothelial", "T cell", "NK cell", "B cell", "Neutrophil", "Activated neutrophil"),
  unique(global_de_results$cell_group)
)

volcano_plots <- list()
for (ct in volcano_celltypes) {
  vp <- make_volcano(global_de_results, ct)
  if (!is.null(vp)) volcano_plots[[ct]] <- vp
}

volcano_mye_subtypes <- intersect(
  c("Inflammatory neutrophil", "Resident Kupffer cell", "Activated neutrophil",
    "IFN-responsive inflammatory neutrophil", "Antigen-presenting macrophage / Mo-Mac",
    "DC-like antigen-presenting myeloid", "Monocyte-derived macrophage"),
  unique(myeloid_de_results$cell_group)
)

for (ct in volcano_mye_subtypes) {
  vp <- make_volcano(myeloid_de_results, ct)
  if (!is.null(vp)) volcano_plots[[paste0("Mye:", ct)]] <- vp
}

if (length(volcano_plots) > 0) {
  ncol_vol <- min(3, length(volcano_plots))
  nrow_vol <- ceiling(length(volcano_plots) / ncol_vol)

  p_volcano <- wrap_plots(volcano_plots, ncol = ncol_vol) +
    plot_annotation(title = "Pseudobulk DE (edgeR QLF): selected cell types")

  pdf(file.path(FIG_DIR, "pseudobulk_volcano_selected_celltypes.pdf"),
      width = 6 * ncol_vol, height = 5 * nrow_vol)
  print(p_volcano)
  dev.off()
  message("Saved: pseudobulk_volcano_selected_celltypes.pdf (",
          length(volcano_plots), " panels)")
}

message("Module 3 complete.\n")

# =========================================================================
# MODULE 4 — Semi-automated evidence table + key findings
# =========================================================================
# Two outputs:
#   1. stage4_interpretation_summary.csv — detailed evidence per source per group
#   2. stage4_key_findings_summary.csv   — one-row-per-group README-ready table
#
# Design rules:
#   - module_score describes functional programs only; no single direction
#   - confidence from composition + pseudobulk DE concordance (not module scores)
#   - analysis_level distinguishes global_celltype vs myeloid_subtype
# =========================================================================

message("===== Module 4: Evidence table + key findings =====")

# --- Helper: compute module score delta per cell group ---
compute_module_deltas <- function(score_by_sample_df, group_col = "celltype_manual") {
  if (is.null(score_by_sample_df) || nrow(score_by_sample_df) == 0) return(NULL)

  score_cols <- grep("_mean$", colnames(score_by_sample_df), value = TRUE)
  if (length(score_cols) == 0) return(NULL)

  deltas <- score_by_sample_df %>%
    group_by(!!sym(group_col)) %>%
    summarise(
      across(all_of(score_cols), ~ mean(.x[group == "APAP"], na.rm = TRUE) -
               mean(.x[group == "Control"], na.rm = TRUE)),
      .groups = "drop"
    )

  deltas_long <- deltas %>%
    pivot_longer(-!!sym(group_col), names_to = "module", values_to = "delta") %>%
    mutate(module = gsub("_mean$", "", module)) %>%
    mutate(module = gsub("1$", "", module)) %>%
    filter(!is.na(delta)) %>%
    # Keep only Stage 4 functional gene sets (exclude Stage 3 scores like "Neutrophil_score")
    filter(module %in% names(functional_gene_sets))

  deltas_long
}

# --- Helper: get top pseudobulk DEGs ---
get_top_degs <- function(de_df, cg, fdr_thresh = 0.05, n_top = 10) {
  sub <- de_df %>% filter(cell_group == cg)

  up <- sub %>% filter(FDR < fdr_thresh, logFC > 0) %>%
    slice_max(n = n_top, order_by = logFC) %>% pull(gene)
  down <- sub %>% filter(FDR < fdr_thresh, logFC < 0) %>%
    slice_min(n = n_top, order_by = logFC) %>% pull(gene)

  n_sig_up <- sub %>% filter(FDR < fdr_thresh, logFC > 0) %>% nrow()
  n_sig_down <- sub %>% filter(FDR < fdr_thresh, logFC < 0) %>% nrow()

  list(
    up = up, down = down,
    n_up = n_sig_up, n_down = n_sig_down,
    n_total_tested = nrow(sub),
    n_sig_total = n_sig_up + n_sig_down
  )
}

# --- Helper: pick best analysis_level for a cell group ---
pick_level <- function(cg) {
  # If cell group exists in myeloid_subtype, use that (more specific)
  if (exists("myeloid_de_results") && cg %in% unique(myeloid_de_results$cell_group)) return("myeloid_subtype")
  if (exists("global_de_results") && cg %in% unique(global_de_results$cell_group)) return("global_celltype")
  # Fall back to composition_delta level
  comp_lvl <- composition_delta$level[composition_delta$cell_group == cg]
  if (length(comp_lvl) > 0) return(comp_lvl[1])
  "unknown"
}

# --- Helper: assign confidence based on comp + DE concordance ---
assign_confidence_v2 <- function(comp_dir, de_dir, de_sig_count) {
  # Module scores are NOT used for direction consistency
  # comp_dir: 1=expanded, -1=contracted, 0=no data
  # de_dir: 1=up, -1=down, 0=no sig DEGs or no DE data
  # de_sig_count: number of significant DEGs (FDR < 0.05)

  if (comp_dir == 0 && de_dir == 0) return("No evidence")

  # Both agree
  if (comp_dir != 0 && comp_dir == de_dir) {
    if (de_sig_count >= 50) return("High")
    if (de_sig_count >= 10) return("Moderate")
    if (de_sig_count > 0) return("Low")
    return("Low")  # comp points but no DE sig
  }

  # Only one evidence type available
  if (comp_dir == 0 || de_dir == 0) {
    if (de_sig_count >= 10) return("Moderate")
    if (de_sig_count > 0) return("Low")
    return("Low")
  }

  # Disagree
  return("Uncertain")
}

# --- Helper: compose semi-auto biological interpretation ---
compose_interp <- function(cg, comp_dir, comp_detail, mod_up, mod_down,
                           de_dir, de_top_up, de_top_down, de_n_up, de_n_down) {
  parts <- c()

  # Composition narrative
  if (comp_dir == 1) parts <- c(parts, sprintf("Expanded in APAP (%.1f%%)", comp_detail * 100))
  if (comp_dir == -1) parts <- c(parts, sprintf("Contracted in APAP (%.1f%%)", comp_detail * 100))

  # Module functional description (no direction)
  if (length(mod_up) > 0) parts <- c(parts, sprintf("UP modules: %s", paste(head(mod_up, 4), collapse = ", ")))
  if (length(mod_down) > 0) parts <- c(parts, sprintf("DOWN modules: %s", paste(head(mod_down, 4), collapse = ", ")))

  # DE summary
  if (de_dir != 0) {
    parts <- c(parts, sprintf("DE: %d up, %d down (FDR<0.05)", de_n_up, de_n_down))
    if (length(de_top_up) > 0) parts <- c(parts, sprintf("Top up: %s", paste(head(de_top_up, 5), collapse = ", ")))
    if (length(de_top_down) > 0) parts <- c(parts, sprintf("Top down: %s", paste(head(de_top_down, 5), collapse = ", ")))
  } else if (de_n_up + de_n_down > 0) {
    parts <- c(parts, sprintf("DE: %d up, %d down (FDR<0.05)", de_n_up, de_n_down))
  } else {
    parts <- c(parts, "DE: no significant DEGs (FDR<0.05)")
  }

  paste(parts, collapse = "; ")
}

# --- Build detailed evidence table (per-source rows) ---
build_evidence_table <- function() {
  evidence <- data.frame(
    analysis_level     = character(0),
    cell_group         = character(0),
    evidence_source    = character(0),
    direction          = character(0),
    evidence_detail    = character(0),
    confidence         = character(0),
    stringsAsFactors   = FALSE
  )

  cg_priority <- intersect(
    c("Inflammatory neutrophil", "Resident Kupffer cell",
      "Monocyte-derived macrophage / DC-like", "Monocyte-derived macrophage",
      "Antigen-presenting macrophage / Mo-Mac", "DC-like antigen-presenting myeloid",
      "Endothelial", "T cell", "NK cell", "B cell",
      "Fibroblast / stellate", "Activated neutrophil", "Mature neutrophil",
      "IFN-responsive inflammatory neutrophil", "Activated inflammatory neutrophil",
      "Kupffer cell", "Neutrophil"),
    unique(c(composition_delta$cell_group,
             if(exists("global_de_results")) unique(global_de_results$cell_group) else NULL,
             if(exists("myeloid_de_results")) unique(myeloid_de_results$cell_group) else NULL))
  )

  global_module_deltas <- if (exists("global_score_by_sample")) {
    compute_module_deltas(global_score_by_sample, "celltype_manual")
  } else NULL

  myeloid_module_deltas <- if (exists("myeloid_score_by_sample")) {
    compute_module_deltas(myeloid_score_by_sample, "myeloid_subtype")
  } else NULL

  for (cg in cg_priority) {
    alevel <- pick_level(cg)
    evidence_lines <- list()
    comp_dir_val <- 0
    de_dir_val <- 0
    de_sig <- 0

    # ---- Evidence 1: Composition ----
    comp_row <- composition_delta %>% filter(cell_group == cg) %>% slice(1)

    if (nrow(comp_row) > 0 && !is.na(comp_row$delta_apap_minus_control)) {
      delta_val <- comp_row$delta_apap_minus_control
      dir <- if (delta_val > 0) "Expanded" else "Contracted"
      comp_dir_val <- sign(delta_val)

      evidence_lines[[length(evidence_lines) + 1]] <- list(
        source = "composition",
        direction = dir,
        detail = sprintf("delta_proportion = %+.4f; mean_APAP = %.4f, mean_Control = %.4f",
                         delta_val, comp_row$mean_apap, comp_row$mean_control)
      )
    }

    # ---- Evidence 2: Module scores (functional description only, no direction) ----
    mod_deltas <- NULL
    mod_level <- NULL
    if (!is.null(global_module_deltas) && cg %in% global_module_deltas$celltype_manual) {
      mod_deltas <- global_module_deltas %>% filter(celltype_manual == cg)
      mod_level <- "global"
    }
    if ((is.null(mod_deltas) || nrow(mod_deltas) == 0) &&
        !is.null(myeloid_module_deltas) && cg %in% myeloid_module_deltas$myeloid_subtype) {
      mod_deltas <- myeloid_module_deltas %>% filter(myeloid_subtype == cg)
      mod_level <- "myeloid"
    }

    if (!is.null(mod_deltas) && nrow(mod_deltas) > 0) {
      top_up <- mod_deltas %>% slice_max(n = 3, order_by = delta) %>%
        filter(delta > 0.01) %>% pull(module)
      top_down <- mod_deltas %>% slice_min(n = 3, order_by = delta) %>%
        filter(delta < -0.01) %>% pull(module)

      detail_parts <- c()
      if (length(top_up) > 0) detail_parts <- c(detail_parts, paste("Up:", paste(top_up, collapse = ", ")))
      if (length(top_down) > 0) detail_parts <- c(detail_parts, paste("Down:", paste(top_down, collapse = ", ")))
      if (length(detail_parts) == 0) detail_parts <- "No module with |delta| > 0.01"

      evidence_lines[[length(evidence_lines) + 1]] <- list(
        source = "module_score",
        direction = "Functional description only",
        detail = paste(detail_parts, collapse = "; ")
      )
    }

    # ---- Evidence 3: Pseudobulk DE (edgeR) ----
    de_df <- NULL
    if (exists("global_de_results") && cg %in% global_de_results$cell_group) {
      de_df <- global_de_results
    }
    if (exists("myeloid_de_results") && cg %in% myeloid_de_results$cell_group) {
      de_df <- myeloid_de_results
    }

    if (!is.null(de_df)) {
      degs <- get_top_degs(de_df, cg, fdr_thresh = FDR_THRESHOLD)
      de_sig <- degs$n_sig_total

      if (degs$n_sig_total > 0) {
        de_dir_val <- sign(degs$n_up - degs$n_down)
        dir <- if (de_dir_val > 0) "Upregulated" else if (de_dir_val < 0) "Downregulated" else "Balanced"

        detail_parts <- c()
        if (length(degs$up) > 0) {
          detail_parts <- c(detail_parts,
                            sprintf("Up (FDR<%s): %s", FDR_THRESHOLD, paste(head(degs$up, 5), collapse = ", ")))
        }
        if (length(degs$down) > 0) {
          detail_parts <- c(detail_parts,
                            sprintf("Down (FDR<%s): %s", FDR_THRESHOLD, paste(head(degs$down, 5), collapse = ", ")))
        }

        evidence_lines[[length(evidence_lines) + 1]] <- list(
          source = "pseudobulk_DE",
          direction = dir,
          detail = sprintf("n_up=%d, n_down=%d, n_total_tested=%d. %s",
                           degs$n_up, degs$n_down, degs$n_total_tested,
                           paste(detail_parts, collapse = "; "))
        )
      } else if (degs$n_total_tested > 0) {
        evidence_lines[[length(evidence_lines) + 1]] <- list(
          source = "pseudobulk_DE",
          direction = "No significant DEGs",
          detail = sprintf("0/%d genes with FDR < %s", degs$n_total_tested, FDR_THRESHOLD)
        )
      }
    }

    # ---- Assemble rows ----
    if (length(evidence_lines) == 0) next

    conf <- assign_confidence_v2(comp_dir_val, de_dir_val, de_sig)

    for (ev in evidence_lines) {
      evidence <- rbind(evidence, data.frame(
        analysis_level   = alevel,
        cell_group       = cg,
        evidence_source  = ev$source,
        direction        = ev$direction,
        evidence_detail  = ev$detail,
        confidence       = conf,
        stringsAsFactors = FALSE
      ))
    }
  }

  evidence
}

# --- Build key findings table (one row per cell group, README-ready) ---
build_key_findings_table <- function() {
  findings <- data.frame(
    analysis_level          = character(0),
    cell_group              = character(0),
    composition_change      = character(0),
    key_up_modules          = character(0),
    key_down_modules        = character(0),
    top_up_DEGs             = character(0),
    top_down_DEGs           = character(0),
    biological_interpretation = character(0),
    confidence              = character(0),
    caveat                  = character(0),
    stringsAsFactors        = FALSE
  )

  cg_priority <- intersect(
    c("Inflammatory neutrophil", "Resident Kupffer cell",
      "Monocyte-derived macrophage / DC-like", "Monocyte-derived macrophage",
      "Antigen-presenting macrophage / Mo-Mac", "DC-like antigen-presenting myeloid",
      "Endothelial", "T cell", "NK cell", "B cell",
      "Fibroblast / stellate", "Activated neutrophil", "Mature neutrophil",
      "IFN-responsive inflammatory neutrophil", "Activated inflammatory neutrophil",
      "Kupffer cell", "Neutrophil"),
    unique(c(composition_delta$cell_group,
             if (exists("global_de_results")) unique(global_de_results$cell_group) else NULL,
             if (exists("myeloid_de_results")) unique(myeloid_de_results$cell_group) else NULL))
  )

  global_module_deltas <- if (exists("global_score_by_sample")) {
    compute_module_deltas(global_score_by_sample, "celltype_manual")
  } else NULL

  myeloid_module_deltas <- if (exists("myeloid_score_by_sample")) {
    compute_module_deltas(myeloid_score_by_sample, "myeloid_subtype")
  } else NULL

  for (cg in cg_priority) {
    alevel    <- pick_level(cg)
    comp_val  <- 0
    de_val    <- 0
    de_sig    <- 0
    comp_d    <- NA_real_
    comp_ctl  <- NA_real_
    comp_apap <- NA_real_

    # ---- Composition ----
    comp_row <- composition_delta %>% filter(cell_group == cg) %>% slice(1)
    if (nrow(comp_row) > 0 && !is.na(comp_row$delta_apap_minus_control)) {
      comp_val  <- sign(comp_row$delta_apap_minus_control)
      comp_d    <- abs(comp_row$delta_apap_minus_control)
      comp_ctl  <- comp_row$mean_control
      comp_apap <- comp_row$mean_apap
    }

    comp_change <- if (comp_val == 1) {
      sprintf("Expanded (%.1f%% -> %.1f%%, delta=+%.4f)", comp_ctl * 100, comp_apap * 100, comp_d)
    } else if (comp_val == -1) {
      sprintf("Contracted (%.1f%% -> %.1f%%, delta=%.4f)", comp_ctl * 100, comp_apap * 100, -comp_d)
    } else {
      ""
    }

    # ---- Module scores (functional context, no direction) ----
    mod_up   <- character(0)
    mod_down <- character(0)
    mod_deltas <- NULL
    if (!is.null(global_module_deltas) && cg %in% global_module_deltas$celltype_manual) {
      mod_deltas <- global_module_deltas %>% filter(celltype_manual == cg)
    }
    if ((is.null(mod_deltas) || nrow(mod_deltas) == 0) &&
        !is.null(myeloid_module_deltas) && cg %in% myeloid_module_deltas$myeloid_subtype) {
      mod_deltas <- myeloid_module_deltas %>% filter(myeloid_subtype == cg)
    }
    if (!is.null(mod_deltas) && nrow(mod_deltas) > 0) {
      mod_up   <- mod_deltas %>% slice_max(n = 5, order_by = delta) %>% filter(delta > 0.01) %>% pull(module)
      mod_down <- mod_deltas %>% slice_min(n = 5, order_by = delta) %>% filter(delta < -0.01) %>% pull(module)
    }

    # ---- Pseudobulk DE ----
    de_top_up   <- character(0)
    de_top_down <- character(0)
    de_n_up <- 0; de_n_down <- 0

    de_df <- NULL
    if (exists("global_de_results") && cg %in% global_de_results$cell_group) de_df <- global_de_results
    if (exists("myeloid_de_results") && cg %in% myeloid_de_results$cell_group) de_df <- myeloid_de_results

    if (!is.null(de_df)) {
      degs <- get_top_degs(de_df, cg, fdr_thresh = FDR_THRESHOLD, n_top = 10)
      de_sig    <- degs$n_sig_total
      de_n_up   <- degs$n_up
      de_n_down <- degs$n_down
      de_val    <- if (degs$n_sig_total == 0) 0 else sign(degs$n_up - degs$n_down)
      de_top_up   <- degs$up
      de_top_down <- degs$down
    }

    # ---- Confidence ----
    conf <- assign_confidence_v2(comp_val, de_val, de_sig)

    # ---- Biological interpretation (semi-auto) ----
    interp <- compose_interp(cg, comp_val, comp_d, mod_up, mod_down,
                             de_val, de_top_up, de_top_down, de_n_up, de_n_down)

    # ---- Caveat ----
    caveat <- switch(conf,
      "High"    = "n = 3 per group; moderate sample size for definitive conclusions.",
      "Moderate" = "n = 3 per group limits statistical power; FDR-significant DEGs detected but effect size may vary.",
      "Low"     = "n = 3 per group; trend-level observation requiring validation.",
      "Uncertain" = "Evidence direction between composition and DE is discordant. n = 3 per group limits resolution.",
      "n = 3 per group; hypothesis-generating."
    )

    findings <- rbind(findings, data.frame(
      analysis_level            = alevel,
      cell_group                = cg,
      composition_change        = comp_change,
      key_up_modules            = if (length(mod_up) > 0) paste(mod_up, collapse = ", ") else "",
      key_down_modules          = if (length(mod_down) > 0) paste(mod_down, collapse = ", ") else "",
      top_up_DEGs               = if (length(de_top_up) > 0) paste(head(de_top_up, 10), collapse = ", ") else "",
      top_down_DEGs             = if (length(de_top_down) > 0) paste(head(de_top_down, 10), collapse = ", ") else "",
      biological_interpretation = interp,
      confidence                = conf,
      caveat                    = caveat,
      stringsAsFactors          = FALSE
    ))
  }

  findings
}

# Build and save
interpretation_summary <- build_evidence_table()
key_findings_summary   <- build_key_findings_table()

write_csv(interpretation_summary,
          file.path(TAB_DIR, "stage4_interpretation_summary.csv"))
message("Saved: stage4_interpretation_summary.csv (",
        nrow(interpretation_summary), " rows)")

write_csv(key_findings_summary,
          file.path(TAB_DIR, "stage4_key_findings_summary.csv"))
message("Saved: stage4_key_findings_summary.csv (",
        nrow(key_findings_summary), " rows)")

# =========================================================================
# Final: save outputs and clean up
# =========================================================================

message("\n===== Saving final outputs =====")

# Rename AddModuleScore columns to clean names (strip the "1" suffix)
rename_score_cols <- function(meta, score_prefixes) {
  for (sp in score_prefixes) {
    old_col <- paste0(sp, "1")
    if (old_col %in% colnames(meta)) {
      colnames(meta)[colnames(meta) == old_col] <- sp
    }
  }
  meta
}

global_prefixes <- names(functional_gene_sets)
merged@meta.data <- rename_score_cols(merged@meta.data, global_prefixes)

myeloid_prefixes <- names(functional_gene_sets)
myeloid@meta.data <- rename_score_cols(myeloid@meta.data, myeloid_prefixes)

# --- Build output list ---
# Contains: the annotated global Seurat object (with Stage 4 module scores),
# the annotated myeloid Seurat object (with Stage 4 module scores),
# and all summary tables for CellChat (Stage 5) and downstream use.

stage4_output <- list(
  seurat_global          = merged,
  seurat_myeloid         = myeloid,
  composition_delta      = composition_delta,
  module_score_summary   = if (exists("global_score_by_sample")) global_score_by_sample else NULL,
  pseudobulk_de_global   = global_de_results,
  pseudobulk_de_myeloid  = if (exists("myeloid_de_results")) myeloid_de_results else NULL,
  pseudobulk_skip_log    = skip_log,
  interpretation_summary = interpretation_summary,
  key_findings_summary   = key_findings_summary
)

saveRDS(stage4_output,
        file.path(PROC_DIR, "seurat_functional_comparison.rds"))
message("Saved: seurat_functional_comparison.rds")
message("  Contains: seurat_global, seurat_myeloid, composition_delta,",
        " module_score_summary, pseudobulk_de_global, pseudobulk_de_myeloid,",
        " pseudobulk_skip_log, interpretation_summary, key_findings_summary")

# --- Session info ---
writeLines(capture.output(sessionInfo()),
           file.path(TAB_DIR, "sessionInfo_04_APAP_Control_Functional_Comparison.txt"))
message("Saved: sessionInfo_04_APAP_Control_Functional_Comparison.txt")

# =========================================================================
# Completion
# =========================================================================

message("")
message("===================================================================")
message("=== 04_APAP_Control_Functional_Comparison.R complete ===")
message("===================================================================")
message("")
message("Output summary:")
message("  Object:     data/processed/04_APAP_Control_Functional_Comparison/seurat_functional_comparison.rds")
message("  Figures:    results/figures/04_APAP_Control_Functional_Comparison/ (",
        length(list.files(FIG_DIR, pattern = "\\.pdf$")), " PDFs)")
message("  Tables:     results/tables/04_APAP_Control_Functional_Comparison/ (",
        length(list.files(TAB_DIR, pattern = "\\.csv$")), " CSVs)")
message("")
message("Key outputs for README:")
message("  composition_delta_summary.csv          — Unified composition delta (global + myeloid)")
message("  stage4_interpretation_summary.csv       — Detailed evidence per source per group")
message("  stage4_key_findings_summary.csv          — One-row-per-group README-ready table")
message("  module_scores_per_cell.csv          — Per-cell module scores (global)")
message("  module_scores_by_sample.csv         — Per-sample aggregated module scores")
message("  pseudobulk_de_by_celltype.csv       — edgeR pseudobulk DE (global)")
message("  pseudobulk_de_by_myeloid_subtype.csv — edgeR pseudobulk DE (myeloid)")
message("  pseudobulk_skip_log.csv             — Skipped populations and reasons")
message("")
message("Limitations to communicate:")
message("  - n = 3 Control vs 3 APAP — all p-values are EXPLORATORY")
message("  - Pseudobulk DE uses edgeR QLF with TMM normalisation (within cell group)")
message("  - Module scores describe gene expression, NOT pathway activation")
message("  - No independent validation cohort")
message("  - All conclusions are hypothesis-generating, not confirmatory")
