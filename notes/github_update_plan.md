# GSE228305 AILI scRNA-seq Reconstruction Plan

本计划用于约束后续脚本实现。项目不是复刻一套标准 Seurat 教程，也不是提前规定结果必须落在某个故事上，而是围绕 APAP 诱导急性肝损伤中非实质细胞免疫微环境的变化，建立一条可解释、可复查、可调整的探索性单细胞分析流程。

由于当前只有一套 3 Control vs 3 APAP 数据，后续结果存在不确定性。本计划中的细胞群、marker 和通路均作为优先检验的候选方向，不作为必须得到的结论。


## 1. 目录与命名规范

```text
data/
  external/                         # 原始 GEO 矩阵，只读
  processed/
    01_QC_and_Integration/
    02_Global_Annotation/
    03_Myeloid_Subclustering/
    04_Pseudotime_Trajectory/
    05_Cell_Communication/

results/
  figures/
    01_QC_and_Integration/
    02_Global_Annotation/
    03_Myeloid_Subclustering/
    04_Pseudotime_Trajectory/
    05_Cell_Communication/
  tables/
    01_QC_and_Integration/
    02_Global_Annotation/
    03_Myeloid_Subclustering/
    04_Pseudotime_Trajectory/
    05_Cell_Communication/

scripts/
  01_QC_and_Integration.R
  02_Global_Annotation.R
  03_Myeloid_Subclustering.R
  04_Pseudotime_Trajectory.R
  05_Cell_Communication.R
```

推荐每个脚本开头固定包含：

- 包加载与版本提示。
- 输入文件存在性检查。
- 输出目录创建。
- 关键参数区，方便人工调整。
- 固定随机种子。

## 2. 分析原则：固定流程，开放结论

### 固定的部分

- 数据流向固定：`data/external/` -> `data/processed/` -> `results/`。
- 五阶段脚本结构固定，便于复现和审阅。
- 每阶段必须保存 `.rds` 中间对象。
- QC、整合、注释、聚类和可视化的基础评估必须完成。
- conda 环境固定为 `singlecell`。

### 开放的部分

- 不预设 APAP 组一定出现某类细胞大规模扩增。
- 不预设髓系细胞一定是最终主线；它只是最优先检查的候选方向。
- 不预设一定能分出清晰的促炎型和修复型巨噬细胞。
- 不预设 Monocle3 一定能形成可信轨迹。
- 不预设 CellChat 一定能给出特定通路故事。

### 分支与止损规则

- 如果全局注释后髓系细胞数量过少、marker 不清楚或组间变化不明显，Stage 3 应转向数据中更稳定的非实质细胞群，例如内皮细胞、T/NK 细胞、B 细胞或其他变化明显群体。
- 如果二次聚类只产生技术噪声或过度碎片化 cluster，应降低 resolution 或停止高分辨率拆分，改为粗粒度状态比较。
- 如果拟时序轨迹对 root、降维或亚群选择高度敏感，Stage 4 应降级为 marker/module score 沿状态变化的描述，不强行解释方向性演化。
- 如果 CellChat 结果主要由低细胞数亚群、少数离群基因或背景通路驱动，Stage 5 应只保留表格结果或作为补充，不画夸张网络图。
- README 只总结数据实际支持的发现；候选假设、阴性结果和不稳定分析可以在 notes 中诚实记录。

## 3. Conda 环境管理

本项目统一使用 conda 环境 `singlecell`。后续所有脚本运行、R 包安装、包版本检查和结果复现，都默认在该环境中完成。

### 执行规则

- 运行脚本前先核对环境：`conda activate singlecell`。
- 不要在 `base` 环境安装项目依赖。
- 不要为了某个脚本临时创建 `singlecell_test`、`scRNA` 等额外环境。
- 如果脚本报包缺失，先记录缺失包、所属阶段和用途，再补进 `singlecell`。
- 同一类功能优先选择稳定常用包，不为了“显得高级”额外堆依赖。

### 预计核心依赖

后续环境至少会围绕以下 R 包整理：

- 单细胞主体分析：`Seurat`、`SeuratObject`、`Matrix`
- 批次整合：`harmony`
- 可视化：`ggplot2`、`patchwork`、`cowplot`、`RColorBrewer`
- 数据整理：`dplyr`、`tidyr`、`readr`、`tibble`
- 富集分析：`clusterProfiler`、`org.Mm.eg.db`、`msigdbr` 或 `GSVA`
- 拟时序：`monocle3`
- 细胞通讯：`CellChat`

具体版本以实际成功运行为准，不在脚本完成前硬写死。

### environment.yml 生成原则

项目跑通并核对主要图表可复现后，在仓库根目录生成：

```text
environment.yml
```

推荐做法是先从真实环境导出，再人工精简掉平台缓存和无关包。最终文件需要能说明：

- 环境名为 `singlecell`。
- conda channels 清楚。
- R 版本和关键 R 包版本明确。
- 必要时把无法通过 conda 稳定安装的 R 包写入 README 的补充安装说明。

## 4. 样本设计

GSE228305 当前使用 6 个样本：

| group | sample_id | input_dir |
|---|---|---|
| Control | con1 | `data/external/GSM7118537_con1_matrix` |
| Control | con2 | `data/external/GSM7118537_con2_matrix` |
| Control | con3 | `data/external/GSM7118537_con3_matrix` |
| APAP | AP-300 | `data/external/GSM7118538_AP-300_matrix` |
| APAP | AP-300-2 | `data/external/GSM7118538_AP-300-2_matrix` |
| APAP | AP-300-3 | `data/external/GSM7118538_AP-300-3_matrix` |

脚本中应显式构建 sample metadata，至少包含：

- `sample_id`
- `group`
- `condition`
- `orig.ident`

不要从文件夹名临时猜组别，以免后续样本命名变化导致错误。

## 5. Stage 1: QC and Integration

**脚本**：`scripts/01_QC_and_Integration.R`

**输入**：

- `data/external/*_matrix/matrix.mtx`
- `data/external/*_matrix/barcodes.tsv`
- `data/external/*_matrix/genes.tsv`

**核心步骤**：

1. 读取 6 个 10x 风格矩阵，构建单样本 Seurat 对象。
2. 写入样本元信息并合并对象。
3. 计算 `percent.mt`、必要时计算 `percent.ribo`。
4. 先输出过滤前 QC 图，再根据图形观察设置过滤阈值。
5. 过滤低质量细胞和潜在 doublet 风险细胞。
6. 标准化、寻找高变基因、ScaleData、PCA。
7. 在整合前先输出 PCA/UMAP 按样本和组别着色图。
8. 运行 Harmony，输出整合后 UMAP，对比样本混合是否合理。

**人工决策点**：

- 根据 QC violin/scatter 判断 `nFeature_RNA`、`nCount_RNA`、`percent.mt` 阈值。
- Harmony 后检查：样本是否过度分离，或 Control/APAP 的真实差异是否被完全抹平。

**输出对象**：

- `data/processed/01_QC_and_Integration/seurat_integrated.rds`

**输出图表**：

- `results/figures/01_QC_and_Integration/qc_violin_before_filter.pdf`
- `results/figures/01_QC_and_Integration/qc_violin_after_filter.pdf`
- `results/figures/01_QC_and_Integration/umap_before_harmony_by_sample.pdf`
- `results/figures/01_QC_and_Integration/umap_after_harmony_by_sample.pdf`
- `results/tables/01_QC_and_Integration/sample_qc_summary.csv`

## 6. Stage 2: Global Annotation

**脚本**：`scripts/02_Global_Annotation.R`

**输入**：

- `data/processed/01_QC_and_Integration/seurat_integrated.rds`

**核心步骤**：

1. `readRDS()` 读取整合后的 Seurat 对象。
2. 基于 Harmony reduction 进行邻居图、聚类和 UMAP。
3. 输出不同 resolution 的聚类评估图，选择一个主分析 resolution。
4. 使用 `FindAllMarkers()` 提取每个 cluster 的 marker。
5. 根据经典 marker 与参考文献进行粗注释。
6. 计算 Control/APAP 组间细胞类型比例。

**人工决策点**：

- 结合 marker 表和参考文献 Figure 2/补充表核对大类细胞命名。
- 检查 APAP 组是否存在清晰的细胞组成变化。髓系细胞是优先关注对象，但如果变化不明显，应记录这一点并寻找更可靠的后续分析对象。

**输出对象**：

- `data/processed/02_Global_Annotation/seurat_annotated_global.rds`

**输出图表**：

- `results/figures/02_Global_Annotation/global_umap_by_cluster.pdf`
- `results/figures/02_Global_Annotation/global_umap_by_celltype.pdf`
- `results/figures/02_Global_Annotation/global_marker_dotplot.pdf`
- `results/figures/02_Global_Annotation/celltype_fraction_by_group.pdf`
- `results/tables/02_Global_Annotation/global_cluster_markers.csv`
- `results/tables/02_Global_Annotation/global_annotation_table.csv`
- `results/tables/02_Global_Annotation/celltype_fraction_by_sample.csv`

## 7. Stage 3: Focused Subclustering

**脚本**：`scripts/03_Myeloid_Subclustering.R`

脚本名暂时保留 `Myeloid`，因为髓系细胞仍是本项目的首选假设方向。但执行时必须先根据 Stage 2 结果判断髓系细胞是否适合作为深入对象；如果不适合，应在脚本参数区和输出文件中明确记录实际选择的目标细胞群。

**输入**：

- `data/processed/02_Global_Annotation/seurat_annotated_global.rds`

**核心步骤**：

1. `readRDS()` 读取全局注释对象。
2. 优先评估髓系细胞的数量、marker 清晰度和组间变化幅度。
3. 如果髓系适合深入，使用明确的 cell type 和 marker 双重逻辑提取髓系细胞，避免只按 cluster ID 生硬抽取。
4. 如果髓系不适合深入，选择 Stage 2 中更稳定、更有组间变化的非实质细胞群，并在表格中说明选择依据。
5. 对目标 subset 重新 Normalize/FindVariableFeatures/Scale/PCA/Harmony/UMAP/Cluster。
6. 输出多个 resolution 的 UMAP，选择高分辨率但不过度碎片化的结果。
7. 提取目标亚群 marker。
8. 如果目标为髓系细胞，可重点检查 Trem2、Mertk、Ly6c、S100a8/S100a9、Lyz2、Adgre1 等 marker；如果目标改变，应改用对应细胞类型的经典 marker。
9. 对核心亚群做 GO/GSVA，描述数据实际支持的功能方向。

**人工决策点**：

- resolution 的目标不是 cluster 越多越好，而是能稳定区分有生物学意义的状态。
- 若髓系细胞可用，再对照文献 Figure 3 和补充 marker，检查是否存在促炎/修复相关状态；如果不可分，应如实保留更粗粒度注释。

**输出对象**：

- `data/processed/03_Myeloid_Subclustering/seurat_myeloid_final.rds`

**输出图表**：

- `results/figures/03_Myeloid_Subclustering/myeloid_umap_by_subcluster.pdf`
- `results/figures/03_Myeloid_Subclustering/myeloid_umap_by_group.pdf`
- `results/figures/03_Myeloid_Subclustering/myeloid_marker_heatmap.pdf`
- `results/figures/03_Myeloid_Subclustering/myeloid_function_enrichment_dotplot.pdf`
- `results/tables/03_Myeloid_Subclustering/myeloid_subcluster_markers.csv`
- `results/tables/03_Myeloid_Subclustering/myeloid_annotation_table.csv`
- `results/tables/03_Myeloid_Subclustering/myeloid_enrichment_results.csv`

如果最终目标细胞群不是髓系，文件名可暂时保留以维持五阶段结构，但表格中必须用 `target_celltype` 或类似字段说明真实分析对象。

## 8. Stage 4: Pseudotime Trajectory

**脚本**：`scripts/04_Pseudotime_Trajectory.R`

**输入**：

- `data/processed/03_Myeloid_Subclustering/seurat_myeloid_final.rds`

**核心步骤**：

1. `readRDS()` 读取髓系细胞对象。
2. 先判断目标细胞群是否存在连续状态变化：UMAP 是否呈连续结构，marker 是否呈梯度，亚群是否能形成合理过渡。
3. 若适合轨迹分析，转换为 Monocle3 `cell_data_set`。
4. 继承 Seurat 中的 UMAP 坐标和亚群注释，减少重复降维带来的解释偏移。
5. `learn_graph()` 学习轨迹结构。
6. 基于生物学假设人工指定 root cells。
7. 绘制伪时间轨迹和关键基因随伪时间变化曲线。
8. 若轨迹结构不稳定，输出状态评分或 marker 趋势图作为替代，不强行解释伪时间方向。

**人工决策点**：

- root 不能由代码自动决定。只有在轨迹结构稳定时才指定 root；候选 root 应来自明确的生物学假设和 marker 证据。
- 如果不同 root 或参数导致完全不同解释，应停止使用强方向性叙事。

**输出对象**：

- `data/processed/04_Pseudotime_Trajectory/monocle_cds.rds`

**输出图表**：

- `results/figures/04_Pseudotime_Trajectory/myeloid_pseudotime_umap.pdf`
- `results/figures/04_Pseudotime_Trajectory/myeloid_trajectory_by_state.pdf`
- `results/figures/04_Pseudotime_Trajectory/key_gene_pseudotime_trends.pdf`
- `results/tables/04_Pseudotime_Trajectory/pseudotime_cell_metadata.csv`
- `results/tables/04_Pseudotime_Trajectory/pseudotime_genes.csv`

## 9. Stage 5: Cell Communication

**脚本**：`scripts/05_Cell_Communication.R`

**输入**：

- `data/processed/03_Myeloid_Subclustering/seurat_myeloid_final.rds`
- 必要时可读取 `data/processed/02_Global_Annotation/seurat_annotated_global.rds` 用于非髓系互作背景。

**核心步骤**：

1. `readRDS()` 读取最终注释对象。
2. 构建 CellChat 对象，按 group 或 condition 分别分析 Control/APAP。
3. 使用小鼠配体-受体数据库。
4. 计算通讯概率、过滤低细胞数亚群导致的不稳定互作。
5. 优先审查炎症、趋化、吞噬、细胞迁移和修复相关通路，但不预设一定出现。
6. 输出精简图，而不是展示所有通路。
7. 若通讯结果不稳定或缺乏清楚解释，只输出筛选表格或作为补充结果。

**人工决策点**：

- CellChat 会产生大量背景信号。最终展示应围绕明确问题：哪些细胞发送信号，哪些细胞接收信号，这些互作是否被样本和表达证据支持。
- 可重点审查 CSF、CCL/CXCL、GAS、TGFb、SPP1、MIF 等通路，但最终保留哪些必须由结果和文献共同决定；没有稳定信号时不要硬画故事图。

**输出对象**：

- `data/processed/05_Cell_Communication/cellchat_myeloid.rds`

**输出图表**：

- `results/figures/05_Cell_Communication/cellchat_interaction_count_circle.pdf`
- `results/figures/05_Cell_Communication/cellchat_selected_pathways.pdf`
- `results/figures/05_Cell_Communication/cellchat_apap_vs_control_signaling.pdf`
- `results/tables/05_Cell_Communication/cellchat_significant_interactions.csv`
- `results/tables/05_Cell_Communication/cellchat_selected_pathways.csv`

## 10. 代码风格要求

- 优先写清楚生物学意图，再写参数；注释要像实验记录，不像自动生成说明。
- 所有关键参数集中放在脚本开头，便于后续根据图重新调整。
- 对输入文件、输出目录、关键 metadata 列做显式检查。
- 图表命名稳定，避免后续 README 引用失效。
- 所有 R 脚本默认在 `singlecell` 环境中运行，不在脚本中自动安装包。
- 代码和注释中避免“证明”“确认”“完美区分”等结论先行词汇，优先使用“评估”“检查”“支持/不支持”。
- 不修改 `data/external/`。
- 不批量删除文件；需要清理时逐个明确路径处理。

## 11. 后续 README 展示逻辑

README 最终建议按科研故事组织，而不是按代码功能罗列：

1. 数据集与科学问题。
2. 分析流程图。
3. QC 与整合质量控制。
4. 全局非实质细胞图谱。
5. 数据实际支持的重点细胞群高分辨率分析。
6. 若成立，展示状态连续性或拟时序证据；若不成立，展示替代的状态评分/通路分析。
7. 若成立，展示细胞通讯支持的微环境调控模型；若不成立，将 CellChat 作为补充或省略。
8. 可复现性说明与运行方式。
