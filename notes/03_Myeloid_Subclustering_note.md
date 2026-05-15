# 第三阶段 Myeloid Subclustering 分析笔记

**日期**：2026-05-13
**脚本**：`scripts/03_Myeloid_Subclustering.R`
**版本**：人工 working annotation（16 clusters，基于 res=0.4）

---

## 1. 本阶段为什么选择髓系细胞

髓系细胞是 APAP 急性肝损伤中最早响应且最核心的效应细胞群体。
文献表明，APAP 过量后 Kupffer cell 来源的促炎因子驱动初始损伤，
随后 monocyte-derived macrophage 浸润并参与炎症消退和组织修复。
中性粒细胞在早期炎症中也发挥重要作用。

Stage 2 全局注释识别了 5 个髓系相关 cluster：Kupffer cell、
Monocyte-derived macrophage/DC-like、Inflammatory neutrophil、
Neutrophil、Activated neutrophil，合计 13,556 cells，
数量充足且 marker 表达清晰，具备 subclustering 条件。

本阶段目标：在髓系内部发现更精细的亚群结构，
评估 Control vs APAP 的组成和功能状态变化趋势。

---

## 2. 输入对象、纳入和排除的细胞类型

**输入**：`data/processed/02_Global_Annotation/seurat_annotated.rds`
- 13556 myeloid cells x 24873 genes
- 6 samples（3 Control + 3 APAP）

**纳入的细胞类型**：
- Inflammatory neutrophil：7608 cells
- Monocyte-derived macrophage / DC-like：3266 cells
- Kupffer cell：1501 cells
- Activated neutrophil：596 cells
- Neutrophil：585 cells

**排除的细胞类型**：
- pDC-like（虽为髓系来源但属于 pDC lineage）
- cDC1（属于 DC lineage）
- Mast / basophil-like（细胞数少，身份需确认）
- Endothelial（肝脏结构细胞）
- Low-quality / mixed（质量控制排除）
- B cell / possible doublet
- Hepatocyte contamination
- Cholangiocyte / epithelial

---

## 3. 重新聚类使用的参数

- NormalizeData：LogNormalize，scale factor = 10,000
- FindVariableFeatures：vst，nfeatures = 2,000
- RunPCA：npcs = 30
- Harmony：group.by.vars = "sample_id"，dims = 1:25（成功）
- RunUMAP：reduction = "harmony"，dims = 1:25
- FindNeighbors：reduction = "harmony"，dims = 1:25
- FindClusters：resolution = 0.2, 0.4, 0.6, 0.8
- presto::wilcoxauc() 用于 marker 检测（padj < 0.05, logFC > 0.25）

---

## 4. Resolution 选择依据

测试了 4 个 resolution：
- 0.2：10 clusters
- 0.4：16 clusters
- 0.6：20 clusters
- 0.8：21 clusters

最终选择 **resolution = 0.4**（16 clusters）。

resolution = 0.4 在髓系内部提供了适中的粒度，
能够区分 Kupffer cell、monocyte-derived macrophage、
neutrophil 的不同激活状态和功能亚群，
同时避免 resolution = 0.6/0.8 的过度碎片化。

---

## 5. 每个亚群的 marker 和人工命名（working annotation）

**重要声明：以下标签为 working annotation，不是最终发表级别的细胞注释。**
注释基于 DotPlot、FeaturePlot、top marker 和文献中已知的髓系 marker 综合判断。

当前 16 个 cluster 的人工注释：
- Cluster 0：Inflammatory neutrophil
- Cluster 1：Resident Kupffer cell
- Cluster 2：Antigen-presenting macrophage / Mo-Mac
- Cluster 3：DC-like antigen-presenting myeloid
- Cluster 4：Inflammatory neutrophil
- Cluster 5：Activated neutrophil
- Cluster 6：Mature neutrophil
- Cluster 7：Monocyte-derived macrophage
- Cluster 8：Activated inflammatory neutrophil
- Cluster 9：Possible doublet / hepatocyte ambient-contaminated neutrophil
- Cluster 10：IFN-responsive inflammatory neutrophil
- Cluster 11：Neutrophil with erythroid/ambient signal
- Cluster 12：Antigen-presenting myeloid with lymphoid-like signal
- Cluster 13：Cycling neutrophil
- Cluster 14：CCR7+ DC-like antigen-presenting myeloid
- Cluster 15：Rare injury-associated macrophage-like cells

**重点人工修正（与自动 marker 评分不同的 cluster）**：
- **Cluster 1**：自动评分判为 Monocyte_derived_Mph，人工复核发现高表达 Clec4f/Vsig4/Timd4，
  修正为 Resident Kupffer cell。Kupffer 标记物在此 cluster 的表达模式与 cluster 7 (真正的 Mo-Mac) 明显不同。
- **Cluster 9**：自动评分判为 Neutrophil，人工复核发现同时携带高 Alb/Apoa1 环境信号 + S100a8/S100a9，
  修正为 Possible doublet / hepatocyte ambient-contaminated neutrophil。
  这些细胞可能并非真实的生物学亚群，需在后续分析中审慎对待。
- **Cluster 13**：自动评分判为 Neutrophil，人工复核发现 Mki67/Top2a/Stmn1 阳性，
  修正为 Cycling neutrophil。代表了处于增殖状态的中性粒细胞前体或活化群体。

详细 marker 和评分见：
- `results/tables/03_Myeloid_Subclustering/myeloid_top_markers_per_cluster.csv`
- `results/tables/03_Myeloid_Subclustering/myeloid_preliminary_subcluster_annotation.csv`
- `results/tables/03_Myeloid_Subclustering/myeloid_cluster_markers_presto.csv`（完整 presto 输出含 AUC）
- `results/figures/03_Myeloid_Subclustering/myeloid_dotplot_canonical_markers.pdf`
- `results/figures/03_Myeloid_Subclustering/myeloid_featureplot_key_markers_batch*.pdf`

---

## 6. 哪些亚群最不确定，需要人工复核

以下 cluster 的身份不确定性较高，建议后续实验验证或谨慎解释：

- **Cluster 9**（Possible doublet / hepatocyte ambient-contaminated neutrophil）：
  高 contamination_score，可能为实验 artifact 或双细胞。不建议用于生物学结论。
- **Cluster 11**（Neutrophil with erythroid/ambient signal）：
  可能携带红细胞或环境 RNA 污染信号，需确认其 marker 表达是否独立于 contamination。
- **Cluster 12**（Antigen-presenting myeloid with lymphoid-like signal）：
  同时表达髓系和淋系 marker，需排除 doublet 可能性。
- **Cluster 15**（Rare injury-associated macrophage-like cells）：
  细胞数极少的稀有群体，可能为真实稀有亚群或 outlier。

一般注意事项：
- 如果一个 cluster 同时高表达多个 lineage marker，可能为 doublet 或过渡状态
- 如果 contamination_score 较高，需检查是否为实验 artifact
- 如果某个 cluster 细胞数很少（< 50 cells），可能为 outlier cluster
- 如果某 cluster 仅在单个 sample 中出现，可能为 sample-specific batch effect

---

## 7. Control vs APAP 的比例变化，谨慎描述

**重要声明**：n = 3 Control vs 3 APAP，Wilcoxon 秩和检验 p 值仅作探索性参考，
**不能作为强统计结论**。需要更大样本量或独立的验证队列。

分别输出了按 myeloid_subtype 和按 myeloid_clusters 的 Wilcoxon 结果：
- `myeloid_subcluster_proportion_wilcox.csv`
- `myeloid_cluster_proportion_wilcox.csv`
- `myeloid_cluster_counts_by_sample.csv`
- `myeloid_cluster_proportions_by_sample.csv`

主要趋势（待人工根据实际输出填写）：
- **Inflammatory neutrophil** 是最清楚的 APAP 上升趋势：Control 均值 0.3114，APAP 均值 0.4036。
- **Mature neutrophil**、**Activated neutrophil**、**IFN-responsive inflammatory neutrophil** 和 **Activated inflammatory neutrophil** 也整体偏向 APAP 上升，但幅度较小。
- **Resident Kupffer cell** 在 APAP 中相对下降：Control 均值 0.1550，APAP 均值 0.1049。
- **Antigen-presenting macrophage / Mo-Mac** 和 **DC-like antigen-presenting myeloid** 在 APAP 中相对下降或持平，提示 APAP 后髓系组成更偏向炎症性中性粒细胞，而不是单纯抗原呈递髓系扩增。
- **Possible doublet / hepatocyte ambient-contaminated neutrophil** 在 APAP 中下降，不应作为生物学主线解释。

这些趋势与 `myeloid_subcluster_proportion_stacked.pdf` 一致：按样本展示已经足够说明组成变化来自多个样本的方向性差异，而不是只来自合并细胞后的假象。但 n = 3 vs 3 的样本量有限，p 值不能作为强统计判断。

---

## 8. Module score 的趋势，谨慎描述

计算了 8 个基因集的 AddModuleScore。
Module score 反映的是基因集在单个细胞中的平均表达水平，
不代表通路激活状态，不能作为机制证明。

输出文件：
- `myeloid_module_scores_per_cell.csv`：每个细胞的 module score（per-cell level）
- `myeloid_module_scores_aggregated.csv`：按 sample_id/group/myeloid_clusters 聚合的均值表

主要观察（待人工根据实际输出填写）：
- 中性粒细胞相关 cluster 的 `Neutrophil_score` 和 `Activated_neutrophil_score` 较高，符合其 S100a8/S100a9、Retnlg、Il1b、Cxcl2 等 marker 特征。
- Resident Kupffer cell 的 `Kupffer_resident_score` 明显较高，并伴随 Clec4f、Vsig4、Cd5l、C1qa/b/c 等经典 Kupffer marker，支持人工修正。
- Antigen-presenting macrophage / Mo-Mac 与 DC-like antigen-presenting myeloid 的 `Antigen_presentation_score` 较高，提示这部分细胞更适合放入“抗原呈递/免疫调节”功能模块，而不是简单归为促炎细胞。
- 当前 module score 更适合作为功能状态描述：炎症募集、抗原呈递、Kupffer 稳态、IFN-response、损伤应答。它不足以直接证明某类细胞“好”或“坏”。

---

## 9. 是否建议第四阶段进入 pseudotime

**当前决策：不将 Monocle3 pseudotime 作为第四阶段主线。**

原因：`myeloid_umap_by_subcluster.pdf` 显示髓系结构主要分成多个离散岛：
左侧为中性粒细胞相关状态，右侧为 Mo-Mac/DC-like 抗原呈递细胞，上方为 Resident Kupffer cell。
这些结构反映的是不同 lineage/功能状态的分离，而不是一条干净的连续发育路径。
中性粒细胞内部存在局部邻接和激活梯度，但不足以支撑全髓系 Monocle3 方向性叙事。

替代方案：
- 第四阶段改为 **APAP-Control functional comparison**：围绕每个主要细胞群/髓系亚群做组成变化、pseudobulk 差异表达、marker/module score 和通路解释。
- 使用 module score（Inflammatory/Activated/IFN-response/Antigen presentation/Kupffer resident）在 cluster 和 group 间比较功能状态趋势。
- 使用 top marker 的表达模式描述激活/损伤应答梯度，但不写成伪时间。

**仅在以下情况下启动 Monocle3**：
- UMAP 上中性粒细胞各亚群呈现明显的连续弧线结构
- 或 Kupffer→Mo-Mac 之间存在清晰的过渡路径
- 选择 root 需要明确的生物学依据（如选择 resident Kupffer 作为起始点）

如果后续补充 pseudotime，建议仅对中性粒细胞 lineage 或 Kupffer/Mo-Mac lineage
分别独立试跑，作为补充分析，而非项目主线。

---

## 10. 声明

本阶段注释基于自动 marker 富集评分作为参考，经人工复核后应用 working annotation。
标签不代表最终发表级别的细胞注释，部分 cluster 的身份需要在独立数据集中验证。
module score 和 Wilcoxon 结果仅为探索性描述，不是机制证明。
Cluster 1/9/13 为重点人工修正，与自动评分结果不同，已在第 5 节详细说明了修正依据。

**此版本为 working annotation，不视为最终发表级别的细胞注释。**
