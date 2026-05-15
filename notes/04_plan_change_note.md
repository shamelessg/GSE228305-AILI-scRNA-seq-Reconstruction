# 计划变更笔记：第四阶段从拟时序改为 APAP-Control 功能比较

**日期**：2026-05-14

**变更结论**：第四阶段不再以 Monocle3 pseudotime 作为主线，改为 `04_APAP_Control_Functional_Comparison`。拟时序保留为可选补充分析，仅在后续出现清楚连续状态结构时再做。

---

## 1. 为什么变更

第三阶段已经完成髓系二次聚类，并产出按样本的亚群比例图和髓系 UMAP：

- `results/figures/03_Myeloid_Subclustering/myeloid_subcluster_proportion_stacked.pdf`
- `results/figures/03_Myeloid_Subclustering/myeloid_umap_by_subcluster.pdf`

比例图已经满足"以样本为单位检查 APAP-Control 组成变化"的需求。APAP 组中 inflammatory neutrophil 比例更高，Resident Kupffer cell 和部分 antigen-presenting myeloid/macrophage 亚群相对更低。这是当前最稳妥的生物学主线。

UMAP 结构不支持把全部髓系细胞放入单一拟时序：中性粒、Kupffer cell、Mo-Mac/DC-like cells 形成分离区域，更多反映 lineage 和功能状态差异，而不是一条连续发育或激活轨迹。中性粒内部存在局部状态邻接，但目前不足以支撑强方向性的 Monocle3 叙事。

---

## 2. 新的科学问题

原问题：

> APAP 后髓系细胞是否存在可解释的连续状态轨迹？

调整为：

> APAP 损伤组织中，哪些非实质细胞或髓系亚群发生组成和功能状态改变？这些改变更像损伤放大、炎症募集、抗原呈递、稳态维持丧失，还是修复/清除反应？

这个问题更贴近当前数据实际支持的证据，也更适合 3 Control vs 3 APAP 的设计。

---

## 3. 第四阶段新目标

第四阶段命名：

```text
04_APAP_Control_Functional_Comparison
```

脚本建议命名：

```text
scripts/04_APAP_Control_Functional_Comparison.R
```

核心输出应回答三件事：

1. **Composition**：APAP vs Control 中哪些细胞类型/髓系亚群比例变化最稳。
2. **Expression**：在主要细胞群内，APAP 相比 Control 上调或下调哪些基因程序。
3. **Function**：这些变化对应炎症、趋化、抗原呈递、Kupffer 稳态、IFN-response、损伤修复等哪类生物学功能。

---

## 4. 建议分析模块

### 4.1 组成变化

沿用并整理 Stage 2/3 已有的按样本比例表：

- 全局 cell type 层面：`celltype_proportions_by_sample.csv`
- 髓系 subtype 层面：`myeloid_subcluster_proportions_by_sample.csv`

输出更适合 README 的 summary 表，包含：

- `celltype_or_subtype`
- `mean_control`
- `mean_apap`
- `delta_apap_minus_control`
- `fold_change_apap_vs_control`
- `direction`
- `interpretation_note`

### 4.2 Pseudobulk 差异表达

优先做 sample-level pseudobulk，而不是直接把所有细胞当独立重复。

建议对象：

- 全局主要细胞类型：Inflammatory neutrophil、Kupffer cell、Monocyte-derived macrophage/DC-like、Endothelial、T cell、B cell、NK cell。
- 髓系主要亚群：Inflammatory neutrophil、Activated neutrophil、Resident Kupffer cell、Antigen-presenting macrophage / Mo-Mac、DC-like antigen-presenting myeloid。

每个细胞群内按 `sample_id` 聚合 counts 或 normalized expression，再比较 APAP vs Control。若某细胞群每组有效样本太少或细胞数太低，应跳过并记录原因。

### 4.3 功能评分和通路解释

优先使用清晰可解释的小型基因集：

- Neutrophil inflammation / activation
- Chemokine recruitment
- IL1/TNF inflammatory response
- IFN response
- Antigen presentation
- Kupffer resident identity
- Phagocytosis / efferocytosis
- Tissue repair / ECM remodeling
- Endothelial activation

输出 per-cell 和 sample-aggregated score，避免只展示合并细胞层面的漂亮图。

---

## 5. 对"好/坏作用"的表述原则

不要直接把某类细胞命名为"好细胞"或"坏细胞"。更稳妥的表达是：

- 可能促损伤：炎症因子、趋化因子、ROS/脱颗粒、组织降解相关程序增强。
- 可能保护或修复：吞噬清除、efferocytosis、Kupffer 稳态、抗炎/修复、ECM remodeling 相关程序增强。
- 可能免疫调节：抗原呈递、IFN-response、T/NK activation 相关程序改变。

最终 README 中只写数据直接支持的功能状态，不把候选机制写成定论。

---

## 6. 第五阶段如何衔接

第五阶段 CellChat 仍保留，但它应基于第四阶段筛出的稳定细胞和功能轴，而不是全量跑完后挑故事。

优先互作方向：

- Inflammatory neutrophil / macrophage 与 endothelial 的趋化和黏附信号。
- Kupffer/macrophage 的补体、吞噬和损伤清除相关信号。
- Antigen-presenting myeloid 与 T/NK cell 的抗原呈递或共刺激信号。
- Stellate/fibroblast/endothelial 的修复、ECM 或血管反应信号。

如果 CellChat 结果主要由低细胞数群体或背景通路驱动，应作为补充表格，不进入主图。
