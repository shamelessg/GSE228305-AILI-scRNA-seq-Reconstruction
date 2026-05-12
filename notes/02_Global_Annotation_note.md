# 第二阶段 Global Annotation 分析笔记

**日期**：2026-05-12
**脚本**：`scripts/02_Global_Annotation.R`
**版本**：第二版（人工注释修正）

---

## 1. 输入与输出

**输入**：
- `data/processed/01_QC_and_Integration/seurat_integrated.rds`
- 88,060 cells × 24,873 genes
- 6 samples (3 Control + 3 APAP)
- Harmony-integrated, PCA + UMAP 已存在

**输出**：
- `data/processed/02_Global_Annotation/seurat_annotated.rds` — 带 celltype_manual 注释的 Seurat 对象
- `results/figures/02_Global_Annotation/` — UMAP、DotPlot、FeaturePlot、比例图
- `results/tables/02_Global_Annotation/` — marker 表、细胞计数/比例表、Wilcoxon 结果、myeloid assessment

---

## 2. Resolution 选择

测试了 4 个 resolution：

| Resolution | Clusters |
|---|---|
| 0.2 | 21 |
| **0.4** | **27** |
| 0.6 | 34 |
| 0.8 | 39 |

选择 **resolution = 0.4**，理由：

- 目标是粗粒度全局注释，不是精细亚群分析
- 0.2 只有 21 个 cluster，明显合并过度（例如将 T / NK / B 混在一起）
- 0.6 有 34 个 cluster，对全局层面来说偏碎片化
- 0.4 在 27 个 cluster 时可以较好地区分主要 lineage：T、B、NK、髓系（Kupffer / monocyte-derived / neutrophil）、内皮、上皮等

**注意**：resolution = 0.4 是人工选择，不是算法最优。如果后续 DotPlot / FeaturePlot 发现某些 cluster 内明显混合了多种细胞类型（over-merging），需要回到 0.6 复核。

---

## 3. 人工注释依据

每个 cluster 的注释基于 presto::wilcoxauc() 输出的 cluster-enriched genes，结合已知小鼠肝脏 NPC marker 文献判断。以下是每个 cluster 的标注和关键 marker 证据：

### Lymphoid

| Cluster | Label | Key markers |
|---|---|---|
| 0 | T cell | Lef1, Tcf7, Il7r, Trac, Cd28 → naive/conventional T |
| 1 | B cell | Pax5, Ms4a1, Cd79a, Cd79b, Fcer2a, Bank1 |
| 3 | NK cell | Ncr1, Klra9, Klrb1b, Gzma, Prf1, Gzmb, Xcl1 |
| 10 | Cytotoxic T / NK-like | Gzmk, Ifit1, Ifit3, Mx1, Cd8a, Cd8b1, Ccl5；IFN-responsive |
| 13 | B cell | 同 cluster 1，Pax5, Ms4a1, Cd79a, Cd79b, Fcer2a, Bank1 |
| 16 | Plasma cell | Jchain, Tnfrsf17, Igha, Ighg, Mzb1/Xbp1 |

### Myeloid

| Cluster | Label | Key markers |
|---|---|---|
| 2 | Inflammatory neutrophil | Cxcr2, Il1f9, Acod1, Ptgs2, Cxcl2, S100a8/a9 |
| 8 | Monocyte-derived Mph/DC-like | Clec4b1, Ms4a6c, Cd209a, Clec4a, Cx3cr1, Lyz2, Fcgr1 |
| 11 | pDC-like | Siglech, Ccr9, Flt3, Cd300c |
| 14 | Kupffer cell | Cd5l, Clec4f, Folr2, Fcna, Cd163, Vsig4, C6 |
| 15 | cDC1 | Xcr1, Clec9a, Flt3, Itgax；注意不是 monocyte |
| 17 | Activated neutrophil | Il1r2, Trem1, Csf3r, Cxcr2, Il1b |
| 18 | Neutrophil | Ltf, Camp, Ngp, Cd177, Ly6g |
| 22 | Mast / basophil-like | Cpa3, Gata2, Kit, Prtn3；不是 monocyte |
| 26 | pDC-like | 同 cluster 11，Siglech, Ccr9, Flt3, Cd300c |

### Stromal / Structural

| Cluster | Label | Key markers |
|---|---|---|
| 4 | Cholangiocyte / epithelial | Epcam, Krt8, Krt18, Krt19, Sox9, Dmbt1 |
| 5 | Endothelial | Robo4, Lyve1, Mmrn2, Ptprb, Kdr, Cdh5, Stab2, Pecam1 |
| 9 | Mesothelial | Msln, Muc16, Upk3b, Rspo1 |
| 19 | VSMC / pericyte | Lmod1, Cnn1, Myh11, Tagln |
| 21 | Fibroblast / stellate | Dpep1, Dpt, Gdf10, Mmp2, Tcf21, Pdgfra, Lum |
| 23 | Cholangiocyte / epithelial | 同 cluster 4 |
| 25 | Endothelial | Thbd, Nos3, Lyve1, Stab2, Cdh5, Kdr, Pecam1 |

### Other

| Cluster | Label | Key markers |
|---|---|---|
| 6 | Hepatocyte contamination | Alb, Apoa1, Ttr, Cyp genes, Hnf4a |
| 7 | Low-quality / mixed | mt-Rnr2, mt-Rnr1, Gm42418, Malat1, Cst3 |
| 12 | Cycling cell | Mki67, Top2a, Stmn1, Hist1h3c, Pclaf, Ccna2 |
| 20 | B cell / possible doublet | Ms4a1, Pax5, Fcer2a + Alb/Apoa1/Ttr background |
| 24 | Hepatocyte contamination | 同 cluster 6 |

---

## 4. 不确定 / 需复核的 cluster

以下 cluster 需要进一步验证，标注不能视为最终结论：

### Cluster 7 — Low-quality / mixed
- 高表达 mt-Rnr2, mt-Rnr1, Gm42418, Malat1, Cst3
- 不是干净的生物学细胞类型，可能包含 stressed/dying cells 或 empty droplets
- 建议检查该 cluster 的 mt% 和 nFeature 分布

### Cluster 10 — Cytotoxic T / NK-like
- 同时表达 CD8 T 和 NK 相关基因（Gzmk, Cd8a, Ccl5 + Ifit1/3 IFN 信号）
- 可能是 CD8+ effector T 或 NK-like CTL，resolution 0.4 下未完全分开
- 后续若需区分 CD8 T vs NK，可对此 cluster 做 subclustering 或回到 0.6

### Cluster 20 — B cell / possible doublet
- B cell markers (Ms4a1, Pax5, Fcer2a, Ighd, Bank1, Cd79a) 明确
- 但同时有 Alb/Apoa1/Ttr 背景，可能是 hepatocyte-B cell doublet 或 ambient RNA 污染
- 建议后续考虑 SoupX 或 CellBender 去 ambient RNA 后再评估

### Cluster 22 — Mast / basophil-like
- Cpa3, Gata2, Kit, Prtn3 指向 mast cell 或 basophil
- 细胞数较少，需确认 Cpa3/Kit 在 FeaturePlot 上的共表达模式
- 如果确实为 mast cell，第三阶段髓系分析时需决定是否纳入

---

## 5. 细胞比例变化（探索性）

n = 3 Control vs 3 APAP，Wilcoxon 秩和检验 p 值仅作探索性参考，**不能作为显著差异结论**。

粗看趋势：
- **Inflammatory neutrophil** (cluster 2)：APAP 组比例均数从 4.2% → 9.8%，p = 0.38，趋势上升但不显著（样本量太小）
- **Activated neutrophil** (cluster 17)：APAP 组略升
- **Hepatocyte contamination** (clusters 6+24)：APAP 组明显降低（8.2% → 2.4%），但这更可能反映 APAP 处理后肝细胞损伤/丢失，而非真实的 NPC 比例变化
- **Endothelial**：APAP 组下降趋势
- T cell 和 B cell 比例变化不大

**再次强调**：n = 3 vs 3，这些 p 值不能支持任何统计推断。需要更大样本量或独立的验证队列。

---

## 6. 第三阶段建议

### 优先进入 subclustering 的 lineage：髓系

理由：
- 髓系在 APAP 肝损伤中是最核心的效应细胞群体
- 本阶段已识别多个髓系亚群：Kupffer cell (cluster 14)、Monocyte-derived macrophage/DC-like (cluster 8)、Inflammatory neutrophil (cluster 2)、Neutrophil (cluster 18)、Activated neutrophil (cluster 17)
- 比例数据提示 neutrophil/monocyte-derived 群体在 APAP 下有上升趋势

### subclustering 前需要先排除的群体：
- **pDC-like** (cluster 11, 26)：虽是髓系来源但功能上属于 pDC lineage，混入髓系对象会引入噪声
- **cDC1** (cluster 15)：同理，属于 DC lineage
- **Mast / basophil-like** (cluster 22)：需先确认身份，若确认则单独处理
- **Endothelial** (cluster 5, 25)：确保不因 UMAP 邻近而被误混入

### 推荐髓系 subclustering 对象包括：
- Kupffer cell (14)
- Monocyte-derived macrophage / DC-like (8)
- Inflammatory neutrophil (2)
- Neutrophil (18)
- Activated neutrophil (17)

共 5 个 cluster，细胞数如下：

| Cluster | Label | Cell count |
|---|---|---|
| 14 | Kupffer cell | 1,501 |
| 8 | Monocyte-derived macrophage / DC-like | 3,266 |
| 2 | Inflammatory neutrophil | 7,608 |
| 18 | Neutrophil | 585 |
| 17 | Activated neutrophil | 596 |
| **Total** | | **13,556** |

合计 13,556 cells，占总细胞数 15.4%，细胞量充足，适合 subclustering。

---

## 7. 声明

本阶段注释是**人工初版**，基于 resolution = 0.4 的聚类结果和已知 marker 基因手工标注。

后续可能在以下情况下修订：
- 更精细的 marker 文献复查
- subclustering 揭示内部异质性
- 去 ambient RNA 后的 cluster 重新分配
- 与公共数据集（如 Tabula Muris、Liver Cell Atlas）的 label transfer 比对

**此版本不作为最终发表级别的细胞注释。**
