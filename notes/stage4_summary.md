# Stage 4: APAP vs Control 功能比较

**目的**：回答 APAP 急性肝损伤 24h 后，肝脏非实质细胞（尤其是髓系）的组成和功能状态变化。

**背景**：原计划 Monocle3 拟时序，但髓系 UMAP 呈分离 lineage 结构（中性粒、Kupffer、Mo-Mac/DC 各自成簇），不支持单一连续轨迹，改为组间比较。

---

## 分析管线

### Module 1 — 组成变化
- 方法：从 Stage 2 (global celltype) 和 Stage 3 (myeloid subtype) 的比例表合并，计算 APAP vs Control 的比例差值和 fold change
- 工具：dplyr / tidyr / ggplot2
- 产出：`composition_delta_summary.csv`（37 行，人+ 髓系合并），火山图

### Module 2 — 功能模块打分
- 12 个功能基因集：Neutrophil_inflammation / Chemokine_recruitment / IL1_TNF_inflammatory / IFN_response / Antigen_presentation / Kupffer_resident_identity / Phagocytosis_efferocytosis / Tissue_repair_ECM / Endothelial_activation / ROS_degranulation / T_NK_activation / Anti_inflammatory_repair
- 方法：Seurat AddModuleScore → 样本级聚合 → APAP vs Control 比较
- 产出：per-cell 分数 × 2（全局 + 髓系），per-sample 聚合分数 × 2，3 个 figure

### Module 3 — edgeR 伪批量差异表达（核心改动）
- 之前用 Welch t-test + 全局 CPM 归一化，现改为 edgeR 标准管线
- 每个 cell group 独立构建 gene × sample count 矩阵
- TMM 归一化在 cell group 内计算，不混入其他类型的 library size
- filterByExpr → calcNormFactors → estimateDisp → glmQLFit → glmQLFTest
- 输出：logFC, logCPM, PValue, **FDR**
- 覆盖 11 个 global cell types + 9 个 myeloid subtypes，跳过 2 个（Control 样本不足）
- 产出：`pseudobulk_de_by_celltype.csv`（93,738 行）, `pseudobulk_de_by_myeloid_subtype.csv`（49,828 行）, 火山图

### Module 4 — 证据汇总
- 每细胞群收集三方向证据：composition / module_score / pseudobulk_DE
- confidence 仅取决于 composition + DE 方向一致性（module score 只做功能描述）
- 两张表：
  - `stage4_interpretation_summary.csv` — 逐证据源明细（51 行）
  - `stage4_key_findings_summary.csv` — 每群一行，README 主结论表（17 行）

---

## 核心结果

| Cell group | 组成变化 | n DEGs (FDR<0.05) | 功能特征 | 置信度 |
|---|---|---|---|---|
| Inflammatory neutrophil | 4.2% → 9.8% | 0 | Chemokine ↑, IFN ↓ | Low |
| Activated neutrophil | 0.4% → 1.0% | 12 | Chemokine ↑, Neutrophil_inflammation ↑ | **Moderate** |
| Resident Kupffer cell | 15.5% → 10.5% | 430 | 炎症模块 ↑, 噬菌/修复模块 ↓ | Uncertain |
| Endothelial | 7.1% → 3.8% | 297 | 炎症模块 ↑, 激活模块高变异 | Uncertain |
| Mo-Mac/DC-like | 3.5% → 6.2%（稳定） | 161 | AP 程序维持, IL1/TNF 低于中性粒 | Uncertain |
| Kupffer cell (global) | 稳定 | **420**（最多） | 全局转录重编程 | High |

---

## 生物学解读

**1. 中性粒细胞扩张是最稳健的观察结果。** 三个中性粒状态（Inflammatory / Activated / Mature）在 APAP 后比例均升高。Inflammatory neutrophil 无显著 DEG——功能变化可能更多体现在组成扩增而非群体内转录重塑，也不排除 n=3 对于该群体检验效力不足。

**2. Kupffer 细胞不能简单下"消失"的结论。** Resident Kupffer cell 比例从 15.5% 降至 10.5%，但 Kupffer cell 全局有 420 个显著 DEG（所有群体中最多）。更准确的表述：Kupffer 发生大规模转录重编程，marker 表达变化可能导致部分细胞在分类中迁移而非真正死亡。

**3. Mo-Mac/DC-like 维持免疫监视/抗原呈递功能背景。** 组成稳定，AP 和 phagocytosis 程序维持，IL1/TNF 炎症水平低于中性粒细胞——定位偏向组织清理和免疫监视，而非主动促炎。

**4. 内皮细胞存在损伤/激活信号，但方向不明确。** 组成下降 + 297 DEGs 但上下调平衡（149 up / 148 down），内皮激活模块样本间差异大，难以区分是 LSEC 损伤还是活化。

**5. Module score 作为功能描述，不参与方向一致性判断。** 最终证据强度只看 composition + DE 的一致性。

---

## 局限性

- n = 3 vs 3，统计效力有限，所有推断为探索性
- "比例下降"不等于绝对细胞损失（可能被浸润细胞稀释）
- Module scores 描述基因表达程序，非直接通路活性测量
- 无独立验证队列

---

## 产出总结

```
results/tables/04_APAP_Control_Functional_Comparison/   — 11 CSV
results/figures/04_APAP_Control_Functional_Comparison/  — 5 PDF
  ├ key output: stage4_key_findings_summary.csv         — README-ready 结论表
  ├ key output: pseudobulk_de_by_celltype.csv            — edgeR 差异表达（全局）
  └ key output: pseudobulk_de_by_myeloid_subtype.csv     — edgeR 差异表达（髓系）
```

---

## 下一步：Stage 5 CellChat

基于 Stage 4 筛出的稳定信号轴优先分析：
1. Inflammatory neutrophil ↔ Endothelial（趋化/黏附）
2. Kupffer cell（补体/吞噬）
3. Antigen-presenting myeloid ↔ T/NK（共刺激/抗原呈递）
4. Stellate / Endothelial（修复/ECM）
