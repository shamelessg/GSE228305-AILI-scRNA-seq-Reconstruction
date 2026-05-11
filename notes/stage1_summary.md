# Stage 1 总结：QC 与 Harmony 整合

## 核心问题

6 个样本（3 Control + 3 APAP）来自不同测序批次，原始数据存在：
- 技术性批次效应（样本间系统偏差）
- 低质量细胞（空液滴、死细胞、双细胞）
- 高维稀疏数据（每个细胞只有少数基因被检测到）

目标是获得一个干净、去批次效应、保留生物学差异的整合数据集。

## 数据概况

| | |
|---|---|
| 合并细胞数 | 92,015 |
| QC 过滤后 | 90,973（去除 1,042, 1.1%） |
| 基因数 | 24,873（小鼠 mm10） |
| 高变基因 | 2,000 |
| 使用 PC 数 | 30 |
| Harmony 迭代 | 5 轮收敛 |

样本细胞数分布：

| 样本 | 组别 | 过滤后细胞数 |
|---|---|---|
| con1 | Control | 4,814 |
| con2 | Control | 20,204 |
| con3 | Control | 17,837 |
| AP-300 | APAP | 5,194 |
| AP-300-2 | APAP | 22,390 |
| AP-300-3 | APAP | 20,534 |

## 生物学思考

### QC 判断

- 保守阈值策略（nFeature 200-6000, nCount < 40000, MT < 25%），过滤率仅 1.1%。这是刻意宽松的——第一阶段目标不是做高纯度数据集，而是先看清数据全貌，避免在未充分了解前误删有生物学意义的细胞群。
- con3 丢失 3%，略高于其他样本。可能是该样本本身 RNA 质量稍差，但也可能是 QC 指标在样本间本身存在系统差异。后续注释阶段如果发现 con3 在某些细胞类型中缺失，需要回溯是否过滤策略对它的影响过大。
- con1 和 AP-300 仅有约 5000 细胞，远低于其他样本（2 万+）。这可能导致这两个样本在聚类时权重不足。若后续发现某些簇仅由这两个少数样本贡献或排除，需记录为潜在批次效应残留。

### Harmony 整合

- 使用 theta=2（默认值），在"消除样本间技术差异"和"保留 APAP vs Control 生物学差异"之间取平衡。
- 5 轮收敛说明批次效应不算强烈，数据内在一致性较好。
- 整合变量为 `sample_id`，意味着 Harmony 被告知"细胞来自不同样本"这个信息，目标是消除与此相关的变异。但没有把 `group` 告诉 Harmony——APAP vs Control 的差异应当被保留。
- 关键评估点：整合后 UMAP 中，同组样本（如 con1/con2/con3）的细胞应混合良好；但 Control 与 APAP 之间不应完全揉成一团。如果在整合后两组完全重叠，说明 theta=2 过强，需要调低（theta=0 即完全不强制多样性）；如果样本仍各自成群，说明批次效应强于预期，需要调高 theta 或检查是否存在未考虑的技术协变量。

### PCA 维度选择

- 使用 30 个 PCs，这是中等偏保守的选择。PCA 肘图用于判断是否需要更多或更少的维度。
- 维度太少 → 丢失稀有的生物学信号（如稀有细胞类型仅在高维 PC 上体现）
- 维度太多 → 引入更多噪音，Harmony 可能去拟合噪音而非生物学结构

## 技术问题与解决

### 1. Rscript 路径检测
`sys.frame(1)$ofile` 在交互式 R 中有效，但 `Rscript` 执行时调用栈为空。改用 `commandArgs(trailingOnly=FALSE)` 解析 `--file=` 参数来定位脚本路径。

### 2. Harmony 参数名不匹配
harmony 1.2.3 的 Seurat S3 方法使用 `reduction.use` 而非 `reduction`，且不存在 `assay.use` 参数。部分匹配导致"参数3有多个与之相对应的正式参数"错误。查阅 `?RunHarmony.Seurat` 确认正确签名后修正。

### 3. macOS 内存 OOM (16GB)
92K 细胞的 Seurat 对象，在 macOS 默认 `plan("multicore")` 下，R 的垃圾回收会打破 fork 的 COW 共享内存，导致物理内存翻倍耗尽。
解决方案：
- `plan("sequential")` —— 不使用多进程，避免内存膨胀
- `ScaleData(features = VariableFeatures(...))` —— 仅缩放 2,000 高变基因而非全 24,873 基因
- 在 ScaleData / RunPCA / RunUMAP / RunHarmony 之间插入 `gc()` 手动触发垃圾回收

### 4. QC summary 表格结构异常
`pre_qc_n` 和 `post_qc_n` 是 R 的 `table` 类型，除法和四舍五入后产生带 `Var1`/`Freq` 列的 data.frame，而非预期的一列数值。用 `as.numeric()` 显式转换后修复。

## 需要人工判断的节点

以下位置在当前保守参数下已自动运行，但代码中已预留注释标记，后续可根据实际数据特征调整：

1. **QC 阈值**（脚本 line ~165）：nFeature_RNA 上下限、nCount_RNA 上限、percent.mt 上限
2. **PCA 维度数**（line ~251）：当前 NPCS=30，需对照肘图决定
3. **Harmony theta**（line ~288）：当前 theta=2，需根据整合前后 UMAP 判断是否过强/过弱
4. **regress 变量**（line ~240）：当前仅回归 percent.mt，若样本间细胞周期信号强可加入 `cc.difference`

## 待确认的观察

- con1（Control）和 AP-300（APAP）细胞量明显偏低（~5000 vs ~20000）。这两个样本分别是各自组内的"少数派"，后续分析中需注意它们是否系统性偏离组内其他样本。
- 当前宽松 QC 留下了 99% 的细胞。如果后续聚类出现来源不明的小簇（可能是双细胞或 ambient RNA 污染），需回退收紧过滤后重跑。


我看到的判断依据：

QC 阈值偏宽
当前 90,973 个细胞只过滤掉约 1%，percent.mt < 25 和 nFeature_RNA > 200 比较松。con2/con3 有较多 percent.mt 在 15-25 的细胞，尤其 con3，这更像样本质量尾部，不太像 APAP 生物学信号。改成 percent.mt < 20 会更稳妥；nFeature_RNA > 300 可以去掉低基因数尾部。

30 PCs 合适
pca_elbow.pdf 里前 10 个 PC 下降很快，15-20 后进入长尾，但 30 PCs 作为后续 UMAP/Harmony 的探索参数是合理的。没有必要因为 elbow 改成很小的 PC 数。

Harmony 整合效果可以接受
Harmony 后各样本基本覆盖同一套 UMAP 结构，说明样本批次没有明显把细胞类型拆开。没有看到“某个样本完全单独成岛”的严重问题。

APAP vs Control 信号还在，但原图有绘图偏差
umap_after_harmony_by_group.pdf 没有 shuffle = TRUE，APAP/Control 叠点会互相遮挡，所以不能直接凭那张判断“APAP 是否铺满全图”。我重画成分面高亮后看，APAP 和 Control 都覆盖主要结构，但部分区域确实有组间富集差异，说明没有被 Harmony 完全抹平。

脚本本身主要问题有两个：

没有检查当前 conda 环境是不是 singlecell。项目规则要求后续运行前核对环境，这个建议补上。
UMAP 图建议加 shuffle = TRUE，最好再输出 sample/group 的 facet 高亮图，否则 9 万细胞叠点容易误判。
所以我的建议是：**让他小幅重跑第一步，不改 PC/Harmony，只收紧 QC 并补图。**重跑后如果细胞数变化不大、UMAP 结构稳定，就可以进入第二步 Global Annotation。