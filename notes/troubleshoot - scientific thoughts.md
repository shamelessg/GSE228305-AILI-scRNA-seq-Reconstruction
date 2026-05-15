# Troubleshooting and Scientific Thoughts

This note records questions that came up during the project. They are not polished conclusions. They are the kinds of practical and scientific questions that decide whether a single-cell project is believable.

## 1. How should clustering resolution be selected?

Resolution is not simply "higher is better".

A low resolution may merge biologically different cell types. A high resolution may split one cell type into many small clusters that are driven by noise, sample bias, cell cycle or mitochondrial stress.

In this project, resolution was selected by asking:

- Do known major lineages separate clearly?
- Are small clusters supported by meaningful marker genes?
- Does a cluster appear across multiple samples, or only in one sample?
- Is the current stage asking for broad cell types or fine subtypes?
- Does increasing resolution add biological information, or only fragment the map?

For global annotation, resolution should be coarse enough to identify major lineages. For myeloid subclustering, a slightly finer resolution is acceptable because the biological question is focused on internal heterogeneity.

## 2. How can cluster marker genes be used to infer cell identity?

A marker table is only a starting point. A cell label should come from a combination of:

- top enriched genes;
- canonical marker genes;
- DotPlot and FeaturePlot expression patterns;
- whether markers are co-expressed in the same cluster;
- whether the label makes sense in tissue and disease context;
- whether the cluster appears consistently across samples.

For example, a Kupffer cell label is stronger when `Clec4f`, `Vsig4`, `Cd5l`, `Folr2` and complement genes appear together. A neutrophil label is stronger when `S100a8`, `S100a9`, `Cxcr2`, `Csf3r`, `Ltf`, `Camp` or `Ly6g` support the same direction.

If a cluster expresses markers from two unrelated lineages, it may be a doublet, ambient RNA contamination or a transitional state. It should not be assigned a confident label too quickly.

## 3. How do I distinguish real biology from doublets or ambient RNA?

Warning signs include:

- one cluster expressing markers from two very different lineages;
- immune cells carrying strong hepatocyte genes such as `Alb`, `Apoa1` or `Ttr`;
- very small clusters appearing mainly in one sample;
- high mitochondrial or ribosomal/stress signals;
- cluster identity depending on only one marker gene.

In this project, some neutrophil-like clusters carried hepatocyte ambient signals. These were not used as main biological evidence. This is important because APAP liver injury can release abundant hepatic RNA, increasing the risk of ambient contamination.

## 4. When is pseudotime appropriate?

Pseudotime should not be used just because it can produce an attractive trajectory plot.

It is more appropriate when:

- cells form a visible continuous structure;
- the biological process is expected to be gradual;
- root selection has a clear biological basis;
- marker genes change smoothly along the inferred path;
- the trajectory is not just connecting unrelated lineages.

In this project, the myeloid UMAP showed separated neutrophil, Kupffer and Mo-Mac/DC-like regions. That structure did not support one global myeloid pseudotime. Therefore, APAP-Control functional comparison was a better main analysis.

## 5. Why use pseudobulk instead of cell-level differential expression?

Single cells from the same mouse are not independent biological replicates. Treating thousands of cells as thousands of independent samples can produce inflated significance.

Pseudobulk analysis aggregates cells within each sample and cell type. This makes the mouse/sample the unit of comparison, which is closer to the real experimental design.

However, pseudobulk is still limited here because there are only 3 Control and 3 APAP samples. It is more reliable than cell-level tests, but still exploratory.

## 6. How should cell proportion changes be interpreted?

Cell proportions are relative. If neutrophils expand strongly, other cell types may appear to decrease even if their absolute numbers did not change.

Therefore:

- "proportion decreased" does not necessarily mean cell loss;
- "proportion increased" does not necessarily mean local proliferation;
- validation would need flow cytometry, tissue imaging or absolute cell counting.

This is why endothelial and Kupffer proportion changes were interpreted cautiously.

## 7. What does module score actually mean?

Module score is a summary of expression for a selected gene set. It does not directly prove pathway activation, protein activity or cell function.

A useful module score should:

- use a small and interpretable gene set;
- be checked at sample level, not only pooled-cell level;
- be used as functional description, not mechanism proof;
- be interpreted together with marker genes and DE results.

In this project, module scores helped describe inflammatory, chemokine, IFN, antigen-presentation and Kupffer-resident programs, but they were not treated as direct evidence of pathway activity.

## 8. What can CellChat or ligand-receptor analysis prove?

Ligand-receptor analysis can suggest possible communication axes, but it cannot prove real signaling.

It does not directly measure:

- protein abundance;
- ligand secretion;
- receptor activation;
- spatial proximity;
- functional consequence.

Therefore, Stage 5 was designed as targeted screening rather than full mechanism discovery. Candidate pairs are useful for wet-lab follow-up, not final conclusions.

## 9. How strong is the evidence in this project?

The strongest signal is APAP-associated neutrophil expansion. It appears consistently and matches APAP liver injury biology.

Kupffer and endothelial changes are interesting but more complex. They may reflect transcriptional remodeling, cell-state shifts, injury response, dilution by infiltrating cells or annotation effects.

The safest project-level conclusion is:

> The data support exploratory evidence of APAP-induced immune microenvironment remodeling, but not a validated molecular mechanism.

This boundary is part of the result, not a weakness to hide.
