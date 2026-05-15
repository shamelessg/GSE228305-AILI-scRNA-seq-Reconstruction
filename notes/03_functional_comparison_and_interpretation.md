# Note 3: Functional Comparison and Interpretation

## Purpose

This note records the major scientific turn in the project: Stage 4 was changed from pseudotime analysis to APAP-Control functional comparison.

The change was made because the data did not support a clean continuous trajectory across the whole myeloid compartment. The project therefore moved toward a more defensible question:

> Which non-parenchymal cell types or myeloid subtypes show composition and functional-state changes after APAP treatment?

## Why Pseudotime Was Not Used as the Main Line

The original plan considered Monocle3 pseudotime analysis. After myeloid subclustering, the UMAP showed several separated structures:

- neutrophil-related states;
- Resident Kupffer cells;
- monocyte-derived macrophage / DC-like antigen-presenting cells.

This structure looked more like separate lineages or functional compartments than one continuous differentiation path.

Using all myeloid cells to force a single pseudotime trajectory would create an attractive figure, but the biological direction would be weak. Root selection would also be difficult to justify. Therefore, pseudotime was not used as the Stage 4 main analysis.

Pseudotime could still be considered later for a restricted lineage, such as neutrophil states only, but only if a clearer local continuum is observed.

## Stage 4 Analysis Design

Stage 4 used three evidence layers:

1. Composition change  
   Sample-level proportions were compared between Control and APAP.

2. Functional module score  
   Small interpretable gene sets were used to describe inflammatory, chemokine, IFN, antigen-presentation, Kupffer-resident, phagocytosis and repair-related programs.

3. Pseudobulk differential expression  
   Counts were aggregated at sample level within each cell group, then APAP vs Control was tested using edgeR.

This design was chosen because n = 3 vs 3 is too small for strong claims, but still useful for exploratory effect-size and direction-based interpretation.

## Key Findings

### Neutrophil Expansion

Inflammatory, activated and mature neutrophil-related populations increased in APAP samples. This was the most stable observation across the project.

The strongest interpretation is composition expansion, not necessarily strong within-cell transcriptional remodeling. Some neutrophil groups had few or no FDR-significant DEGs, which may reflect limited sample size, low power or a real composition-dominant effect.

### Kupffer Cell Reprogramming

Resident Kupffer cell proportion decreased in the myeloid compartment, while global Kupffer cells showed many differentially expressed genes.

The safer interpretation is not "Kupffer cells disappear", but:

> APAP is associated with Kupffer-related transcriptional remodeling, and apparent proportion changes may partly reflect marker shifts, tissue injury or dilution by infiltrating neutrophils.

### Endothelial Change

Endothelial cells showed a decrease in proportion and clear pseudobulk DE signals. However, the direction was not simple: up- and down-regulated genes were balanced, and functional modules were variable across samples.

This supports endothelial injury or state disturbance as a candidate signal, but not a confident claim of endothelial activation or loss.

### Antigen-Presenting Myeloid and Mo-Mac/DC-Like Cells

These groups showed antigen-presentation and immune-regulatory features, but APAP-Control direction was not fully consistent across composition, module score and DE.

They are better treated as exploratory immune-state findings rather than the main mechanistic axis.

## Why CellChat Was Downgraded

Stage 5 was not performed as a full CellChat discovery analysis. Instead, it was treated as targeted ligand-receptor screening.

The reason is simple: ligand-receptor inference from scRNA-seq is hypothesis-generating. It does not directly measure physical interaction, protein abundance or functional signaling.

The project therefore screened only a few axes suggested by Stage 4:

- neutrophil or activated neutrophil to endothelial;
- Kupffer/macrophage to neutrophil or endothelial;
- antigen-presenting myeloid to T/NK;
- fibroblast/stellate to endothelial or myeloid.

Candidate pairs such as `Cxcl2-Cxcr2`, `Ccl6-Ccr1`, `Thbs1-Cd47` and `App-Cd74` were kept as weak supporting clues, not independent discoveries.

## Interpretation Rules

The project follows these boundaries:

- n = 3 Control vs 3 APAP limits statistical power.
- Cell proportion changes are relative, not absolute cell counts.
- Module scores describe expression programs, not direct pathway activity.
- Pseudobulk DE is more appropriate than treating cells as independent replicates, but it is still limited by sample number.
- Cell communication analysis is candidate screening, not mechanism proof.

## Final Position

The project supports the following cautious conclusion:

> APAP treatment is associated with remodeling of the liver non-parenchymal immune microenvironment, especially neutrophil expansion and Kupffer/endothelial/myeloid state changes. These findings are exploratory and should be treated as candidate directions for validation rather than established mechanisms.

This interpretation is less dramatic than a full mechanistic story, but it is closer to what the data can actually support.
