# Note 2: Annotation and Myeloid Subclustering

## Purpose

This note records how the global cell annotation was made, why myeloid cells were selected for focused subclustering, and which labels should be treated as working annotations rather than final cell identities.

The main idea was to move from a global atlas to a focused APAP-relevant compartment without pretending that marker-based annotation is absolute.

## Global Clustering Decision

Four resolutions were compared during global clustering:

| Resolution | Cluster number |
|---|---:|
| 0.2 | 21 |
| 0.4 | 27 |
| 0.6 | 34 |
| 0.8 | 39 |

Resolution 0.4 was selected for global annotation. It gave enough separation between major lineages while avoiding excessive fragmentation. At this stage, the goal was not to define every fine subtype, but to obtain stable broad cell classes for downstream comparison.

The selected resolution separated major liver non-parenchymal and immune populations, including T cells, B cells, NK cells, Kupffer cells, monocyte-derived macrophage/DC-like cells, neutrophil-related populations, endothelial cells, cholangiocyte/epithelial cells and stromal cells.

## Marker-Based Annotation Logic

Cell labels were assigned by combining cluster markers, canonical marker expression, DotPlot/FeaturePlot patterns and biological context.

Representative marker logic:

| Cell type | Marker evidence |
|---|---|
| T cell | `Trac`, `Cd3d`, `Tcf7`, `Il7r` |
| B cell | `Ms4a1`, `Cd79a`, `Cd79b`, `Pax5` |
| NK cell | `Ncr1`, `Klra` genes, `Gzma`, `Prf1` |
| Kupffer cell | `Clec4f`, `Vsig4`, `Cd5l`, `Folr2`, `Cd163` |
| Monocyte-derived macrophage / DC-like | `Lyz2`, `Cx3cr1`, `Fcgr1`, `Cd209a`, `Clec4` genes |
| Neutrophil | `S100a8`, `S100a9`, `Cxcr2`, `Csf3r`, `Ltf`, `Camp`, `Ly6g` |
| Endothelial | `Pecam1`, `Cdh5`, `Kdr`, `Lyve1`, `Stab2` |
| Cholangiocyte / epithelial | `Epcam`, `Krt8`, `Krt18`, `Krt19` |
| Fibroblast / stellate | `Dpt`, `Lum`, `Pdgfra`, `Mmp2`, `Tcf21` |

Annotation was treated as a biological judgement step, not an automatic output from a marker table.

## Uncertain Global Clusters

Several clusters were marked as uncertain:

- Low-quality or mixed clusters with high mitochondrial/ribosomal/stress signals.
- B cell clusters carrying hepatocyte ambient RNA signals such as `Alb`, `Apoa1` and `Ttr`.
- Cytotoxic T/NK-like clusters where CD8 T and NK markers overlapped.
- Rare mast/basophil-like clusters with limited cell number.

These clusters were not used as strong biological anchors.

## Why Myeloid Cells Were Selected

Myeloid cells were selected for focused subclustering because they are central to APAP-induced acute liver injury. The global annotation identified several relevant populations:

- inflammatory neutrophils;
- activated neutrophils;
- mature neutrophils;
- Kupffer cells;
- monocyte-derived macrophage / DC-like cells.

Together, these cells formed a large enough compartment for subclustering and showed APAP-related composition trends. This made myeloid cells a stronger focus than lymphoid or stromal populations for this dataset.

## Myeloid Subclustering Decision

In the myeloid object, four resolutions were compared:

| Resolution | Cluster number |
|---|---:|
| 0.2 | 10 |
| 0.4 | 16 |
| 0.6 | 20 |
| 0.8 | 21 |

Resolution 0.4 was selected again, but for a different reason than in the global object. Here it provided enough granularity to distinguish Kupffer, monocyte-derived macrophage/DC-like and several neutrophil states, while avoiding excessive splitting of small or noisy clusters.

## Working Myeloid Labels

The final myeloid labels were working annotations. The most important groups were:

- Inflammatory neutrophil
- Activated neutrophil
- Mature neutrophil
- IFN-responsive inflammatory neutrophil
- Resident Kupffer cell
- Monocyte-derived macrophage
- Antigen-presenting macrophage / Mo-Mac
- DC-like antigen-presenting myeloid
- Cycling neutrophil

Several manual corrections were important:

- A cluster initially resembling monocyte-derived macrophage was revised to Resident Kupffer cell because it strongly expressed `Clec4f`, `Vsig4` and `Timd4`.
- A neutrophil-like cluster carrying high hepatocyte ambient markers was labeled as possible doublet or ambient-contaminated neutrophil.
- A proliferative neutrophil cluster was identified by `Mki67`, `Top2a` and `Stmn1`.

## Interpretation Boundary

The myeloid annotation supports exploratory analysis, but not publication-level cell identity claims by itself.

The safest conclusion from this stage is:

> APAP samples show a shift toward neutrophil-related myeloid states, while Kupffer and antigen-presenting myeloid compartments show composition or transcriptional changes that require cautious interpretation.

This stage justified moving into functional comparison rather than treating the subclusters as final mechanistic discoveries.
