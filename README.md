# GAP43 as a synaptic marker in brain-derived extracellular vesicles in Alzheimer's disease

Brain-derived extracellular vesicle (BD-EV) proteomics analysis for:

> Lan G., Li B., Wang M., **Barve A.**, Colin M., Buée L., … Richetin K., … Guo T. **Plasma growth-associated protein 43 correlates with synaptic loss in Alzheimer's disease.** *Cell Reports Medicine* (2025). [doi:10.1016/j.xcrm.2025.102508](https://doi.org/10.1016/j.xcrm.2025.102508) · [PubMed 41406948](https://pubmed.ncbi.nlm.nih.gov/41406948/)

The paper shows that plasma GAP43 tracks CSF GAP43, synaptic-density PET, tau aggregation, brain atrophy and cognitive decline in two Alzheimer's disease (AD) cohorts (n = 785). This repository holds the post-mortem part I contributed: proteomics of BD-EVs from human prefrontal cortex (n = 21), showing that GAP43 is a hub in synaptic protein networks that are depleted in AD.

## Analyses

| Analysis | What it shows | Files |
|---|---|---|
| **Proteome-wide GAP43 correlation** | Correlation of GAP43 with each of 2,085 BD-EV proteins; top 20 positive and negative partners | `GAP43_correlate.xlsx` (sheet 1), `GAP43_correlation.xlsx` |
| **Synaptic and neurodegeneration panels** | GAP43 correlation with a synaptic protein panel and with tau, APP and other neurodegeneration proteins | `GAP43_correlate.xlsx` (sheets `Synaptic`, `NDD`) |
| **Pathway enrichment** | GSEA (GO, KEGG, Reactome, WikiPathways) and extraction of the enriched pathways whose core genes include GAP43 | `GAP43_GSEA.R` |
| **Interaction network** | SIGNOR protein–protein interaction network of GAP43 and its correlated proteins | `gap43_correlate_signor.cys` (Cytoscape) |
| **Group comparisons** | GAP43 levels by amyloid/tau (A/T) status | `Classification All Patient A+T+.prism` (GraphPad Prism) |

## How to run

`GAP43_GSEA.R` requires R ≥ 4.2 with `clusterProfiler`, `ReactomePA`, `org.Hs.eg.db`, `enrichplot`, `DOSE`, `ComplexHeatmap`, `tidyverse`, `ggrepel` and `patchwork` (CRAN and Bioconductor).

1. Replace the `setwd(...)` path at the top of the script with your local folder.
2. Provide a ranked fold-change table (gene symbol, log2 fold change) as the input file named in the script.
3. Run the script top to bottom. The last section filters the results for pathways that include GAP43.

The `.cys` file opens in Cytoscape 3.10 or later; the `.prism` files open in GraphPad Prism 10.

## Citation

If you use this code, please cite the paper above.

## Contact

Aatmika Barve · [ORCID 0000-0003-2331-5825](https://orcid.org/0000-0003-2331-5825) · [LinkedIn](https://www.linkedin.com/in/aatmika-b-504aaa212)
