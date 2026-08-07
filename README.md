# Simultaneous False Discovery Rate Methods for Genomics

This repository contains the supplementary material for my master's thesis.

## Terms and symbols

The layout below uses some terms and symbols unique to my thesis. Their definitions appear below as a quick reference.

| Term or symbol    | Meaning |
| ----------------- | ------- |
| $\alpha_{\textnormal{min}}$    | The smallest significance level for which a closed method rejects a pathway (or any discovery set). |
| Decision boundary | A boundary based on volcano plot filters that determines which discovery sets belong to the set of valid discovery sets for a given method. |

## Layout

- `figures/`:               Figures in Sections 5.2 and 5.3 for all ten real data sets. <hr>
  - `pathway_comp/`:        Pairwise comparisons of each method's $\alpha_{\textnormal{min}}$/p-values.
  - `pathway_size/`:        Relationships between $\alpha_{\textnormal{min}}$/p-values and pathway size for each method.
  - `volcano/`:             Comparison of each method's decision boundary. <hr>
- `notebooks/`:             Files that explain how various figures/results were obtained.<br>To view a notebook, download the file to your computer and open it locally.<hr>
  - `fdp_comp.nb.html`:     Comparison of the mean and median false discovery proportion on simulated data.
  - `pathway.nb.html`:      Performing a pathway analysis on example data.
  - `venn_diagram.nb.html`: Demonstration of how to plot a Venn diagram to compare enriched pathways between methods.
  - `volcano.nb.html`:      Demonstration of how to plot each method's decision boundary.
  - `example_data.csv`:     Example data for the pathway analysis.
  - `functions.R`:          Implementations of algorithms, data generation tools, and helper functions for closed e-BH.<hr>
- `output/`:                Output files for KEGG pathway analyses performed on all ten real data sets.

## Data availability

The real data sets referenced in the layout are not publicly available and are therefore not included in this repository. However, the analyses of the real data sets were performed in the same manner as demonstrated in the notebooks `volcano.nb.html` and `pathway.nb.html`.
