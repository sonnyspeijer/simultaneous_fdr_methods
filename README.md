# Simultaneous False Discovery Rate Methods for Genomics

This repository contains the accompanying code for my master's thesis.

## Terms and symbols

The layout below uses some terms and symbols unique to my thesis. Their definitions appear below as a quick reference.

| Term or symbol    | Meaning |
| ----------------- | ------- |
| $\alpha_{min}$    | The smallest significance level for which a closed method rejects a pathway (or any discovery set). |
| Decision boundary | A boundary based on volcano plot filters that determines which discovery sets belong to the set of valid discovery sets for a given method. |

## Layout

- `/figures`:               Figures in Sections 5.2 and 5.3 for all ten real data sets. <hr>
  - `/pathway_comp`:        Pairwise comparisons of each method's $\alpha_{min}$/p-values.
  - `/pathway_size`:        Relationships between $\alpha_{min}$/p-values and pathway size for each method.
  - `/volcano`:             Comparison of each method's decision boundary. <hr>
- `/notebooks`:             Files that explain how various figures/results were obtained.<br>To view a notebook, download the file to your computer and open it locally.<hr>
  - `fdp_comp.nb.html`:     Comparison of the mean and median false discovery proportion on simulated data.
  - `pathway.nb.html`:      Performing a pathway analysis on example data.
  - `venn_diagram.nb.html`: Demonstration of how to plot a Venn diagram to compare enriched pathways between methods.
  - `volcano.nb.html`:      Demonstration of how to plot each method's decision boundary.
  - `functions.R`:          Implementations of algorithms, data generation, and helper functions for closed e-BH.
  <br><b>Note</b>:          The decision boundary algorithm is not optimised and may run very slowly for large data sets. <hr>
- `/output`:                Output files for all ten real data sets of pathway analyses performed on KEGG pathways.
