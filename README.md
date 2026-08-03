# Simultaneous False Discovery Rate Methods for Genomics

This repository contains the accompanying code for my master's thesis.

## Terms and symbols

| Term or symbol    | Meaning |
| ----------------- | ------- |
| $\alpha_{min}$    | The smallest significance level for which a method rejects a pathway (or any discovery set). |
| Decision boundary | The boundary of volcano plot filters for which a method maintains FDR control over the discovery set. |

## Layout

- `/figures`:               Figures in Sections 5.2 and 5.3 for all ten real data sets.
  - `/pathway_comp`:        Pairwise comparisons of each method's $\alpha_{min}$/p-values.
  - `/pathway_size`:        Relationships between $\alpha_{min}$/p-values and pathway size for each method.
  - `/volcano`:             Comparison of each method's decision boundary.
- `/notebooks`:             ...
  - `fdp_comp.nb.html`:     Comparison of the mean and median false discovery proportion on simulated data.
  - `pathway.nb.html`:      Performing a pathway analysis on real example data.
  - `venn_diagram.nb.html`: Plotting a Venn diagram to compare which pathways are deemed enriched between methods.
  - `volcano.nb.html`:      Plotting each method's decision boundary based on simulated data.
- `/output`:                Output files for all ten real data sets of pathway analyses performed on KEGG pathways.
- `functions.R`:            Implementations of algorithms, data generation, and helper functions for closed e-BH.
  <br><b>Note</b>:          The decision boundary algorithm is not optimised and may run very slowly.