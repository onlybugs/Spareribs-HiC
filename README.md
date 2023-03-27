# Spareribs-HiC: a One-stop, User-friendly, and Scalable Hi-C Data Analysis Pipeline


## Summary
Spareribs-HiC is implemented by snakemake, the most popular process writing language, which enables one-stop sequencing of downstream files `(.fq)` with various data preprocessing operations to obtain a classical representation of chromatin 3D structure data `(.cool/.mcool)`. In addition, Spareribs-HiC can automatically call embedded scripts `(.py/.R/.sh)` for chromatin compartment **(A/B)** analysis and topologically associating domain **(TAD)** analysis. More importantly, users can write their own personalized analysis scripts to be embedded in Spareribs-HiC to perform a variety of downstream analysis tasks in a synergistic manner.


## Dependency
