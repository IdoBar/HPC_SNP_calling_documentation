[![DOI](https://zenodo.org/badge/927221841.svg)](https://doi.org/10.5281/zenodo.16901697)

## Background

*Ascochyta rabiei* is a necrotrophic fungal disease of chickpea (*Cicer
arietnum*) with the potential to cause total yield loss under conducive
conditions. The Australian *A. rabiei* population presents a challenge
to produce reproducible SNPs due to its clonal nature and repetitive
genome. To be able to study *A. rabiei* in the field using an
experimental population a workflow needs to be developed that will
produce both high-resolution genome wide-SNP haplotypes that are
applicable for downstream genetic analysis and re-identify clonal
lineages.

An experimental population of six *A. rabiei* genotypes were defined as
genetically distinct by principal component analysis (PCA). Low-medium
coverage whole genome sequencing was used on field population of 35 *A.
rabiei* genotypes that had been re-sampled two years after the
experimental population was released in the field and allowed to
proliferate. A variant calling pipeline using Freebayes, followed by
variant filtering for polymorphic informative loci with sequencing depth
&gt; 30 and genotype call rate &gt; 99% resulted in a VCF calling
pipeline with high confidence SNP variants. High-quality haplotypes and
clonal lineages of the experimental population were successfully defined
using this variant calling method. This strategy utilised discriminatory
analysis of principal components (DAPC) where 8 clusters were defined
and separated the experimental population into separate groups, which
was confirmed by phylogenetic tree constructed using maximum-likelihood
analysis.

The workflow demonstrated by this study presents the opportunity to
trace an experimental population while still capturing sufficient
resolution for downstream genetic analysis and taking into consideration
the handling of clonal data.

### Aims

This study aims to develop a workflow utilising low to medium coverage
whole genome sequencing to determine a reliable genome-wide SNP-based
haplotype for use in experimental evolution of fungal clonal pathosystem
for which there are two main outcomes. The first is the obtainment of
high-resolution haplotype clones within a cluster for the analysis of
selection and divergence. The second, is the reidentification of an
experimental population in the field for fitness purposed and assignment
of haplotype-based clusters. This repository contains the bioinformatics
pipeline used to perform the following research aims:

-   Compare variant detection pipelines that allow accurate
    identification of Multi-Locus-Haplotypes (MLH)  
-   Identify MLHs descending from an original population and measure
    their frequencies  
-   Asses the impact of host genotype on MLHs frequencies  
-   Identify genomic loci and regions that drive isolate
    adaptation/evolution under different selection pressures

Detailed information on the experimental and analysis approaches and
methods, including detailed bioinformatics pipelines and code can be
found at
<https://github.com/IdoBar/HPC_SNP_calling_documentation/index.html>.

This repository was created by Ido Bar (<i.bar@griffith.edu.au>) and
Hayley Wilson (<hayley.wilson7@griffithuni.edu.au>).  
Please contact us for any additional information or collaboration
opportunities.
