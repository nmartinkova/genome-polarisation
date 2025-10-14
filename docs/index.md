---
title: "Genome polarisation with diemr"
author: "Natália Martínková"
date: "Last updated: 2025-10-14"
site: bookdown::bookdown_site
documentclass: book
bibliography: [refs.bib]
biblio-style: apalike
link-citations: true
output:
  bookdown::gitbook:
    css: style.css
---



# Introduction

Genome polarisation is a genome-painting approach based on the likelihood-based diagnostic index expectation maximisation (_diem_) algorithm [@Baird2023]. It identifies which alleles of single-nucleotide variants (SNV) belong to either side of a barrier to gene flow, co-estimating both the assignment of individuals to a barrier side and the diagnosticity of each marker, meaning how consistently individuals on one side are homozygous for the allele associated with that side.

By inferring which parts of the genome correspond to each parental lineage, genome polarisation provides a direct view of **how barriers to gene flow shape genomic architecture**. It can detect and quantify hybridisation, distinguish introgressed from non-introgressed regions, and reveal how species boundaries evolve during speciation and diversification. Compared with methods such as STRUCTURE, ADMIXTURE, or PCA, which summarise population structure statistically, genome polarisation explicitly identifies the genomic segments that define the divergence between taxa or lineages.

The diagnostic index computed by _diem_ **highlights markers that are most informative** for the primary axis of genetic differentiation. These diagnostic loci can then be used to describe patterns of hybridisation, assess barrier strength, or visualise the genomic distribution of ancestry.

This book provides a step-by-step guide to performing genome polarisation analyses in `R` using [the `diemr` package](https://github.com/nmartinkova/diemr), also available from [CRAN](https://CRAN.R-project.org/package=diemr). The package includes functions for input validation, file-format conversion, visualisation, and diagnostic summaries, with examples based on typical genomic datasets such as variant call format (VCF) files and SNP matrices. By the end of the book, users will be able to run the complete analysis workflow, interpret its outputs, and generate graphical representations of genome polarisation.

The _diem_ algorithm itself is also implemented in `Mathematica` and `Python`, available [here](https://github.com/StuartJEBaird/diem), allowing reproducibility and interoperability across analytical environments.

## How to cite

If you use the *diemr* package or this documentation in your work, please cite:

Baird, S. J. E., J. Petružela, I. Jaroň, P. Škrabánek, and N. Martínková. 2023. *Genome polarisation for detecting barriers to geneflow.* Methods in Ecology and Evolution 14: 512–28. https://doi.org/10.1111/2041-210X.14010.

and this online documentation for the `diemr` package:

Martínková, N. 2025. *Genome polarisation with diemr.* Bookdown online documentation.  
Available at: <https://nmartinkova.github.io/genome-polarisation/>.

