This repository provides the R scripts and supporting materials for the transcriptomic analyses conducted in the study Explainable Transcriptomic Framework in Sepsis.
The analytical workflow includes transcriptomic data preprocessing, biological module and axis scoring, sepsis endotype characterization, cross-cohort comparisons, longitudinal analyses, and sensitivity analyses.

The repository is intended to support transparency and reproducibility of the results presented in the manuscript.

Repository Structure
scripts/: The 13 main R scripts used for data processing, statistical analyses, and figure generation.
lib/: Supporting R functions.
data/: Reference data, gene/module definitions, annotation files, and other analysis inputs.
README.md: Documentation of the analytical workflow.

Analysis Workflow
1 Definition of analytical parameters and dictionaries
2 Download of publicly available transcriptomic datasets
3 Transcriptomic preprocessing and computation of axis/module scores
4 Consensus Transcriptomic Subtype (CTS) classification
5 Construction of analysis datasets
6 Figure 1 and Supplemental Tables 1–2
7 Figure 2: Gene coverage and inter-axis correlations
8 Figure 3: Endotype-associated biological profiles
9 Figure 4: MARS–CTS correspondence
10 Supplemental Figure 1: Gene-deletion sensitivity analysis
11 Supplemental Table 3 and Supplemental Figure 2: Longitudinal analyses
12 Supplemental Table 5: ssGSEA sensitivity analysis
13 Supplemental Table 4 and Supplemental Figure 3: Transcriptomic neutrophil-to-lymphocyte ratio analysis

Data Sources
The analyses use publicly available transcriptomic datasets, including data deposited in the Gene Expression Omnibus (GEO) and other public transcriptomic resources.
Dataset identifiers and source information are documented in the analysis scripts and reference data files.

Software Requirements
The analyses were developed in R using CRAN and Bioconductor packages.
The workflow uses packages including survival, ggplot2, data.table, limma, edgeR, GSVA, sva, and randomForest, among others.
Exact R and package versions, external dependencies, and installation instructions should be documented before the final reproducibility release.

Reproducibility
The scripts are organized according to the analytical workflow.
Reproduction of the analyses requires the appropriate source datasets, reference files, supporting functions, and R package dependencies.
The supplied scripts currently include machine-specific file paths that must be configured for a different computing environment. End-to-end execution and correspondence with manuscript results remain to be verified.

Data Availability and Restrictions
Public transcriptomic data should be obtained from their original repositories, subject to the applicable access and reuse conditions.
Any restricted clinical information or data not authorized for public redistribution must be excluded from the publicly released repository.

Citation
Please cite the accompanying manuscript when using this analytical workflow.
Manuscript citation: To be updated upon publication.
