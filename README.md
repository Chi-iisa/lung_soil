This repository contains the custom R scripts used to process, analyze, and visualize transcriptomic data associated with the manuscript:

draft lung_soil

---

# Overview
The code provided here reproduces the main computational analyses and figures presented in the manuscript, including:

- Data preprocessing and filtering
- Integration with clinical metadata
- Z-score transformation
- Heatmap generation
- Boxplot visualization
- Statistical testing

---

# Repository Structure
- The data folder contains metadata used for downstream analyses.
- The data-raw folder contains the transcriptomic data matrices for each dataset used in this study.
- The functions folder contains R files with custom functions and some supporting variables.
- The scripts folder contains R files used to generate the figures and supplementary figures.

---

# Data Availability
The transcriptomic datasets analyzed in this study are publicly available from:
- García-Mulero et al, J Immunother Cancer, 2020 (https://github.com/odap-ubs/mets-immunecluster/)
- Pancreatic ductal adenocarcinoma (PDAC) data, with GEO accession number GSE205154 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205154)

---

# Requirements

## Software
R ≥ 4.2

## R Packages
Main packages used:

tidyverse
pheatmap
ggplot2
gridExtra
stringr
ConsensusTME

## Code Availability Statement
All custom scripts used to generate the results reported in the manuscript are available in this repository.
