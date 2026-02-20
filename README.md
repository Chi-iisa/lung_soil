This repository contains the custom R scripts used to process, analyze, and visualize transcriptomic data associated with the manuscript:

**"A signature of NK cell effector function characterizes lung metastases and reveals therapeutic antimetastatic potential"**
Moreno-Manuel et al.

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

├── data/ # Metadata
│ ├── NKgenes_list.csv
│ ├── IPS_genes.txt
│ └── datasets_lung_soil.RData
|
├── data-raw/ # Numerical data
│ ├── ex_pancan.txt
│ ├── ex_crc.txt
│ └── ex_mel.txt
│ └── ex_all.txt
|
├── functions/ # Custom functions and supporting variables
│ ├── functions.R
│ ├── variables.R
│ └── variables_all.R
│ └── variables_crc.R
│ └── variables_pancan.R
│ └── variables_mel.R
│
├── scripts/ # Figures
│ ├── 01_figure1.R
│ ├── 02_figure2.R
│ ├── 03_suppl.fig1.R
│ ├── 04_suppl.fig2.R
│ └── 05_suppl.fig3.R
│ └── 06_suppl.fig4.R
│ └── 07_suppl.fig6.R
│
└── README.md


---

# Data Availability

The transcriptomic datasets analyzed in this study are publicly available from:

- García-Mulero et al, J Immunother Cancer, 2020 (https://github.com/odap-ubs/mets-immunecluster/)
- Pancreatic ductal adenocarcinoma (PDAC) data, with GEO accession number GSE205154 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205154)

---

# Requirements

## Software

- R ≥ 4.2

# R Packages

Main packages used:

tidyverse
pheatmap
ggplot2
gridExtra
gtable
stringr


# Code Availability Statement

All custom scripts used to generate the results reported in the manuscript are available in this repository.
