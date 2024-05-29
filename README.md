# Analysis of Squamous Lung Cell Carcinoma Progression

## Overview

This repository contains code and data for the project titled "Analysis of Squamous Lung Cell Carcinoma Progression Using an Integrative Bioinformatics Approach." This study integrates single-cell RNA sequencing and spatial transcriptomics to provide a detailed understanding of the heterogeneity of squamous lung cell carcinomas.

### Authors
- Archana Balasubramanian
- Shivaramakrishna Srinivasan
- Aaryan Sharma
- Aathira Sarath Chandran
- Hersh Gupta

## Repository Structure

- `data/`: Contains the datasets used for the analysis, including single-cell RNA sequencing data and spatial transcriptomics data. Access the data files via the provided [Google Drive link](https://drive.google.com/drive/folders/1DLSwo194McYE81v2Zl5NL9aBTAZiGFsP?usp=sharing)
- `Analysis of Squamous Lung Cell Carcinoma Progression Using an Integrative Bioinformatics Approach.pdf`: Detailed PDF report encompassing in-depth analysis, findings, and results from the project. The report provides comprehensive insights into the analysis process and outcomes.
- Run the individual Jupyter Notebooks for analysis
  
## Data

### Single-cell RNA Sequencing Data

The single-cell RNA sequencing data for this project were obtained from the Human Cell Atlas database, including the following sources:

- Barbry_unpubl: Nasal and tracheobronchial samples from Nice University Hospital (Dec 2020). Processed with Cell Ranger 6.0.0 and Scanpy.
- Schiller_2021: Tumor-free lung samples from Asklepios Fachkliniken München-Gauting. Processed with Cell Ranger 3.1.0.
- Duong_lungMAP_unpubl: Postmortem lung samples from University of Rochester. Processed with Cell Ranger 3, analyzed with Pagoda2 and Seurat.
- Jain_Misharin_2021: Nasal epithelial samples from Northwestern Medicine. Processed with Cell Ranger 3.1.0 and Seurat.
- Schultze_unpubl: Normal lung tissues from Hannover Medical School. Processed with Seq-Well and Seurat v3.2.
- These datasets provide a comprehensive view of lung tissue characteristics and were processed using state-of-the-art bioinformatics pipelines.

### Spatial Transcriptomics Data

The spatial transcriptomics dataset originates from human lung cancer tissue samples diagnosed as squamous cell carcinoma, acquired from Avaden Biosciences by 10x Genomics. 

- Preservation: FFPE method for long-term storage.
- Preparation: Followed Visium CytAssist Spatial Gene Expression for FFPE protocol (CG000518).
- Sectioning: Tissue sliced into 5 µm sections.
- Staining: H&E staining after deparaffinization.
- Imaging: Olympus VS200 Slide Scanning Microscope.
- Processing: Included coverslipping, decoverslipping, dehydration, and decrosslinking.
- Analysis: Analytes transferred to a Visium CytAssist Spatial Gene Expression slide.
- Sequencing: Libraries prepared per Visium CytAssist Reagent Kits for FFPE User Guide (CG000495) and sequenced on an Illumina NovaSeq platform.

This dataset reveals the spatial distribution of gene expression within lung cancer tissues, providing insights into squamous cell carcinoma progression and informing therapeutic strategies. 
Key metrics include the number of spots detected, median genes per spot, and median UMI counts per spot.

### Prerequisites

Ensure you have the following software and libraries installed:
- Python (>=3.6)
- Scanpy
- Squidpy
- LIANA
- Other dependencies 

### Running the Analysis

1. **Clone the repository:**
   ```bash
   git clone https://github.com/archana-b24/sc-spatial-analysis.git
   cd sc-spatial-analysis

2. **Run the Jupyter Notebooks for analysis:**

- Data Filtering: 01_dataFiltering.ipynb
Filters raw single-cell RNA sequencing data to remove low-quality cells and genes.

- Doublet Removal: 02_doublet_removal.ipynb
Identifies and removes doublets (artificially merged cells) from the dataset.

- Main Analysis: 03_analysis.ipynb
Conducts primary analysis including normalization, clustering, and initial exploratory data analysis of single-cell data.

- Trajectory Analysis: 04_trajectory.ipynb
Investigates cellular trajectories and differentiation paths.

- Data Integration: 05_integration.ipynb
Integrates single-cell datasets with spatial data to uncover intricate cellular interactions and spatial organization within the tissues.

- Ligand-Receptor Interaction Analysis: 06_liana.ipynb
Analyzes cell-cell communication through ligand-receptor interactions using the LIANA package.

- Gene Set Enrichment Analysis: 07_GSEA.ipynb
Performs gene set enrichment analysis to identify significantly enriched biological pathways.

## Note:
Feel free to explore the repository, access the data, and review the analysis report to gain a thorough understanding of the "Analysis of Squamous Lung Cell Carcinoma Progression Using an Integrative Bioinformatics Approach." 
If you have any questions or require further information, don't hesitate to reach out.
