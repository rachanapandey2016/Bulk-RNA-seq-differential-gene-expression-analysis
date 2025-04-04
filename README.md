# Bulk RNA-seq-differential-gene-expression-analysis
---

**Project Title**  
Transcriptomic Signatures in Stage IV Melanoma Patients: Gene Expression Patterns in Long Responders vs. Short Progressors to Immunotherapy  

**Project Description**  
This repository contains the complete analysis pipeline for identifying differentially expressed genes (DEGs) in peripheral **CD4+ T cells** from **Stage IV melanoma patients** who have undergone **immune checkpoint inhibitor (ICI)** therapy. I performed this analysis as part of a class project for the *Computational Genomics* course at the University of Minnesota.
The project aims to distinguish **long responders** (relapse-free survival >12 months) from **short progressors** (<12 months) based on gene expression signatures using publicly available bulk RNA-seq data.  

**Reference paper**  
*Transcriptomic signatures in peripheral CD4+T-lymphocytes may reflect melanoma staging and immunotherapy responsiveness prior to ICI initiation*  

**Research Objectives**  
General Objective: 
To identify transcriptomic markers that distinguish Long Responders from Short Progressors in Stage IV melanoma patients after ICI therapy.  
Specific Objectives: 
1. To identify differentially expressed genes (DEGs) in CD4+ T cells between S4L and S4S. 
2. To visualize expression patterns using PCA, volcano plots, and heatmaps. 
3. To interpret the biological relevance of DEGs in the context of immune response and melanoma progression.

**Data information**  
- Bulk RNA seq data from pheripheral CD4+ T cells, single end, 75bp
- Immunotherapy provided- Nivolumab, Pembrolizumab, or Nivolumab/Ipilimumab
- Groups:
  - 6 long responders (Relapse free survival >12 months)
  - 6 short progressors (Relapse free survival <12 months)
  - GEO Accession number-  GSE292798
---   
## 🧪 Methods Overview
- **Preprocessing**: FastQC, MultiQC, fastp
- **Alignment**: STAR
- **Quantification**: Salmon
- **Gene-level summarization**: tximport
- **Differential Expression Analysis**: DESeq2 with LFC shrinkage (`ashr`)
- **Visualization**: PCA, MA Plot, Volcano Plot, Heatmap

## 📊 Key Findings
- Identified **14 DEGs** (adj. p < 0.05, |log2FC| > 1)
- Notable genes:
  - **Upregulated**:- 5 upregulated genes: AMPD2, NPIPA5, SAMD5, RNU1-3, UPF3AP1
  - **Downregulated**:9 downregulated genes PROSER3, BCL2L1, CARNS1, JAKMIP2, ZNF808, IST1, AL022328.2, AC134407.2, MEG3


