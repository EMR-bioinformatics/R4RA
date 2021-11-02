# R4RA Molecular Analysis

This repository contains scripts to support the publication: 'Deep Molecular Pathology Profiling of Synovial Tissue Biopsies in the R4RA Randomised Clinical Trial Identifies Predictive Signatures of Response/Resistance to Biologic Therapies in Rheumatoid Arthritis'


## Repository Structure


```
.
├── 1_R4RA_power_calculation
├── 2_R4RA_crossover_analysis
├── 3_Subset_Module_Scores
├── 4_WTA_GeoMx_QC_Normalization_DEG
├── 5_R4RA_longitudinal_pathway_analysis
└── 6_CDAI50_response_prediction
    ├── 1_Data_Exploration
    ├── 2_Machine_Learning_binary
    └── 3_Plotting_and_summary_scripts
```

- **1\_R4RA\_power\_calculation**: statistical power calculations to ensure that RNA-Seq gene expression studies were conducted with enough sequencing depth and sample size. 
- **2\_R4RA\_crossover\_analysis (Figure 2g-j)**: 3-way differential gene expression analysis on baseline synovial biopsies of patients that switched drug over time. Based on their response to different treatments, patients were classified as Pro-RTX, Pro-TOC, and Refractory.
- **3\_Subset\_Module\_Scores (Figure 3a)**: cell subset module scores
- **4\_WTA\_GeoMx\_QC\_Normalization\_DEG (Figure 3d-h)**: Digital Spatial Profiling (DSP) data analysis using DESeq2 for preprocessing and DEG analysis. 
- **5\_R4RA\_longitudinal\_pathway\_analysis (Figure 5g-i)**: longitudinal differential gene expression analysis using mixed-effects model. Genes that showed expression change over time were enriched via pathway analysis. 
- **6\_CDAI50\_response\_prediction (Figure 6)**: analysis to predict CDAI50% response to rituximab and tocilizumab using clinical variables and gene expression.

## Additional software packages

For the purposes of this publication an R package, _**glmmSeq**_, was designed to model gene expression with a general linear mixed model (glmm). This is publicly available on [CRAN](https://cloud.r-project.org/web/packages/glmmSeq/index.html). The source code can be found [here](https://github.com/KatrionaGoldmann/glmmSeq). 
  
## Data exploration website

The supporting website is available at: [https://r4ra.hpc.qmul.ac.uk/](https://r4ra.hpc.qmul.ac.uk/)
