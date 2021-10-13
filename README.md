# R4RA Molecular Analysis

This repository contains scripts to support the publication: 'Deep Molecular Pathology Profiling of Synovial Tissue Biopsies in the R4RA Randomised Clinical Trial Identifies Predictive Signatures of Response/Resistance to Biologic Therapies in Rheumatoid Arthritis'


## Repository Structure


```
.
├── 1_Power_Calculations
├── 2_Crossover Analysis
├── 3_Longitudinal Analysis
└── 4_CDAI50_response_prediction
    ├── 1_Data_Exploration
    ├── 2_Machine_Learning_binary
    └── 3_Plotting_and_summary_scripts
```

- 1_Power_Calculations: statistical power calculations to ensure that RNA-Seq gene expression studies were conducted with enough sequencing depth and sample size. 
- 2_Crossover Analysis: 3-way differential gene expression analysis on baseline synovial biopsies of patients that switched drug over time. Based on their response to different treatments, patients were classified as Pro-RTX, Pro-TOC, and Refractory (see Figure 2g of the paper).
- 3_Longitudinal Analysis: longitudinal differential gene expression analysis using mixed-effects model. Genes that showed expression change over time were enriched via pathway analysis. 
- 4_CDAI50_response_prediction: analysis to predict CDAI50% response to rituximab and tocilizumab using clinical variables and gene expression.
  
  
## Data exploration website

The supporting website is available at: https://r4ra.hpc.qmul.ac.uk/
