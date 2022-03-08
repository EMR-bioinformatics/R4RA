# R4RA machine learning


This repo contains code for a nested machine learning pipeline used to predict binary (responder vs non-responder) [CDAI50](https://pubmed.ncbi.nlm.nih.gov/22454398/) response in the R4RA cohort. 

---

### Directory Structure

```
.
├── 1_Data_Exploration
├── 2_Machine_Learning_binary
└── 3_Plotting_and_summary_scripts

```

* **1\_Data\_Exploration** contains scripts for curating, analysing, and exploring the features and response variables.  
* **2\_Machine\_Learning_binary** contains scripts for performing machine learning predictions on binary/discrete response metrics. This is split into classic cross-validation and nested cross-validation pipelines. 
* **3\_Plotting\_and\_summary\_scripts** contains scripts for downstream visualisations such as ROC curves and variable importance plots. 



---

### Model Features 

Model features should include: 

* gene expression
* histology
* clinical parameters

We implement a number of feature selection/filtering methods including: 

* remove co-linear features
* subsetting to differentially expressed genes (DEG)
* univariate filtering see [here](https://topepo.github.io/caret/feature-selection-using-univariate-filters.html)
* recursive feature elimintaion (RFE) see [here](https://topepo.github.io/caret/recursive-feature-elimination.html)

---

### Model Predictors

Models were set up to predict two response variables: 

1. Final/Secondary CDAI50 response (at visit 9, post treatment cross-over)
1. Initial/Primary CDAI50 response (at visit 7, prior to treatment cross-over)

As well as response to different treatments: 

1. Rituximab
1. Tocilizumab
1. Refractory (response to either treatment, at final response only)

---

### Objectives

To predict binary CDAI50 response with R's caret package using _**nested**_ cross-validation. 
    
For more information on nested cross-validation with caret see:

* [SO discussion](https://stats.stackexchange.com/questions/125843/outer-crossvalidation-cycle-in-caret-package-r/126131)
* [Applied Predictive Modelling in R](http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd)  
 