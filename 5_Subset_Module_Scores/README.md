# R4RAcell subset module scores


This repo contains code for estimating relative abundance of cell subsets using module scores. 
The following [Zhang et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602051/) and [Alivernini et al.](https://www.nature.com/articles/s41591-020-0939-8) subsets were analyzed:

---

*  SC-F1 = CD34+ sublining
*  SC-F2 = HLA+ sublining
*  SC-F3 = DKK3+ sublining
*  SC-F4 = CD55+ lining

*  SC-M1 = IL1B+ pro-inflammatory
*  SC-M2 = NUPR1+
*  SC-M3 = C1QA+
*  SC-M4 = IFN-activated

*  SC-T1 = CCR7+ CD4+
*  SC-T2 = FOXP3+ Tregs
*  SC-T3 = PD-1+ Tph/Tfh
*  SC-T4 = GZMK+ CD8+
*  SC-T5 = GNLY+ GZMB+
*  SC-T6 = GZMK+/GZMB+

*  SC-B1 = IGHD+ CD270 naive
*  SC-B2 = IGHG3+ CD27- memory
*  SC-B3 = Autoimmune associated
*  SC-B4 = Plasmablasts

---

### Directory Structure

```
.
├── 1_calculate_subset_modules_bulk.R

```

* **1\_Data\_Exploration** contains scripts for module score calculation, statistical comparison of module scores between responders and refractory patients and correlation of module scores with clinical variables


---
