# Sipes-2025-Spatial-Transcriptome-Fallopian-Tube

## Summary

This R Project was created for the paper "Spatial Transcriptomic Profiling of the Human Fallopian Tube Epithelium Reveals Region-specific Gene Expression Patterns" to be published in _Communications Biology_ in 2025. 

This project takes data generated on the GeoMx DSP platform and performs analysis using code produced by Nanostring. 
This analysis is based on the following Bioconductor vignette:
[Analyzing GeoMx-NGS RNA Expression Data with GeomxTools (bioconductor.org)](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#5_Normalization)


## Overview of Main Folders and Files
A brief summary of the contents of this project:

Files:
* **data_analysis.Rmd** - contains preliminary QC, Filtering, and Normalization of the GeoMx Dataset
* **main_figures.Rmd** - contains the code used to generate the main figures for this paper
* **supplementary_figures.Rmd** - contains the code used to generate the supplemental figures for this paper
* **functions.R** - contains helper functions for generating graphs that were taking up too much space in the main markdown documents. 

Folders:
* **/data_input** - contains all input files for the GeoMx platform, including:
  * _/dccs_ - contains .DCC files with probe counts (one .dcc per segment)
  * _/pkcs_ - contains .pkc file with probe information
  * _/annotations_ - contains an excel document with annotation information
* **/stat_analysis** - a file containing excel files with results of statistical analyses performed by a collaborator
* **/pictures** - figures generated on Biorender needed for creating some figures
* **/IHC_analysis** - contains excel document with the results of IHC analysis performed by Pathologist Dr. Rashna Madan, needed for generating some figures
* **/GO_analysis** - contains files needed for generating GO plots
  * /Source Files - contains source files used to generate GO files
  * /All_GO - contains results of GO analysis when considering all data
  * /Discovery_GO - results of GO analysis using only the discovery cohort
  * /Validation_Go - results of GO analysis using only the validation cohort
* **cell_markers** - contains markers of cell types identified in Ulrich 2022 paper



