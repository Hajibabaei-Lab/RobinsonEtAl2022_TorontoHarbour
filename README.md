# README

This repository contains the dataflow and scripts used to process the COI and 18S metabarcode reads in the paper Robinson et al., 2022 (Scientific Reports, Accepted).

Infiles, R scripts, and supplementary data files can be downloaded from https://github.com/Hajibabaei-Lab/RobinsonEtAl2022_TorontoHarbour/releases .

## Infiles

SitesCSOs.csv  
metadata.csv  
TH_extra_data.txt  
diversity.metrics.csv  
results_F230R_10Aug20.csv  
results_ml-jg_10Aug20.csv  
morph_new.csv  
results_18S_v4.1.csv  
EPA_families.csv  
Freshwaterecol_Moog.csv  
Freshwater_Tachet.csv  

## R Scripts

To try to improve reproducibility, script libraries are loaded using groundhog but libraries can be loaded the usual way with only minor changes.

Fig1_mapped_gradients.R uses SitesCSOs.csv, metadata.csv, TH_extra_data.txt, and diversity.metrics.csv (from Fig5_hierarchical_partitioning.R) to produce Fig1_mapped_gradients.pdf .

Fig2_richness.R uses results_F230R_10Aug20.csv, results_ml-jg_10Aug20.csv, morph_new.csv, and results_18S_v4.1.csv to produce Fig2_richness.pdf .

Fig3_PCA_NMDS.R uses results_F230R_10Aug20.csv, results_ml-jg_10Aug20.csv, momrph_new.csv, SitesCSOs.csv, metadata.csv, TH_extra_data.txt, and results_18S_v4.1.csv to produce Fig3_PCA_NMDS.pdf as well as a suite of supporting files.

Fig4_relative_abundance.R uses results_F230R_10Aug20.csv, results_ml-jg_10Aug20.csv, morph_new.csv, results_18S_v4.1.csv, EPA_families.csv, Freshwaterecol_Moog.csv, and Freshwater_Tachet.csv to produce Fig4_relative_abundance.pdf.

Fig5_hierarchical_partitioning.R uses results_F230R_10Aug20.csv, results_ml-jg_10Aug20.csv, morph_new.csv, results_18S_v4.1.csv, metadata.csv, and TH_extra_data.txt to produce Fig5_hierarchical_partitioning.pdf, the results summarized in Table S5 and Table S6, and a couple other supporting files.

## Supplementary data files

FigS1_resolution_rarefaction.R uses results_F230R_10Aug20.csv, results_ml-jg_10Aug20.csv, results_18S_v4.1.csv, and morph_new.csv to produce FigS1_resolution_rarefaction.pdf.  

FigS2_predictor_correlations.R uses metadata.csv and TH_extra_data.txt to produce FigS2_predictor_correlations.pdf.

TableS3_checklist.R uses morph_new.csv, results_F230R_10Aug20.csv, results_ml-jg_10Aug20.csv, and results_18S_v4.1.csv to produce checklists that were compiled into Table S3 (excel spreadsheet).

## References

Robinson, C.V., Porter, T.M., McGee, K.M., McCusker, M., Wright, M.T.G., Hajibabaei, M. 2022.  Multi-marker DNA metabarcoding detects suites of environmental gradients from an urban harbour.  Scientific Reports, Accepted.


Last updated: May 12, 2022
