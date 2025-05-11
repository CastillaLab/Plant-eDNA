# eDNA Metabarcoding for Monitoring Plant Biodiversity in Michigan Inland Lakes

## Overview

This repository contains the data, code, and documentation for the manuscript titled:

"Seasonal and landscape influences on plant diversity in freshwater ecosystems: Insights from eDNA metabarcoding."

Freshwater ecosystems are essential to global biodiversity and human well-being, yet they face escalating threats from anthropogenic activities. Environmental DNA (eDNA) metabarcoding offers a non-invasive method to monitor plant diversity at the catchment-level, gathering information of freshwater and upland plant communities. This study evaluates the use of eDNA to quantify plant community composition and diversity across 22 lakes inland lakes in Michigan, including aquatic, wetland and terrestrial species.       

## Folder description

### DATA

- **FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv**: Plant community matrix including the number of sequences per sample for each taxon across the 22 study lakes. Data also include for each sample the lake and whether the sample was collected from the surface or benthic area.
- **Alpha Diversity Comparison Across Plant Communities**: This dataset was used to generate Figure S4 and contains species richness and Shannon-Wiener diversity index values for each plant community type (terrestrial, wetland, and aquatic) across 22 study lakes.
- **prop.reads.richness.groups.csv**: For each lake and plant community (i.e., terrestrial, wetland, and aquatic), the data include the proportion of sequence reads and species richness.

### SCRIPTS

#### ALPHA DIVERSITY

- **ALPHA DIVERSITY COMPARISON ACROSS PLANT COMMUNITIES**: This folder contains the script used to compare species richness and Shannon-Wiener diversity indices across terrestrial, wetland, and aquatic plant communities from 22 study lakes.
- **ALPHA DIVERSITY MODEL SELECTION**:This folder contains the script used to evaluate the effects of lake and landscape characteristics on species richness and Shannon-Wiener diversity indices for terrestrial, wetland, and aquatic plant communities, using AIC-based model selection.

#### BETA DIVERSITY

- **HOMOGENEITY OF VARIANCES**: This folder contains the script used to analyze the multivariate homogeneity of group dispersions (variances) across lakes.
- **PERMANOVA_dbRDA_BETAPART**: This folder contains the script used to perform PERMANOVA, distance-based redundancy analysis (dbRDA), and to partition beta diversity into nestedness and turnover components for terrestrial, wetland, and aquatic plant communities.

#### HEAT MAPS

This folder contains two R scripts:
- *CountSpecies.R*: This script generates species richness counts used for the heatmaps.
- *IDWplants_22PanelPlots_2025.R*: This script creates heatmaps for each plant community across the 22 study lakes.

#### MINIMUM NUMBER OF READS PER SAMPLE

This folder contains the script used to reproduce Figure S1, which compares the number of sequence reads per sample both combining all lakes and by individual lake. 

#### PROPORTION OF SEQUENCE READS ACROSS PLANT COMMUNITIES

- **PLOT RICHNESS AND PROPORTION OF SEQUENCE READS PER LAKE**: This folder contains the script for generating stacked bar plots that compare species richness and the proportion of sequence reads for each plant community across all study lakes.
- **MODEL SELECTION PLANT COMMUNITY PROPORTIONS**: This folder contains the script used to evaluate the effects of lake and landscape characteristics on the proportion of sequence reads of terrestrial, wetland, and aquatic plant communities using AIC-based model selection.



