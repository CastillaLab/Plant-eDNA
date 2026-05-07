# eDNA Metabarcoding for Monitoring Plant Biodiversity in Michigan Inland Lakes

## Overview

This repository contains the data, code, and documentation for the manuscript titled:

"Seasonal and landscape influences on plant diversity in freshwater ecosystems: Insights from eDNA metabarcoding."

Freshwater ecosystems are essential to global biodiversity and human well-being, yet they face escalating threats from anthropogenic activities. Environmental DNA (eDNA) metabarcoding offers a non-invasive method to monitor plant diversity at the catchment-level, gathering information of freshwater and upland plant communities. This study evaluates the use of eDNA to quantify plant community composition and diversity across 22 lakes inland lakes in Michigan, including aquatic, wetland and terrestrial species.       

## Folder description

### DATA

- **FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv**: Plant community matrix including the number of sequences per sample for each taxon across the 22 study lakes. For each sample, the dataset includes lake name, sample name, sample type (surface or benthic), latitude, longitude, and the shortest distance to the shoreline.
- **alpha diversity lakes.csv**: This dataset was used to generate Figure S4 and contains species richness and Shannon-Wiener diversity index values for each plant community type (terrestrial, wetland, and aquatic) across 22 study lakes.
- **lakes traits.csv**: Dataset containing variables describing sampling timing (year, Julian day), lake characteristics (lake area and maximum depth), and surrounding landscape features, including latitude, hydrological connectivity (IWS stream density) and the percentage of different land cover types (agriculture, developed land, forest).
- **prop reads richness groups.csv**: For each lake and plant community (i.e., terrestrial, wetland, and aquatic), the data include the proportion of sequence reads and species richness. This dataset was used to generate Figure 2.
- **spinfo.csv**: Dataset containing species-level information for all species detected through eDNA metabarcoding analysis, including habitat association (terrestrial, wetland, or aquatic), native status (native or alien), and growth habit (e.g., forb, grass, vine, tree). 

### SCRIPTS

#### WITHIN-LAKE PLANT COMMUNITY DIVERSITY
This folder contains two R scripts:
- *Within_lake_diversity.R*: This script evaluates the effects of distance to shoreline and sample type (surface vs. benthic) on alpha diversity metrics, including species richness and the Shannon–Wiener diversity index, and reproduces Figure 1.
- *Rarefaction curves.R*: This script creates the species accumulation curves for each of the 22 study lakes based on presence-absence data derived from filtered eDNA samples (≥ 1,000 reads) and reproduces Figure S2. 

#### ALPHA DIVERSITY

- **ALPHA DIVERSITY COMPARISON ACROSS PLANT COMMUNITIES**: This folder contains the script used to compare species richness and Shannon-Wiener diversity indices across terrestrial, wetland, and aquatic plant communities from 22 study lakes and reproduce Figure S4.
- **ALPHA DIVERSITY MODEL SELECTION**: This folder contains the scripts used to evaluate the effects of sampling timing, lake characteristics, and landscape variables on species richness and Shannon–Wiener diversity indices of terrestrial, wetland, and aquatic plant communities using AIC-based model selection. This script reproduces Figure 3.

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



