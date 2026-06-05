# eDNA Metabarcoding for Monitoring Plant Biodiversity in Michigan Inland Lakes

## Overview

This repository contains the data, code, and documentation for the manuscript titled:

"Seasonal and landscape influences on plant diversity in freshwater ecosystems: Insights from eDNA metabarcoding."

Freshwater ecosystems are essential to global biodiversity and human well-being, yet they face escalating threats from anthropogenic activities. Environmental DNA (eDNA) metabarcoding offers a non-invasive method to monitor plant diversity at the catchment-level, gathering information of freshwater and upland plant communities. This study evaluates the use of eDNA to quantify plant community composition and diversity across 22 lakes inland lakes in Michigan, including aquatic, wetland and terrestrial species.       

## Folder description

### Bioinformatic Analysis
This folder contains scripts used for processing plant *rcbL* metabarcoding data, including sequencing cleaning, alignment, duplicate removal, GenBank retrieval, taxonomy assignment, and OUT clustering using Mothur.

#### Data
- **Plant_rbcL_align_noprimers_092419.fas**: Curated reference alignment of plant *rbcL* sequences used to construct the metabarcoding reference database. Sequences were retrieved from GenBank, aligned with DECIPHER, manually curated, trimmed to the target barcode region, deduplicated, and supplemented with newly generated Sanger sequences.
- **Plant_rbcL_taxonomy_file_092419.txt**: Reference taxonomic file associated with the curated *rbcL* sequence database. Each sequence is linked to a hierarchical classification derived from taxonomic databases and used for Mothur-based taxonomic assignment.

#### Scripts

#### 1. Remove duplicate seqs
This folder contains one script:
- *Rem_Dups_Gaps_By_Species.R*: This R script processes a FASTA alignment to removes= duplicate sequences within each species while retaining the longest representative sequence. The script calculates pairwise distance within species, identify identical sequences (distance = 0), remove redundant sequences, keeps the longest sequence among duplicates, and outputs a filtered FASTA alignment.

#### 2. Sequence alignment with Decipher
This folder contains one script:
- *DECIPHER_AlignTranslation.R*: This R script performs codon-aware multiple sequence alignment of protein-coding DNA sequences using the DECIPHER package. It translates DNA sequences into amino acids, aligns protein sequences to improve accuracy, back-translates aligned sequences to nucleotides, and uses NCBI genetic code 11 (plant plastid code).

#### 3. Sequence retrieval from GenBank
This folder contains R scripts used to query NCBI GenBank for DNA sequences based on species lists and target loci. The scripts automate sequence retrieval, formatting, and compilation for downstream phylogenetic or metabarcoding analyses.
- *Find_NCBI.R*: This script defines the function *FindNCBI()*, which retrieves DNA sequences from NCBI GenBank using the rentrez package. The function searches GenBank using species names and target locus, retrieves accessions IDs and FASTA sequences, writes one FASTA file per species, and generates a summary table of sequence counts per species.
- *RunFunctions_NCBI_051819.R*: This script is a wrapper used to execute *FindNCBI()* for a given species list and retrieve sequences from GenBank for a specific locus. It reads a formatted species list, constructs scientific names, queries GenBank for the target locus, and runs batch sequence retrieval.

#### 4. Taxonomy file
This folder contains an R script used to generate a Mothur-compatible taxonomy file by retrieving taxonomic classifications from the ITIS database using the taxize package.
- *Taxize_taxonomy_file_prep_plant.R*: This script extracts taxonomic classifications for a list of plant sequences and constructs a formatted taxonomy file suitable for downstream use in Mothur. The script parses sequence headers to extract species names, queries ITIS for taxonomic hierarchy, extract ranks (Kingdom to Species), handles missing classifications with NA values, and outputs a Mothur-compatible taxonomy file.

#### 5. Mothur analysis
This folder contains one script:
- *Plant.rbcL.interactive.mothur.script.0diff_LP*: This script defines a complete Mothur-based workflow for processing paired-end plant metabarcoding data (i.e., *rbcL*). It includes read assembly, quality filtering, alignment, chimera removal, taxonomic classification, and OTU clustering.


### Plant Community Analyses 

#### Data
- **FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv**: Plant community matrix including the number of sequences per sample for each taxon across the 22 study lakes. For each sample, the dataset includes lake name, sample name, sample type (surface or benthic), latitude, longitude, and the shortest distance to the shoreline.
- **alpha diversity lakes.csv**: This dataset was used to generate Figure S4 and contains species richness and Shannon-Wiener diversity index values for each plant community type (terrestrial, wetland, and aquatic) across 22 study lakes.
- **lakes traits.csv**: Dataset containing variables describing sampling timing (year, Julian day), lake characteristics (lake area and maximum depth), and surrounding landscape features, including latitude, hydrological connectivity (IWS stream density) and the percentage of different land cover types (agriculture, developed land, forest).
- **prop reads richness groups.csv**: For each lake and plant community (i.e., terrestrial, wetland, and aquatic), the data include the proportion of sequence reads and species richness. This dataset was used to generate Figure 2.
- **spinfo.csv**: Dataset containing species-level information for all species detected through eDNA metabarcoding analysis, including habitat association (terrestrial, wetland, or aquatic), native status (native or alien), and growth habit (e.g., forb, grass, vine, tree). 

### Scripts

#### alpha diversity
- **alpha diversity comparison across plant communities**: This folder contains the script used to compare species richness and Shannon-Wiener diversity indices across terrestrial, wetland, and aquatic plant communities from 22 study lakes and reproduce Figure S4.
- **alpha diversity model selections**: This folder contains the scripts used to evaluate the effects of sampling timing, lake characteristics, and landscape variables on species richness and Shannon–Wiener diversity indices of terrestrial, wetland, and aquatic plant communities using AIC-based model selection. This script reproduces Figure 3.

#### beta diversity
- **homogeneity of variances**: This folder contains the script used to analyze the multivariate homogeneity of group dispersions (variances) across lakes.
- **permanova dbRAD betapart**: This folder contains the script used to perform PERMANOVA, distance-based redundancy analysis (dbRDA), and to partition beta diversity into nestedness and turnover components for terrestrial, wetland, and aquatic plant communities.

#### heat maps
This folder contains two R scripts:
- *CountSpecies.R*: This script generates species richness counts used for the heatmaps.
- *IDWplants_22PanelPlots_2025.R*: This script creates heatmaps for each plant community across the 22 study lakes.

#### minimum number of reads per sample
This folder contains the script used to reproduce Figure S1, which compares the number of sequence reads per sample both combining all lakes and by individual lake. 

#### proportion of sequence reads across plant communities
- **Model selection plant community proportions**: This folder contains the script used to evaluate the effects of lake and landscape characteristics on the proportion of sequence reads of terrestrial, wetland, and aquatic plant communities using AIC-based model selection.
- **Plot richness and proportion of sequence reads per lake**: This folder contains the script for generating stacked bar plots that compare species richness and the proportion of sequence reads for each plant community across all study lakes.

#### within-lake plant community diversity
This folder contains two R scripts:
- *Within_lake_diversity.R*: This script evaluates the effects of distance to shoreline and sample type (surface vs. benthic) on alpha diversity metrics, including species richness and the Shannon–Wiener diversity index, and reproduces Figure 1.
- *Rarefaction curves.R*: This script creates the species accumulation curves for each of the 22 study lakes based on presence-absence data derived from filtered eDNA samples (≥ 1,000 reads) and reproduces Figure S2. 





