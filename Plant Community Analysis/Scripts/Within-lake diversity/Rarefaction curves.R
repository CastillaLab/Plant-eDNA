################################################################################
# Script: Species accumulation curves for plant metabarcoding data
#
# Description:
# This script calculates and visualizes species accumulation curves for
# metabarcoding samples collected across multiple lakes.
#
# The workflow:
#   1. Imports metabarcoding read-count data
#   2. Filters out samples with fewer than 1000 sequencing reads
#   3. Converts read counts to presence/absence data
#   4. Calculates species accumulation curves for each lake using random
#      permutations
#   5. Generates a multi-panel figure containing one accumulation curve
#      per lake
#
#
# Required input file:
#   - FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv
#
################################################################################



#######################################
# Load required packages
#######################################
library(vegan)
library(tidyverse)

#######################################
# Import plant community matrix
#######################################
df <- read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv")

######################################################
#  FILTERING OUT SAMPLES WITH FEWER THAN 1000 READS  #
######################################################

#######################################
# Calculate total reads per sample
#######################################

# Sum sequencing reads across all species columns
# Species data are assumed to begin at column 7
df$total <- rowSums(df[,7:396])

#######################################
# Remove low-read samples
#######################################

# Retain only samples with at least 1000 reads
# Remove temporary "total" column afterward
df.filtered <- df %>% filter(total >= 1000) %>% select(-total)

#######################################
# Convert abundance data to presence/absence
#######################################

# Create a copy of filtered dataset
df.filtered_pa <- df.filtered

# Convert all species read counts:
#   > 0  -> TRUE (species detected)
#   = 0  -> FALSE (species absent)
df.filtered_pa[, 7:ncol(df.filtered_pa)] <- df.filtered_pa[, 7:ncol(df.filtered_pa)] > 0


#######################################
# Split data by lake
#######################################

# Create a list where each element contains
# presence/absence data for a single lake
lake_list <- split(df.filtered_pa[, 7:ncol(df.filtered_pa)], df.filtered_pa$lake_name)



#######################################
# Calculate species accumulation curves
#######################################

# For each lake:
#   - Calculate species accumulation curves
#   - Use random sample ordering
#   - Perform 100 permutations
acc_data <- lapply(names(lake_list), function(lake_name) {
  acc <- specaccum(lake_list[[lake_name]], method = "random", permutations = 100)
  data.frame(
    Lake = lake_name,
    Samples = acc$sites,
    Richness = acc$richness,
    SD = acc$sd
  )
})

#######################################
# Combine all lake results
#######################################

# Merge all accumulation data frames
# into a single table for plotting
acc_df <- bind_rows(acc_data)


################################################################################
# GENERATE FIGURE WITH MULTIPLE PANELS
################################################################################

#######################################
# Create species accumulation figure
#######################################

p_all <- ggplot(acc_df, aes(x = Samples, y = Richness)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD),
              alpha = 0.2, fill = "grey70") +
  facet_wrap(~ Lake, scales = "free") +
  theme_classic() +
  labs(
    title = "Species Accumulation Curves per Lake",
    x = "Number of Samples",
    y = "Species Richness"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(1, "lines")
  )


#######################################
# Save figure
#######################################
ggsave("Figure S2 - Accumulation_curves_22_lakes.png",
       plot = p_all,
       width = 18, height = 12, dpi = 300)
