# Load packages
library(vegan)
library(tidyverse)

# Read your data
df <- read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv")

######################################################
#  FILTERING OUT SAMPLES WITH FEWER THAN 1000 READS  #
######################################################

# Total reads
df$total <- rowSums(df[,7:396])

# Filter samples with >= 1000 reads
df.filtered <- df %>% filter(total >= 1000) %>% select(-total)

# Convert to presence/absence
df.filtered_pa <- df.filtered
df.filtered_pa[, 7:ncol(df.filtered_pa)] <- df.filtered_pa[, 7:ncol(df.filtered_pa)] > 0

# Split data by lake
lake_list <- split(df.filtered_pa[, 7:ncol(df.filtered_pa)], df.filtered_pa$lake_name)

# Build accumulation-curve data frame for each lake
acc_data <- lapply(names(lake_list), function(lake_name) {
  acc <- specaccum(lake_list[[lake_name]], method = "random", permutations = 100)
  data.frame(
    Lake = lake_name,
    Samples = acc$sites,
    Richness = acc$richness,
    SD = acc$sd
  )
})

acc_df <- bind_rows(acc_data)

##################################################
#            ONE FIGURE WITH 22 PANELS
##################################################

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

# Save the figure
ggsave("Figure S2 - Accumulation_curves_22_lakes.png",
       plot = p_all,
       width = 18, height = 12, dpi = 300)
