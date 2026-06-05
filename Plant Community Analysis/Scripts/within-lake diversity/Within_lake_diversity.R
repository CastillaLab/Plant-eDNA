#############################################################
# LOAD REQUIRED LIBRARIES
#############################################################

library(dplyr)
library(vegan)
library(MASS)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(ggplot2)
library(cowplot)


#############################################################
# LOAD DATA
#############################################################
# Plant community matrix (samples x species counts)

CM<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv",
             header=TRUE)

# Ensure distance to shoreline is numeric
CM$shortest.dist.shoreline<-as.numeric(CM$shortest.dist.shoreline)

#############################################################
# FILTER LOW-SEQUENCING-DEPTH SAMPLES
#############################################################
# Total reads per sample (columns 7–396 = species counts)
CM$total=rowSums(CM[,7:396])

# Inspect sequencing depth distribution
hist(CM$total)
range(CM$total)

# Retain samples with >= 1000 reads
CM.filtered<-filter(CM,total >= 1000)

# Check filtered distribution
hist(CM.filtered$total)
range(CM.filtered$total)

# Remove total reads column
CM.filtered<-CM.filtered %>% dplyr::select(-total)


#############################################################
# CALCULATE DIVERSITY METRICS
#############################################################
# Identify species columns (numeric columns excluding metadata)
species_cols <- CM.filtered %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-latitude) %>%
  dplyr::select(-longitude) %>%
  dplyr::select(-shortest.dist.shoreline) %>%
  colnames()

# Compute species richness and Shannon diversity per sample
CM.filtered <- CM.filtered %>%
  mutate(
    Richness = specnumber(dplyr::select(., all_of(species_cols))),
    Shannon  = diversity(dplyr::select(., all_of(species_cols)), index = "shannon")
  )

#############################################################
# PREPARE VARIABLES FOR MODELING
#############################################################
CM.filtered$lake_name<-as.factor(CM.filtered$lake_name)
CM.filtered$sample<-as.factor(CM.filtered$sample)
CM.filtered$sample_type<-as.factor(CM.filtered$sample_type)

# Log-transform distance (avoid log(0))
CM.filtered$log.shortest.dist.shoreline<-log(CM.filtered$shortest.dist.shoreline+1)


#############################################################
# SPECIES RICHNESS MODELS
#############################################################

# Poisson GLMMs
mod.poisson.a <- glmer(Richness ~ sample_type * log.shortest.dist.shoreline + (1 | lake_name),
            family = poisson(link = "log"),
            data = CM.filtered)

mod.poisson.b <- glmer(Richness ~ sample_type * shortest.dist.shoreline + (1 | lake_name),
             family = poisson(link = "log"),
             data = CM.filtered)

summary(mod.poisson.a)
summary(mod.poisson.b)

# Function to test overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}

# Check overdispersion
overdisp_fun(mod.poisson.a)
overdisp_fun(mod.poisson.b)


#############################################################
# NEGATIVE BINOMIAL GLMMs (HANDLE OVERDISPERSION)
#############################################################
mod.neg.bin.a <- glmmTMB(Richness ~ sample_type * log.shortest.dist.shoreline + (1 | lake_name),
                family = nbinom2,
                data = CM.filtered)

mod.neg.bin.b <- glmmTMB(Richness ~ sample_type * shortest.dist.shoreline + (1 | lake_name),
              family = nbinom2,
              data = CM.filtered)

summary(mod.neg.bin.a)
summary(mod.neg.bin.b)

# Model comparison
AIC(mod.neg.bin.a)
AIC(mod.neg.bin.b)

#############################################################
# PLOT: SPECIES RICHNESS
#############################################################

# Prediction grid (population-level predictions)
newdat <- expand.grid(
  log.shortest.dist.shoreline = seq(min(CM.filtered$log.shortest.dist.shoreline),
                                max(CM.filtered$log.shortest.dist.shoreline),
                                length.out = 100),
  sample_type = levels(CM.filtered$sample_type),
  lake_name = NA
)

# Predictions + SE
pred <- predict(mod.neg.bin.a, newdata = newdat, type = "response", se.fit = TRUE)

# Confidence intervals
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lwr <- newdat$fit - 1.96 * newdat$se
newdat$upr <- newdat$fit + 1.96 * newdat$se

# Plot
p1<-ggplot(newdat, aes(x = log.shortest.dist.shoreline, y = fit, 
                   color = sample_type, fill = sample_type)) +
  geom_jitter(data = CM.filtered, 
              aes(x = log.shortest.dist.shoreline, y = Richness, color = sample_type),
              width = 0.1, alpha = 0.4, size = 1.5, inherit.aes = FALSE) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("benthic" = "darkblue", "surface" = "skyblue")) +
  scale_fill_manual(values = c("benthic" = "darkblue", "surface" = "skyblue")) +
  labs(x = "log(distance to shoreline)",
       y = "Species richness",
       color = "Sample type",
       fill  = "Sample type") +
  theme_classic()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14))+
  theme(legend.position = "none")


#############################################################
# SHANNON DIVERSITY MODELS
#############################################################

mod.shannon.a <- lmer(Shannon ~ sample_type * log.shortest.dist.shoreline 
            + (1 | lake_name), data = CM.filtered)

mod.shannon.b <- lmer(Shannon ~ sample_type * shortest.dist.shoreline 
            + (1 | lake_name), data = CM.filtered)

summary(mod.shannon.a)
summary(mod.shannon.b)

AIC(mod.shannon.a)
AIC(mod.shannon.b)

#############################################################
# PLOT: SHANNON DIVERSITY
#############################################################

# Prediction grid
pred_grid <- expand.grid(
  log.shortest.dist.shoreline = seq(min(CM.filtered$log.shortest.dist.shoreline, na.rm = TRUE),
                     max(CM.filtered$log.shortest.dist.shoreline, na.rm = TRUE),
                     length.out = 100),
  sample_type = levels(CM.filtered$sample_type)
)

# Predictions (fixed effects only)
pred_grid$Shannon_pred <- predict(mod.shannon.a,
                                  newdata = pred_grid,
                                  re.form = NA) # fixed effects only

# Confidence intervals
pred_grid$se <- predict(mod.shannon.a, newdata = pred_grid, re.form = NA, se.fit = TRUE)$se.fit
pred_grid <- pred_grid %>%
  mutate(lower = Shannon_pred - 1.96 * se,
         upper = Shannon_pred + 1.96 * se)

# Plot
p2<-ggplot(CM.filtered, aes(x = log.shortest.dist.shoreline, y = Shannon, color = sample_type)) +
  geom_point(size=2, alpha = 0.5) +  # raw data
  geom_line(data = pred_grid, aes(x = log.shortest.dist.shoreline, y = Shannon_pred, color = sample_type), linewidth = 1) +
  geom_ribbon(data = pred_grid,
              aes(x = log.shortest.dist.shoreline, ymin = lower, ymax = upper, fill = sample_type),
              alpha = 0.2, inherit.aes = FALSE)+
  scale_color_manual(values = c("benthic" = "darkblue", "surface" = "skyblue")) +
  scale_fill_manual(values = c("benthic" = "darkblue", "surface" = "skyblue")) +
  labs(x = "log(distance to shoreline)",
       y = "Shannon-Wiener \ndiversity index",
       color = "Sample Type",
       fill = "Sample Type") +
  theme_classic()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        legend.text=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)))


combinedPlot <- plot_grid(p1, p2,ncol=2, align='hv',
                          labels='auto',label_y= 1,label_x=-0.01,
                          label_size=25, 
                          rel_widths=1, rel_heights= 1)

#############################################################
# COMBINE AND EXPORT FIGURE
#############################################################
png("Figure 1.png",width = 3500, height = 1000, units = "px", res=300)
combinedPlot
dev.off()




