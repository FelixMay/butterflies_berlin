# Load packages
library(tidyverse)
library(vegan)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(GGally)
library(multcomp)
library(multcompView)
library(car)
library(patchwork)
library(permute)
library(indicspecies)
library(devtools)
library(dplyr)
library(ggplot2)

# Installing pairwiseAdonis package from GitHub if not already installed
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Read the data ----------------------------------------------------------------
dat1 <- read_csv("dat_clean_final.csv") 
names(dat1)
str(dat1)

# 1. GLMM analyses of abundance and biodiversity indices -----------------------

#Calculate species richness and diversity indices
dat_abund <- dat1 |> dplyr::select("Pieris_sp.":"Pararge_aegeria")
dat1$N <- rowSums(dat_abund)
dat1$S <- specnumber(dat_abund)
dat1$shannon <- diversity(dat_abund, index = "shannon")
dat1$even <- dat1$shannon / log(dat1$S)

# Statistical modelling

# Preparations of the dataset

# Center and scale numeric predictors
dat2 <- dat1 |> 
  mutate(across(c(cloud_cover, temperature, wind), ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))

# Select dataset without the NAs of wind
dat3 <- dat2 |> filter(!is.na(wind))

# Check correlations among predictors
dat3 |> 
  dplyr::select(cloud_cover, temperature, wind) |> 
  ggpairs()
ggsave("FigA1.png", width = 12, height = 12, units = "cm")

# The correlation is relatively high for wind and cloud cover. Since wind
# will influence butterflies directly, but not cloud cover, we exclude cloud cover
# in the following analyses.

# Model for N > 0 -- detection of any butterflies

# Create binary variable
dat3$N_pos <- as.numeric(dat3$N > 0)

# Check number of sites with and without butterflies
dat3 |> count(N_pos, method)
dat3 |> filter(N_pos == 0)

## 1.1 Fit model for probability of zero vs. non-zero ---------------------------
glmm1_npos <- glmmTMB(N_pos ~ method * green_area_type + wind + 
                      temperature + 
                      (1 | site_id/replicate), family = binomial, data = dat3)


# Check the model
sim_res <- simulateResiduals(glmm1_npos)
plot(sim_res)

# All good :o)

# Hypothesis tests
Anova(glmm1_npos, type = "II")

# Remove non-significant interactions:
glmm2_npos <- glmmTMB(N_pos ~ method + green_area_type + wind + 
                      temperature + 
                      (1 | site_id/replicate), family = binomial, data = dat3)

Anova(glmm2_npos, type = "II")

# Effect sizes and directions
summary(glmm2_npos)

# The probability of detecting any butterflies decreases with increasing wind speed,
# but there are no significant differences between methods or green area types.

# Graphics

# Model predictions for wind speed
new_data <- tibble(wind = seq(min(dat3$wind), max(dat3$wind), length.out = 100),
                   method = "TM",
                   green_area_type = "BC",
                   temperature = 0,
                   site_id = 0) # centered and scaled)
new_data$N_pos_pred <- predict(glmm2_npos, newdata = new_data, type = "response", re.form = NA)

ggplot(dat3, aes(x = wind, y = N_pos)) +
  geom_jitter(height = 0.1, width = 0.1) +
  geom_line(data = new_data, aes(x = wind, y = N_pos_pred), color = "black") +
  theme_bw()

# Well, the effect is there detectable but really weak.

# Data set with positive abundance values only
dat_npos <- dat3 |> filter(N > 0)

## 1.2 Abundance ----

# Fit most complex model
glmm1_n <- glmmTMB(N ~ method * green_area_type + wind + temperature +
               (1 | site_id/replicate), family = poisson, data = dat_npos)

# Check model diagnostics
sim_res <- simulateResiduals(glmm1_n)
plot(sim_res)

# Not perfect, but also OKish. Try with negative binomial

glmm1_n_nb <- glmmTMB(N ~ method * green_area_type +  wind + temperature +
               (1 | site_id/replicate), family = truncated_nbinom2, data = dat_npos)
sim_res <- simulateResiduals(glmm1_n_nb)
plot(sim_res)

# This looks better! So we use the negative binomials distribution  

#Hypothesis tests
Anova(glmm1_n_nb, type = "II")

# Remove insignificant interaction
glmm2_n_nb <- glmmTMB(N ~ method + green_area_type +  wind + temperature +
               (1 | site_id/replicate), family = truncated_nbinom2, data = dat_npos)

Anova(glmm2_n_nb, type = "II")

# We keep the non-significant main effect of wind. Potential justification:
# The interaction between sampling method and green area type was not supported
# and was removed from the final model. 
# Weather covariates were retained to account for variation in butterfly activity and detectability.

# Effect sizes and directions
summary(glmm2_n_nb)

# Post-Hoc Tests

# Pairwise comparisons for green area type
emm_hab <- emmeans(glmm2_n_nb, ~ green_area_type)
pairs(emm_hab, adjust = "tukey", type = "response")

# Pairwise comparisons for method
emm_met <- emmeans(glmm2_n_nb, ~ method)
pairs(emm_met, adjust = "tukey", type = "response")

# Get compact letter display (CLD)
cld_hab <- cld(emm_hab, Letters = letters, decreasing = TRUE)
cld_hab_df <- as.data.frame(cld_hab)

# Figure
ymin <- 0
ymax <- max(dat_npos$N)
yrange <- ymax - ymin

# Relative annotation positions
letters_y_n <- ymin + 1.15 * yrange
method_y  <- ymin + 1.3 * yrange

fig_n <- ggplot(dat_npos,
       aes(x = green_area_type,
           y = N,
           fill = method)) +

  geom_boxplot(position = position_dodge(width = 0.8),
               outliers = FALSE) +

  geom_jitter(aes(color = method),
              position = position_jitterdodge(
                jitter.width = 0.3,
                dodge.width = 0.8),
              show.legend = FALSE,
              alpha = 0.5) +
  geom_text(data = cld_hab_df,
            aes(x = green_area_type,
                y = letters_y_n,
                label = .group),
            inherit.aes = FALSE,
            size = 6) +
  annotate("text",
         x = 2,
         y = method_y,
         label = "Method: p < 0.0001",
         size = 6) +
  labs(
    x = "Green area type",
    y = "Abundance (N > 0)",
    fill = "Method",
    color = "Method"
  ) +
  theme_classic(base_size = 14) +
  coord_cartesian(ylim = c(ymin, method_y))
 
fig_n

## 1.3 Species richness ---------------------------------------------------------

# Fit most complex model
glmm1_sr <- glmmTMB(S ~ method * green_area_type + wind + temperature +
               (1 | site_id/replicate), family = truncated_poisson, data = dat_npos)

# Check model diagnostics
sim_res <- simulateResiduals(glmm1_sr)
plot(sim_res)

# Not perfect, but also OKish. Try with negative binomial
glmm1_sr_nb <- glmmTMB(S ~ method * green_area_type + cloud_cover + wind + temperature +
               (1 | site_id/replicate), family = truncated_nbinom2, data = dat_npos)
sim_res <- simulateResiduals(glmm1_sr_nb)
plot(sim_res)

# This does not help, so we stick with Poisson.  

# Hypothesis tests
Anova(glmm1_sr, type = "II")

# Remove insignificant interaction
glmm2_sr <- glmmTMB(S ~ method + green_area_type + wind + temperature +
               (1 | site_id/replicate), family = truncated_poisson, data = dat_npos)

Anova(glmm2_sr, type = "II")

# Effect sizes and directions

summary(glmm2_sr)
# Post-Hoc Tests

# Pairwise comparisons for green area type
emm_hab <- emmeans(glmm2_sr, ~ green_area_type)
pairs(emm_hab, adjust = "tukey", type = "response")

# Pairwise comparisons for method
emm_met <- emmeans(glmm2_sr, ~ method)
pairs(emm_met, adjust = "tukey", type = "response")

# Get compact letter display (CLD)
cld_hab <- cld(emm_hab, Letters = letters, decreasing = TRUE)
cld_hab_df <- as.data.frame(cld_hab)

# Figure
ymin <- 0
ymax <- max(dat_npos$S) 

yrange <- ymax - ymin

## Relative annotation positions
letters_y_s <- ymin + 1.15 * yrange
method_y  <- ymin + 1.3 * yrange

fig_s <- ggplot(dat_npos,
                aes(x = green_area_type, y = S, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outliers = FALSE,
               alpha = 0.5) +

  geom_jitter(aes(color = method),
              position = position_jitterdodge(
                jitter.width = 0.3,
                dodge.width = 0.8
              ),
              show.legend = FALSE,
              alpha = 0.5) +
  geom_text(data = cld_hab_df,
            aes(x = green_area_type,
                y = letters_y_s,
                label = .group),
            inherit.aes = FALSE,
            size = 6) +
  annotate("text",
         x = 2,
         y = method_y,
         label = "Method: p = 0.07",
         size = 6) +
  labs(
    x = "Green area type",
    y = "Species richness (S > 0)",
    fill = "Method"
  ) +
  guides(color = "none", fill = "none") +
  theme_classic(base_size = 14) +
  coord_cartesian(ylim = c(ymin, method_y))

fig_s

## 1.4 Shannon-diversity -------------------------------------------------------

# Fit most complex model
lmm1_shannon <- glmmTMB(shannon ~ method * green_area_type + wind + temperature +
                        (1 | site_id/replicate), family = gaussian, data = dat_npos)

# Check model diagnostics
sim_res <- simulateResiduals(lmm1_shannon)
plot(sim_res)

# All good :o)

# Hypothesis tests
Anova(lmm1_shannon, type = "II")

# Remove insignificant interaction
lmm2_shannon <- glmmTMB(shannon ~ method + green_area_type + wind + temperature +
                        (1 | site_id/replicate), family = gaussian, data = dat_npos)

Anova(lmm2_shannon, type = "II")

# Effect sizes and directions
summary(lmm2_shannon)

# Post-Hoc Tests
# Pairwise comparisons for green area type

emm_hab <- emmeans(lmm2_shannon, ~ green_area_type)
pairs(emm_hab, adjust = "tukey", type = "response")

# Get compact letter display (CLD)
cld_hab <- cld(emm_hab, Letters = letters, decreasing = TRUE)
cld_hab_df <- as.data.frame(cld_hab)

# Figure
ymin <- 0
ymax <- max(dat_npos$shannon)

yrange <- ymax - ymin

# Relative annotation positions
letters_y_shannon <- ymin + 1.15 * yrange
method_y  <- ymin + 1.3 * yrange

fig_shannon <- ggplot(dat_npos,
       aes(x = green_area_type,
           y = shannon)) +

  geom_boxplot(position = position_dodge(width = 0.8),
               outliers = FALSE) +

  geom_jitter(position = position_jitterdodge(
                jitter.width = 0.3,
                dodge.width = 0.8
              ),
              show.legend = FALSE,
              alpha = 0.5) +
  geom_text(data = cld_hab_df,
            aes(x = green_area_type,
                y = letters_y_shannon,
                label = .group),
            inherit.aes = FALSE,
            size = 6) +
  annotate("text",
         x = 2,
         y = method_y,
         label = "Method: p = 0.31",
         size = 6) +
  labs(
    x = "Green area type",
    y = "Shannon-diversity (S > 0)",
  ) +
  theme_classic(base_size = 14) +
  coord_cartesian(ylim = c(ymin, method_y))

fig_shannon


## 1.4 Evenness ---------------------------------------------------------------

### Select finite evenness only
dat_s2 <- dat_npos |> filter(S > 1)

# Transform  1 values

# This is necessary because thee beta-distribution cannot handle 1 values. 
# We set them to 0.999, which is a common approach.

dat_s2$even_beta <- pmin(dat_s2$even, 0.999)
range(dat_s2$even_beta)

# Fit most complex model
glmm1_even <- glmmTMB(even_beta ~ method * green_area_type +  wind + temperature +
                        (1 | site_id/replicate), family = beta_family, data = dat_s2)

#Check model diagnostics
sim_res <- simulateResiduals(glmm1_even)
plot(sim_res)

# OKish. I have no clue how to improve this further.

# Hypothesis tests
Anova(glmm1_even, type = "II")

# Remove insignificant interaction
glmm2_even <- glmmTMB(even_beta ~ method + green_area_type +  wind + temperature +
                             (1 | site_id/replicate), family = beta_family, data = dat_s2)

Anova(glmm2_even, type = "II")

# Effect sizes and directions
summary(glmm2_even)

# Post-Hoc Tests

# Pairwise comparisons for green area type
emm_hab <- emmeans(glmm2_even, ~ green_area_type)
pairs(emm_hab, adjust = "tukey", type = "response")

# Pairwise comparisons for method
emm_met <- emmeans(glmm2_even, ~ method)
pairs(emm_met, adjust = "tukey", type = "response")

# Get compact letter display (CLD)
cld_hab <- cld(emm_hab, Letters = letters, decreasing = TRUE)
cld_hab_df <- as.data.frame(cld_hab)

# Figure
ymin <- 0.5
ymax <- max(dat_s2$even_beta)

yrange <- ymax - ymin

# Relative annotation positions
letters_y_even <- ymin + 1.15 * yrange
method_y  <- ymin + 1.3 * yrange

fig_even <- ggplot(dat_s2,
       aes(x = green_area_type,
           y = even,
           fill = method)) +

  geom_boxplot(position = position_dodge(width = 0.8),
               outliers = FALSE) +

  geom_jitter(aes(color = method),
              position = position_jitterdodge(
                jitter.width = 0.3,
                dodge.width = 0.8
              ),
              show.legend = FALSE,
              alpha = 0.5) +
  geom_text(data = cld_hab_df,
            aes(x = green_area_type,
                y = letters_y_even,
                label = .group),
            inherit.aes = FALSE,
            size = 6) +
  annotate("text",
         x = 2,
         y = method_y,
         label = "Method: p < 0.0001",
         size = 6) +
  labs(
    x = "Green area type",
    y = "Evenness (S > 1)",
    fill = "Method"
  ) +
  theme_classic(base_size = 14) +
  coord_cartesian(ylim = c(ymin, method_y))

fig_even

## 1.5 Multipanel Figure 2 ------------------------------
combined_plot <-
  ((fig_n | fig_s) /
   (fig_shannon | fig_even)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  plot_layout(guides = "collect")

combined_plot &
  
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("model_results_jitter.png", width = 21, height = 18, units = "cm")
ggsave("model_results_jitter.pdf", width = 21, height = 18, units = "cm")

# 2. Multivariate analysis of species composition --------------------------------

# Aggregate by replicate ID and method
dat_agg <- dat1 |> 
  group_by(replicate, method, green_area_type) |> 
  summarise(across(`Pieris_sp.`:Pararge_aegeria, ~ mean(.x, na.rm = TRUE)), .groups = "drop")

comm <- dat_agg |> dplyr::select(`Pieris_sp.`:Pararge_aegeria)
env <- dat_agg |> dplyr::select(replicate:green_area_type)

# Hellinger transformation
comm_hel <- decostand(
  comm,
  method = "hellinger"
)

## 2.1 NMDS plot ---------------------------------------------------------------
metadata <- dat_agg[, 1:3]

metadata$green_area_type_method <- paste0(
  metadata$green_area_type,
  metadata$method
)

butterfly.NMDS <- metaMDS(comm_hel) # choose ESV level -> no high stress above 0.2 (20%)

zNMDS.df = data.frame(NMDS1 = butterfly.NMDS$points[,1], 
                      NMDS2 = butterfly.NMDS$points[,2],
                      TREAT = metadata$replicate)

# first overview plot
plot(butterfly.NMDS) 
data.scores = as.data.frame(scores(butterfly.NMDS)$sites)

#add columns to data frame 
data.scores$MGAT = metadata$green_area_type_method
data.scores$TF = metadata$replicate

# second dataframe to extract and plot centroid
GROUP <- metadata$green_area_type_method
data.scores$GROUP <- GROUP # group for computing the group centroid

scrs <- scores(butterfly.NMDS, display = 'sites') # extract sample nmds scores
scrs <- cbind(as.data.frame(scrs), GROUP = GROUP)

cent <- aggregate(cbind(NMDS1, NMDS2) ~ GROUP, data = scrs, FUN = mean) # centroid
segs <- merge(scrs, setNames(cent, c('GROUP','oNMDS1','oNMDS2')),
              by = 'GROUP', sort = FALSE) # data for segment between centroid and samples

# give colors "#505050"
zpal <- c("#F18309","#E9DC65", "#89B975", "#95AF2E","#537799","#94BCEF")

# plot NMDS
ggplot(data = scrs, aes(NMDS1, NMDS2, color = GROUP)) +
  geom_point(size = 1) +                                    # sample scores
  geom_segment(data = segs,mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent, size = 3) +    # centroids
  coord_fixed() +
  scale_color_manual(name = "method_green_area_type", values = zpal) +
  theme_bw()

## 2.2 PERMANOVA with permutation within sites only ----------------------------

# Testing differences in species composition using PERMANOVA

# Check dispersion within green area types
bd <- betadisper(
  vegdist(comm_hel, "bray"),
  env$green_area_type
)

anova(bd) # This test does indicate that there is no problem with the dispersion 
# of the groups, so we can proceed with PERMANOVA

# Homogeneity of multivariate dispersion among green area types was assessed 
# using betadisper and was not violated (F = ..., p = 0.27).

# use seed to have stable results in the permutation-based analysis
set.seed(567)

# Table 3
adonis2(
  comm_hel ~ green_area_type + method,
  data = env,
  method = "bray",
  permutations = 999,
  by = "margin"
)

# Table 4
pairwise.adonis(
  comm_hel,
  factors = env$green_area_type,
  perm = 9999, #,
  p.adjust.m = "holm"
)

# pairwise.adonis2(
#   comm_hel ~ green_area_type,
#   data = env,
#   perm = 9999,
#   p.adjust.m = "holm"
# )

# Species abundance data were Hellinger-transformed prior to community analyses 
# to emphasize differences in relative species composition. Community composition
# was visualized using non-metric multidimensional scaling (NMDS) based on
# Bray-Curtis dissimilarities and tested using PERMANOVA with Bray-Curtis dissimilarities.

# 3. Indicator species analysis -------------------------------------------------
set.seed(567)
indic <- multipatt(comm, env$green_area_type, duleg = T) # !!! consider factor combination with duleg = T
indic.tab <- indic$sign
indic.tab.signif <- indic.tab[indic.tab$p.value<=0.05,]
summary(indic)
