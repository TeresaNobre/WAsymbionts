require(tidyverse)
require(magrittr)
require(lme4)

# odds = p/(1-p)
# log(odds)
# link = log(p/(1-p))
# inverse link or logistic = 1/(1 + exp(-x))




# 1. Data ####
# 1.1 Load data ####
data <- read_csv("~/Desktop/WAsymbionts/presence_absence.csv") %T>%
  print()

# Add proportion column
data %<>%
  mutate(proportion = obs/n) %T>%
  print()

# Sanity check
data %>% filter(proportion > n)

# Exploratory GLMs
# mod <- glm(obs/n ~ lat * clade, weights = n,
#           family = binomial(link = "logit"), 
#          data = data)
# mod
# 
# mod <- glm(cbind(obs, n - obs) ~ lat * clade,
#            family = binomial(link = "logit"), 
#            data = data)
# mod

# Expand data to long format: one row per individual sample (presence/absence)
data_long <- data %>%
  rowwise() %>%
  mutate(presence = list( c( rep(1, obs) , rep(0, n - obs) ) )) %>%
  unnest(presence) %T>%
  print()

# Center latitude at -22°S (Ningaloo) so model intercept is interpretable
data_long %<>%
  mutate(lat_c = lat - (-22)) %T>% 
  print() 

# 2. Models ####
# 2.1  Clade-level GLM ####
# Models presence/absence of each symbiont clade across latitude, allowing
# each clade to have its own baseline presence and latitudinal response

# Fit GLM with interaction between latitude and clade
mod_clade <- glm(presence ~ lat_c * clade, 
                 family = binomial(link = "logit"), 
                 data = data_long %>%
                   mutate(clade = clade %>% fct_relevel("C")))

mod_clade_simple <- glm(presence ~ clade, 
                 family = binomial(link = "logit"), 
                 data = data_long %>%
                   mutate(clade = clade %>% fct_relevel("C")))

mod_clade_simpleA <- glm(presence ~ clade, 
                        family = binomial(link = "logit"), 
                        data = data_long %>%
                        mutate(clade = clade %>% fct_relevel("A")))

summary(mod_clade_simpleA)

mod_clade_simpleB <- glm(presence ~ clade, 
                         family = binomial(link = "logit"), 
                         data = data_long %>%
                           mutate(clade = clade %>% fct_relevel("B")))

summary(mod_clade_simpleB)

mod_clade_simpleD <- glm(presence ~ clade, 
                         family = binomial(link = "logit"), 
                         data = data_long %>%
                           mutate(clade = clade %>% fct_relevel("D")))

summary(mod_clade_simpleD)

# 1/(1 + exp(--1.7772)) = probability (-1.7772 is intercept estimate)
# Inv(alpha) = prevalence at baseline as a proportion/probability

# OUTCOME: 
# - Clade C (reference): 51% presence at -22°S, increases 15% per degree north (p < 0.001)
# - Clade D: Lower baseline than C (p < 0.001), no different latitudinal trend
# - Clade B: Similar baseline to C, similar latitudinal trend
# - Clades A, F, G: Too rare for reliable slope estimates (NA for interactions)
# - Model explains 7% of deviance (Null: 7456 → Residual: 6954)

# Data structure check
# View Clade C rows without haplotype-level data (n = 181)
data_long %>%
  filter(clade == "C" & is.na(haplotype))
  
# View Clade C rows with haplotype data (n = 4,727)
data_long %>%
  filter(clade == "C") %>%
  drop_na(haplotype)

# MODEL SUMMARY
# with lat: 37.381
# with lat_c: -1.77803 inv(-1.77803)= 0.1445466; 14.5% likelihood of finding C in Ningaloo

inv(0.043235)

summary(mod_clade)

# KEY RESULTS:
# Intercept (Clade C at -22°S): 0.043 → 51% presence probability
# lat_c (Clade C slope): 0.141 → 15% odds increase per degree north (highly significant ***)
# cladeD: -1.918 → Much lower presence than C (highly significant ***)
# cladeB: 0.402 → Slightly higher presence than C (marginal)
# cladeG: -1.505 → Lower presence than C (significant *)
# lat_c:cladeF & lat_c:cladeG: NA due to singularities (too few observations)

# GENERATE PREDICTIONS ACROSS LATITUDE FOR EACH CLADE
new_clade <- expand_grid(
  lat_c = data_long %$% seq(min(lat_c), max(lat_c), 0.1), # generate new data
  clade = data_long %>% distinct(clade) %>% pull(clade)
) %T>%
  print()
# RESULT: 1,284 rows (214 latitude values × 6 clades)

# Get inverse logit function to convert log-odds to probabilities
inv <- family(mod_clade)$linkinv # calculate inverse of link function

# Generate predictions with standard errors
predicted <- predict(mod_clade, type = "link", se.fit = TRUE, newdata = new_clade) # predict
# Warning: "rank-deficient fit" because Clades F & G have insufficient data for interaction terms

# Convert predictions to probability scale
new_clade$fit <- inv( predicted$fit )
new_clade$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new_clade$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new

# 2.2 Global clade-level GLMM ####
# GLOBAL CLADE-LEVEL MIXED EFFECTS MODEL
# Models presence/absence of all clades simultaneously, allowing each clade
# to have its own baseline presence and its own latitudinal response
# (random intercept + random slope per clade)

mod_global <- glmer(presence ~ lat + (lat|clade), 
                  family = binomial(link = "logit"), 
                  data = data_long)
# WARNING: "boundary (singular) fit" - the model is overfitted

# Fixed effects:
  #   Intercept: -0.013 → ~50% presence at equator(not centered lat)
  #   lat: 0.061 → Small positive trend northward on average
# Random effects:
  #   Intercept SD: 3.11 → HUGE variation in baseline presence between clades
  #   Slope SD: 0.081 → Moderate variation in latitudinal responses between clades

mod_global

# Exploratory
# mmod <- glmer(presence ~ lat * clade + (lat|ref), # intercept only + (1|ref)
#               family = binomial(link = "logit"), 
#               data = data_long)
# mmod

# mod <- glm(presence ~ lat * haplotype, 
#            family = binomial(link = "logit"), 
#            data = data_long %>% drop_na(haplotype))
# 
# mod

# Prediction grid: all combinations of latitude and clade for smooth prediction curves
inv <- family(mod_global)$linkinv  # calculate inverse of link function
new <- expand_grid(
  lat = data_long %$% seq(min(lat), max(lat), 0.1), # generate new data
  clade = data_long %>% distinct(clade) %>% pull(clade)
) %T>%
  print()

predicted <- predict(mod_global, type = "link", se.fit = TRUE, newdata = new) # predict
new$fit <- inv( predicted$fit )
new$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new

# 2.2.1 Visualisation ####
new %>%
  # filter(clade %in% c("C", "B", "D")) %>%
  ggplot() +
    geom_line(aes(lat, fit, colour = clade)) +
    geom_ribbon(aes(lat, ymin = lower, ymax = upper, fill = clade),
                alpha = 0.2) +
    geom_point(data = data, # %>% filter(clade %in% c("C", "B", "D")),
               aes(lat, proportion, colour = clade), shape = 16, alpha = 0.2) +
    facet_grid(~ clade) +
    coord_flip() +
    theme_minimal() +
    theme(legend.position = "none")

# some indication of increase of prevalence of c and b with latitude and no hange in d but most observations at centre latitudes and uneven data spread, 
# and too few A,B, F, G observation to make a reliable predition

# 2.3 Haplotype-level GLM (Clade C only) ####
# Models presence/absence of each Clade C haplotype across latitude, allowing
# each haplotype to have its own baseline presence AND latitudinal response
# Only includes rows with haplotype-level identification (drops clade-only studies)

mod_hap <- glm(presence ~ lat_c * haplotype,         # Latitude centered at 22ºS
                 family = binomial(link = "logit"), 
                 data = data_long %>% 
                 filter(clade == "C") %>%
                 drop_na(haplotype)) 

# Prediction grid: 
# Only for Clade C haplotypes with identification data
inv <- family(mod_hap)$linkinv  # inverse logit function (converts log-odds → probabilities)
new <- expand_grid(
  lat_c = data_long %>% 
    filter(clade == "C") %>%
    drop_na(haplotype) %$% 
    seq(min(lat_c), max(lat_c), 0.1),
  haplotype = data_long %>% 
    filter(clade == "C") %>%
    drop_na(haplotype) %>% 
    distinct(haplotype) %>% 
    pull(haplotype) # All 56 Clade C haplotypes
)%T>%
  print()

# Predict on link (log-odds) scale with standard errors
predicted <- predict(mod_hap, type = "link", se.fit = TRUE, newdata = new) 
new$fit <- inv( predicted$fit )
new$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

# Back-calculate actual latitude from centered latitude
new %<>%
  mutate(lat = lat_c - 22) %T>%
  print()

# 2.3.1 Visualisation ####

# All haplotypes on single panel (exploratory - 56 lines, no legend)
new %>%
  ggplot() +
  geom_line(aes(lat, fit, colour = haplotype), alpha = 0.5) +
  geom_point(data = data,
             aes(lat, proportion, colour = haplotype), shape = 16, alpha = 0.2) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none")


# Plot each individual haplotype
ggplot() +
  geom_line(data = new, aes(x = lat, y = fit, color = haplotype)) +
  geom_ribbon(data = new, aes(x = lat, ymin = lower, ymax = upper, fill = haplotype), 
              alpha = 0.2) +
  facet_wrap(~haplotype) +
  theme_minimal() +
  labs(
    title = "Haplotype presence probability across latitude",
    x = "Latitude (°S)",
    y = "Probability of presence"
  )


# 2.4 Clade C GLMM - population-level predictions  ####
# GLMM with random intercept + slope per haplotype
# Allows each haplotype its own baseline prevalence and latitudinal response
# Includes 181 clade-only rows (no haplotype ID) unlike mod_hap
# Best balance of fit and stability (AIC: 4616, ΔAIC = 1913 vs simple GLM)

mod_global_c <- glmer(presence ~ lat_c + (lat_c|haplotype),
                      family = binomial(link = "logit"),
                      data = data_long %>% 
                      filter(clade == "C"))
mod_global_c

# OUTCOME:
# Fixed effects (average across all 56 haplotypes):
#   Intercept: -0.211 → 45% average presence at -22°S because exp(-0.21) / (1 + exp(-0.21)) = 0.45 (45%)
#   lat_c: -0.095 → Slight decrease northward on average
# Random effects:
#   Intercept SD: 2.006 → Huge variation in baseline presence between haplotypes
#   Slope SD: 0.767 → Large variation in latitudinal responses between haplotypes
#   Correlation: 0.47 → Common haplotypes tend to increase more northward
# AIC: 4616 | n = 4,727 obs across 56 haplotypes

# Population-level prediction grid (no haplotype column - predicts average haplotype)
inv <- family(mod_global_c)$linkinv
new_mod_global_c <- expand_grid(
  lat_c    = data_long %>%
    filter(clade == "C") %$%
    seq(min(lat_c), max(lat_c), 0.1),
  presence = 0   # Placeholder required by model.matrix below
) %T>%
  print()
# RESULT: 214 rows (one per 0.1° latitude step)

# Predict on log-odds scale, ignoring random effects (population-level average)
# re.form = NA means "predict for an average/new haplotype not in the data"
new_mod_global_c$fit <- predict(mod_global_c, 
                                type = "link", re.form = NA, 
                                newdata = new_mod_global_c) # fit on link scale

# Propagate BOTH sources of uncertainty into confidence intervals:
# 1. Fixed effect uncertainty (precision of average trend estimate)
# 2. Random effect variance (spread of haplotypes around average)
mm   <- model.matrix(terms(mod_global_c), new_mod_global_c)  # Fixed effect model matrix
fvar <- diag(mm %*% tcrossprod(vcov(mod_global_c), mm))       # Fixed effect variance
tvar <- fvar + VarCorr(mod_global_c)$haplotype[1]             # Total variance (fixed + random intercept)

# Convert to probability scale with 95% prediction intervals
new_mod_global_c$upper <- inv(new_mod_global_c$fit + qnorm(0.975)*sqrt(tvar))
new_mod_global_c$lower <- inv(new_mod_global_c$fit - qnorm(0.975)*sqrt(tvar))
new_mod_global_c$fit <- inv(new_mod_global_c$fit)

# Back-calculate actual latitude from centered latitude
new_mod_global_c %<>%
  mutate(lat = lat_c - 22) %T>%
  print()

# To isolate predictions for particular coordinates (e.g. -15º lat)
new_mod_global_c %>%
  filter(lat == -15)

# OUTCOME: 214 rows with population-level predictions:
# At -33.6°S: fit = 71%, CI = (2.4%, 99.6%) → wide CI reflects haplotype variation
# At -22°S:   fit = 45% → average presence at reference latitude
# Wide CIs are expected - they represent the full range across all haplotypes

# 2.4.1 Visualisation ####

degE <- function(x) paste0(abs(x), "°E")
degS <- function(y) paste0(abs(y), "°S")

mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2), "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 10, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 10),
                 legend.text.align = 0,
                 legend.title = element_text(size = 10),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 10, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

# Exploratory plot
new_mod_global_c %>%
  ggplot() +
  geom_line(data = new, aes(x = lat, y = fit, group = haplotype), 
            alpha = 0.1, colour = "#71874E") +
  geom_line(aes(lat, fit), linewidth = 1, colour = "#71874E") +
  # geom_line(data = new_lat, aes(lat, fit), linewidth = 1, linetype = 2) + # uncomment to see WRONG prediction
  geom_ribbon(aes(lat, ymin = lower, ymax = upper),
              alpha = 0.2) +
  geom_point(data = data %>% filter(clade == "C") %>% drop_na(haplotype),
             aes(lat, proportion, size = n), 
             shape = 21,
             fill = "#71874E50",    
             colour = "#71874E") +
  coord_flip() +
  theme_minimal()

# FINAL FIGURE PLOT (used in combined Figure 1)
# Thin lines = individual haplotype predictions (from mod_hap, stored in new)
# Bold line = population-level average trend (from mod_global_c)
# Points = observed prevalence, sized by sample size
model_plot <- new_mod_global_c %>%
  ggplot() +
  geom_line(data = new,
            aes(x = lat, y = fit, group = haplotype),
            alpha  = 0.1,
            colour = "#71874E") +
  geom_line(aes(lat, fit),
            linewidth = 1,
            colour    = "#71874E") +
  geom_point(data   = data %>% filter(clade == "C") %>% drop_na(haplotype),
             aes(lat, proportion, size = n),
             shape  = 21,
             fill   = "#71874E50",
             colour = "#71874E") +
  scale_size_continuous(
    range  = c(2, 6),
    breaks = c(5, 10, 20, 30, 60),
    guide  = "none"
  ) +
  scale_x_continuous(
    limits = c(-35, -10),
    breaks = seq(-35, -10, 5),
    labels = degS # Keep this to match map
  ) +
  labs(y = "Prevalence") +
  scale_y_continuous(labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 1))) +
  coord_flip(expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.x = element_text(hjust = 0.5))

model_plot

# Exploration
# 2.5 Model latitude #
# mod_lat <- glm(presence ~ lat_c,
#                   family = binomial(link = "logit"),
#                   data = data_long %>% 
#                   filter(clade == "C"))
# mod_lat
# # (Intercept)        lat_c  
# # 0.04323      0.14066  
# 
# # At the reference latitude (-22°S):
# # Log-odds of presence = 0.04323
# # Probability = exp(0.04323) / (1 + exp(0.04323)) = 0.51 (51% chance)
# # At -22°S, Clade C haplotypes are present in about half the samples
# 
# new_lat <- expand_grid(
#   lat_c = data_long %>% 
#     filter(clade == "C") %$% 
#     seq(min(lat_c), max(lat_c), 0.1)
# )%T>%
#   print()
# 
# inv <- family(mod_lat)$linkinv # calculate inverse of link function
# predicted <- predict(mod_lat, type = "link", se.fit = TRUE, newdata = new_lat) # predict
# 
# new_lat$fit <- inv( predicted$fit )
# new_lat$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
# new_lat$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )
# 
# new_lat %<>%
#   mutate(lat = lat_c - 22) %T>%
#   print()


# Visualise sampling effort across latitude
data_long %>%
  filter(clade == "C") %>%
  group_by(lat, ref) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  ggplot(aes(x = n_obs, y = reorder(factor(round(lat, 1)), lat), fill = ref)) +
  geom_col() +
  scale_y_discrete(labels = function(y) paste0(abs(as.numeric(y)), "°S")) +
  scale_fill_manual(values = c(
    "#A8B5A2", "#C9A96E", "#7B9E87", "#D4956A",
    "#8FA8BF", "#C4876A", "#6B8C6B"
  )) +
  labs(x = "Number of observations",
       y = "Latitude (°S)",
       fill = "Study",
       title = "Sampling effort across latitude") +
  theme_minimal()

# Check how many unique sites per latitude band
data_long %>%
  filter(clade == "C") %>%
  mutate(lat_band = round(lat, 0)) %>%
  group_by(lat_band) %>%
  summarise(
    n_sites = n_distinct(site),
    n_obs = n(),
    n_studies = n_distinct(ref)
  ) %>%
  arrange(lat_band)


# 3. Base map ####
library(maps)
library(sf)
library(patchwork)

# One row per unique site (Clade C only, haplotype-level data only)
site_data <- data %>%
  filter(clade == "C") %>%
  drop_na(lat, lon, haplotype) %>%  # removes sites w/o coordinates and clade-only studies
  distinct(lat, lon, n)             # one row per unique site

australia_map <- map_data("world", region = "Australia")
reef <- read_sf("/Users/teresa/Desktop/WAsymbionts/WA_reef_extent/Reef-Extent/reefextent.shp")

WA <- ggplot() +
  geom_polygon(
    data = australia_map,
    aes(long, lat, group = group),
    colour    = "#000000",
    fill      = "#A19D6A",
    linewidth = 0.2
  ) +
  coord_fixed(
    ratio = 1,
    xlim  = c(110, 130),
    ylim  = c(-35, -10),
    expand = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(
    labels = degE
  ) +
  scale_y_continuous(
    breaks = seq(-35, -10, 5),
    labels = degS
  ) +
  mytheme

WA

# 3.1 Map reef + sites #####
WA_sites <- WA +
  geom_sf(
    data        = reef,
    mapping     = aes(),
    colour      = NA,
    fill        = "#043743",
    inherit.aes = FALSE
  ) +
  geom_point(
    data   = site_data,
    aes(x = lon, y = lat, size = n),
    shape  = 21,
    fill   = "#71874E50",
    colour = "#71874E",
    stroke = 0.5
  ) +
  scale_size_continuous(
    name   = NULL,
    range  = c(2, 6),
    breaks = c(5, 10, 20, 30, 60)
  ) +
  coord_sf(
    xlim   = c(110, 130),
    ylim   = c(-35, -10),
    expand = FALSE
  ) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.margin   = margin(5, 5, 5, 5),
    axis.title.y    = element_blank(),
    axis.title.x    = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.line       = element_blank(),
    panel.border    = element_rect(fill = NA, linewidth = 1)
  )

WA_sites

# 4. Figure 1 ####
combined <- (
  WA_sites |
    model_plot +
    theme(
      axis.text.y  = element_blank(),
      axis.line.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
) +
  plot_layout(widths = c(1, 1))

combined

# Save Figure 1
combined %>%
  ggsave(
    filename = "Fig_1c.pdf",
    path     = "figures",
    device   = cairo_pdf,
    units    = "cm",
    width    = 16,
    height   = 11
  )


# Clean
# rm(list=ls())

