require(tidyverse)
require(magrittr)

# 1. Data ####
# 1.1 Load data ####
data <- read_csv("~/Desktop/WAsymbionts/presence_absence.csv") %T>%
  print()

data %<>%
  mutate(proportion = obs/n) %T>%
  print()

data %>% filter(proportion > n)

# mod <- glm(obs/n ~ lat * clade, weights = n,
#           family = binomial(link = "logit"), 
#          data = data)
# 
# mod
# 
# mod <- glm(cbind(obs, n - obs) ~ lat * clade,
#            family = binomial(link = "logit"), 
#            data = data)
# 
# mod

# data_long <- data %>%
#   rename(ones = obs) %>%
#   mutate(zeros = n - ones) %>%
#   pivot_longer(c(ones, zeros),
#                values_to = "count") %>%
#   mutate(presence = if_else(name == "ones", 1, 0)) %>%
#   uncount(count) %>%
#   select(-name) %T>%
#   print()
  
data_long <- data %>%
  rowwise() %>%
  mutate(presence = list( c( rep(1, obs) , rep(0, n - obs) ) )) %>%
  unnest(presence) %T>%
  print()

# data_long %<>% 
#   mutate(lat_c = lat - mean(lat)) to center

data_long %<>%
  mutate(lat_c = lat - (-22)) %T>%
  print() #Ningaloo

# 2. Model ####
# 2.1 Clade model ####

mod_clade <- glm(presence ~ lat_c * clade, 
                 family = binomial(link = "logit"), 
                 data = data_long %>%
                   mutate(clade = clade %>% fct_relevel("C")))

data_long %>%
  filter(clade == "C" & is.na(haplotype))
  

data_long %>%
  filter(clade == "C") %>%
  drop_na(haplotype)

data_long %>%
  filter(clade == "C" & !is.na(haplotype))

# with lat: 37.381
# with lat_c: -1.77803 inv(-1.77803)= 0.1445466; 14.5% likelihood of finding C in Ningaloo

inv(0.043235)

summary(mod_clade)

summary(mod_clade)

new_clade <- expand_grid(
  lat_c = data_long %$% seq(min(lat_c), max(lat_c), 0.1), # generate new data
  clade = data_long %>% distinct(clade) %>% pull(clade)
) %T>%
  print()

inv <- family(mod_clade)$linkinv # calculate inverse of link function
predicted <- predict(mod_clade, type = "link", se.fit = TRUE, newdata = new_clade) # predict

new_clade$fit <- inv( predicted$fit )
new_clade$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new_clade$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new

inv(-1.77803)*100
inv(-1.77803 - 1.45787 * qnorm(0.975))*100
inv(-1.77803 + 1.45787 * qnorm(0.975))*100

# 2.2 Model global ####

mod_global <- glmer(presence ~ lat + (lat|clade), 
                  family = binomial(link = "logit"), 
                  data = data_long)

mod_global


# 
# require(lme4)
# 
# mmod <- glmer(presence ~ lat * clade + (lat|ref), # intercept only + (1|ref)
#               family = binomial(link = "logit"), 
#               data = data_long)
# mmod
# 
# 
# 
# mod <- glm(presence ~ lat * haplotype, 
#            family = binomial(link = "logit"), 
#            data = data_long %>% drop_na(haplotype))
# 
# mod
# 
# inv <- family(mod)$linkinv
# 
# inv


new <- expand_grid(
  lat = data_long %$% seq(min(lat), max(lat), 0.1), # generate new data
  clade = data_long %>% distinct(clade) %>% pull(clade)
) %T>%
  print()

inv <- family(mod_clade)$linkinv # calculate inverse of link function
predicted <- predict(mod_clade, type = "link", se.fit = TRUE, newdata = new) # predict

new$fit <- inv( predicted$fit )
new$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new


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

# 2.3 Haplotype model ####

mod_hap <- glm(presence ~ lat_c * haplotype,         # Latitude centered at 22ºS
                 family = binomial(link = "logit"), 
                 data = data_long %>% 
                 filter(clade == "C") %>%
                 drop_na(haplotype)) 

mod_hap
# (Intercept)     1.80345

# Create new data with CENTERED latitude
new <- expand_grid(
  lat_c = data_long %>% 
    filter(clade == "C") %>%
    drop_na(haplotype) %$% 
    seq(min(lat_c), max(lat_c), 0.1),
  haplotype = data_long %>% 
    filter(clade == "C") %>%
    drop_na(haplotype) %>% 
    distinct(haplotype) %>% 
    pull(haplotype)
)%T>%
  print()

# Predictions
inv <- family(mod_hap)$linkinv # calculate inverse of link function
predicted <- predict(mod_hap, type = "link", se.fit = TRUE, newdata = new) # predict

new$fit <- inv( predicted$fit )
new$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new %<>%
  mutate(lat = lat_c - 22) %T>%
  print()


#2.3.1 Visualisation ####

new %>%
  ggplot() +
  geom_line(aes(lat, fit, colour = haplotype), alpha = 0.5) +
  # geom_ribbon(aes(lat, ymin = lower, ymax = upper, fill = clade),
  #             alpha = 0.2) +
  geom_point(data = data,
             aes(lat, proportion, colour = haplotype), shape = 16, alpha = 0.2) +
  #facet_grid(~ clade) +
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


# 2.4 Global pred for C ####

mod_global_c <- glmer(presence ~ lat_c + (lat_c|haplotype),
                      family = binomial(link = "logit"),
                      data = data_long %>% 
                      filter(clade == "C"))
mod_global_c
# (Intercept)  lat  
# -2.30386     -0.09514  
# At the reference latitude (-22°S):
# Log-odds of presence = -0.21
# Probability = exp(-0.21) / (1 + exp(-0.21)) = 0.45 (45%)
# On average across all haplotypes, there's a 45% chance of presence at -22°S

# Create new data with CENTERED latitude
new_mod_global_c <- expand_grid(
  lat_c = data_long %>% 
    filter(clade == "C") %$% 
    seq(min(lat_c), max(lat_c), 0.1),
  presence = 0) %T>%
  print()

# Predictions
inv <- family(mod_global_c)$linkinv # calculate inverse of link function
new_mod_global_c$fit <- predict(mod_global_c, type = "link", re.form = NA, newdata = new_mod_global_c) # fit on link scale
mm <- model.matrix(terms(mod_global_c), new_mod_global_c)

fvar <- diag(mm %*% tcrossprod(vcov(mod_global_c), mm)) # fixed uncertainty
tvar <- fvar + VarCorr(mod_global_c)$haplotype[1] # total (fixed + random) uncertainty 2.0060^2

new_mod_global_c$upper <- inv(new_mod_global_c$fit + qnorm(0.975)*sqrt(tvar))
new_mod_global_c$lower <- inv(new_mod_global_c$fit - qnorm(0.975)*sqrt(tvar))
new_mod_global_c$fit <- inv(new_mod_global_c$fit)


new_mod_global_c %<>%
  mutate(lat = lat_c - 22) %T>%
  print()


# 2.4.1 Visualisation

new_mod_global_c %>%
  ggplot() +
  geom_line(data = new, aes(x = lat, y = fit, group = haplotype), 
            alpha = 0.1, colour = "forestgreen") +
  geom_line(aes(lat, fit), linewidth = 1, colour = "forestgreen") +
  # geom_line(data = new_lat, aes(lat, fit), linewidth = 1, linetype = 2) + # uncomment to see WRONG prediction
  # geom_ribbon(aes(lat, ymin = lower, ymax = upper),
  #             alpha = 0.2) +
  geom_point(data = data %>% filter(clade == "C") %>% drop_na(haplotype),
             aes(lat, proportion, size = n), 
             shape = 16, alpha = 0.2, colour = "forestgreen") +
  coord_flip() +
  theme_minimal()
  # theme(legend.position = "none")


summary(mod_global_c)

new_mod_global_c %>%
  filter(lat == -15)


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


# 2.5 Model latitude ####
mod_lat <- glm(presence ~ lat_c,
                  family = binomial(link = "logit"),
                  data = data_long %>% 
                  filter(clade == "C"))
mod_lat
# (Intercept)        lat_c  
# 0.04323      0.14066  

# At the reference latitude (-22°S):
# Log-odds of presence = 0.04323
# Probability = exp(0.04323) / (1 + exp(0.04323)) = 0.51 (51% chance)
# At -22°S, Clade C haplotypes are present in about half the samples

new_lat <- expand_grid(
  lat_c = data_long %>% 
    filter(clade == "C") %$% 
    seq(min(lat_c), max(lat_c), 0.1)
)%T>%
  print()

inv <- family(mod_lat)$linkinv # calculate inverse of link function
predicted <- predict(mod_lat, type = "link", se.fit = TRUE, newdata = new_lat) # predict

new_lat$fit <- inv( predicted$fit )
new_lat$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new_lat$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new_lat %<>%
  mutate(lat = lat_c - 22) %T>%
  print()



# 3. Map ####

library(tidyverse)
library(maps)
library(scales)

mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2), "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_text(size = 12),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 15, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

# Load your data
data <- read_csv("~/Desktop/WAsymbionts/presence_absence.csv")

# Prepare site data - one point per site with max sample size
site_data <- data %>%
  drop_na(lat, lon) %>%  # Remove sites without coordinates
  group_by(lat, lon) %>%
  summarise(
    n = max(n),  # Use maximum sample size at this site
    .groups = 'drop'
  )

site_data <- data %>%
  filter(clade == "C") %>%
  drop_na(lat, lon, haplotype) %>%
  distinct(lat, lon, n)


# WA map
australia_map <- map_data("world", region = "Australia")

degE <- function(x) paste0(abs(x), "°E")
degS <- function(y) paste0(abs(y), "°S")

WA <- ggplot() +
  geom_polygon(
    data = australia_map,
    aes(long, lat, group = group),
    colour    = "#898b8e",
    fill      = "#DCDEE0",
    linewidth = 0.2
  ) +
  coord_fixed(
    ratio = 1,
    xlim  = c(110, 130),
    ylim  = c(-36, -9),
    expand = FALSE
  ) +
  scale_x_continuous(
    labels = degE
  ) +
  scale_y_continuous(
    labels = degS
  ) +
  mytheme

WA

# Add sites to your existing WA map
WA_sites <- WA +
  geom_point(
    data = site_data, 
    aes(x = lon, y = lat, size = n),
    shape = 16,              # Circle with border
    colour = "forestgreen",        # Border color
           # Fill color (your golden color)
    alpha = 0.2,
    stroke = 0.5             # Border thickness
  ) +
  scale_size_continuous(
    name = "Sample size",
    range = c(2, 10),        # Adjust size range as needed
    breaks = pretty_breaks(4) # Automatic nice breaks
  ) +
  theme(
    legend.position = c(0.85, 0.80),  # Top right corner
    legend.background = element_rect(
      fill = "white", 
      colour = "black", 
      linewidth = 0.3
    ),
    legend.margin = margin(5, 5, 5, 5)
  )

WA_sites

# Clean

rm(list=ls())

