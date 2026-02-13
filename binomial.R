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


# 2. Model ####
# 2.1 Clade model ####

mod_clade <- glm(presence ~ lat * clade, 
                 family = binomial(link = "logit"), 
                 data = data_long)

mod_clade



# 2. Model ####
# 2. Model global ####



mod_global <- glmer(presence ~ lat + (lat|clade), 
                  family = binomial(link = "logit"), 
                  data = data_long)

mod_global

# mod <- glm(obs/n ~ lat * clade, weights = n,
#            family = binomial(link = "logit"), 
#            data = data)
# 
# mod
# 
# mod <- glm(cbind(obs, n - obs) ~ lat * clade,
#            family = binomial(link = "logit"), 
#            data = data)
# 
# mod
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


mod_hap <- glm(presence ~ lat * haplotype, 
                 family = binomial(link = "logit"), 
                 data = data_long %>% 
                 filter(clade == "C")) 

mod_hap

mod_global <- glmer(presence ~ lat + (lat|haplotype),
              family = binomial(link = "logit"),
              data = data_long %>% 
                filter(clade == "C"))
mod_global

mod_global <- glm(presence ~ lat,
                    family = binomial(link = "logit"),
                    data = data_long %>% 
                      filter(clade == "C"))
mod_global

new <- expand_grid(
  lat = data_long %>% 
    filter(clade == "C") %$% 
    seq(min(lat), max(lat), 0.1), # generate new data
  haplotype = data_long %>% 
    filter(clade == "C") %>% 
    distinct(haplotype) %>% 
    pull(haplotype)
) %T>%
  print()

inv <- family(mod_hap)$linkinv # calculate inverse of link function
predicted <- predict(mod_hap, type = "link", se.fit = TRUE, newdata = new) # predict

new$fit <- inv( predicted$fit )
new$upper <- inv( predicted$fit + predicted$se.fit * qnorm(0.975) )
new$lower <- inv( predicted$fit - predicted$se.fit * qnorm(0.975) )

new


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


str_extract()

data %>%
  mutate(sample2 = interaction(site, sample),
         sample3 = str_c(site, sample, sep = "_"))

data %>%
  distinct(site)



mod

# logit is log odds, meaning log(p/(1-p))
# the inverse logits is just the standard logistic 1/(1+exp(-x))
1/(1+exp(-(-1.63900+0.01977)))

predict(mod)


