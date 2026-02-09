require(tidyverse)
require(magrittr)
data <- read_csv("~/Desktop/presence_absence.csv") %T>%
  print()

data %>%
  mutate(sample2 = interaction(site, sample),
         sample3 = str_c(site, sample, sep = "_"))

data %>%
  distinct(site)

mod <- glm(presence ~ site, family = binomial(link = "logit"), data = data)

mod

# logit is log odds, meaning log(p/(1-p))
# the inverse logits is just the standard logistic 1/(1+exp(-x))
1/(1+exp(-(-1.63900+0.01977)))

predict(mod)
