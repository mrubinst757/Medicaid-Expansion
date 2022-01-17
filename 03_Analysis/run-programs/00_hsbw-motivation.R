# program 00_hsbw-motivation.R
# purpose: simple simulation to show solutions of SBW v H-SBW
# date modified: november 25, 2020
# author: max rubinstein

# load libraries and set simulation parameters -------------------------------------------------------
library(tidyverse)
library(optweight)

expit = function(x) exp(x)/(1 + exp(x))

set.seed(100)
n = 800 # number of units
M = 40 # number of blocks
n_j = 20 # number of units per block 
mu_j = rbinom(M, 1, 0.5) # mean X
n_j <- round(rnorm(M, 20, 3))
n <- sum(n_j)

# generate data ---------------------------------------------------------------------------------------

data = tibble(
  id = 1:n,
  M = unlist(map2(1:M, n_j, ~rep(.x, each = .y))), 
) %>%
  group_by(M) %>%
  mutate(
    mu_j = rbinom(1, 1, 0.5)
  ) %>%
  ungroup() %>%
  mutate(X = rnorm(n, mu_j, 1)) %>%
  group_by(M) %>%
  mutate(A = rbinom(1, 1, expit(mean(X))))

cdat = subset(data, A == 0)
tdat = subset(data, A == 1)
M_c = length(unique(cdat$M))
targets = mean(tdat$X)

# estimate models -------------------------------------------------------------------------------
gen_models = function(tols) {
  group_n <- as.numeric(table(cdat$M))
  models = map(c(0, 0.5, 0.99), ~ optweight.svy(~ X, data = cdat, tols = tols, targets = targets, std.cont = FALSE,
                                               sigma2.y = 1 - .x, re = .x, group_n = group_n, target_n = 100))
  models
}

# models for sequence of error tolerance
models_sequence = gen_models(0)  

# plot results ----------------------------------------------------------------------------------
grays <- gray.colors(4)[1:3]

plot.dat <- cdat %>%
  ungroup() %>%
  mutate(`rho = 0 (SBW)` = models_sequence[[1]]$weights,
         `rho = 0.5` = models_sequence[[2]]$weights,
         `rho = 0.99` = models_sequence[[3]]$weights) %>%
  gather(key, value, `rho = 0 (SBW)`, `rho = 0.5`, `rho = 0.99`) %>%
  group_by(key) %>%
  mutate(var = var(value)) %>%
  group_by(key, M) %>%
  summarize(value = sum(value), `Variance of \n weights` = as.numeric(mean(var))) %>%
  ungroup() 

values <- rev(unique(plot.dat$`Variance of \n weights`))
names(values) <- grays

plot.dat %>%
  ggplot(aes(x = factor(M), y = value, fill = factor(`Variance of \n weights`))) +
  geom_bar(stat = 'identity') +
  facet_wrap(~key) +
  theme_minimal() +
  scale_fill_manual(values = rev(grays)) +
  xlab('Group') + ylab("Sum of Weights by Group") +
  theme(axis.text.x = element_blank()) +
  labs(subtitle = "Note : Darker shading indicates higher variance weights") +
  theme(legend.position = "none") + 
  ggsave('../../02_Paper/01_Plots/proofofconcept.png')

