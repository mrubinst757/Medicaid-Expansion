source("03_Analysis/06_simulations.R")

################################################################################
########################## process and plot results ############################
################################################################################

# set parameter space and run simulations --------------------------------------
simulation_parameters <- expand.grid(
  pop_size = 5000,
  hetero = c(TRUE, FALSE),
  SS_cor = c(0, 0.25, 0.5),
  X_cor = 0.25,
  SNR.xw = c(0.95, 0.9, 0.85),
  rho = 0.25,
  sample_states = 25,
  num_sims = 1000
)

all_sims <- pmap(simulation_parameters, RunSims)
saveRDS(all_sims, "../04_Output/simulations-05-06-22.rds")

# create plots summarizing results ---------------------------------------------
all_sims <- readRDS("../04_Output/simulations-05-06-22.rds")

sim_performance <- all_sims %>%
  map(~invoke(rbind, .)) %>%
  map(~group_by(.x, Xset, re, jackest) %>%
        dplyr::summarize(bias = mean(bias), 
                         covered = mean(covered),
                         N = mean(N),
                         ci.width = mean(abs(uci - lci)),
                         se.est = mean(se),
                         var.obs = var(est),
                         mse = var.obs + bias^2))

sim_performance <- map2(1:nrow(simulation_parameters), sim_performance, 
                        ~cbind(simulation_parameters[.x,], .y)) %>%
  map2(1:nrow(simulation_parameters), ~mutate(.x, id = .y)) %>%
  invoke(rbind, .) %>%
  rename(`Balance variables` = Xset) %>%
  mutate_at("re", factor) %>%
  mutate_at("SS_cor", ~paste0("Rho_x = ", .)) %>%
  mutate_at("SNR.xw", ~paste0("Tau = ", .)) %>%
  rename(`H-SBW Rho` = re) 

bias.plot <- sim_performance %>%
  filter(`Balance variables` %in% c("X", "Xhat.hom", "Xhat.het", "Xhat.cor", "W"), hetero == TRUE) %>%
  ggplot(aes(y = bias, x = `Balance variables`, fill = `H-SBW Rho`)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(SS_cor~SNR.xw) +
  theme_bw() +
  xlab("Covariate set") +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  ylab("Bias \n (equivalently: total imbalance across all covariates)") +
  scale_fill_brewer(palette = "Dark2")

ggsave("../../02_Paper/01_Plots/bias-plot.png", bias.plot)

var.plot <- sim_performance %>%
  filter(`Balance variables` %in% c("X", "Xhat.hom", "Xhat.het", "Xhat.cor", "W"), hetero == TRUE) %>%
  ggplot(aes(y = var.obs, x = `Balance variables`, fill = `H-SBW Rho`)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(SS_cor~SNR.xw) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  xlab("Covariate set") +
  ylab("Variance") +
  scale_fill_brewer(palette = "Dark2")

ggsave("../../02_Paper/01_Plots/var-plot.png", var.plot)

mse.plot <- sim_performance %>%
  filter(`Balance variables` %in% c("X", "Xhat.hom", "Xhat.het", "Xhat.cor", "W"), hetero == TRUE) %>%
  ggplot(aes(y = mse, x = `Balance variables`, fill = `H-SBW Rho`)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(SS_cor~SNR.xw) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  xlab("Covariate set") +
  ylab("Mean Square Error") +
  scale_fill_brewer(palette = "Dark2")

ggsave("../../02_Paper/01_Plots/mse-plot.png", mse.plot)

coverage.plot.1 <- sim_performance %>%
  filter(`Balance variables` %in% c("X", "Xhat.hom", "Xhat.het", "Xhat.cor", "W"), hetero == TRUE) %>%
  ggplot(aes(y = covered, x = `Balance variables`, fill = `H-SBW Rho`)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim=c(0.5, 1)) +
  facet_grid(SS_cor~SNR.xw) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  xlab("Covariate set") +
  ylab("Coverage rate \n (across 1000 simulations for each specification)") +
  geom_hline(yintercept = 0.95, color = "firebrick", linetype = "dashed") +
  scale_fill_brewer(palette = "Dark2")

ggsave("../../02_Paper/01_Plots/coverage-plot-1.png", coverage.plot.1)

coverage.plot.2 <- sim_performance %>%
  filter(`Balance variables` %in% c("Xhat.hom", "Xhat.het", "Xhat.cor"), hetero == TRUE) %>%
  filter(rho == 0.25) %>%
  mutate(`H-SBW Rho \nVariance Est` = paste0("H-SBW Rho: ", `H-SBW Rho`, ",\nXhat Jackknife: ", jackest)) %>%
  ggplot(aes(y = covered, x = `Balance variables`, fill = `H-SBW Rho \nVariance Est`)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim=c(0.5, 1)) +
  facet_grid(SS_cor~SNR.xw) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  xlab("Covariate set") +
  ylab("Coverage rate \n (across 1000 simulations for each specification)") +
  geom_hline(yintercept = 0.95, color = "firebrick", linetype = "dashed") +
  scale_fill_brewer(palette = "Paired")

ggsave("../../02_Paper/01_Plots/coverage-plot-2.png", coverage.plot.2)

ci.length.plot <- sim_performance %>%
  filter(`Balance variables` %in% c("X", "Xhat.hom", "Xhat.het", "Xhat.cor", "W"), hetero == TRUE) %>%
  ggplot(aes(y = ci.width, x = `Balance variables`, fill = `H-SBW Rho`)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(SS_cor~SNR.xw) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  xlab("Covariate set") +
  ylab("Average CI Length \n (across 1000 simulations for each specification)") +
  scale_fill_brewer(palette = "Dark2")

ggsave("../../02_Paper/01_Plots/ci-length-plot.png", ci.length.plot)

################################################################################
#####################additional simulations: x observed ########################
################################################################################
simulation_parameters_x <- expand.grid(
  pop_size = 5000,
  sample_states = 25,
  hetero = FALSE,
  SS_cor = 0.25,
  X_cor = 0.25,
  SNR.xw = 1,
  rho = c(0, 0.25, 0.5, 0.75, 0.99),
  num_sims = 1000
)

all_sims_X <- pmap(simulation_parameters_x, RunSimsX)

saveRDS(all_sims_X, "../04_Output/simulations-x-12-27-21.rds")

all_sims_X <- readRDS("../04_Output/simulations-x-12-27-21.rds") %>%
  map(~invoke(rbind, .)) %>%
  map2(1:nrow(simulation_parameters_x), ~cbind(.x, simulation_parameters_x[.y,])) %>%
  invoke(rbind, .) %>%
  group_by(re, SS_cor, rho) %>%
  summarize(var = var(est), ciw = mean(abs(uci - lci)), covered = mean(covered),
            var.est = mean((25/24)*se^2)) %>%
  ungroup() %>%
  mutate(truth = ifelse(rho == re, 1, 0)) 

coverage.x.plot <- all_sims_X %>%
  mutate_at("rho", ~paste0("Rho* = ", .)) %>%
  ggplot(aes(x = factor(re), y = 100*covered, fill = factor(truth))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap( ~ rho) +
  coord_cartesian(ylim=c(90, 100)) +
  theme_bw() +
  scale_fill_manual(values = c("darkgray", "firebrick")) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "firebrick") +
  ylab("Coverage rate \n (averaged over 1000 simulations)") +
  xlab("H-SBW Rho") +
  labs(fill = "Optimal Rho")

ggsave("../../02_Paper/01_Plots/coverage-x-plot.png", coverage.x.plot)

variance.x.plot <- all_sims_X %>%
  mutate_at("rho", ~paste0("Rho* = ", .)) %>%
  ggplot(aes(x = re, y = var, fill = factor(truth))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap( ~ rho) +
  theme_bw() +
  scale_fill_manual(values = c("darkgray", "firebrick")) +
  ylab("Variance of estimator \n (averaged across 1000 simulations)") +
  xlab("H-SBW rho") +
  labs(fill = "Optimal Rho")

ggsave("../../02_Paper/01_Plots/variance-x-plot.png", variance.x.plot)

variance.x.bias.plot <- all_sims_X %>%
  mutate(bias.var = var.est - var) %>%
  ggplot(aes(x = re, y = bias.var, fill = factor(truth))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap( ~ rho) +
  theme_bw() +
  scale_fill_manual(values = c("darkgray", "firebrick")) +
  ylab("Jaccknife bias \n (averaged across 1000 simulations)") +
  xlab("H-SBW rho") +
  labs(fill = "Optimal Rho")

all_sims_X %>%
  mutate(bias.var = var.est - var) %>%
  arrange(bias.var) %>%
  mutate(bias.abs = abs(bias.var)) %>%
  arrange(bias.abs) 

ggsave("../../02_Paper/01_Plots/variance-x-bias-plot.png", variance.x.bias.plot)

################################################################################
###############additional simulations: xhat cor consistency ####################
################################################################################
cor.dat <- GenerateHierPopulation(number_states = 10000, SS_cor = 0.5, SNR.xw = 0.85,
                                  hetero = TRUE, rho = 0.25, X_cor = 0.25)

cor.sims <- map(c(25, 50, 100, 200), ~RunCorSims(cor.dat, .x, nsims = 500))

cor.sims <- readRDS("../04_Output/xhat-cor.rds")

map(cor.sims, ~invoke(rbind, .)) %>%
  map2(c(25, 50, 100, 200), ~mutate(.x, M = .y)) %>%
  invoke(rbind, .) %>%
  filter(Xset == "Xhat.cor") %>%
  group_by(M, re, Xset) %>%
  summarize(bias = mean(bias), mse = mean(mse), var = var(ests)) %>%
  arrange(re) 

