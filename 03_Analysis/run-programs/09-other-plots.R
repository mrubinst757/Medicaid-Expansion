# program: other-plots.R
# purpose: create maps/plots pertaining to medicaid expansion & political ideology
# author: max rubinstein
# date modified: december 28, 2020

library(tidyverse)
library(usmap)

data <- read_csv("../../00_Data/Ideology/stateideology_v2018.csv")
expansion <- read_csv("../../00_Data/Medicaid/raw_data.csv", skip = 2)

plot.data <- data %>%
  filter(year == 2013) %>%
  left_join(expansion, by = c("statename" = "Location")) %>%
  select(state = statename, fips = state, citi, inst, expansion = `Jan-14`) %>%
  mutate(treatment = if_else(expansion >=1.38 | state %in% c("Michigan", "New Hampshire"), "Expansion", "Nonexpansion")) %>%
  arrange(inst) %>%
  mutate(state = factor(state, levels = rev(.$state))) 

yint1 <- mean(plot.data$inst[plot.data$treatment == "Expansion"])
yint0 <- mean(plot.data$inst[plot.data$treatment == "Nonexpansion"])

plot.data %>%
  ggplot(aes(x = state, y = inst, color = factor(treatment))) +
  geom_point() +
  ggthemes::theme_hc() +
  theme(axis.text.x = element_text(angle = 90),
        axis.ticks = element_blank()) +
  scale_color_manual(values = c("firebrick", "gray")) +
  geom_hline(yintercept = yint1, color = "firebrick", linetype = "longdash") +
  geom_hline(yintercept = yint0, color = "gray", linetype = "longdash") +
  ylab("Institutional Ideology Score") +
  xlab("") +
  guides(color = guide_legend(title = "")) +
  labs(caption = "Source: Richard C. Fording (2018) and Kaiser Family Foundation") +
  ylim(c(0, 80)) +
  ggsave("../../02_Paper/01_Plots/political-expansion-plot.png")

plot.data %>% 
  select(state, treatment) %>%
  bind_rows(tibble(state = "District of Columbia", treatment = "Expansion")) %>%
  plot_usmap(data = ., values = "treatment") + 
  theme(legend.position = "right", title = element_text(size = 14)) +
  ggtitle("Medicaid Expansion Decisions (2014)") +
  labs(caption = "Source: Kaiser Family Foundation") +
  scale_fill_manual(values = c("firebrick", "lightgray"), name = "") +
  ggsave("../../02_Paper/01_Plots/expansion-map.png")

plot.data %>% 
  bind_rows(tibble(state = "District of Columbia", treatment = "Expansion")) %>%
  mutate(treatment = if_else(treatment == "Expansion", "Treatment", "Control")) %>%
  mutate(treatment = if_else(state %in% c("District of Columbia", "Vermont", "New York", "Massachusetts", "Delaware"), 
                             "Excluded", treatment)) %>%
  select(state, treatment) %>%
  plot_usmap(data = ., values = "treatment") + 
  theme(legend.position = "right", title = element_text(size = 14)) +
  scale_fill_manual(values = c("navyblue", "darkgray", "firebrick"), name = "") +
  ggsave("../../02_Paper/01_Plots/political-expansion-plot-preferred.png")

kstner <- read_csv("../../00_Data/Medicaid/kaestnerclassification.csv") %>%
  select(state = statename, expansion = category)

kplot_data <- kstner %>%
  left_join(plot.data %>% select(-expansion) %>% mutate(state = as.character(state)), by = c("state")) %>%
  mutate(treatment = if_else(state == "District of Columbia", "Expansion", treatment)) %>%
  mutate(value = 1) %>%
  arrange(expansion) 

kplot_data %>%
  mutate(expansion = factor(expansion, levels = c("No prior expansion", "Prior limited", "Prior expansion (Kaestner)",
                                                  "Prior expansion (Frean)", "Prior full"))) %>%
  ggplot(aes(x = state, y = value, fill = expansion, group = expansion)) +
  geom_bar(stat = "identity") +
  ggthemes::theme_hc() +
  theme(axis.text.x = element_text(angle = 90),
        axis.ticks = element_blank()) +
  ylab("") +
  xlab("") +
  theme(axis.text.y = element_blank()) +
  guides(fill = guide_legend(title = "")) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~treatment, scales = "free") +
  labs(caption = "Source: Kaestner et al. (2016), Frean et al. (2017)") +
  ggsave("../../02_Paper/01_Plots/expansion-heterogeneity.png")

kplot_data %>%
  filter(treatment == "Expansion") %>%
  mutate_at("expansion", ~gsub("Prior expansion \\(Kaestner\\)", "No prior expansion", .)) %>%
  mutate_at("expansion", ~gsub("\\(Frean\\)", "", .)) %>%
  mutate_at("expansion", ~trimws(.)) %>%
  mutate(expansion = factor(expansion, levels = c("No prior expansion",
                                                  "Prior expansion", "Prior full"))) %>%
  arrange(expansion) %>%
  mutate(state = factor(state, levels = .$state))  %>%
  ggplot(aes(x = state, y = value, fill = factor(expansion))) +
  geom_bar(stat = "identity") +
  labs(caption = "Source: Frean et al (2017); Kaestner et al. (2016)") +
  theme_minimal() +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_blank()) +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(title = "")) +
  labs(caption = "Source: Kaestner et al. (2016), Frean et al. (2017)") +
  ggsave("../../02_Paper/01_Plots/expansion-heterogeneity-txonly.png")

kplot_data %>%
  filter(treatment == "Nonexpansion") %>%
  mutate_at("expansion", ~gsub("Prior expansion \\(Kaestner\\)", "No prior expansion", .)) %>%
  mutate_at("expansion", ~gsub("\\(Frean\\)", "", .)) %>%
  mutate(expansion = factor(expansion, levels = c("No prior expansion",
                                                  "Prior limited"))) %>%
  arrange(expansion) %>%
  mutate(state = factor(state, levels = .$state))  %>%
  ggplot(aes(x = state, y = value, fill = factor(expansion))) +
  geom_bar(stat = "identity") +
  labs(caption = "Source: Frean et al (2017); Kaestner et al. (2016)") +
  theme_minimal() +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_blank()) +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = brewer.pal(5, "Dark2")[c(1, 5)]) +
  labs(caption = "Source: Kaestner et al. (2016), Frean et al. (2017)") +
  ggsave("../../02_Paper/01_Plots/expansion-heterogeneity-ctonly.png")

kstner %>%
  left_join(plot.data %>% select(-expansion) %>% mutate(state = as.character(state)), by = c("state")) %>%
  mutate(treatment = if_else(state == "District of Columbia", "Expansion", treatment)) %>%
  mutate(expansion = factor(expansion, levels = c("No prior expansion", "Prior limited", "Prior expansion (Kaestner)",
                                                  "Prior expansion (Frean)", "Prior full"))) %>%
  filter(state != "New Hampshire", expansion != "Prior full") %>%
  arrange(expansion) %>%
  mutate(state = factor(state, levels = .$state)) %>%
  ggplot(aes(x = state, y = value, fill = expansion, group = expansion)) +
  geom_bar(stat = "identity") +
  ggthemes::theme_hc() +
  theme(axis.text.x = element_text(angle = 90),
        axis.ticks = element_blank()) +
  ylab("") +
  xlab("") +
  theme(axis.text.y = element_blank()) +
  guides(fill = guide_legend(title = "")) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~treatment, scales = "free") +
  labs(caption = "Source: Kaestner et al. (2016), Frean et al. (2017)") +
  ggsave("../../02_Paper/01_Plots/expansion-heterogeneity-classification.png")

# etu v ett plot
set.seed(34)
tibble(
  y0.p  = rep(seq(1, 10, 1), 2),
  Y0    = y0.p + rnorm(20, 0, 0.15),
  x     = rep(seq(1, 10, 1), 2),
  grp   = c(rep("Treatment", 10), rep("Control", 10)),
  a     = if_else(x >= 8, 1, 0),
  Y1    = if_else(a == 1 & grp == "Treatment", Y0 + 2, 
                  if_else(a == 1 & grp == "Control", Y0 - 2, Y0)),
  Y     = if_else(a == 1 & grp == "Treatment", Y1, Y0)
) %>%
  gather(key, value, Y0, Y1, Y) %>%
  unite(variable, grp, key, remove = FALSE) %>%
  filter(!variable %in% c("Treatment Y", "Control Y")) %>%
  mutate(panel = case_when(
    key == "Y0" & grp == "Control" ~ "Factual outcomes",
    key == "Y0" & grp == "Treatment" ~ "Unobserved counterfactuals",
    key == "Y1" & grp == "Treatment" ~ "Factual outcomes",
    key == "Y1" & grp == "Control" ~ "Unobserved counterfactuals"
  )) %>%
  mutate_at("variable", funs(gsub("_", " ", .))) %>%
  mutate(panel = factor(panel, levels = c("Factual outcomes", "Unobserved counterfactuals"))) %>%
  mutate(variable = factor(variable, levels = c("Control Y0", "Control Y1", "Treatment Y0", "Treatment Y1"))) %>%
  filter(!is.na(panel)) %>%
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~panel) +
  theme_minimal() +
  theme(axis.ticks = element_blank()) +
  xlab("Time") +
  ylab("Outcome") +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = 7, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_color_brewer(palette = "BrBG") +
  ggsave("../../02_Paper/01_Plots/ett-v-etu-motivation.png")

