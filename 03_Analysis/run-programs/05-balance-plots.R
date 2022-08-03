# program: 05-balance-plots.R
# purpose: generate balance plots
# author: max rubinstein
# date modified: december 14, 2020

source("03_Analysis/04_balance-plots.R")

# read data and variable names -----------------------------------------------------------------------
variable_names = read_csv('../02_Specs/tol_specs.csv') %>%
  filter(`Reference Variable` == 0) %>%
  arrange(Variable) %>%
  select(Variable, `Raw Variable`, Denominator)

variables <- variable_names$Variable

varnames <- read_csv("../02_Specs/codebook.csv") %>%
  filter(Pct == 1)
var_names <- varnames$`Variable Name`
names(var_names) <- varnames$Variable

codebook <- read_csv("../02_Specs/codebook.csv")

merged_data_c1 <- readRDS("../01_ProcessedData/calibrated-data-all.rds") %>%
  filter(set == "true") %>%
  unnest() %>%
  nest(-key) 

merged_data_c2 <- readRDS("../01_ProcessedData/calibrated-data-c2.rds") %>%
  filter(set == "true") %>%
  unnest() %>%
  nest(-key) 

# extract preferred weights ---------------------------------------------------------------------------
etc_c1_results <- readRDS("../04_Output/c1-results.rds")
etc_c2_results <- readRDS("../04_Output/c2-results.rds")
weights_etu_c1 <- etc_c1_results$Preferred$None$`H-SBW`$sigma_uu_avg$weights
weights_etu_c2 <- etc_c2_results$`Early expansion`$None$`H-SBW`$sigma_uu_avg$weights
weights_etu_c1_unadj <- etc_c1_results$Preferred$None$`H-SBW`$sigma_zero$weights

balance_comp1 <- hsbw_baltab(merged_data_c1$data[[3]], weights_etu_c1_unadj, variables) %>%
  filter(grepl("unins|unemp", Variables)) %>%
  select(Variables, `Unweighted Diff`, `Weighted Diff (none)` = `Weighted Diff`)

balance_comp2 <- hsbw_baltab(merged_data_c1$data[[2]], weights_etu_c1_unadj, variables) %>%
  filter(grepl("unins|unemp", Variables)) %>%
  select(Variables, `Weighted Diff (homogeneous)` = `Weighted Diff`)

balance_comp3 <- hsbw_baltab(merged_data_c1$data[[1]], weights_etu_c1_unadj, variables) %>%
  filter(grepl("unins|unemp", Variables)) %>%
  select(Variables, `Weighted Diff (heterogeneous)` = `Weighted Diff`)

final_bcomp <- balance_comp1 %>%
  left_join(balance_comp2, by = "Variables") %>%
  left_join(balance_comp3, by = "Variables") %>%
  mutate_at("Variables", ~stringr::str_replace_all(., var_names))

balance_etc_c1 <- hsbw_baltab(merged_data_c1$data[[2]], weights_etu_c1, variables)
balance_etc_c2 <- hsbw_baltab(merged_data_c2$data[[2]], weights_etu_c2, variables)

total_etc <- balance_etc_c1 %>%
  mutate(`Preferred (Pct Diff)` = paste0("(", round(`Unweighted Diff`, 2), ", ", round(`Weighted Diff`, 2), ")")) %>%
  mutate(`Preferred (SMD Diff)` = paste0("(", round(`Unweighted SMD`, 2), ", ", round(`Weighted SMD`, 2), ")")) %>%
  select(-contains("weighted")) %>%
  left_join(balance_etc_c2 %>%
              mutate(`Early excluded (Pct Diff)` = paste0("(", round(`Unweighted Diff`, 2), ", ", round(`Weighted Diff`, 2), ")")) %>%
              mutate(`Early excluded (SMD Diff)` = paste0("(", round(`Unweighted SMD`, 2), ", ", round(`Weighted SMD`, 2), ")")) %>%
              select(-contains("weighted"))) %>%
  mutate_at("Variables", ~stringr::str_replace_all(., var_names))

balance_etc_c1 %>%
  mutate_if(is.numeric, ~(ifelse(abs(.) > 0.1, 1, 0))) %>%
  summarize_if(is.numeric, sum)

sd(merged_data_c1$data[[3]]$female_pct)

# output plots and tables ---------------------------------------------------
print(xtable::xtable(final_bcomp), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(total_etc), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

hsbw_balplot_smd <- function(results, file_extension) {
  plot <- results %>%
    nest(`Standardized mean difference` = c("Unweighted SMD", "Weighted SMD"), 
         `Percentage point difference` = c("Unweighted Diff", "Weighted Diff"))  %>%
    gather(key, value, -Variables) %>%
    mutate_at("value", ~map(., ~set_names(.x, c("Unweighted", "Weighted")))) %>%
    unnest(cols = c(value)) %>%
    mutate_at("Variables", ~stringr::str_replace_all(., var_names)) %>%
    gather(weighted, difference, Unweighted, Weighted) %>%
    mutate_at("weighted", ~factor(., levels = c("Weighted", "Unweighted"))) %>%
    group_by(Variables) %>%
    mutate(pctdiff = difference[key == "Standardized mean difference" & weighted == "Unweighted"]) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(Variables, abs(pctdiff)), y = difference, fill = weighted)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~key, scales = "free_x") +
    ggthemes::theme_pander() +
    scale_fill_grey(name = "Difference", breaks = c("Unweighted", "Weighted")) +
    coord_flip() +
    theme(panel.grid.major = element_line(colour = "black")) +
    ylab("Mean Difference (treatment mean minus control mean)") +
    xlab("Variables \n (ordered by absolute magnitude \n of unweighted percentage point imbalance)")
  return(plot)  
}

p.c1 <- hsbw_balplot_smd(balance_etc_c1, "c1")
p.c2 <- hsbw_balplot_smd(balance_etc_c2, "c2")

ggsave(paste0("../../02_Paper/01_Plots/balance-plot-all-etu-c1.png"), p.c1, width = 12, height = 7, type = "cairo")
ggsave(paste0("../../02_Paper/01_Plots/balance-plot-all-etu-c2.png"), p.c2, width = 12, height = 7, type = "cairo")
