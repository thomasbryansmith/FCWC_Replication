#==============================================================================#

# load packages (simulations) ====
library(dplyr)
library(brms)
library(tidyr)
library(ggplot2)
library(ggridges)
library(posterior)
library(bayestestR)
library(bayesplot)
library(purrr)
library(truncdist)
library(flextable)
library(broom.mixed)

#==============================================================================#
# Read Data                                                                ====#
#==============================================================================#

exp_dat <- readxl::read_xlsx("./data/Data_2026_02_25.xlsx", sheet = 2)

exp_dat <- exp_dat %>% pivot_longer(
  cols = matches("^(Acc|IG|Con)"),
  names_to = c("Tactic", "Outcome"),
  names_pattern = "^(Acc|IG|Con)(.*)",
  values_to = "Count")

exp_dat <- exp_dat[which(!is.na(exp_dat$Count)),]

exp_dat <- exp_dat %>% pivot_wider(
  names_from = Outcome,
  values_from = Count
)

exp_dat <- exp_dat[which(!is.na(exp_dat$FC)),]

exp_dat <- exp_dat %>%
  mutate(N = coalesce(FC, 0) + coalesce(FD, 0) + coalesce(TC, 0) + coalesce(TD, 0),
         Y = FC,
         Y_fd = FD,
         Location = factor(Location),
         Paradigm = factor(Paradigm),
         Tactic = factor(Tactic,
                         levels = c("Con",
                                    "Acc",
                                    "IG")),
         Ref_complete = Ref,
         Ref = as.numeric(factor(Ref)))

#==============================================================================#
# Fit models                                                               ====#
#==============================================================================#

# False Confession (FC|Tactic) Model

## Check feasibility of interaction between tactic and paradigm
with(exp_dat, tapply(Y, list(Tactic, Paradigm), sum))

  ### Not feasible because of the sparseness of the observations
  ### appearing under information gathering tactics. We need:
  ### 1. Comparison of information gathering to "Other" paradigms.
  ### 2. More alt-key studies of information gathering.
  ### 3. More information gathering studies, in general.
  ### This is, in part, because information gathering is
  ### so unlikely to elicit a false confession.

## Fit the most advanced model possible given data limitations
### Cannot nest within both location and study because we do not
### have sufficient variability in location (we have 4 countries,
### we would need ~10+ with sufficient within and between location
### variance on the observed variables).

### Given that we already nest by study, and there is no within-study
### variation by location, we can safely ignore location and allow
### the random effect for study capture the between-location variability.

### However, this also applies to Paradigm, where we only have 3 levels.
### Problematically, each study is only ever 1 paradigm. This means there
### is no within-study variability on Paradigm for the grouping variable.
### Less problematically, grouping with study will also capture
### between-paradigm variance (because each study can only be one paradigm,
### in the same way that each study is only based in one country).

### We could maybe accommodate these variables by specifying heavily 
### informative priors, but we are mostly interested in differences
### between tactics while controlling the study (and, by extension,
### paradigm and location).

### The most parsimonious model simply nests tactic within study:

#### Find empirical log-odds for each tactic
exp_dat %>% 
  group_by(Tactic) %>% 
  summarize(emp_logodds = log((sum(Y)+0.5)/(sum(N-Y)+0.5)))
##### We observe between -0.874 and -2.6, so we set the sd
##### of the weakly informative prior to 1.5. This roughly covers
##### +/- 3 SDs, including plausible extremes but avoiding implausible
##### unobserved values. Given that we are setting "control" as our 
##### reference tactic, and the average empirical log odds
##### of this no-tactic control condition is -1.73, we
##### set the prior for the intercept to match this observed log-odds.

#### Find empirical log-odds for each tactic across studies
exp_dat %>%
  group_by(Ref) %>%
  summarize(emp_logodds = log((sum(Y)+0.5)/(sum(N-Y)+0.5))) %>%
  summarise(sd(emp_logodds))
##### We observe a standard deviation of 1.87 between studies,
##### which is greater than 1. Exponential(0.5) places more prior
##### mass on larger SDs. Decreases regularization to reflect heterogeneity 
##### between studies, allowing the modelling of more extreme study effects.

### Estimate prelim. model with weakly informative priors
m1 <- brm(
  formula = Y | trials(N) ~ Tactic + (1 | Ref),
  data = exp_dat,
  family = binomial(),
  prior = c(
    prior(normal(-1.73, 0.8), class = "Intercept"), # intercept
    prior(normal(0, 1.5), class = "b"),               # slope
    prior(exponential(0.5), class = "sd")             # re
  ),
  cores = 4,
  chains = 4,
  iter = 4000,
  #warmup = 1000,
  seed = 38655,
  #control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE)
)

## Posterior predictive checks
pp_check(m1, type = "hist", ndraws = 11)
pp_check(m1, type = "bars_grouped", group = "Tactic")

## Check shrinkage magnitude of studies with 1 ob
## Should be ~0 for studies #4, 7, & 11
ranef(m1)$Ref

## Check influence of observations
loo_res <- loo(m1)
pareto_k_ids(loo_res, threshold = 0.7)

## K-fold cross-validation
kfold(m1, k = 10,
      chains = 4,
      cores = 4)
### Low SE on the elpd kfold estimate (along with the
### fact k-fold cross-validation converged) implies that
### while there are a large number of highly influential
### observations in the data, incrementally removing 10%
### of them and re-estimating does not dramatically influence
### the parameter estimates.

### Re-estimate the final model with M&A iterations
m1 <- brm(
  formula = Y | trials(N) ~ Tactic + (1 | Ref),
  data = exp_dat,
  family = binomial(),
  prior = c(
    prior(normal(-1.73, 0.8), class = "Intercept"), # intercept
    prior(normal(0, 1.5), class = "b"),               # slope
    prior(exponential(0.5), class = "sd")             # re
  ),
  cores = 4,
  chains = 4,
  iter = 251000,
  warmup = 1000,
  seed = 38655,
  #control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE)
)

## Visualize results (posterior contrasts)
newdata <- expand.grid(
  Tactic = c("Con","Acc","IG"),
  Ref = NA,   # for population-level estimates
  N = 1
)

post <- posterior_epred(m1, 
                        newdata = newdata,
                        allow_new_levels = TRUE)

post_acc <- post[, newdata$Tactic == "Acc"]
post_info <- post[, newdata$Tactic == "IG"]
post_ctrl <- post[, newdata$Tactic == "Con"]

diff_acc_ctrl <- post_acc - post_ctrl
diff_info_ctrl <- post_info - post_ctrl
diff_acc_info <- post_acc - post_info

posterior_contrasts <- data.frame(
  contrast = c("Acc vs Ctrl", "Info vs Ctrl", "Acc vs Info"),
  mean = c(mean(diff_acc_ctrl), mean(diff_info_ctrl), mean(diff_acc_info)),
  lower = c(quantile(diff_acc_ctrl, 0.025),
            quantile(diff_info_ctrl, 0.025),
            quantile(diff_acc_info, 0.025)),
  upper = c(quantile(diff_acc_ctrl, 0.975),
            quantile(diff_info_ctrl, 0.975),
            quantile(diff_acc_info, 0.975))
)

posterior_contrasts

df_post <- tibble(
  `Acc. vs Ctrl` = as.vector(diff_acc_ctrl),
  `Info. vs Ctrl` = as.vector(diff_info_ctrl),
  `Acc. vs Info.` = as.vector(diff_acc_info)
) %>% pivot_longer(cols = everything(), names_to = "contrast", values_to = "prob_diff")

df_post$contrast <- factor(df_post$contrast, 
                           levels = c("Info. vs Ctrl",
                                      "Acc. vs Ctrl",
                                      "Acc. vs Info."))

ggplot(df_post, aes(x = prob_diff, fill = contrast)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Posterior Difference in Probability of FC", 
       y = "Density",
       fill = "Contrast") +
  scale_fill_manual(values = c("grey5", "grey80", "red")) + 
  theme_classic() + 
  theme(text = element_text(size = 27, family = "serif"),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
        plot.title = element_text(size = 26,
                                  margin = margin(b = 5, unit = "pt")),
        plot.subtitle = element_text(size = 23, 
                                     margin = margin(b = 20, unit = "pt")),
        plot.caption = element_text(size = 17, hjust = 0, 
                                    margin = margin(t = 20, unit = "pt")),
        strip.text = element_text(size = 23,
                                  margin = margin(t = 10, b = 10, 
                                                  unit = "pt")))

#-----------------------------------------------------------------------------#

# False Denial (FD|Tactic) Model

### Subset data so that it includes complete FD cases
exp_dat_fd <- exp_dat %>%
  filter(!is.na(FD))

#### Find empirical log-odds for each tactic
exp_dat_fd %>% 
  group_by(Tactic) %>% 
  summarize(emp_logodds = log((sum(Y)+0.5)/(sum(N-Y)+0.5)))
##### We observe between -1.52 and -3.09, so we set the sd
##### of the weakly informative prior to 1.5. This roughly covers
##### +/- 3 SDs, including plausible extremes but avoiding implausible
##### unobserved values. Given that we are setting "control" as our 
##### reference tactic, and the average empirical log odds
##### of this no-tactic control condition is -2.58, we
##### set the prior for the intercept to match this observed log-odds.

#### Find empirical log-odds for each tactic across studies
exp_dat_fd %>%
  group_by(Ref) %>%
  summarize(emp_logodds = log((sum(Y)+0.5)/(sum(N-Y)+0.5))) %>%
  summarise(sd(emp_logodds))
##### We observe a standard deviation of 1.34 between studies,
##### which is greater than 1. Exponential(0.5) places more prior
##### mass on larger SDs. Decreases regularization to reflect heterogeneity 
##### between studies, allowing the modelling of more extreme study effects.

### Estimate the model with M&A iterations
m2 <- brm(
  formula = Y_fd | trials(N) ~ Tactic + (1 | Ref),
  data = exp_dat_fd,
  family = binomial(),
  prior = c(
    prior(normal(-2.58, 0.8), class = "Intercept"),   # intercept
    prior(normal(0, 1.5), class = "b"),               # slope
    prior(exponential(0.5), class = "sd")             # re
  ),
  cores = 4,
  chains = 4,
  iter = 251000,
  warmup = 1000,
  seed = 38655,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)

## Visualize results (posterior contrasts)
newdata_fd <- expand.grid(
  Tactic = c("Con","Acc","IG"),
  Ref = NA,   # for population-level estimates
  N = 1
)

post_fd <- posterior_epred(m2, 
                           newdata = newdata_fd,
                           allow_new_levels = TRUE)

post_acc_fd <- post_fd[, newdata_fd$Tactic == "Acc"]
post_info_fd <- post_fd[, newdata_fd$Tactic == "IG"]
post_ctrl_fd <- post_fd[, newdata_fd$Tactic == "Con"]

diff_acc_ctrl_fd <- post_acc_fd - post_ctrl_fd
diff_info_ctrl_fd <- post_info_fd - post_ctrl_fd
diff_acc_info_fd <- post_acc_fd - post_info_fd

posterior_contrasts_fd <- data.frame(
  contrast = c("Acc vs Ctrl", "Info vs Ctrl", "Acc vs Info"),
  mean = c(mean(diff_acc_ctrl_fd), 
           mean(diff_info_ctrl_fd), 
           mean(diff_acc_info_fd)),
  lower = c(quantile(diff_acc_ctrl_fd, 0.025),
            quantile(diff_info_ctrl_fd, 0.025),
            quantile(diff_acc_info_fd, 0.025)),
  upper = c(quantile(diff_acc_ctrl_fd, 0.975),
            quantile(diff_info_ctrl_fd, 0.975),
            quantile(diff_acc_info_fd, 0.975))
)

posterior_contrasts_fd

df_post_fd <- tibble(
  `Acc. vs Ctrl` = as.vector(diff_acc_ctrl_fd),
  `Info. vs Ctrl` = as.vector(diff_info_ctrl_fd),
  `Acc. vs Info.` = as.vector(diff_acc_info_fd)
) %>% pivot_longer(cols = everything(), 
                   names_to = "contrast", 
                   values_to = "prob_diff")

df_post_fd$contrast <- factor(df_post_fd$contrast, 
                           levels = c("Info. vs Ctrl",
                                      "Acc. vs Ctrl",
                                      "Acc. vs Info."))

ggplot(df_post_fd, aes(x = prob_diff, fill = contrast)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Posterior Difference in Probability of FD", 
       y = "Density",
       fill = "Contrast") +
  scale_fill_manual(values = c("grey5", "grey80", "red")) + 
  theme_classic() + 
  theme(text = element_text(size = 27, family = "serif"),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
        plot.title = element_text(size = 26,
                                  margin = margin(b = 5, unit = "pt")),
        plot.subtitle = element_text(size = 23, 
                                     margin = margin(b = 20, unit = "pt")),
        plot.caption = element_text(size = 17, hjust = 0, 
                                    margin = margin(t = 20, unit = "pt")),
        strip.text = element_text(size = 23,
                                  margin = margin(t = 10, b = 10, 
                                                  unit = "pt")))

#-----------------------------------------------------------------------------#

# Combined figure

df_post$y <- "False Confession"
df_post_fd$y <- "False Denial"

png("Fig4_posterior_contrasts_experi.png", width = 1000, height = 500, type = "cairo")
ggplot(rbind(df_post, df_post_fd), 
       aes(x = prob_diff, 
           fill = contrast)) +
  geom_density(alpha = 0.5,
               color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Posterior Difference in Probability", 
       y = "Density",
       fill = NULL,
       title = "Figure 4.",
       subtitle = "Posterior contrast plot comparing accusatory, information gathering, and direct questioning\n(control) interrogation tactics under varying experimental paradigms",
       caption = "Note that only the 'cheating' experimental paradigm includes estimates of false denial, hence the reduced distributional\nheterogeneity in false confession compared to false denial.") +
  scale_fill_manual(values = c("grey5", "grey80", "red")) + 
  facet_wrap(~y, scales = "free") +
  theme_classic() + 
  theme(text = element_text(size = 27, family = "serif"),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
        plot.title = element_text(size = 26,
                                  margin = margin(b = 5, unit = "pt")),
        plot.subtitle = element_text(size = 23, 
                                     margin = margin(b = 20, unit = "pt")),
        plot.caption = element_text(size = 17, hjust = 0, 
                                    margin = margin(t = 20, unit = "pt")),
        strip.text = element_text(size = 23,
                                  margin = margin(t = 10, b = 10, 
                                                  unit = "pt")),
        legend.key.spacing.x = unit(1.5, "cm"),
        legend.position = "top")
dev.off()

#-----------------------------------------------------------------------------#

# Output results tables

library(flextable)
library(broom.mixed)

# Model 1: Pr(FC|Tactic)
## Fixed effects
fixed_tbl <- broom.mixed::tidy(m1, effects = "fixed", conf.level = 0.95) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  rename(
    Parameter = term,
    `Posterior Mean` = estimate,
    `95% CI Lower` = conf.low,
    `95% CI Upper` = conf.high
  )

## Random effects (SDs)
ranef_tbl <- broom.mixed::tidy(m1, effects = "ran_pars", conf.level = 0.95) %>%
  select(term, estimate) %>%
  mutate(
    Parameter = term,
    `Posterior Mean` = round(estimate, 2),
    `95% CI Lower` = NA_real_,
    `95% CI Upper` = NA_real_
  ) %>%
  select(Parameter, `Posterior Mean`, `95% CI Lower`, `95% CI Upper`)

## Model diagnostics
bayes_r2_val <- round(bayes_R2(m1)[1, "Estimate"], 2)

loo_val <- loo(m1)
loo_elpd <- round(loo_val$estimates["elpd_loo", "Estimate"], 2)
loo_se   <- round(loo_val$estimates["elpd_loo", "SE"], 2)

diag_tbl <- tibble(
  Parameter = c("Bayes R²", "LOO elpd (SE)"),
  `Posterior Mean` = c(bayes_r2_val, loo_elpd),
  `95% CI Lower` = c(NA_real_, loo_se),
  `95% CI Upper` = NA_real_
)

## Section headers + separator
header_row <- function(label) {
  tibble(
    Parameter = label,
    `Posterior Mean` = NA_real_,
    `95% CI Lower` = NA_real_,
    `95% CI Upper` = NA_real_
  )
}

separator <- header_row("---")

## Combine everything
full_tbl <- bind_rows(
  header_row("Fixed Effects"),
  fixed_tbl,
  separator,
  header_row("Random Effects"),
  ranef_tbl,
  separator,
  diag_tbl
)

## Spot fixing
full_tbl$Parameter[2] <- "Intercept"
full_tbl$Parameter[3] <- "Accusatory"
full_tbl$Parameter[4] <- "Info. Gather."
full_tbl$Parameter[7] <- "Study"
full_tbl$Parameter[10] <- "ELPD (SE)"

## Render table
set_flextable_defaults(font.family = "Times New Roman")
flextable(full_tbl) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("Table 1. Bayesian generalized linear mixed model (binomial link) estimating the posterior probability of false confession (FC) conditional on interrogation tactic") %>%
  save_as_docx(path = "Table1.docx")


# Model 2: Pr(FD|Tactic)
## Extract fixed effects
## Fixed effects
fixed_tbl <- broom.mixed::tidy(m2, effects = "fixed", conf.level = 0.95) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  rename(
    Parameter = term,
    `Posterior Mean` = estimate,
    `95% CI Lower` = conf.low,
    `95% CI Upper` = conf.high
  )

## Random effects (SDs)
ranef_tbl <- broom.mixed::tidy(m2, effects = "ran_pars", conf.level = 0.95) %>%
  select(term, estimate) %>%
  mutate(
    Parameter = term,
    `Posterior Mean` = round(estimate, 2),
    `95% CI Lower` = NA_real_,
    `95% CI Upper` = NA_real_
  ) %>%
  select(Parameter, `Posterior Mean`, `95% CI Lower`, `95% CI Upper`)

## Model diagnostics
bayes_r2_val <- round(bayes_R2(m2)[1, "Estimate"], 2)

loo_val <- loo(m2)
loo_elpd <- round(loo_val$estimates["elpd_loo", "Estimate"], 2)
loo_se   <- round(loo_val$estimates["elpd_loo", "SE"], 2)

diag_tbl <- tibble(
  Parameter = c("Bayes R²", "LOO elpd (SE)"),
  `Posterior Mean` = c(bayes_r2_val, loo_elpd),
  `95% CI Lower` = c(NA_real_, loo_se),
  `95% CI Upper` = NA_real_
)

## Section headers + separator
header_row <- function(label) {
  tibble(
    Parameter = label,
    `Posterior Mean` = NA_real_,
    `95% CI Lower` = NA_real_,
    `95% CI Upper` = NA_real_
  )
}

separator <- header_row("---")

## Combine everything
full_tbl <- bind_rows(
  header_row("Fixed Effects"),
  fixed_tbl,
  separator,
  header_row("Random Effects"),
  ranef_tbl,
  separator,
  diag_tbl
)

## Spot fixing
full_tbl$Parameter[2] <- "Intercept"
full_tbl$Parameter[3] <- "Accusatory"
full_tbl$Parameter[4] <- "Info. Gather."
full_tbl$Parameter[7] <- "Study"
full_tbl$Parameter[10] <- "ELPD (SE)"

## Render table
set_flextable_defaults(font.family = "Times New Roman")

flextable(full_tbl) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("Table 2. Bayesian generalized linear mixed model (binomial link) estimating the posterior probability of false denial (FD) conditional on interrogation tactic") %>%
  save_as_docx(path = "Table2.docx")