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

#==============================================================================#
# Read Data                                                                ====#
#==============================================================================#

# Mourtgos & Adams (2026)
## Structure (WC)
### StudyID  : study identifier (j)
### Bin      : bin i in study j
### Period   : publication year (used for recency weight)
### N        : total respondents in study j
### n_bin    : respondents in bin i
### p_est    : midpoint or conservative lower bound of bin i in study j,
###            or point estimate(s) for study j
### Type     : evidence-type label (survey = 0.9, inc_pers = 0.8, obs = 1)

wc_dat <- read.csv("./data/error_data.csv")

## Structure (FC)
### StudyID  : study identifier (j)
### Period   : publication year (used for recency weight)
### N        : total respondents in study j
### n_fc     : respondents who falsely confessed
### p_est    : point estimate(s) for study j
### Type     : evidence-type label (survey = 0.9, inc_pers = 0.8, obs = 1)

fc_dat <- read.csv("./data/fc_data.csv")    

#==============================================================================#
# 3.1.1 Estimating False Confession - Wrongful Conviction (FCWC)           ====#
#==============================================================================#

# 3.1.1.1. Conviction Error Rate ====

  ## Define function(s)
  min_max_scale <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  ## Rescale p_est, create weight variable(s)
  epsilon <- 1e-6
  wc_dat <- wc_dat %>%
    mutate(
      recency = min_max_scale(Period),
      pij = pmin(pmax(p_est, epsilon), 1 - epsilon),
      type_wt = ifelse(Type == "Survey", 0.9, 
                       ifelse(Type == "Observational", 1, 
                              0.8)),
      wij = recency * type_wt * (n_bin / N),  # total weight
    )
  
  ## formula
  wc_formula <- brmsformula(
    pij ~ 1 + (1 | StudyID),
    phi ~ 1 + wij + (1 | StudyID)
  )
  
  ## priors
  wc_priors <- c(
    ### Grand mean prior
    prior(normal(logit(0.03), 0.75), class = "Intercept"),
    
    ### Precision fixed effects
    prior(normal(0, 1), class = "Intercept", dpar = "phi"),
    prior(normal(0, 1), class = "b", dpar = "phi"),
    
    ### Random-effect SDs
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sd", dpar = "phi")
  )
  
  ## family
  wc_family <- Beta(link = "logit",
                 link_phi = "log")
  
    ### prior predictive check
    prior_check <- brm(
      data    = wc_dat,
      formula = wc_formula,
      prior   = wc_priors,
      family  = wc_family,
      sample_prior = "only",
      chains  = 4,
      iter    = 2000,
      warmup  = 0,
      seed    = 32608
    )
    
    pp_prior <- posterior_predict(prior_check)
    
      ### visualize
      stack(as.data.frame(t(pp_prior[1:200,]))) %>%
        ggplot(aes(x = values, fill = ind)) +
        geom_histogram() + 
        theme(legend.position = "none")
      
      ### Note: 
      ### Mourtgos & Adams (2026) report "that the induced distribution spans
      ### the full empirical range while concentrating probability mass between
      ### 0 and 0.25, where most estimates lie". However, what they fail to mention
      ### is the somewhat extreme zero-inflation (based on our replication).
      ### In our replication, this zero inflation is apparent in later posterior
      ### estimates, whereas they are not apparent in the publication we are
      ### replicating.It does not affect many (or any) of our estimates, but
      ### it does not pass a "face validity" check. 
  
  ## posterior sampling
  wc_model <- brm(
    data    = wc_dat,
    formula = wc_formula,
    prior   = wc_priors,
    family  = wc_family,
    chains  = 4,
    iter    = 251000,
    warmup = 1000,
    cores   = 8,
    seed    = 32608,
    control = list(adapt_delta = 0.995)
  )
  
    ### posterior predictive checks
    pp_post <- posterior_predict(wc_model)
    
      #### Overlay densities
      ppc_dens_overlay(
        y = wc_dat$pij,
        yrep = pp_post[1:200, ]
      )
  
  ## extract posterior distribution for conviction error rate
  post_mu_wc <- posterior_linpred(
    wc_model,
    transform = TRUE,
    re_formula = NULL
  )
  post_mu_wc <- as.vector(post_mu_wc)
  
    ### descriptives
    quantile(post_mu_wc, c(.01, .05, .5, .95, .99))
    mean(post_mu_wc)

#==============================================================================#

# 3.1.1.2. False confession prevalence within wrongful convictions ====

  ## Rescale p_est, introduce weight variable(s)
  fc_dat <- fc_dat %>%
    mutate(
      recency = min_max_scale(Period),
      pij = pmin(pmax(p_est, epsilon), 1 - epsilon),
      type_wt = ifelse(Type == "Survey", 0.9, 
                       ifelse(Type == "Observational", 1, 
                              0.8)),
      wj = recency * type_wt * (N / max(N)),  # total weight
    )
      
  ## Formula
  fc_formula <- brmsformula(
    pij ~ 1 + (1 | StudyID),
    phi ~ 1 + wj + (1 | StudyID)
  )
    
  ## Priors
  fc_priors <- c(
    ### Grand mean prior
    prior(normal(logit(0.15), 0.50), class = "Intercept"),
    
    ### Precision fixed effects
    prior(normal(0, 1), class = "Intercept", dpar = "phi"),
    prior(normal(0, 1), class = "b", dpar = "phi"),
    
    ### Random-effect SDs
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sd", dpar = "phi")
  )
  
  ## Family
  fc_family <- Beta(link = "logit",
                    link_phi = "log")
  
  ## Modeling false confession (abbreviated code)
  fc_model <- brm(
    formula = fc_formula,
    data = fc_dat,
    family = fc_family,
    prior = fc_priors,
    chains = 4,
    iter = 251000,
    warmup = 1000,
    cores = 8,
    control = list(adapt_delta = 0.995)
  )
  
  ## posterior predictive checks
  pp_post <- posterior_predict(fc_model)
  
    ### Overlay densities.
    ppc_dens_overlay(
      y = fc_dat$pij,
      yrep = pp_post[1:200, ]
    )

  ## extract posterior distribution for conviction error rate
  post_mu_fc <- posterior_linpred(
    fc_model,
    transform = TRUE,
    re_formula = NULL
  )
  post_mu_fc <- as.vector(post_mu_fc)
  
  quantile(post_mu_fc, c(.01, .05, .5, .95, .99))
  mean(post_mu_fc)

#==============================================================================#
  
# 3.1.1.3 Joint probability base rate
  
  draws_wc <- as_draws_df(wc_model)$b_Intercept %>%
    plogis() %>%
    sample(size = 1e6, replace = TRUE)
  
  draws_fc <- as_draws_df(fc_model)$b_Intercept %>%
    plogis() %>%
    sample(size = 1e6, replace = TRUE)
  
  ## Calculate FCWC base rate
  fcwc_base_rate <- draws_wc * draws_fc
  
  ## Descriptives for FCWC base rate
  mean(fcwc_base_rate)
  median(fcwc_base_rate)
  quantile(fcwc_base_rate, c(0.025, 0.975))
  
#==============================================================================#    
  
# 3.2 FCWC probabilities given a single interrogation technique 
  
  ## Set plausible levels for interrogation
  sens_ranges <- list(
    low      = c(0.1, 0.3),
    moderate = c(0.4, 0.6),
    high     = c(0.7, 0.9)
  )
  
  spec_ranges <- sens_ranges
    
    ### This represents the probability that a single interrogation
    ### technique within the class of potentially problematic
    ### tactics coincides with FCWC.
  
  ## Define sensitivity and specificity
  sample_from_tiers <- function(ranges, n) {
    tier <- sample(names(ranges), n, replace = TRUE)
    sapply(tier, function(t)
      runif(1, min = ranges[[t]][1], max = ranges[[t]][2])
    )
  }
  
  sensitivity  <- sample_from_tiers(sens_ranges, 1000000)
  specificity  <- sample_from_tiers(spec_ranges, 1000000)
  
  ## Set false positive rate
  false_positive_rate <- 1 - specificity
  
  ## Calculate the conditional probability using Bayes' Theorem
  p_fcwc_given_t <- (sensitivity * fcwc_base_rate) /
    ((sensitivity * fcwc_base_rate) +
       (false_positive_rate * (1 - fcwc_base_rate)))
  
  ## Describe posterior conditional probability distribution
  summary(p_fcwc_given_t)
  
  median(p_fcwc_given_t)
  quantile(p_fcwc_given_t, c(0.025, 0.5, 0.975))
  
  
#==============================================================================#
  
# 3.3 Upper-bound estimates under alternative attribution assumptions
  
  ## Note: The figure is not the "point" of this section, it is simply where
  ## the figure appears, and I am organizing this script to mirror the paper.
  
  ## Figure 1. Posterior distribution of FCWC across tactic scenarios
  
    ### Attribution assumptions
    rho_values <- c(
      baseline = 1.00,
      police_75 = 0.75,
      police_50 = 0.50
    )
    
    ### Scale base rate under alternative assumptions
    fcwc_scaled <- lapply(rho_values, function(rho) {
      rho * fcwc_base_rate
    })
    
    ### Re-estimate the conditional probability for each level of attribution
    p_fcwc_given_t_scaled <- lapply(fcwc_scaled, function(base_rate) {
      (sensitivity * base_rate) /
        ((sensitivity * base_rate) +
           (false_positive_rate * (1 - base_rate)))
    })
    
    ### Describe the posterior conditional probability by level of attribution
    summaries <- lapply(p_fcwc_given_t_scaled, function(x) {
      c(
        mean   = mean(x),
        median = median(x),
        q025   = quantile(x, 0.025),
        q975   = quantile(x, 0.975)
      )
    })
    
    do.call(rbind, summaries)
    
    ### Generate graph data
    gdat <- bind_rows(
      lapply(names(p_fcwc_given_t_scaled), function(name) {
        data.frame(
          fcwc_given_t = p_fcwc_given_t_scaled[[name]] * 100,
          attribution = factor(name,
                               levels = c("police_50",
                                          "police_75",
                                          "baseline"),
                               labels = c("p = 0.5 (50%)",
                                          "p = 0.75 (75%)",
                                          "p = 1 (Full)")))}))
    
    ### Calculate group medians
    group_medians <- gdat %>%
      group_by(attribution) %>%
      summarise(fcwc_given_t = median(fcwc_given_t, na.rm = TRUE))
    
    ### Calculate density by group
    group_density <- gdat %>%
      group_by(attribution) %>%
      summarise(
        dens = max(density(fcwc_given_t)$y)/2)
    
    ### Calculate HDI by group
    group_HDI <- gdat %>%
      group_by(attribution) %>%
      summarise(
        lower_hdi = hdi(fcwc_given_t, ci = 0.95)[[2]],
        upper_hdi = hdi(fcwc_given_t, ci = 0.95)[[3]])
    
    ### Merge labelling data
    labs_dat <- left_join(group_medians, group_density, by = "attribution")
    labs_dat <- left_join(labs_dat, group_HDI, by = "attribution")
    
    ### Replicate figure 1
    png("FigS1_replication.png", width = 1000, height = 750, type = "cairo")
    ggplot(gdat, aes(x=fcwc_given_t, fill = attribution)) +
      geom_density(color = NA,
                   alpha = 0.8) +
      scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                         labels = c("0%", "1%", "2%", "3%", "4%", "5%",
                                    "6%", "7%", "8%", "9%", "10%"),
                         limits = c(0, 10)) +
      scale_fill_manual(values = c("#5653ab",
                                   "#de7ba0",
                                   "#f6fa66")) + 
      labs(x = "Posterior Probability (%)",
           y = "Density",
           title = "Figure S1.",
           subtitle = "Posterior Distribution of FCWC Risk Across All Tactic Scenarios\n(Unadjusted replication of Mourtgos and Adams, 2026)",
           caption = "Three Overall Posteriors: p = 0.5, 0.75, 1\nRed dashed = median; annotation = 95% HDI") +
      facet_wrap(~attribution,
                 ncol = 1, 
                 scales = "free_y") + 
      geom_vline(data = labs_dat, 
                 aes(xintercept = fcwc_given_t), 
                 color = "red", linetype = "dashed", size = 1) +
      geom_text(data = labs_dat,
                aes(label = paste0("median: ", round(fcwc_given_t, 2),"%\n",
                                   "95% HDI: [",round(lower_hdi, 2),"%, ",
                                   round(upper_hdi, 2),"%]"),
                    y = dens),
                x = 9, hjust = 1, nudge_x = -0.1, family = "serif", size = 7) +
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
            legend.position = "none")
    dev.off()