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

# 3.1.1.2. False confession prevalence within wrongful convictions ====

  ## data
    ### id        : factor identifying study j
    ### p_ij           : reported prevalence estimate (e.g., 0.14, 0.179, 0.25)
    ### N               : sample size underlying the estimate
    ### recency_weight  : linear recency weight scaled to [0, 1]
    ### type_weight     : evidence-type weight (observational = 1.0, survey = 0.9, convicted person = 0.8)
    ### Data
    
    data_fc <- data.frame("id" = c("BedauRadelet1987", "Connors1996", "nre2024", "InnocenceProjectND"),
                          "pij" = c(0.14, 0.179, 0.15, 0.29),
                          "N" = c(350, 28, 147, 4045),
                          "year" = c(1987, 1996, 2025, 2025),
                          "type" = c(1, 1, 1, 1)
                          )
    
    ### Recency weight variable
    data_fc$recency_weight <- data_fc$year / 2025
    
    ### Weight variable
    data_fc$wj <- with(data_fc,
                    recency_weight *
                      (N / max(N)) *
                      type)
    
    ### Adjust p_ij so that beta distribution is viable (does not accept 0 or 1)
    epsilon <- 1e-6
    data_fc$pij <- pmin(pmax(data_fc$pij, epsilon), 1 - epsilon)
  
  ## Modeling false confession (abbreviated code)
  fc_model <- brm(
    bf(
      
      pij ~ 1 + (1 | id),
      phi ~ 1 + wj + (1 | id)
      
    ),
    data = data_fc,
    family = Beta(),
    prior = c(
      
      # Grand mean prior centered at 15%
      prior(normal(logit(0.15), 0.50), class = "Intercept"),
      
      # Random-effect SDs
      prior(exponential(2), class = "sd"),
      prior(exponential(2), class = "sd", dpar = "phi"),
      
      # Precision model priors
      # prior(normal(0, 1), class = "Intercept", dpar = "phi"),
      prior(normal(0, 1), class = "b", dpar = "phi")
      
    ),
    chains = 4,
    iter = 251000,
    warmup = 1000,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )
  
  ## posterior predictive checks
  pp_post <- posterior_predict(fc_model)
  
    ### Overlay densities.
      #### This looks a little chaotic given the small number of data points.
      #### Do M&A report more that I am not aware of?
    ppc_dens_overlay(
      y = data_fc$pij,
      yrep = pp_post[1:200, ]
    )
  
  ## extract posterior distribution for conviction error rate
  post_mu <- posterior_linpred(
    fc_model,
    transform = TRUE,
    re_formula = NULL
  )
  conviction_error_rate <- as.vector(post_mu)
  
  quantile(conviction_error_rate, c(.01, .05, .5, .95, .99))
  mean(conviction_error_rate)

#==============================================================================#
  
  # 3.1.1.3 Joint probability base rate

    draws_fc <- as_draws_df(fc_model)
    draws_fc <- plogis(draws_fc$b_Intercept)
    
#==============================================================================#

  # 3.2 FCWC probabilities given a single interrogation technique
  
    ## Set plausible levels for interrogation
      ### This represents the probability that a single interrogation
      ### technique within the class of potentially problematic
      ### tactics coincides with FCWC.
      sens_ranges <- list(
        low      = c(0.1, 0.3),
        moderate = c(0.4, 0.6),
        high     = c(0.7, 0.9)
      )
      
      spec_ranges <- sens_ranges
    
    ## Define sensitivity and specificity
    sample_from_tiers <- function(ranges, n) {
      tier <- sample(names(ranges), n, replace = TRUE)
      sapply(tier, function(t)
        runif(1, min = ranges[[t]][1], max = ranges[[t]][2])
      )
    }

    sensitivity  <- sample_from_tiers(sens_ranges, 1000000)
    specificity  <- sample_from_tiers(spec_ranges, 1000000)

#==============================================================================#
    
    # Narrowing focus to wrongful convictions, P(FC|WC) rather than P(FCWC) = P(WC) * P(FC|WC)
    
    ## Figure 2. Posterior distribution of P(FC|WC, T) with attribution
    
    ## Draw posterior distribution for P(FC|WC)
    p_FC_given_WC <- draws_fc 
    
    ### Attribution assumptions
    rho_values <- c(
      baseline = 1.00,
      police_75 = 0.75,
      police_50 = 0.50
    )
    
    ### Scale base rate under alternative assumptions
    fc_wc_t_scaled <- lapply(rho_values, function(rho) {
      rho * p_FC_given_WC
    })
    
    ### Re-estimate the conditional probability for each level of attribution
    p_fcwc_given_t_scaled <- lapply(fc_wc_t_scaled, function(base_rate) {
      (sensitivity * base_rate) /
        ((sensitivity * base_rate) +
           ((1 - specificity) * (1 - base_rate)))
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
    
    ### Figure 2
    png("Fig3_FC_cond_WC_Risk.png", width = 1000, height = 750, type = "cairo")
    ggplot(gdat, aes(x=fcwc_given_t, fill = attribution)) +
      geom_density(color = NA,
                   alpha = 0.8) +
      scale_x_continuous(breaks = seq(0, 75, by = 15),
                         labels = paste0(seq(0, 75, by = 15), "%")) +
      scale_fill_manual(values = c("#5653ab",
                                   "#de7ba0",
                                   "#f6fa66")) + 
      labs(x = "Posterior Probability (%)",
           y = "Density",
           title = "Figure 3.",
           subtitle = "Posterior distribution of FC|WC risk across all tactic scenarios",
           caption = "Posteriors represent the probability of false confession conditional on wrongful conviction and interrogation tactic, P(FC|WC, T),\nas opposed to the probability of false confession-wrongful conviction conditional upon interrogation tactic, or P(FCWC|T).\nIntuitively, P(FCWC|T) represents an average treatment effect (ATE), where P(FC|WC, T) represents the treatment on treated (ATT).\nThree Overall Posteriors: p = 0.5, 0.75, 1. Red dashed = median; annotation = 95% HDI.") +
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
                x = 70, hjust = 1, nudge_x = -0.1, family = "serif", size = 7) +
      theme_classic() + 
      theme(text = element_text(size = 27, family = "serif"),
            axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
            axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
            legend.position = "none",
            plot.title = element_text(size = 27),
            plot.subtitle = element_text(size = 23),
            strip.text = element_text(size = 23, 
                                      margin = margin(t = 10, b = 10, 
                                                      unit = "pt")),
            plot.caption = element_text(size = 17, hjust = 0, 
                                        margin = margin(t = 20, unit = "pt")))
    dev.off()