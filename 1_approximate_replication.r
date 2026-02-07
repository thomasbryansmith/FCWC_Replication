# Install packages ====
install.packages(c("dplyr", "brms", "tidyr", "ggplot2", "posterior", "bayestestR"))

# load packages (simulations) ====
library(dplyr)
library(brms)
library(tidyr)
library(ggplot2)
library(posterior)
library(bayestestR)

#==============================================================================#
# 3.1.1 Estimating False Confession - Wrongful Conviction (FCWC)           ====#
#==============================================================================#

# 3.1.1.1. Conviction Error Rate ====

  ## data
    ### id : study identifier (j)
    ### p_ij     : midpoint or conservative lower bound of bin i in study j
    ### n_ij     : respondents in bin i
    ### N_j      : total respondents in study j
    ### recency  : scaled recency score in [0,1]
    ### type_wt  : evidence-type weight (1.0, 0.9, 0.8)
  
  
    ### Ramsey & Frank, 2007
    RF2007_national <- data_frame("id" = "RF2007_national",
                                  "year" = 2007,
                                  "pij_lower" = rep(c(0, 0, 0.005, 0.01, 0.04, 0.06, 0.11, 0.16, 0.21, 0.25), 4),
                                  "pij_upper" = rep(c(0, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.20, 0.25, 0.5), 4),
                                  "nij" = c(3, 57, 63, 79, 38, 14, 6, 6, 1, 0,
                                            1, 29, 30, 23, 10, 1, 2, 0, 0, 0,
                                            1, 4, 13, 31, 56, 42, 19, 28, 5, 17,
                                            0, 26, 31, 43, 28, 16, 6, 2, 0, 1),
                                  "Nj" = c(rep(267, 10),
                                           rep(96, 10),
                                           rep(216, 10),
                                           rep(153, 10)),
                                  "type_wt" = 0.9,
                                  "resptype" = c(rep("Police", 10),
                                                  rep("Prosecutors", 10),
                                                  rep("Defense Attorney", 10),
                                                  rep("Judge", 10)))
    
    RF2007_local <- data_frame("id" = "RF2007_local",
                               "year" = 2007,
                               "pij_lower" = rep(c(0, 0, 0.005, 0.01, 0.04, 0.06, 0.11, 0.16, 0.21, 0.25), 4),
                               "pij_upper" = rep(c(0, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.20, 0.25, 0.5), 4),
                               "nij" = c(91, 119, 37, 17, 9, 1, 0, 0, 0, 0,
                                         29, 49, 13, 6, 1, 2, 0, 0, 0, 0,
                                         4, 5, 21, 60, 42, 39, 21, 16, 12, 6,
                                         26, 52, 36, 32, 10, 9, 2, 1, 0, 0
                                             ),
                               "Nj" = c(rep(274, 10),
                                        rep(100, 10),
                                        rep(226, 10),
                                        rep(168, 10)),
                               "type_wt" = 0.9,
                               "resptype" = c(rep("Police", 10),
                                               rep("Prosecutors", 10),
                                               rep("Defense Attorney", 10),
                                               rep("Judge", 10)))
    
    
    ### Zalman et al., 2008
    ZSK2008_national <- data_frame("id" = "ZSK2008_national",
                                   "year" = 2008,
                                   "pij_lower" = rep(c(0, 0, 0.005, 0.01, 0.04, 0.06, 0.11, 0.16, 0.21, 0.25), 4),
                                   "pij_upper" = rep(c(0, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.20, 0.25, 0.5), 4),
                                   "nij" = c(1, 41, 46, 51, 28, 4, 1, 1, 0, 1,
                                             1, 10, 3, 5, 0, 3, 0, 0 ,0, 0,
                                             0, 2, 5, 16, 30, 26, 24, 15, 14, 13,
                                             0, 10, 7, 39, 21, 15, 7, 6, 0, 1),
                                   "Nj" = c(rep(174, 10),
                                            rep(22, 10),
                                            rep(145, 10),
                                            rep(106, 10)),
                                   "type_wt" = 0.9,
                                   "resptype" = c(rep("Police", 10),
                                                   rep("Prosecutors", 10),
                                                   rep("Defense Attorney", 10),
                                                   rep("Judge", 10)))
    
    ZSK2008_local <- data_frame("id" = "ZSK2008_local",
                                "year" = 2008,
                                "pij_lower" = rep(c(0, 0, 0.005, 0.01, 0.04, 0.06, 0.11, 0.16, 0.21, 0.25), 4),
                                "pij_upper" = rep(c(0, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.20, 0.25, 0.5), 4),
                                "nij" = c(73, 73, 19, 9, 3, 1, 0, 0, 0, 0,
                                          11, 10, 1, 1, 0, 0, 0, 0, 0, 0,
                                          0, 4, 6, 24, 35, 32, 14, 9, 10, 9,
                                          6, 32, 23, 30, 11, 7, 2, 1, 0, 1),
                                "Nj" = c(rep(178, 10),
                                         rep(23, 10),
                                         rep(143, 10),
                                         rep(113, 10)),
                                "type_wt" = 0.9,
                                "resptype" = c(rep("Police", 10),
                                                rep("Prosecutors", 10),
                                                rep("Defense Attorney", 10),
                                                rep("Judge", 10)))
    
    
    ### Compile and tweak data
    
    data <- rbind(RF2007_national,
                  RF2007_local,
                  ZSK2008_national,
                  ZSK2008_local)
    
    data$recency <- data$year / 2025
    
    epsilon <- 1e-6
    data <- data %>%
      mutate(
        pij = pmin(pmax((pij_lower + pij_upper) / 2, epsilon), 1 - epsilon),
        wij = recency * type_wt * (nij / Nj),  # total weight
      )
  
  ## brms model formula
    formula <- brmsformula(
      pij ~ 1 + (1 | id),
      phi ~ 1 + wij + (1 | id)
    )
    
    
  ## priors
    ### grand mean
    weighted.mean(data$pij, data$nij)
    ### observed maximum
    max(data$pij)
    
    priors <- c(
      # Grand mean prior
      prior(normal(logit(0.03), 0.75), 
            class = "Intercept"),
      
      # Precision fixed effects
      prior(normal(0, 1), class = "Intercept", dpar = "phi"),
      prior(normal(0, 1), class = "b", dpar = "phi"),
      
      # Random-effect SDs
      prior(exponential(2), class = "sd"),
      prior(exponential(2), class = "sd", dpar = "phi")
    )
    
  ## prior predictive checks (skipping this step later)
    priors_fc <- brm(
      formula = formula,
      data    = data,
      family  = Beta(),
      prior   = priors,
      sample_prior = "only",
      chains  = 4,
      iter    = 2000,
      warmup = 0,
      seed    = 123
    )
    
    pp_prior <- posterior_predict(priors_fc)
    
  ## visualize
    hist(pp_prior[1, ], breaks = 50,
         main = "Prior Predictive Distribution",
         xlab = "Wrongful Conviction Rate")
    
  ## posterior sampling
    wc_model <- brm(
      formula = formula,
      data    = data,
      family  = Beta(),
      prior   = priors,
      chains  = 4,
      iter    = 251000,
      warmup = 1000,
      cores   = 4,
      seed    = 123,
      control = list(adapt_delta = 0.99)
    )
    
  ## posterior predictive checks
    pp_post <- posterior_predict(wc_model)
    
    ### Overlay densities
    ppc_dens_overlay(
      y = data$pij,
      yrep = pp_post[1:200, ]
    )
    
    ### Grouped densities
    ppc_dens_overlay_grouped(y = data$pij, 
                             yrep = pp_post[1:200, ], 
                             group = data$resptype)
    
  ## extract posterior distribution for conviction error rate
    post_mu <- posterior_linpred(
      wc_model,
      transform = TRUE,
      re_formula = NULL
    )
    conviction_error_rate <- as.vector(post_mu)
    
    quantile(conviction_error_rate, c(.01, .05, .5, .95, .99))
    mean(conviction_error_rate)
  
    
#==============================================================================#
  
# 3.1.1.2. False confession prevalence within wrongful convictions ====
  ## data
    ### id        : factor identifying study j
    ### p_ij           : reported prevalence estimate (e.g., 0.14, 0.179, 0.25)
    ### N               : sample size underlying the estimate
    ### recency_weight  : linear recency weight scaled to [0, 1]
    ### type_weight     : evidence-type weight (observational = 1.0, survey = 0.9, convicted person = 0.8)
    
    data_fc <- data.frame("id" = c("BedauRadelet1987", "Connors1996", 
                                   "nre2024", "InnocenceProjectND"),
                          "pij" = c(0.14, 0.179,
                                    0.15, 0.29),
                          "N" = c(350, 28,
                                  147, 4045),
                          "year" = c(1987, 1996,
                                     2025, 2025),
                          "type" = c(1, 1, 
                                     1, 1)
                          )
    ### Recency weight variable
    data_fc$recency_weight <- data_fc$year / 2025
    
    ### Weight variable
    data_fc$wj <- with(data_fc,
                    recency_weight *
                      (N / max(N)) *
                      type
    )
    
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
      
      # Precision model priors
      prior(normal(0, 1), class = "Intercept", dpar = "phi"),
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
  
    ## Draw posteriors from FC and WC models
    draws_wc <- as_draws_df(wc_model)
    draws_wc <- plogis(draws_wc$b_Intercept)
      
    draws_fc <- as_draws_df(fc_model)
    draws_fc <- plogis(draws_fc$b_Intercept)
    
    ## Calculate FCWC base rate
    fcwc_base_rate <- draws_wc * draws_fc
  
    ## Descriptives for FCWC base rate
    mean(fcwc_base_rate)
    median(fcwc_base_rate)
    quantile(fcwc_base_rate, c(0.025, 0.975))
  
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
    
    ## Set false positive rate
    false_positive_rate <- 1 - specificity
    
    ## Calculate the conditional probability using Bayes' Theorem
    p_fcwc_given_t <- (sensitivity * fcwc_base_rate) /
      ((sensitivity * fcwc_base_rate) +
         (false_positive_rate * (1 - fcwc_base_rate)))
    
    ## Describe posterior conditional probabilty distribution
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
      png("Fig1_replication.png", width = 1000, height = 750, type = "cairo")
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
               title = "Figure 1.",
               subtitle = "Posterior Distribution of FCWC Risk Across All Tactic Scenarios\n(Replication of Mourtgos and Adams, 2026)",
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
      
