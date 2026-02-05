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
  
    ### Example structure
    #conviction_df <- tibble(
    #  study_id = factor(c(1,1,2,2,3,4,5)),
    #  p_mid    = c(.01, .03, .015, .025, .00, .04, .06),   # bin midpoints / lower bounds
    #  n_bin    = c(40, 60, 30, 70, 120, 80, 50),
    #  N_study  = c(100,100,100,100,120,80,50),
    #  recency  = c(.2,.2,.4,.4,.6,.8,1.0),                # linear rescaling
    #  type_wt  = c(1,1,.9,.9,1,1,.8)                       # evidence type weight
    #) %>%
    #  mutate(
    #    w_ij = recency * type_wt * (n_bin / N_study)
    #  )



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
      # prior(normal(0, 1), class = "Intercept", dpar = "phi"),
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
    
    ### Example data structure
    #conviction_df <- tibble(
    #  study_id = factor(c(1,1,2,2,3,4,5)),
    #  p_mid    = c(.01, .03, .015, .025, .00, .04, .06),   # bin midpoints / lower bounds
    #  n_bin    = c(40, 60, 30, 70, 120, 80, 50),
    #  N_study  = c(100,100,100,100,120,80,50),
    #  recency  = c(.2,.2,.4,.4,.6,.8,1.0),                # linear rescaling
    #  type_wt  = c(1,1,.9,.9,1,1,.8)                       # evidence type weight
    #) %>%
    #  mutate(
    #    w_ij = recency * type_wt * (n_bin / N_study)
    #  )
    
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

  # Accounting for the very real likelihood that the likelihood of a single
  # interrogation technique is inconsequential, and we should be modelling
  # the intensity of the interrogation.
    
    ## Kaplan & Cutler (2021) Co-occurrences among interrogation tactics in actual
    ## criminal investigations. Psychology, Crime & Law, 28(1): 1-19.
    
    ## 25 unique tactics encoded by Kaplan & Cutler, though they only observe 18.
    ## So we will set the theoretical maximum to 18.
    ## They report of 5 or 6 tactics per observation, but these cases were all
    ## severe (homicide, robbery and fraud, sexual assault, sex crimes, terrorism).
    
    ## Simulating plausible range of latent intensity
      
      ### We have already estimated the FCWC base rate
      pi_s <- fcwc_base_rate
    
      ### plausible interrogation intensity parameters
      mu_1 <- log(0.95 * 6)   # FCWC mean interrogation tactic count is equal
                              # to the average number of interrogation tactics 
                              # multiplied by the probability of a FCWC case
                              # involving an interrogation, ≈ 6 * 0.95.
                              # This is effectively a revised sensitivity.
      
      mu_0 <- log(0.14 * 6)   # non-FCWC mean interrogation tactic count is equal
                              # to the average number of interrogation tactics
                              # multiplied by the probability of any case
                              # involving an interrogation, ≈ 6 * 0.14.
                              # While it should be noted that "non-FCWC" is
                              # not the same as "any case" (it refers to true 
                              # confession and conviction), voluntary admissions
                              # presumably require no interrogation and make up
                              # a non-trivial proportion of confession-convictions.
                              # This will "drive down" the number of tactics in
                              # non-FCWC cases. This is revised specificity.
      
      sigma_I <- 1            # heterogeneity in intensity (standard deviation),
                              # we use the same value as Moutgos & Adams (2026).
      
      #### Note: we are assuming that FCWC interrogations are ~1.5x as intense.
      #### This will need to be revised so that it is in line with the literature.
      
      ### Monte Carlo to approximate likelihoods
      #### Define function to generate average likelihoods for each level of K
      #### Poisson-lognormal mixture:
      likelihood_K <- function(k, mu, sigma, n_mc = 50000) {
        I <- rnorm(n_mc, mu, sigma) # simulated (plausible) distribution of interrogation intensity
        lambda <- exp(I)
        mean(dpois(k, lambda)) # average on the Poisson
      }
      
      likelihood_K_full <- function(k, mu, sigma, n_mc = 50000) {
        I <- rnorm(n_mc, mu, sigma) # simulated (plausible) distribution of interrogation intensity
        lambda <- exp(I)
        dpois(k, lambda) # average on the Poisson
      }
      
      #### Calculate posterior FCWC probability given k instances of a tactic
      posterior_fcwc <- function(k, pi, mu1, mu0, sigma) {
        
        L1 <- likelihood_K(k, mu1, sigma) # P(K = k|FCWC)
        L0 <- likelihood_K(k, mu0, sigma) # P(K = k|¬FCWC) 
        
        # L1 and L0 use difference probability distributions / point estimates
        # because we have reason to believe that FCWC and non-FCWC are subject
        # to different interrogation intensities.
        # Pi is P(FCWC), as estimated via Mourtgos & Adams' (2026) joint probability.
        
        (L1 * pi) / (L1 * pi + L0 * (1 - pi)) # Bayes' Theorem
      }
      
      
      #### Visualize likelihood distributions for FCWC v. non-FCWC
      ##### Ridge Plot
      n_tactics <- seq(1, 18, by = 1)
      
      L1_dist <- lapply(n_tactics, function(x) likelihood_K_full(x, mu_1, sigma_I))
      L0_dist <- lapply(n_tactics, function(x) likelihood_K_full(x, mu_0, sigma_I))
      
      data.frame("FCWC" = c(rep("FCWC", length(unlist(L1_dist))),
                            rep("Non-FCWC", length(unlist(L0_dist)))),
                 "NTact" = c(rep(n_tactics, each = length(L1_dist[[1]])),
                             rep(n_tactics, each = length(L0_dist[[1]]))),
                 "likelihood" = c(unlist(L1_dist),
                                  unlist(L0_dist))) %>%
        ggplot(aes(x = likelihood, y = factor(NTact), fill = factor(NTact))) +
        geom_density_ridges(scale = 5) +
        facet_wrap(~FCWC)
      theme_ridges()
      
      ##### Line Graph
      n_tactics <- seq(1, 18, by = 1)
             
      L1_dist <- unlist(lapply(n_tactics, function(x) likelihood_K(x, mu_1, sigma_I)))
      L0_dist <- unlist(lapply(n_tactics, function(x) likelihood_K(x, mu_0, sigma_I)))
      
      png("Fig1a_Likelihood_function_latent_intensity.png", 1000, 500, type = "cairo")
      data.frame("FCWC" = c(rep("FCWC", length(L1_dist)),
                            rep("Non-FCWC", length(L0_dist))),
                 "NTact" = c(n_tactics,
                             n_tactics),
                 "likelihood" = c(L1_dist,
                                  L0_dist)) %>%
        ggplot(aes(x = NTact, y = likelihood, color = FCWC)) +
        geom_point(size = 2) +
        geom_line(size = 1) +
        labs(title = "Figure 1a.",
             subtitle = "Simulated likelihood function for latent interrogation intensity",
             color = NULL,
             y = "P(K=k|FCWC)",
             x = "Number of Interrogation Tactics") +
        theme_classic() + 
        theme(text = element_text(size = 27, family = "serif"),
              axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
              axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
              plot.title = element_text(size = 27,
                                        margin = margin(b = 5, unit = "pt")),
              plot.subtitle = element_text(size = 23, 
                                           margin = margin(b = 20, unit = "pt")),
              plot.caption = element_text(size = 17, hjust = 0, 
                                          margin = margin(t = 20, unit = "pt")),
              legend.key.spacing.y = unit(1.0, "cm"))
      dev.off()
      
      #### Visualize posterior distributions for each
      png("Fig1b_Interrogation_intensity_figure.png", 1000, 500, type = "cairo")
      data.frame("NTact" = factor(rep(1:6, each = length(pi_s))),
                 "posterior" = unlist(lapply(1:6, 
                                             function(x) 
                                               posterior_fcwc(x, pi_s, 
                                                              mu_1, mu_0, 
                                                              sigma_I)))) %>%
        ggplot(aes(x = posterior, fill = NTact)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_x_continuous(breaks = seq(0, 0.15, by = 0.01),
                           labels = paste0(0:15,"%"),
                           limits = c(0, 0.15)) +
        scale_fill_viridis_d() +
        labs(y = "Density",
             x = "P(FCWC)",
             title = "Figure 1b.",
             subtitle = "Posterior distribution of FCWC risk by the number of tactics employed during interrogation",
             fill = "Tactics",
             caption = "") + 
        theme_classic() + 
        theme(text = element_text(size = 27, family = "serif"),
              axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
              axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
              plot.title = element_text(size = 27,
                                        margin = margin(b = 5, unit = "pt")),
              plot.subtitle = element_text(size = 23, 
                                           margin = margin(b = 20, unit = "pt")),
              plot.caption = element_text(size = 17, hjust = 0, 
                                          margin = margin(t = 20, unit = "pt")),
              legend.title = element_text(margin = margin(b = 0.3, unit = "cm")),
              legend.key.spacing.y = unit(0.3, "cm"))
      dev.off()