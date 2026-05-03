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
library(ggh4x)

#==============================================================================#
# Read Data                                                                ====#
#==============================================================================#

# Catlin (2026)
## Structure
### MAref         : ID in Mourtgos & Adams' data
### Ref           : ID
### StudyYear     : publication year
### DataYear      : range of years that data were collected
### Included      : included in Mourtgos & Adams' data
### Estimate      : FC ("FCshare") v. WC ("error")
### Location      : data geolocation
### Participants  : respondent type
### N             : sample size
### bin           : bin i in study j
### nbin          : respondents in bin i
### %error        : % WC
### nerror        : N WC
### weight        : evidence-type for weighting
### independent   : ???
### shareddata    : data intersections between studies
### notes         : notes about data
### comments      : additional observations

mar_dat <- readxl::read_xlsx("./data/Data_2026_02_25_cleaned.xlsx")
mar_dat <- mar_dat[which(mar_dat$include),]

wc_dat <- mar_dat[which(mar_dat$Estimate == "Error"),]

fc_dat <- mar_dat[which(mar_dat$Estimate == "FCshare"),]
fc_dat <- fc_dat %>%
  group_by(Ref, Location) %>%
  mutate(p_est_ungr = p_est,
         n_ungr = n,
         N_ungr = N,
         n = sum(n),
         N = sum(N),
         p_est = weighted.mean(p_est, N)) %>%
  ungroup() %>%
  filter(US == 0)

#==============================================================================#
# 3.1.1 Estimating False Confession - Wrongful Conviction (FCWC)           ====#
#==============================================================================#

# 3.1.1.1. Conviction Error Rate ====

  ## Define function(s)
  min_max_scale <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  ## Generate Period variable
  wc_dat$Period <- as.numeric(gsub("[0-9]{4}-","", wc_dat$DataYear))
  wc_dat$Period[is.na(wc_dat$Period)] <- wc_dat$StudyYear[is.na(wc_dat$Period)]
  
  ## Rescale p_est, create weight variable(s)
  epsilon <- 1e-6
  wc_dat <- wc_dat %>%
    mutate(
      recency = min_max_scale(Period),
      pij = pmin(pmax(p_est, epsilon), 1 - epsilon),
      type_wt = ifelse(weight == "secondary", 1, 0.9),
      wij = recency * type_wt * (n / N),  # total weight
    )
  
  ## formula
  wc_formula <- brmsformula(
    pij ~ 1 + (1 | Ref),
    phi ~ 1 + wij + (1 | Ref)
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
    
  ## Generate Period variable
  fc_dat$Period <- as.numeric(gsub("[0-9]{4}-","", fc_dat$DataYear))
  fc_dat$Period[is.na(fc_dat$Period)] <- fc_dat$StudyYear[is.na(fc_dat$Period)]
  
  ## Rescale p_est, introduce weight variable(s)
  fc_dat <- fc_dat %>%
    mutate(
      recency = min_max_scale(Period),
      pij = pmin(pmax(p_est, epsilon), 1 - epsilon),
      type_wt = ifelse(weight == "secondary", 1, 0.9),
      wj = recency * type_wt * (N / max(N)),  # total weight
    )
  
  ## Formula
  fc_formula <- brmsformula(
    pij ~ 1 + (1 | Ref),
    phi ~ 1 + wj + (1 | Ref)
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

# Critique 2 
  ## We could argue that the probability of a single interrogation tactic 
  ## is inconsequential. Miranda v. Arizona (1966) set the precedent for
  ## interrogations being "inherently coercive".
  
  ## Mourtgos & Adams (2026) argue that their specificity / sensitivity
  ## Likelihood statistics (which are uniformly sampled for both FCWC and
  ## non-FCWC cases) are the probability of a "problematic" interrogation
  ## tactic, but if this were true they would have attempted to model
  ## the fact that these more problematic are more likely to be present in
  ## FCWC cases. Even if divorced from the intent of the law enforcement
  ## officers, the existing literature clearly shows that FCWC interrogations
  ## are, on average, much longer than non-FCWC interrogations (see: 
  ## Leo, 2004; Kassin et al., 2007). Longer interrogations necessarily
  ## provide greater opportunity (and arguably an impetus) for problematic
  ## interrogation tactics.

  ## This could be approached by modelling the number of interrogation tactics:
  ## Kaplan & Cutler (2021) Co-occurrences among interrogation tactics in actual
  ## criminal investigations. Psychology, Crime & Law, 28(1): 1-19.
    
  ## However, we instead opt to simulate the plausible range of latent 
  ## interrogation intensities for FCWC and non-FCWC separately based on
  ## previous research on the average interrogation length for each.
      
    ### We have already estimated the FCWC base rate
    pi_fcwc <- fcwc_base_rate
    
    ### plausible interrogation intensity parameters
    mu_1 <- log(16.23)  # Average length of an FCWC interrogation is 16.3 hours.
    sigma_1 <- 15.05    # Standard deviation for this estimate is 15.05 hours.
                        # This is effectively a revised sensitivity. (N = 44)
                        # Leo, Richard A. (2004). The Problem of False 
                        # Confessions in the Post-DNA World. N.C.L.Rev 891.
      
    mu_0 <- log(1.6)   # Average length of a non-FCWC interrogation is 1.6 hours.
    sigma_0 <- 0.89    # Standard deviation for this estimate is 0.89 hours.
                       # This is effectively specificity. (N = 601)
                       # Kassin, Saul M., Leo, Richard, A., Meissner, 
                       # Christian A., et al.(2007). Police Interviewing and
                       # Interrogation: A Self-Report Survey of Police Practices
                       # and Beliefs. Law and Human Behavior, 31: 381 - 400.
        
    ### Monte Carlo to approximate likelihoods
      #### Goal: Define function to generate average likelihoods for each level of K
      #### Truncated lognormal distribution:
      rlnorm_trunc <- function(n, meanlog, sdlog, lower = 0, upper = 48) {
          
          if (lower < 0) stop("Lower bound must be >= 0 for lognormal")
            
          # CDF bounds
          p_lower <- if (lower == 0) 0 else plnorm(lower, meanlog, sdlog)
          p_upper <- plnorm(upper, meanlog, sdlog)
            
          if (p_upper <= p_lower)
            stop("Invalid bounds: upper CDF <= lower CDF")
            
          # Draw from truncated uniform
          u <- runif(n, min = p_lower, max = p_upper)
            
          # Inverse CDF
          return(qlnorm(u, meanlog, sdlog))
          
        }
        
      #### Monte Carlo likelihood functions:
      mc_likelihood_H <- function(h, n_mc, meanlog, sdlog, upper = 48){
            
          sims <- rlnorm_trunc(n = n_mc,
                               meanlog = meanlog,
                               sdlog = sdlog,
                               upper = upper)
          dens <- density(sims, from = 0, to = upper, n = 2048)
          return(approx(dens$x, dens$y, xout = h, rule = 2)$y)
            
        }
          
      mc_likelihood_H_distribution <- function(h, n_rep = 1000, n_mc = 5000,
                                               meanlog, sdlog, upper = 48,
                                               bw = "nrd0"){
          
        # n_rep = number of likelihood draws
        # n_mc = MC size per draw
          
        mc_draw <- function(h, n_mc, meanlog, sdlog, upper, bw){
            
            sims <- rlnorm_trunc(n = n_mc, meanlog = meanlog, sdlog = sdlog, 
                                 upper = upper)
            dens <- density(sims, from = 0, to = upper, n = 2048, bw = bw)
            return(approx(dens$x, dens$y, xout = h, rule = 2)$y)
              
          }
            
          return(replicate(n_rep,
                           mc_draw(h = h, n_mc = n_mc, 
                                   meanlog = meanlog, sdlog = sdlog,
                                   upper = upper, bw = bw)))
          
        }
              
      #### Define Bayes theorem with stochastic likelihoods:
      posterior_fcwc_given_H_mc <- function(h, pi, mu1, sigma1, mu0, sigma0,
                                            n_mc = 20000, upper = 48){
            
          L1 <- mc_likelihood_H(h, n_mc, mu1, sigma1, upper)
          L0 <- mc_likelihood_H(h, n_mc, mu0, sigma0, upper)
            
          (L1 * pi) / (L1 * pi + L0 * (1 - pi))
            
      }
      
      #### Visualize likelihood distributions for FCWC v. non-FCWC
        ##### Ridge Plot
        hour_ranges_r <- seq(6, 48, by = 2)
        
        L1_dist <- lapply(hour_ranges_r, function(h) mc_likelihood_H_distribution(h, 1000, 5000, mu_1, sigma_1, 48))
        L0_dist <- lapply(hour_ranges_r, function(h) mc_likelihood_H_distribution(h, 1000, 5000, mu_0, sigma_0, 48))
        
        png("FigS2_Likelihood_function_latent_intensity_w_stochasticity.png", 1000, 500, type = "cairo")
        data.frame("FCWC" = c(rep("FCWC", length(unlist(L1_dist))),
                              rep("Non-FCWC", length(unlist(L0_dist)))),
                   "NTact" = c(rep(hour_ranges_r, each = length(L1_dist[[1]])),
                               rep(hour_ranges_r, each = length(L0_dist[[1]]))),
                   "likelihood" = c(unlist(L1_dist),
                                    unlist(L0_dist))) %>%
          ggplot(aes(x = likelihood, y = factor(NTact), fill = factor(NTact))) +
          geom_density_ridges(scale = 5, color = NA) +
          scale_fill_viridis_d(option = "D") +
          labs(title = "Figure S2.",
               subtitle = "Simulated likelihood function for interrogation length including random distributions\nat each hourly interval",
               y = "Interrogation Hours",
               x = "Likelihood, P(h|θ)",
               caption = "We define separate likelihood functions conditional upon the mean, μ, and standard deviation, σ, for hours interrogated for FCWC\nand non-FCWC respectively. This is more appropriate given that prior literature has observed different average interrogation lengths\ndependent on FCWC. Per Kassin et al. (2007), the average length of an interrogation is approximately 1.6 hours (σ = 0.89).Whereas,\nper Leo (2004), the average length of a FCWC interrogation is 16.23 hours (σ = 15.05). We retain the uncertainty of these estimates\nvia Monte Carlo marginalization, treating the interrogation hours as stochastic. Hours 2 and 4 are omitted from this visualization\nfor clarity (the probability of a 1-to-5 hour interrogation is inordinately higher for non-FCWC than FCWC), a more complete\ndistributional comparison is reported in the main text.") +
          coord_flip() +
          facet_wrap(~FCWC) + 
          theme_classic() + 
          theme(legend.position = "none",
                text = element_text(size = 27, family = "serif"),
                axis.text.x = element_text(size = 15),
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
      
      
        ##### Line Graph
        hour_ranges <- seq(0, 48, by = 0.05)
        
        L1_dist_pe <- unlist(lapply(hour_ranges, function(x) mc_likelihood_H(x, 20000, mu_1, sigma_I)))
        L0_dist_pe <- unlist(lapply(hour_ranges, function(x) mc_likelihood_H(x, 20000, mu_0, sigma_I)))
        
        png("Fig2a_Likelihood_function_interrogation_hours.png", 1000, 500, type = "cairo")
        data.frame("FCWC" = c(rep("FCWC", length(L1_dist_pe)),
                              rep("Non-FCWC", length(L0_dist_pe))),
                   "NTact" = c(hour_ranges,
                               hour_ranges),
                   "likelihood" = c(L1_dist_pe,
                                    L0_dist_pe)) %>%
          ggplot(aes(x = NTact, y = likelihood, color = FCWC)) +
          #geom_point(size = 0.25) +
          geom_line(linewidth = 0.75) +
          scale_x_continuous(breaks = seq(0, 48, by = 2)) +
          labs(title = "Figure 2a.",
               subtitle = "Simulated likelihood function for interrogation length",
               caption = "We define separate likelihood functions conditional upon the mean, μ, and standard deviation, σ, for hours\ninterrogated for FCWC and non-FCWC respectively. This is more appropriate given that prior literature has\nobserved different average interrogation lengths dependent on FCWC. Per Kassin et al. (2007), the average\nlength of an interrogation is approximately 1.6 hours (σ = 0.89).Whereas, per Leo (2004), the average length\nof an FCWC interrogation is 16.23 hours (σ = 15.05). We retain the uncertainty of these estimates via Monte\nCarlo marginalization, treating the interrogation hours as stochastic. See Figure S1 for hourly distributions.",
               color = NULL,
               y = "Likelihood, P(h|θ)",
               x = "Interrogation Hours") +
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
        
        
#==============================================================================#
        
        
    ### Visualize posterior distributions for each hour of interrogation
      #### Draw posteriors for P(FCWC|H)
      hour_ranges <- c(1.6, 16.23)
      
      posterior_fcwc_given_H_mc_fi <- function(h, pi, mu1, sigma1, mu0, sigma0,
                                               n_mc = 5000, upper = 48){
        
        L1 <- mc_likelihood_H(h, n_mc, mu1, sigma1, upper)
        L0 <- mc_likelihood_H(h, n_mc, mu0, sigma0, upper)
        
        return(c("h" = h,
                 "pi" = pi,
                 "mu1" = mu1,
                 "sigma1" = sigma1,
                 "mu0" = mu0,
                 "sigma0" = sigma0,
                 "L1" = L1,
                 "L0" = L0,
                 "posterior" = (L1 * pi) / (L1 * pi + L0 * (1 - pi))))
      }
      
      posterior_draws <- lapply(hour_ranges, function(h){
        lapply(sample(pi_fcwc, 5000), function(pi_s) {
          posterior_fcwc_given_H_mc_fi(
            h = h,
            pi = pi_s,
            mu1 = mu_1,  sigma1 = sigma_1,
            mu0 = mu_0,  sigma0 = sigma_0
          )
        })})
      posterior_draws <- lapply(posterior_draws, function(x) as.data.frame(do.call(rbind, x)))
      posterior_draws <- do.call("rbind", posterior_draws)
      
      #### Visualization
      gdat <- posterior_draws
      gdat$NTact <- paste(gdat$h, "hours")

      ### Calculate group medians
      group_medians <- gdat %>%
        group_by(NTact) %>%
        summarise(median = median(posterior, na.rm = TRUE))
      
      ### Calculate density by group
      group_density <- gdat %>%
        group_by(NTact) %>%
        summarise(dens_y = max(density(posterior)$y) * 0.8,
                  dens_x = max(density(posterior)$x) * 0.95)
      
      ### Calculate HDI by group
      group_HDI <- gdat %>%
        group_by(NTact) %>%
        summarise(
          lower_hdi = hdi(posterior, ci = 0.95)[[2]],
          upper_hdi = hdi(posterior, ci = 0.95)[[3]])
      
      ### Merge labelling data
      labs_dat <- left_join(group_medians, group_density, by = "NTact")
      labs_dat <- left_join(labs_dat, group_HDI, by = "NTact")
      labs_dat$dens_x[which(labs_dat$NTact=="16.23 hours")] <- 0.5
      
      png("Fig2b_Interrogation_intensity_figure.png", 1000, 500, type = "cairo")
        ggplot(gdat, aes(x = posterior)) +
        geom_density(alpha = 0.8, color = "black", fill = "grey20") + 
        scale_fill_viridis_d() +
        labs(y = "Density",
             x = "P(FCWC)",
             title = "Figure 2b.",
             subtitle = "Posterior distribution of FCWC risk by length of interrogation",
             fill = "Tactics",
             caption = "") + 
        facet_wrap2(~NTact, scale = "free") +
        facetted_pos_scales(
            x = list(
              NTact == "16.23 hours" ~ scale_x_continuous(limits = c(0, 0.5)) # Second plot
            )
          ) +
        geom_vline(data = labs_dat, 
                     aes(xintercept = median), 
                     color = "red", linetype = "dashed", size = 1) +
        geom_text(data = labs_dat,
                  aes(label = paste0("median: ", round(median, 3),"\n",
                                     "95% HDI: [",round(lower_hdi, 3),", ",
                                     round(upper_hdi, 3),"]"),
                      y = dens_y, x = dens_x),
                  hjust = 1, family = "serif", size = 7) +
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

      
#==============================================================================#
      
      
      ## GAM distribution visualization
      hour_ranges <- 1:48
      posterior_draws <- lapply(hour_ranges, function(h){
        sapply(sample(pi_fcwc, 1000), function(pi_s) {
          posterior_fcwc_given_H_mc(
            h = h,
            pi = pi_s,
            mu1 = mu_1,  sigma1 = sigma_1,
            mu0 = mu_0,  sigma0 = sigma_0
          )
        })})
      
      png("Fig2c_Posterior_distributions_gam.png", 1000, 500, type = "cairo")
      data.frame("FCWC" = rep(hour_ranges, each = length(posterior_draws[[1]])),
                 "likelihood" = unlist(posterior_draws)) %>%
        ggplot(aes(y = likelihood, x = FCWC)) +
        geom_point(alpha = 0.05) +
        geom_smooth(color = "black") +
        scale_fill_viridis_d(option = "D") +
        labs(title = "Figure 2c.",
             subtitle = "Posterior probability distributions for each hour of interrogation with fitted Generalized\nAdditive Model",
             x = "Interrogation Hours",
             y = "P(FCWC|H)") +
        xlim(0, 36) +
        theme_classic() + 
        theme(legend.position = "none",
              text = element_text(size = 27, family = "serif"),
              axis.text.x = element_text(size = 15),
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
      
      
      
#==============================================================================#
      
      
      ### Visualize posterior distributions at 1.6 and 16 hours across sensitivity and specificity
      hour_ranges <- c(1.6, 6, 12, 16.23)
      
      posterior_fcwc_given_H_mc_fi <- function(h, pi, mu1, sigma1, mu0, sigma0,
                                            n_mc = 5000, upper = 48){
        
        L1 <- mc_likelihood_H(h, n_mc, mu1, sigma1, upper)
        L0 <- mc_likelihood_H(h, n_mc, mu0, sigma0, upper)
        
        return(c("h" = h,
                 "pi" = pi,
                 "mu1" = mu1,
                 "sigma1" = sigma1,
                 "mu0" = mu0,
                 "sigma0" = sigma0,
                 "L1" = L1,
                 "L0" = L0,
                 "posterior" = (L1 * pi) / (L1 * pi + L0 * (1 - pi))))
      }
      
      posterior_draws <- lapply(hour_ranges, function(h){
        lapply(sample(pi_fcwc, 5000), function(pi_s) {
          posterior_fcwc_given_H_mc_fi(
            h = h,
            pi = pi_s,
            mu1 = mu_1,  sigma1 = sigma_1,
            mu0 = mu_0,  sigma0 = sigma_0
          )
        })})
      posterior_draws <- lapply(posterior_draws, function(x) as.data.frame(do.call(rbind, x)))
      posterior_draws <- do.call("rbind", posterior_draws)
      
      library(interp)
      
      interp_data <- posterior_draws %>%
        mutate(
          posterior = (L1 * pi) / (L1 * pi + L0 * (1 - pi))
        )
      
      interp_data <- split(interp_data, interp_data$h)
      
      interp_data <- lapply(interp_data, function(x) with(x,
                                                        interp::interp(
                                                          x = L0,
                                                          y = L1,
                                                          z = posterior,
                                                          duplicate = "mean",
                                                          nx = 100, 
                                                          ny = 100
                                                        )))
      
      grid_df <- list()
      for(i in seq_along(interp_data)){
        grid_df[[i]] <- expand.grid(sp = interp_data[[i]]$x,
                                    se = interp_data[[i]]$y)
        grid_df[[i]]$posterior <- as.vector(interp_data[[i]]$z)
        grid_df[[i]]$posterior_odds <- grid_df[[i]]$posterior / (1 - grid_df[[i]]$posterior)
        grid_df[[i]]$posterior_log_odds <- log(grid_df[[i]]$posterior / (1 - grid_df[[i]]$posterior))
        grid_df[[i]]$h <- names(interp_data[i])
        grid_df[[i]]$f_lab <- paste0(names(interp_data[i]), " hours")
      }
      
      l_lots <- c(-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1)
      u_lots <- c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2)
      
      grid_df <- do.call("rbind", grid_df)
      grid_df$h <- factor(grid_df$h, levels = c(1.6, 6, 12, 16.23))
      grid_df$f_lab <- factor(grid_df$f_lab, levels = c("1.6 hours",
                                                        "6 hours",
                                                        "12 hours",
                                                        "16.23 hours"))
      
      png("Fig4c_posterior_contours.png", width = 1100, height = 800, type = "cairo")
      ggplot(grid_df,
             aes(x = sp, y = se, z = posterior_log_odds)) +
        geom_contour_filled() +
        labs(
          x = expression(Specificity ~~~~~ italic(f)(h ~ "|" ~ theta["\u00ACFCWC"])),
          y = expression(Sensitivity ~~~~~ italic(f)(h ~ "|" ~ theta[FCWC])),
          fill = expression(Pr(FCWC ~ "|" ~ H == h))
        ) +
        scale_fill_viridis_d(option = "turbo",
                             labels = paste0("(", signif(exp(l_lots) / (1 + exp(l_lots)), 2), 
                                             ",", signif(exp(u_lots) / (1 + exp(u_lots)), 2), "]")
                             ) +
        facet_wrap(~ f_lab, scale = "free") +
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
              plot.margin = margin(30, 30, 30, 30, "pt"),
              legend.position = "bottom",
              legend.title = element_text(margin = margin(r = 20)))
      dev.off()