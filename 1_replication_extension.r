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
library(interp)

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
#wc_dat_vis <- wc_dat %>%
#  group_by(Ref, Participants) %>%
#  mutate(p_est_ungr = p_est,
#         n_ungr = n,
#         N_ungr = N,
#         n = sum(n),
#         N = sum(N),
#         p_est = ifelse(!is.na(nbin), weighted.mean(p_est, nbin), p_est),
#         collapsed = ifelse(!is.na(nbin), TRUE, FALSE))
#wc_dat_vis[,c("Ref", "StudyYear", "DataYear", "weight","Location", "Participants", "N", "n", "p_est", "collapsed")] %>%
#  unique() %>% 
#  flextable() %>%
#  theme_booktabs() %>%
#  autofit() %>%
#  set_caption("Table 1. Erroneous conviction rate (WC) estimation, revised data") %>%
#  save_as_docx(path = "table1.docx")

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

#fc_dat[,c("Ref", "StudyYear", "DataYear", "weight","Location", "Participants", "N", "n", "p_est")] %>%
#  unique() %>% 
#  flextable() %>%
#  theme_booktabs() %>%
#  autofit() %>%
#  set_caption("Table 2. False confession rate (WC) estimation, revised data") %>%
#  save_as_docx(path = "table2.docx")



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
    pp_check(wc_model, type = "hist", ndraws = 11)

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
  
  ## Confirmed method for constructing base rate
  n_sim <- 1000000
  
  wc_mu_draws <- posterior_epred(wc_model, ndraws = n_sim)
  fc_mu_draws <- posterior_epred(fc_model, ndraws = n_sim)
  
  w_wc <- model.frame(wc_model)$wij
  w_wc <- w_wc / sum(w_wc)
  
  w_fc  <- model.frame(fc_model)$wj
  w_fc  <- w_fc / sum(w_fc)
  
  wc_hat <- as.numeric(wc_mu_draws %*% w_wc)
  fc_hat <- as.numeric(fc_mu_draws %*% w_fc)
  
  combined_df <- tibble(
    ErrorRate     = wc_hat,
    FalseConfPrev = fc_hat,
    FCWC_BaseRate = wc_hat * fc_hat
  )
  
  fcwc_base_rate <- combined_df$FCWC_BaseRate
  
  ## Descriptives for FCWC base rate
  mean(fcwc_base_rate)
  median(fcwc_base_rate)
  quantile(fcwc_base_rate, c(0.025, 0.975))
  
  
#==============================================================================#    
  
# 3.2 FCWC probabilities given a single interrogation technique 
  
  ## Set plausible levels for interrogation
  ll <- 0.1
  ul <- 0.9
  sens_ranges <- list(
    low      = c(ll,                     ll + (((ul-ll)/3) * 1)),
    moderate = c(ll + (((ul-ll)/3) * 1), ll + (((ul-ll)/3) * 2)),
    high     = c(ll + (((ul-ll)/3) * 2), ll + (((ul-ll)/3) * 3))
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
  
  
    ## M&A original binned, discrete
    #bucket_vals <- list(low = seq(.1,.3, by = .1),
    #                    moderate = seq(.4,.6, by = .1),
    #                    high = seq(.7,.9, by = .1))
    
    #sens_lvl <- sample(names(bucket_vals), n_sim, TRUE)
    #spec_lvl <- sample(names(bucket_vals), n_sim, TRUE)
    
    #sensitivity <- map_dbl(sens_lvl, ~ sample(bucket_vals[[.x]], 1))
    #specificity <- map_dbl(spec_lvl, ~ sample(bucket_vals[[.x]], 1))
    
    #names(sensitivity) <- case_when(sensitivity < 0.33 ~ "Low",
    #                                sensitivity < 0.66 ~ "Moderate",
    #                                TRUE               ~ "High")
    #names(specificity) <- case_when(specificity < 0.33 ~ "Low",
    #                                specificity < 0.66 ~ "Moderate",
    #                                TRUE               ~ "High")
  
    
  ## Calculate the conditional probability using Bayes' Theorem
  p_fcwc_given_t <- (sensitivity * fcwc_base_rate) /
    ((sensitivity * fcwc_base_rate) +
       ((1 - specificity) * (1 - fcwc_base_rate)))
  
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
    p_fcwc_given_t_scaled <- lapply(fcwc_scaled, function(fcwc_base_rate) {
      (sensitivity * fcwc_base_rate) /
        ((sensitivity * fcwc_base_rate) +
           ((1 - specificity) * (1 - fcwc_base_rate)))
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
        lower_hdi = bayestestR::hdi(fcwc_given_t, ci = 0.95)[[2]],
        upper_hdi = bayestestR::hdi(fcwc_given_t, ci = 0.95)[[3]])
    
    ### Merge labelling data
    labs_dat <- left_join(group_medians, group_density, by = "attribution")
    labs_dat <- left_join(labs_dat, group_HDI, by = "attribution")
    
    ### Replicate figure 1
    #png("Fig1_replication_w_mary_data.png", width = 1000, height = 750, type = "cairo")
    ggplot(gdat, aes(x=fcwc_given_t, fill = attribution)) +
      geom_density(color = NA,
                   alpha = 0.8) +
      #geom_histogram(color = NA,
      #               alpha = 0.8,
      #               bins = 200) +
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
           subtitle = "Posterior Distribution of FCWC Risk Across All Tactic Scenarios\n(replication and extension of Mourtgos and Adams, 2026 to include neglected publications,\nand omit inappropriate publications)",
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
    #dev.off()
  
    
    ### Scatter plot
    png("Fig3_posterior_by_loglik.png", width = 1000, height = 450, type = "cairo")
    data.frame(prior = fcwc_base_rate,
               se = sensitivity,
               sp = specificity) %>%
      slice_sample(n = 10000) %>%
      mutate(
        sens_bin = cut(se, breaks = c(0, .33, .66, 1)),
        posterior = (se * prior) /
            (se * prior + (1 - sp) * (1 - prior)),
        lr_pos = se / (1 - sp)
      ) %>%
      ggplot(aes(x = lr_pos, y = posterior)) +
      geom_point(aes(color = se,
                     size = 1 - sp),
                 alpha = 0.3,
                 stroke = NA) +
      labs(
        x = expression(log~bgroup("(", over( "Pr(T|FCWC)", "Pr⁡(T│¬FCWC"), ")")),
        y = "Pr(FCWC|T)",
        color = "Pr(T|FCWC)",
        size = "Pr⁡(T│¬FCWC)"
      ) +
      scale_color_gradient(low = "skyblue",
                          high = "red") +
      scale_size(transform = "reverse",
                 range = c(0.1, 3)) +
      scale_x_log10() +
      guides(colour = guide_colourbar(barwidth = 1, barheight = 10)) +
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
            legend.position = "right")
    dev.off()
    
    
    ### Contour plot
    interp_data <- data.frame(prior = fcwc_base_rate,
               se = sensitivity,
               sp = 1 - specificity) %>%
      mutate(
        posterior = (se * prior) /
          (se * prior + sp * (1 - prior))
      ) %>%
      slice_sample(n = 100000) 
    
    interp_data <- with(interp_data,
      interp::interp(
        x = sp,
        y = se,
        z = posterior,
        duplicate = "mean",
        nx = 100, 
        ny = 100
        ))

    grid_df <- expand.grid(
      sp = interp_data$x,
      se = interp_data$y
    )
    
    grid_df$posterior <- as.vector(interp_data$z)
    
    png("Fig3_posterior_contours.png", width = 1000, height = 800, type = "cairo")
      ggplot(grid_df,
             aes(x = sp, y = se, z = posterior)) +
        geom_contour_filled() +
        labs(
          x = "Pr⁡(T│¬FCWC)",
          y = "Pr(T|FCWC)",
          fill = "Pr(FCWC|T)"
        ) +
        scale_fill_viridis_d(option = "magma") +
        scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                           limits = c(0.09, 0.91),
                           expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                           limits = c(0.09, 0.91),
                           expand = c(0, 0)) +
        theme_minimal() + 
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
              legend.position = "right")
     dev.off()