library(tidyverse)
library(reshape2)

# Number of simulations
N <- 1000000

# Simulate prior distributions, sensitivity, and specificity to examine
dat <- data.frame("fcwc_sim" = fcwc_base_rate,
                  "Uniform" = rbeta(N, 1, 1),
                  "Center" = rbeta(N, 4, 4),
                  "EdgeAverse" = rbeta(N, 2, 2),
                  "SkewLeft" = rbeta(N, 5, 2),
                  "SkewRight" = rbeta(N, 2, 5),

                  "Se" = runif(N, 0.1, 0.9),
                  "Sp" = runif(N, 0.1, 0.9))

# Calculate posterior distributions
## Simplified version of Bayes' Theorem
dat <- dat %>%
  mutate(Diff = Se / (1 - Sp),
         fcwc_sim_Posterior = 1 / (1 + ((1 - fcwc_sim) / fcwc_sim) * ((1 - Sp) / Se)),
         Uniform_Posterior = 1 / (1 + ((1 - Uniform) / Uniform) * ((1 - Sp) / Se)),
         Center_Posterior = 1 / (1 + ((1 - Center) / Center) * ((1 - Sp) / Se)),
         EdgeAverse_Posterior = 1 / (1 + ((1 - EdgeAverse) / EdgeAverse) * ((1 - Sp) / Se)),
         SkewLeft_Posterior = 1 / (1 + ((1 - SkewLeft) / SkewLeft) * ((1 - Sp) / Se)),
         SkewRight_Posterior = 1 / (1 + ((1 - SkewRight) / SkewRight) * ((1 - Sp) / Se)))

## Traditional specification
#dat <- dat %>%
#  mutate(Diff =  Se / (1 - Sp),
#         Uniform_Posterior = (Se * Uniform) / (Se * Uniform + (1 - Sp) * (1 - Uniform)),
#         Center_Posterior = (Se * Center) / (Se * Center + (1 - Sp) * (1 - Center)),
#         EdgeAverse_Posterior = (Se * EdgeAverse) / (Se * EdgeAverse + (1 - Sp) * (1 - EdgeAverse)),
#         SkewLeft_Posterior = (Se * SkewLeft) / (Se * SkewLeft + (1 - Sp) * (1 - SkewLeft)),
#         SkewRight_Posterior = (Se * SkewRight) / (Se * SkewRight + (1 - Sp) * (1 - SkewRight)))



# Prep data for graph
dat <- melt(dat)
dat$graph <- NA
dat$graph[grep("^Se$", dat$variable)] <- "Sensitivity"
dat$graph[grep("^Sp$", dat$variable)] <- "Specificity"
dat$graph[grep("^Diff$", dat$variable)] <- "Sens. / (1 - Spec.)"
dat$graph[grep("^fcwc_sim", dat$variable)] <- "FCWC"
dat$graph[grep("^Uniform", dat$variable)] <- "Uniform, Beta(1, 1)"
dat$graph[grep("^Center", dat$variable)] <- "Center, Beta(4, 4)"
dat$graph[grep("^EdgeAverse" ,dat$variable)] <- "Edge Averse, Beta(2, 2)"
dat$graph[grep("^SkewLeft" ,dat$variable)] <- "Left Skew, Beta(5, 2)"
dat$graph[grep("^SkewRight" ,dat$variable)] <- "Right Skew, Beta(2, 5)"

dat$col <- "Prior"
dat$col[grep("_Posterior", dat$variable)] <- "Posterior"

dat$graph <- factor(dat$graph, 
                    levels = c("Sensitivity", 
                               "Specificity", 
                               "Sens. / (1 - Spec.)",
                               "FCWC",
                               "Uniform, Beta(1, 1)", 
                               "Center, Beta(4, 4)", 
                               "Edge Averse, Beta(2, 2)", 
                               "Left Skew, Beta(5, 2)", 
                               "Right Skew, Beta(2, 5)"))

dat$col <- factor(dat$col,
                  levels = c("Prior", "Posterior"))

# Plot posterior distributions
library(ggh4x)

png("prior_posterior_comparisons.png", 1000, 500, type = "cairo")
dat[which(!dat$variable %in% c("Se", 
                               "Sp", 
                               "Diff")), ] %>%
  ggplot(aes(x = value, group = variable, fill = col)) +
    geom_density(alpha = 0.5, color = NA) +
    scale_fill_manual(values = c("black", "red")) + 
    facet_wrap(~ graph, scales = "free") +
    labs(
      title = "Posterior Distributions with Uniformly Sampled Sensitivity & Specificity",
      x = "Probability",
      y = "Density"
    ) +
    facetted_pos_scales(
      x = list(
        graph == "FCWC" ~ scale_x_continuous(limits = c(0, 0.1))
      )) +
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
          legend.position = "bottom")
dev.off()