
# Packages
library(ggplot2)
library(future.apply)
plan(multiprocess, workers = 16)

# Functions ---------------------------------------------------------------

#' simulate data from lognormal dist
#'
#' @param n sample size per group
#' @param meanlog mean of lognormal dist
#' @param sdlog sd of lognormal dist
#' @param delta treatment effect on log scale
gen_data_ln <- function(n, meanlog, sdlog, delta) {
  
  tx <- rlnorm(n, meanlog + delta, sdlog)
  cc <- rlnorm(n, meanlog, sdlog)
  
  data.frame(y = c(cc, tx), 
             treatment = rep(c(0,1), each = n))
  
}

#' perform bootstrap replication
#' 
#' @param d a data.frame
#' @param func summary function, e.g. mean or median
do_boot <- function(d, func = mean) {
  tx <- d[d$treatment == 1, ]
  cc <- d[d$treatment == 0, ]
  
  ## sample with replacement
  tx <- tx[sample(1:nrow(tx), nrow(tx), replace = TRUE), ]
  cc <- cc[sample(1:nrow(cc), nrow(cc), replace = TRUE), ]
  
  func(tx$y) - func(cc$y)
}

#' Bootstrap lognormal on data scale
#' 
#' @param d a data.frame
do_boot_ln <- function(d) {
  tx <- d[d$treatment == 1, ]
  cc <- d[d$treatment == 0, ]
  
  ## sample with replacement
  tx <- tx[sample(1:nrow(tx), nrow(tx), replace = TRUE), ]
  cc <- cc[sample(1:nrow(cc), nrow(cc), replace = TRUE), ]
  d_b <- rbind(tx, cc)
  fit <- lm(log(y) ~ treatment, data = d_b)
  
  sdlog <- sigma(fit)
  b <- coef(fit)
  
  # backtransform tx diff
  exp(b[1] + b[2] + sdlog^2/2) - exp(b[1] + sdlog^2/2) 
  
}

#' Simulate linear regression
#' 
#' @param i iteration
#' @param n sample size per group
#' @param mulog see gen_data_ln()
#' @param sdlog ...
#' @param delta ...
do_sim <- function(i, n, mulog, sdlog, delta) {
  d <- gen_data_ln(n, mulog, sdlog, delta)
  fit <- lm(y ~ treatment, data = d)
  
  CI <- confint(fit)
  
  data.frame(delta = coef(fit)["treatment"],
             CI_lwr = CI["treatment", 1],
             CI_upr = CI["treatment", 2])
  
}

#' Simulate with log-transformed outcome
do_sim_ln <- function(i, n, mulog, sdlog, delta) {
  d <- gen_data_ln(n, mulog, sdlog, delta)
  fit <- lm(log(y) ~ treatment, data = d)
  
  b <- coef(fit)
  CI <- confint(fit)
  
  delta <- coef(fit)["treatment"]

  data.frame(delta = exp(delta),
             CI_lwr = exp(CI["treatment", 1]),
             CI_upr = exp(CI["treatment", 2]))
  
}

#' Simulate backtransformed log-outcome using bootstrap
#'
#' @param R the number of bootstrap replications
do_sim_boot_ln <- function(i, n, mulog, sdlog, delta, R = 500) {
  d <- gen_data_ln(n, mulog, sdlog, delta)

  b <- replicate(n = R, 
                 do_boot_ln(d))
  CI <- quantile(b, c(0.025, 0.975))
  
  data.frame(delta = mean(b), 
             CI_lwr = CI[1],
             CI_upr = CI[2])
  
}

#' Simulate bootstrap of difference
do_sim_boot <- function(i, n, mulog, sdlog, delta, func = mean,
                        R = 500) {
  
  d <- gen_data_ln(n, 
                   mulog, 
                   sdlog, 
                   delta)

  b <- replicate(n = R, 
                 do_boot(d, func = func))
  CI <- quantile(b, c(0.025, 0.975))
  
  data.frame(delta = mean(b), 
             CI_lwr = CI[1],
             CI_upr = CI[2])
  
}

#' Summarise simulation results
#' 
#' @param x a data.farme with the results from the simulation
#' @param theta the true parameter value
summarize_sim <- function(x, theta) {
  est <- mean(x$delta)
  data.frame("est" = est,
             "theta" = theta,
             "rel_bias" = (est - theta)/theta,
             "CI_coverage" = mean(with(x, CI_lwr < theta & CI_upr > theta)),
             "power" = mean(with(x, sign(CI_lwr) == sign(CI_upr))) 
             )

}



# Parameters --------------------------------------------------------------
mulog <- 5
sdlog <- 1.5
delta <- log(2)
nsim <- 5000
R <- 5000 # bootstrap resamples


# take a big sample
d <- gen_data_ln(1e5, mulog, sdlog, delta)

# plot
ggplot(d, aes(y, color = treatment, group = treatment)) + 
  geom_density() +
  scale_x_continuous(labels = scales::comma)

# plot on log scale
ggplot(d, aes(y, color = treatment, group = treatment)) + 
  geom_density() +
  scale_x_log10()

# Tx effects
mean(d[d$treatment == 1, "y"])/mean(d[d$treatment == 0, "y"])
mean(d[d$treatment == 1, "y"]) - mean(d[d$treatment == 0, "y"])
median(d[d$treatment == 1, "y"])/median(d[d$treatment == 0, "y"])
mean(d[d$treatment == 1, "y"]  > mean(d[d$treatment == 1, "y"]))

# fit OLS
d <- gen_data_ln(5e2, mulog, sdlog, delta)

fit <- lm(y ~ treatment, data = d)
plot(fit)

# boostrap treatment effect
boot_samples <- replicate(n = 5000, do_boot(d))
quantile(boot_samples, c(0.025, 0.975))

# compare with OLS
confint(fit)


# Simulations -------------------------------------------------------------
#' Here we'll some diffeent simulation to compare the bootstrap to OLS

# simulate OLS
sim_res_lm <- future_lapply(1:nsim, 
                         do_sim,
                  n = 50, 
                  mulog = mulog,
                  sdlog = sdlog, 
                  delta = delta)
sim_res_lm <- do.call(rbind, sim_res_lm)

# summary
# theta is the true treatment effect
theta <- exp(mulog + delta + sdlog^2/2) - exp(mulog + sdlog^2/2) 
res_lm <- summarize_sim(sim_res_lm, theta)


# Sim bootstrap -----------------------------------------------------------
sim_res_boot <- future_lapply(1:nsim, 
                       do_sim_boot,
                  n = 50, 
                  mulog = mulog,
                  sdlog = sdlog, 
                  delta = delta,
                  R = R)
sim_res_boot <- do.call(rbind, sim_res_boot)

res_boot <- summarize_sim(sim_res_boot, theta)

# Bootstrap median
sim_res_boot_median <- future_lapply(1:nsim, 
                              do_sim_boot,
                              n = 50, 
                              mulog = mulog,
                              sdlog = sdlog, 
                              delta = delta,
                              func = median,
                              R = R)
sim_res_boot_median <- do.call(rbind, sim_res_boot_median)

# diff medians
theta_med <- exp(mulog + delta) - exp(mulog) 
res_boot_median <- summarize_sim(sim_res_boot_median, theta_med)


# Sim correct lognormal model ---------------------------------------------
sim_res_ln <- future_lapply(1:nsim, 
                              do_sim_ln,
                              n = 50, 
                              mulog = mulog,
                              sdlog = sdlog, 
                              delta = delta)
sim_res_ln <- do.call(rbind, sim_res_ln)

res_ln <- summarize_sim(sim_res_ln, theta = exp(delta))

## boot backtransformed 
## runtime approx. 3 min
sim_res_ln_boot <- future_lapply(1:nsim, 
                            do_sim_boot_ln,
                            n = 50, 
                            mulog = mulog,
                            sdlog = sdlog, 
                            delta = delta,
                            R = R)
sim_res_ln_boot <- do.call(rbind, sim_res_ln_boot)

res_ln_boot <- summarize_sim(sim_res_ln_boot, theta = theta)



# Combine all results -----------------------------------------------------

res <- cbind(data.frame("model" = c("lm", 
                                    "boot",
                                    "boot median",
                                    "lognormal",
                                    "lognormal (boot, diff)")),
             rbind(res_lm, 
                   res_boot,
                   res_boot_median, 
                   res_ln,
                   res_ln_boot)
)
format(res, digits = 2)

# save
save(sim_res_lm, sim_res_boot, sim_res_boot_median, sim_res_ln, sim_res_ln_boot,
     res_lm, res_boot,res_boot_median, res_ln,res_ln_boot, file = "res.Rdata")