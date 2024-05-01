library(ggplot2)
library(patchwork)
library(extraDistr)
library(dplyr)
library(mvtnorm)
library(atathin)

set.seed(2023)

#############################################################
# Plotting functions

plot_distribution = function(x, density_fn) {

  p = ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(bins = 30, aes(y = stat(density)), fill = "gray") +
    theme_bw() +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0))
  if (!is.null(density_fn)) {
    x_range = seq(min(x), max(x), length.out = 1000)
    p = p + geom_line(data = data.frame(x = x_range, y = density_fn(x_range)), aes(x = x, y = y)) +
      scale_y_continuous(expand=c(0,0), limits = layer_scales(p)$y$range$range)
  }
  p

}

plot_heatmap = function(x1, x2, cdf1_fn, cdf2_fn = cdf1_fn) {

  ggplot(data.frame(x1 = cdf1_fn(x1), x2 = cdf2_fn(x2)), aes(x = x1, y = x2)) +
    geom_bin2d(aes(fill = stat(density)), binwidth = 0.05) +
    scale_y_continuous(expand=c(0,0), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(expand=c(0,0), breaks = c(0, 0.5, 1)) +
    theme_bw() +
    theme(legend.position = "none") +
    coord_fixed() +
    labs(y = expression({Q[theta]^{(2)}}(X^{(2)})), x = expression({Q[theta]^{(1)}}(X^{(1)})))

}

#############################################################
# Setup

n = 100000

#############################################################
# Uniform(0, theta)

theta = 3
uniform_x = theta * runif(n)

uniform_thinned = datathin(uniform_x, "scaled-uniform")[,1,] 

p1 = plot_distribution(uniform_x, function(x) dunif(x, 0, theta)) + xlab(expression(X)) + ylab("Density") + labs(subtitle = expression(Uniform(0, theta == 3)))
p2 = plot_distribution(uniform_thinned[, 1], function(x) dbeta(x / theta, 1/2, 1) / theta) + xlab(expression(X^{(1)})) + ylab(NULL) + labs(subtitle = bquote(paste(theta, "Beta", "(1/2, 1)")))
p3 = plot_distribution(uniform_thinned[, 2], function(x) dbeta(x / theta, 1/2, 1) / theta) + xlab(expression(X^{(2)})) + ylab(NULL) + labs(subtitle = bquote(paste(theta, "Beta", "(1/2, 1)")))
p4 = plot_heatmap(uniform_thinned[, 1], uniform_thinned[, 2], function(x) pbeta(x / theta, 1/2, 1))
uniform_plots = list(p1, p2, p3, p4)

#############################################################
# Gamma(2, 6^-4) Weibull decomposition

alpha = 2
theta = 6
nu = 4
gamma_1_x = rgamma(n, shape = alpha, rate = theta^(-nu))

gamma_1_thinned = datathin(gamma_1_x, family="gamma-weibull", arg=4)[,1,]

p1 = plot_distribution(gamma_1_x, function(x) dgamma(x, alpha, theta^(-nu))) + ylab("Density") + labs(subtitle = bquote(paste("Gamma", "(", alpha == 2, ", ", theta^{-nu} == 6^{-4}, ")")))
p2 = plot_distribution(gamma_1_thinned[, 1], function(x) dweibull(x, shape=nu, scale=theta)) + labs(subtitle = expression(Weibull(theta, nu)))
p3 = plot_distribution(gamma_1_thinned[, 2], function(x) dweibull(x, shape=nu, scale=theta)) + labs(subtitle = expression(Weibull(theta, nu)))
p4 = plot_heatmap(gamma_1_thinned[, 1], gamma_1_thinned[, 2], function(x) pweibull(x, shape=nu, scale=theta))
gamma_1_plots = list(p1, p2, p3, p4)

#############################################################
# Gamma(1, 3) Normal decomposition

alpha = 1
theta = 3
gamma_2_x = rgamma(n, shape = alpha, rate = theta)

gamma_2_thinned = datathin(gamma_2_x, family="chi-squared")[,1,]

p1 = plot_distribution(gamma_2_x, function(x) dgamma(x, alpha, theta)) + ylab("Density") + labs(subtitle = bquote(paste("Gamma", "(", alpha == 1, ", ", theta == 3, ")")))
p2 = plot_distribution(gamma_2_thinned[, 1], function(x) dnorm(x, 0, sqrt(1 / (2 * theta)))) + labs(subtitle = expression(N(0, (2*theta)^{-1})))
p3 = plot_distribution(gamma_2_thinned[, 2], function(x) dnorm(x, 0, sqrt(1 / (2 * theta)))) + labs(subtitle = expression(N(0, (2*theta)^{-1})))
p4 = plot_heatmap(gamma_2_thinned[, 1], gamma_2_thinned[, 2], function(x) pnorm(x, 0, sqrt(1 / (2 * theta))), function(x) pnorm(x, 0, sqrt(1 / (2 * theta))))
gamma_2_plots = list(p1, p2, p3, p4)

#############################################################
# Gamma(7, 3) Gamma decomposition

alpha = 7
theta = 3
gamma_3_x = rgamma(n, shape = alpha, rate = theta)

gamma_3_thinned = datathin(gamma_3_x, family="gamma", arg=7)[,1,]

p1 = plot_distribution(gamma_3_x, function(x) dgamma(x, alpha, theta)) + ylab("Density") + labs(subtitle = bquote(paste("Gamma", "(", alpha == 7, ", ", theta == 3, ")")))
p2 = plot_distribution(gamma_3_thinned[, 1], function(x) dgamma(x, alpha/2, theta)) + labs(subtitle = bquote(paste("Gamma", "(", alpha / 2, ", ", theta, ")")))
p3 = plot_distribution(gamma_3_thinned[, 2], function(x) dgamma(x, alpha/2, theta)) + labs(subtitle = bquote(paste("Gamma", "(", alpha / 2, ", ", theta, ")")))
p4 = plot_heatmap(gamma_3_thinned[, 1], gamma_3_thinned[, 2], function(x) pgamma(x, alpha/2, theta), function(x) pgamma(x, alpha/2, theta))
gamma_3_plots = list(p1, p2, p3, p4)

#############################################################
# Beta(11, 6)

theta = 11
beta = 6
beta_x = rbeta(n, theta, beta)

# Conditional distribution of the beta decomposition up to the normalising constant
#  x is the vector of k-1 folds being sampled
#  z is the observed sufficient statistic
#  k is the number of folds
#  b is the known value of beta
betaGt <- function(x, z, k, b) {
  if (prod(x) <= z^k || sum(x >= 1) > 0) {
    return(0)
  } else {
    return( (prod(x^(((c(1:(k-1))-k)/k)-1))) * (prod(1-x) * (1-((z^k)/prod(x))))^(b/k - 1) )
  }
}

# Sampler for the beta decomposition
betadecomp <- function(z, k, b, M=100, B=100) {
  n <- length(z)
  samps <- matrix(0, nrow=M+B, ncol=k-1)
  
  x <- matrix(0, nrow=n, ncol=k-1)
  
  for (j in 1:n) {
    samps[1,] <- runif(k-1, min=z[j], max=1)
    for (i in 2:(M+B)) {
      prop <- runif(k-1, min=z[j]^k, max=1)
      U <- runif(1)
      alpha <- betaGt(prop, z[j], k, b)/betaGt(samps[i-1,], z[j], k, b)
      
      if (alpha > U) {
        samps[i,] <- prop
      } else {
        samps[i,] <- samps[i-1,]
      }
    }
    x[j,] <- samps[M+B,]
  }
  
  results <- data.frame(X=z, as.data.frame(x)) %>%
    rowwise() %>% 
    mutate("V{k}" := X^k/prod(c_across(starts_with("V")))) %>%
    ungroup()
  
  return(results)
}

thin_beta = function(x, beta) {

  z = betadecomp(x, 2, beta)[, 2:3]
  z = as.matrix(z)
  colnames(z) = paste0("x", 1:ncol(z))

  return(z)

}

beta_thinned = thin_beta(beta_x, beta)

p1 = plot_distribution(beta_x, function(x) dbeta(x, theta, beta)) + ylab("Density") + labs(subtitle = bquote(paste("Beta", "(", theta == 11, ", ", beta == 6, ")")))
p2 = plot_distribution(beta_thinned[, 1], function(x) dbeta(x, theta/2, beta/2)) + labs(subtitle = bquote(paste("Beta", "(", theta/2, ", ", beta/2, ")")))
p3 = plot_distribution(beta_thinned[, 2], function(x) dbeta(x, theta/2 + 1/2, beta/2)) +  labs(subtitle = bquote(paste("Beta", "(", theta/2 + 1/2, ", ", beta/2, ")")))
p4 = plot_heatmap(beta_thinned[, 1], beta_thinned[, 2], function(x)  pbeta(x, theta/2, beta/2), function(x) pbeta(x, theta/2 + 1/2, beta/2))
beta_plots = list(p1, p2, p3, p4)

#############################################################
# Plot (Figure S1 in the supplement)

wrap_plots(c(gamma_3_plots, gamma_2_plots, gamma_1_plots, beta_plots, uniform_plots), ncol = 4, byrow = TRUE)
ggsave("numerical_example_figure.pdf", height = 8, width = 8)
