library(ggplot2)
library(patchwork)
library(dplyr)
library(extraDistr)
library(datathin)

estimate_changepoints = function(x, pen.value = 10, minseglen = 10) {

  c(0, changepoint.np::cpt.np(x, minseglen = minseglen, penalty = "BIC")@cpts)

}

# testing for difference in means between (left+1):middle and (middle+1):right indices
gamma_glm = function(x, left, middle, right) {

  z = numeric(length(x))
  z[(middle + 1):right] = 1

  x = x[(left + 1):right]
  z = z[(left + 1):right]

  summary(glm(x ~ z, family = Gamma()))$coef[2, 4]

}

run_analysis = function(sx, filename, penalty = 10, minseglength = 10, stability = TRUE, zoom = NULL) {
  
  naive_p_values = c()
  naive_changepoints = estimate_changepoints(sx, pen.value = penalty, minseglen = minseglength)
  for (i in 1:(length(naive_changepoints) - 2)) {
    naive_p_values = c(naive_p_values, gamma_glm(sx, naive_changepoints[i], naive_changepoints[i + 1], naive_changepoints[i + 2]))
  }
  
  split_p_values = c()
  split_changepoints = estimate_changepoints(sx[seq(1, 2000, 2)], pen.value = penalty, minseglen = minseglength / 2)
  split_changepoints_even = split_changepoints - 1
  split_changepoints_even[1] = 0
  split_changepoints_even[length(split_changepoints_even)] = length(sx) / 2
  for (i in 1:(length(split_changepoints) - 2)) {
    split_p_values = c(split_p_values, gamma_glm(sx[seq(2, 2000, 2)], split_changepoints_even[i], split_changepoints_even[i + 1], split_changepoints_even[i + 2]))
  }
  split_changepoints = split_changepoints * 2 - 1
  split_changepoints[1] = 0
  split_changepoints[length(split_changepoints)] = length(sx)
  
  x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
  x1 = x_thinned[, 1]
  x2 = x_thinned[, 2]
  
  thinned_p_values = c()
  thinned_changepoints = estimate_changepoints(x1, pen.value = penalty, minseglen = minseglength)
  for (i in 1:(length(thinned_changepoints) - 2)) {
    thinned_p_values = c(thinned_p_values, gamma_glm(x2, thinned_changepoints[i], thinned_changepoints[i + 1], thinned_changepoints[i + 2]))
  }
  
  naive_p_values <- p.adjust(naive_p_values, method="bonferroni")
  split_p_values <- p.adjust(split_p_values, method="bonferroni")
  thinned_p_values <- p.adjust(thinned_p_values, method="bonferroni")
  
  for (xlim in unique(list(NULL, zoom))) {
    
    if (!is.null(xlim)) {
      stability = FALSE
      filename = paste0(filename, "_zoomed")
    } else {
      xlim = c(0, length(x))
    }
    
    nrow = ifelse(stability, 5, 3)
    
    pdf(paste0(filename, ".pdf"), height = 1.5 * nrow, width = 8)
    par(mfcol = c(nrow, 1), mar = rep(2.5, 4))
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Naive (Algorithm 2)":.(length(naive_changepoints) - 2)~"changepoints,"~.(sum(naive_p_values < 0.05))~"with p-value < 0.05 /"~.(length(naive_changepoints) - 2)))
    for (i in naive_changepoints[2:(length(naive_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "red", lwd = 1)
    }
    for (i in naive_changepoints[2:(length(naive_changepoints) - 1)][which(naive_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "red", lwd = 3, pch = "*", cex = 2)
    }
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Order-preserved sample splitting (Algorithm 3)":.(length(split_changepoints) - 2)~"changepoints,"~.(sum(split_p_values < 0.05))~"with p-value < 0.05 /"~.(length(split_changepoints) - 2)))
    for (i in split_changepoints[2:(length(split_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "purple", lwd = 1)
    }
    for (i in split_changepoints[2:(length(split_changepoints) - 1)][which(split_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "purple", lwd = 3, pch = "*", cex = 2)
    }
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Generalized data thinning (Algorithm 4)":.(length(thinned_changepoints) - 2)~"changepoints,"~.(sum(thinned_p_values < 0.05))~"with p-value < 0.05 /"~.(length(thinned_changepoints) - 2)))
    for (i in thinned_changepoints[2:(length(thinned_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "blue", lwd = 1)
    }
    for (i in thinned_changepoints[2:(length(thinned_changepoints) - 1)][which(thinned_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "blue", lwd = 3, pch = "*", cex = 2)
    }
    
    repeated_thinned_changepoints = list()
    repeated_rejected_thinned_changepoints = list()
    for (rep in 1:100) {
      x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
      x1 = x_thinned[, 1]
      x2 = x_thinned[, 2]
      cpts = estimate_changepoints(x1, pen.value = penalty, minseglen = minseglength)
      repeated_thinned_changepoints = c(repeated_thinned_changepoints, list(cpts[2:(length(cpts) - 1)]))
      pvals = c()
      for (i in 1:(length(cpts) - 2)) {
        pvals = c(pvals, gamma_glm(x2, cpts[i], cpts[i + 1], cpts[i + 2]))
      }
      pvals <- p.adjust(pvals, method="bonferroni")
      repeated_rejected_thinned_changepoints = c(repeated_rejected_thinned_changepoints, list(cpts[2:(length(cpts) - 1)][which(pvals < 0.05)]))
    }
    if (stability) {
      hist(unlist(repeated_thinned_changepoints), breaks = length(x) / 10, 
           main = bquote("Percentage of data thinning replicates with changepoint"~""), xlab = NULL, freq = T, xlim = xlim)
      for (i in 1:5) {
        abline(v = i * 365, lty = "dashed")
      }
      hist(unlist(repeated_rejected_thinned_changepoints), breaks = length(x) / 10, 
           main = bquote("Percentage of data thinning replicates with changepoint p-value < 0.05 / (number of detected changepoints)"~""), xlab = NULL, freq = T, xlim = xlim)
      for (i in 1:5) {
        abline(v = i * 365, lty = "dashed")
      }
    }
    
    dev.off()
    
  }
  
}

### Application

set.seed(2023)
data("wind", package = "gstat")
x = diff(wind[, 11])[1:2000]
noise = rnorm(length(x), 0, 1e-12)
x = x + noise
sx = x^2

run_analysis(sx, "application")

### Type 1 error simulation

set.seed(1)
n = 2000
rep = 1000
methods = c("naive", "thinned", "split")
result = list()

for (method in methods) {

  p_values = c()

  for (r in 1:rep) {

    x = rnorm(n)
    sx = x^2

    x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
    x1 = x_thinned[, 1]
    x2 = x_thinned[, 2]

    if (method == "thinned") {
      cpt_data = x1
      test_data = x2
      minseglength = 10
    } else if (method == "naive") {
      cpt_data = sx
      test_data = sx
      minseglength = 10
    } else if (method == "split") {
      cpt_data = sx[seq(1, 2000, 2)]
      test_data = sx[seq(2, 2000, 2)]
      minseglength = 5
    }

    changepoints = estimate_changepoints(cpt_data, minseglen = minseglength)
    
    if (method == "split") {
      changepoints = changepoints - 1
      changepoints[1] = 0
      changepoints[length(changepoints)] = length(test_data)
    }
    
    if (length(changepoints) > 2) { 
      for (i in sample(1:(length(changepoints) - 2), 1)) {
        p_values = c(p_values, gamma_glm(test_data, changepoints[i], changepoints[i + 1], changepoints[i + 2]))
      }
    }

  }

  result[[method]] = p_values

}

pdf("simulation_T1E.pdf", height = 3, width = 8)
par(mfrow = c(1, 3))
plot(sort(result$naive), qunif(ppoints(length(result$naive))), 
     main = bquote("Naive"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
plot(sort(result$split), qunif(ppoints(length(result$split))), 
     main = bquote("Order-preserved sample splitting"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
plot(sort(result$thinned), qunif(ppoints(length(result$thinned))), 
     main = bquote("Generalized data thinning"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
dev.off()

### Visualization of simulation under null and alternative

for (type in c("null", "alt")) {

  if (type == "null") {
    x = rnorm(n)
  } else {
    x = c(rnorm(500, sd = 2), rnorm(1000, sd = 5), rnorm(500, sd = 1))
  }
  
  sx = x^2

  run_analysis(sx, paste0("simulation_", type), stability = FALSE)

}

