# 3.4.26
# (a)
st_norm_monte_carlo <- function(x, N) {
  if (x < 0) stop('x should be greate than or equal to zero.')
  X <- rnorm(N, 0, 1)
  count <- length(X[X >= x | X <= -x])
  return(count/N)
}

st_norm_monte_carlo(2, 20000)
pnorm(-2, 0, 1) + pnorm(2, 0, 1, lower.tail = FALSE)


#(b)
cont_normal_generate <- function(N, eps, sigma) {
  # Check if the epsilon is in the right range
  if (eps > 1 | eps < 0) stop("Epsilon has to be in between (0, 1)")

  I <- rbinom(N, 1, 1 - eps)
  Z <- rnorm(N)
  W <- Z * I + sigma * Z * (1 - I)
  return(W)
}

cont_normal_pdf <- function(x, eps, sigma) {
  # Check if the epsilon is in the right range
  if (eps > 1 | eps < 0) stop("Epsilon has to be in between (0, 1)")
  
  return(dnorm(x) * (1 - eps) + dnorm(x / sigma) * eps / sigma)
}

cont_normal_cdf <- function(q, eps, sigma, lower.tail = TRUE) {
  # Check if the epsilon is in the right range
  if (eps > 1 | eps < 0) stop("Epsilon has to be in between (0, 1)")

  if (lower.tail == TRUE) {
    return(pnorm(q) * (1 - eps) + pnorm(q / sigma) * eps)
  } else {
    return(1 - (pnorm(q) * (1 - eps) + pnorm(q / sigma) * eps))
  }
}

cont_normal_cdf(-2, 0.15, 10) + cont_normal_cdf(2, 0.15, 10, lower.tail = FALSE)

# sample
st norm sample <− rnorm(100)
cont norm sample <− cont normal generate(100, 0.15, 10)


# (c)
cont_normal_cdf(-2, 0.15, 20) + cont_normal_cdf(2, 0.15, 20, lower.tail = FALSE)

# (d)
cont_normal_cdf(-2, 0.25, 20) + cont_normal_cdf(2, 0.25, 20, lower.tail = FALSE)

cont_normal_cdf(-2, 0.15, 10) + cont_normal_cdf(2, 0.15, 10, lower.tail = FALSE)
x <− seq(-6, 6, 0.1)
par(mfrow = c(2, 2))
plot(x, dnorm(x), type = 'l', ylab = '', main = 'Standard Normal pdf')
plot(x, cont normal_pdf(x, 0.15, 10), type = 'l', ylab = '', main = expression(paste(epsilon, '=0.15 ', sigma, '=10')))


# 3.4.27
x <- seq(-6, 6, 0.1)
par(mfrow = c(2, 2))
plot(x, cont_normal_pdf(x, 0.15, 10), type = 'l', ylab = '', main = expression(paste(epsilon, '=0.15 ', sigma, '=10')))
plot(x, cont_normal_pdf(x, 0.15, 20), type = 'l', ylab = '', main = expression(paste(epsilon, '=0.15 ', sigma, '=20')))
plot(x, cont_normal_pdf(x, 0.25, 20), type = 'l', ylab = '', main = expression(paste(epsilon, '=0.25 ', sigma, '=20')))
plot(x, dnorm(x), type = 'n', ylab = '', main = 'Overlapped plot')
lines(x, dnorm(x), type = 'l', col = 'blue', ylab = '', main = 'Overlapped plot')
lines(x, cont_normal_pdf(x, 0.15, 10), type = 'l', col = 'red' )
legend('topright', lty = c(1, 1), col = c('blue', 'red'), legend = c('Standard', 'Contaminated'))

# 4.8.19
AR_t_generate <- function(N, r) {
  # N is the maximum number of iteration and r is the degrees of freedom
  if (r < 2) stop("r should be greater than or equal to 2")

  # Define the objective function: t-distribution / Cauchy-distribution
  objective_function <- function(x) gamma((r + 1) / 2) / (sqrt(r * pi) * gamma(r / 2)) * (1 + x^2 / r)^(-(r + 1) / 2) * (pi * (1 + x^2))

  # Find the upper bound
  M <- optimize(objective_function, c(-10, 10), maximum = TRUE)$objective
  temp <- c()

  while (length(temp) < (N)) {
    # CDF of Cauchy is F = 1/pi * arctan(x) + 1/2
    # Inverse of Cauchy is F^(-1) = tan(pi * (y - 1/2))
    # Apply probabiliy integral transform
    U1 <- runif(1)
    cauchy_rand <- tan(pi * (U1 - 1/2))

    # Accept-Reject algorithm
    U2 <- runif(1)
    if (M * U2 <= objective_function(cauchy_rand)) {
      temp <- c(temp, cauchy_rand)
    }
  }
  return(temp)
}
t_rand <- AR_t_generate(1000, 30)

# Conduct Kolmogorov-Smirnov test to check if this is actually the targeted distribution
true_t <- rt(length(t_rand), 30)
library(stats)
ks.test(t_rand, true_t)

# Plot the histogram and the true pdf to check
par(mfrow = c(1, 1))
hist(t_rand, probability = TRUE, ylab = '', main = 't-distribution with r = 30')
curve(dt(x, 30), col = 'blue', lwd = 2, add = TRUE)
legend('topright', lty = 1, col = 'blue', legend = 'true t')

# 4.8.20
beta_generate <- function(N, alpha, beta) {
  if (alpha <= 0 | beta <= 0) stop('alpha and beta have to be greater than 0')

  X <- c()
  while (length(X) < N) {
    U1 <- runif(1); U2 <- runif(1)
    V1 <- U1^(1/alpha); V2 <- U2^(1/beta)
    W <- V1 + V2

    if (W <= 1) {
      X <- c(X, V1/W)
    }
  }
  return(X)
}

beta_rand <- beta_generate(1000, 3, 2)

# Conduct Kolmogorov-Smirnov test to check
true_beta <- rbeta(length(beta_rand), 3, 2)
ks.test(beta_rand, true_beta)

# Plot the histogram and the true pdf to check
par(mfrow = c(1, 1))
hist(beta_rand, probability = TRUE, ylab = '', ylim = c(0, 2), main = expression(paste('Beta distribution ', alpha, '=3, ', beta, '=2')))
curve(dbeta(x, 3, 2), col = 'blue', lwd = 2, add = TRUE)
legend('topright', lty = 1, col = 'blue', legend = 'true beta')


# 4.8.21
MB_normal_generate <- function(N) {
  # Marsaglia and Bray algorithm
  X <- matrix(0, nrow = 1, ncol = 2)
  while (nrow(X) < N+1) {
    U <- runif(1, -1, 1); V <- runif(1, -1, 1)
    W <- U^2 + V^2
    if (W <= 1) {
      Z <- sqrt(-2 * log(W) / W)
      X <- rbind(X, c(U*Z, V*Z))
    }
  }
  return(X[2:(N+1), ])
}

AR_normal_generate <- function(N) {
  # Rejection sampling algorithm

  # Define the objective function: Normal-distribution / Cauchy-distribution
  objective_function <- function(x) 1 / sqrt(2 * pi) * exp(-x^2 / 2) * pi * (1 + x^2)
  # Find the upper bound
  M <- optimize(objective_function, c(-10, 10), maximum = TRUE)$objective
  temp <- c()
  while (length(temp) < N) {
    # CDF of Cauchy is F = 1/pi * arctan(x) + 1/2
    # Inverse of Cauchy is F^(-1) = tan(pi * (y - 1/2))
    # Apply probabiliy integral transform
    U1 <- runif(1)
    cauchy_rand <- tan(pi * (U1 - 1/2))

    # Accept-Reject algorithm
    U2 <- runif(1)
    if (M * U2 <= objective_function(cauchy_rand)) {
      temp <- c(temp, cauchy_rand)
    }
  }
  return(temp)
}

BM_normal_generate <- function(N) {
  # Box-Muller transform
  X <- matrix(0, nrow = N, ncol = 2)
  for (i in 1:N) {
    U1 <- runif(1)
    U2 <- runif(1)
    X[i, ] <- c(sqrt(-2*log(U1))*cos(2*pi*U2), sqrt(-2*log(U1))*sin(2*pi*U2))
  }
  return(X)
}

# Sample
MB_normal <- MB_normal_generate(1000)
AR_normal <- AR_normal_generate(1000)
BM_normal <- BM_normal_generate(1000)
true_normal <- rnorm(1000)

# Plot
par(mfrow = c(2, 2))
hist(MB_normal, probability = TRUE, main = 'Marsaglia and Bray')
curve(dnorm(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col = 'blue', lty = 1, legend = 'true Gaussian')

hist(AR_normal, probability = TRUE, main = 'Accept-Reject')
curve(dnorm(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col ='blue', lty = 1, legend = 'true Gaussian')

hist(BM_normal, probability = TRUE, main = 'Box-Muller')
curve(dnorm(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col = 'blue', lty = 1, legend = 'true Gaussian')

hist(true_normal, probability = TRUE, main = 'Built-in')
curve(dnorm(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col = 'blue', lty = 1, legend = 'true Gaussian')


# The algorithm of the built-in function 'pnorm' can be found here:
# https://github.com/wch/r-source/blob/776708efe6003e36f02587ad47b2eaaaa19e2f69/src/nmath/snorm.c
# There are several algorithms you can choose from which default to 'INVERSION'.
# INVERSION method seems to be referring to inverse transform but the CDF of Normal distribution
# does not have a closed inverse form..... Failed to crack the code.


# <10>
# (a)
mixture_DE_generate <- function(N) {
  X <- c()
  while (length(X) < N) {
    s <- rexp(1, 1/2)
    Z <- rnorm(1, 0, sqrt(s))
    X <- c(X, Z)
  }
  return(X)
}
mix_DE_rand <- mixture_DE_generate(1000)

# (b)
Y1 <- rexp(1000, 1)
Y2 <- rexp(1000, 1)
laplace_rand <- Y1 - Y2

# command 'require' returns TRUE/FALSE depending on whether the package is installed
# check if {VGAM} is installed and install if there isn't
if (!require(VGAM)) {
  # Get user permission
  user_response <- readline(prompt = '{VGAM} will be downloaded. Press y to continue or n to stop. (y/n)')
  if (user_response == 'y') {
    # Choose South Korea mirror, index = 62
    chooseCRANmirror(ind = 62)

    install.packages('VGAM')
    library(VGAM)
  } else if (user_response == 'n') {
    stop('No permission to proceed.')
  } else {
    stop('Wrong input. Please answer with y/n')
  }
} else {
  library(VGAM)  
}

true_laplace <- rlaplace(1000, 0, 1)

par(mfrow = c(2, 2))
hist(mix_DE_rand, probability = TRUE, ylim = c(0, 0.5), main = 'Double Exponential by Mixture')
curve(dlaplace(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col = 'blue', lty = 1, legend = 'true laplace')

hist(laplace_rand, probability = TRUE, ylim = c(0, 0.5), main = 'Dist. of the difference of two Exp')
curve(dlaplace(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col = 'blue', lty = 1, legend = 'true laplace')

hist(true_laplace, probability = TRUE, ylim = c(0, 0.5), main = 'Package Double exponential')
curve(dlaplace(x), col = 'blue', lwd = 2, add = TRUE)
legend('topright', col = 'blue', lty = 1, legend = 'true laplace')
