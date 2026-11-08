# STAT 544: Bayesian Inference executable checks
# Generated from the R blocks in Bayesian_Inference.md, in document order.
# Base R only. Run: Rscript --vanilla Bayesian_Inference_checks.R

# Document R block 1
options(digits = 10)
checks_run <- 0L
check <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(
    as.numeric(actual), as.numeric(expected), tolerance = tolerance
  )))
  checks_run <<- checks_run + 1L
}
quad <- function(f, lower, upper) {
  integrate(f, lower, upper, rel.tol = 1e-10, abs.tol = 1e-10)$value
}
beta_mean <- function(A, B) A / (A + B)
beta_var <- function(A, B) A * B / ((A + B)^2 * (A + B + 1))
beta_moment <- function(r, A, B) {
  stopifnot(A > 0, B > 0, A + r > 0)
  exp(lbeta(A + r, B) - lbeta(A, B))
}
bb_pmf <- function(k, size, A, B) {
  # Domain validation: these examples use a nonnegative integer size.
  stopifnot(length(size) == 1L, is.finite(size), size >= 0,
            size == floor(size), A > 0, B > 0,
            all(is.finite(k)), all(k == floor(k)))
  ans <- numeric(length(k))
  valid <- k >= 0 & k <= size
  ans[valid] <- exp(
    lchoose(size, k[valid]) +
      lbeta(A + k[valid], B + size - k[valid]) - lbeta(A, B)
  )
  ans
}

# A: manuscript authorship.
prior_author <- c(H = 3/4, M = 1/4)
lik_author <- dpois(6, lambda = c(H = 3, M = 6))
joint_author <- prior_author * lik_author
evidence_author <- sum(joint_author)
post_author <- joint_author / evidence_author
bf_HM <- exp(3) / 64
check(sum(post_author), 1)
check(evidence_author, 0.077962840303841, tolerance = 1e-10)
check(post_author["H"], 0.484936860875, tolerance = 1e-10)
check(lik_author["H"] / lik_author["M"], bf_HM)
check(post_author["H"] / post_author["M"], 3 * bf_HM)
check((3 * bf_HM) / (1 + 3 * bf_HM), post_author["H"])

# B: uniform prior, six successes in eight trials.
A <- 7
B <- 3
evidence_binom <- quad(function(t) dbinom(6, 8, t), 0, 1)
check(evidence_binom, 1/9)
check(exp(-lbeta(A, B)), 252)
check(quad(function(t) dbeta(t, A, B), 0, 1), 1)
check(quad(function(t) 252 * t^6 * (1 - t)^2, 0, 0.5), 23/256)
check(pbeta(0.5, A, B), 23/256)
check(quad(function(t) t * dbeta(t, A, B), 0, 1), 7/10)
check(quad(function(t) t^2 * dbeta(t, A, B), 0, 1), 28/55)
check(quad(function(t) (t - 0.7)^2 * dbeta(t, A, B), 0, 1), 21/1100)
ci_beta <- qbeta(c(0.025, 0.975), A, B)
check(diff(pbeta(ci_beta, A, B)), 0.95)
check(ci_beta, c(0.399906426284, 0.925145368580), tolerance = 1e-10)
theta_test <- seq(0.05, 0.95, length.out = 19)
check(dbinom(6, 8, theta_test) / evidence_binom,
      dbeta(theta_test, A, B))

# Document R block 2
# C: informative prior and sequential data.
alpha <- 2
beta <- 8
n <- 10
x <- 6
A2 <- alpha + x
B2 <- beta + n - x
check(c(A2, B2), c(8, 12))
check(beta_mean(A2, B2), 0.4)
check(beta_var(A2, B2), 2/175)
weight_prior <- (alpha + beta) / (alpha + beta + n)
check(weight_prior * alpha/(alpha + beta) +
        (1 - weight_prior) * x/n, beta_mean(A2, B2))
after_first <- c(alpha + 3, beta + 4 - 3)
after_second <- after_first + c(3, 6 - 3)
check(after_second, c(A2, B2))

sensitivity <- data.frame(alpha = c(1, 2, 20), beta = c(1, 8, 80))
sensitivity$A <- sensitivity$alpha + x
sensitivity$B <- sensitivity$beta + n - x
sensitivity$posterior_mean <- with(sensitivity, A/(A+B))
sensitivity$prob_above_half <- with(
  sensitivity, pbeta(0.5, A, B, lower.tail = FALSE)
)
check(sensitivity$posterior_mean, c(7/12, 2/5, 13/55))
print(sensitivity)

# D: loss minimizers under Beta(7,3).
A <- 7
B <- 3
mu <- beta_mean(A, B)
v <- beta_var(A, B)
sq_risk <- function(d) {
  quad(function(t) (t - d)^2 * dbeta(t, A, B), 0, 1)
}
abs_risk <- function(d) {
  # Split at the kink to make numerical integration more accurate.
  quad(function(t) (d - t) * dbeta(t, A, B), 0, d) +
    quad(function(t) (t - d) * dbeta(t, A, B), d, 1)
}
check(sq_risk(mu), v)
check(sq_risk(0.75) - sq_risk(mu), 0.0025)
check(optimize(sq_risk, c(0, 1), tol = 1e-10)$minimum,
      mu, tolerance = 1e-7)
posterior_median <- qbeta(0.5, A, B)
check(optimize(abs_risk, c(0, 1), tol = 1e-10)$minimum,
      posterior_median, tolerance = 1e-7)
check(beta_moment(2, A, B), 28/55)
check(beta_moment(2, A, B) - mu^2, v)
q_H <- pbeta(0.5, A, B, lower.tail = FALSE)
check(q_H, 233/256)
check(9 * (1 - q_H), 207/256)
check(1 * q_H, 233/256)
check(as.numeric(q_H > 9/10), 1)
check(as.numeric(q_H > 19/20), 0)

# E: four future trials, not a plug-in binomial.
k <- 0:4
pmf4 <- bb_pmf(k, size = 4, A = A, B = B)
check(pmf4, c(15, 70, 168, 252, 210)/715)
check(sum(pmf4), 1)
check(sum(k * pmf4), 2.8)
check(sum((k - 2.8)^2 * pmf4), 294/275)
check(pmf4[5], 42/143)
check(quad(function(t) t^4 * dbeta(t, A, B), 0, 1), pmf4[5])
check(4 * mu * (1 - mu) + 4 * 3 * v, 294/275)
check(bb_pmf(c(-1, 5), size = 4, A = A, B = B), c(0, 0))
check(bb_pmf(0, size = 0, A = A, B = B), 1)
for (kk in 0:4) {
  direct <- quad(function(t) dbinom(kk, 4, t) * dbeta(t, A, B), 0, 1)
  check(direct, pmf4[kk + 1L])
}

# Prior predictive uniformity for several sample sizes.
for (nn in c(1L, 3L, 8L, 12L)) {
  check(bb_pmf(0:nn, size = nn, A = 1, B = 1),
        rep(1/(nn + 1), nn + 1L))
}

# Document R block 3
# F: nine arrivals over three hours with Gamma(2,1) prior.
a_star <- 11
b_star <- 4
gamma_post <- function(lam) dgamma(lam, shape = a_star, rate = b_star)
check(quad(gamma_post, 0, Inf), 1)
check(quad(function(lam) lam * gamma_post(lam), 0, Inf), 11/4)
check(quad(function(lam) (lam - 11/4)^2 * gamma_post(lam), 0, Inf), 11/16)
gamma_ci <- qgamma(c(0.025, 0.975), shape = a_star, rate = b_star)
check(diff(pgamma(gamma_ci, shape = a_star, rate = b_star)), 0.95)
check(pgamma(3, shape = 11, rate = 4, lower.tail = FALSE),
      0.347229417554, tolerance = 1e-10)
u <- 2
prob_nb <- b_star / (b_star + u)
for (kk in 0:8) {
  mixed_pmf <- quad(function(lam) {
    dpois(kk, lambda = u * lam) * gamma_post(lam)
  }, 0, Inf)
  check(mixed_pmf, dnbinom(kk, size = a_star, prob = prob_nb))
}
check(dnbinom(0, size = 11, prob = prob_nb), (2/3)^11)
check(11 * (1 - prob_nb)/prob_nb, 5.5)
check(11 * (1 - prob_nb)/prob_nb^2, 8.25)
check(quad(function(lam) exp(-lam * 0.5) * gamma_post(lam), 0, Inf),
      (4/4.5)^11)
check(quad(function(w) (4/(4+w))^11, 0, Inf), 0.4)

# G: normal mean, known sampling variance.
observations <- c(2, 3, 3, 4)
m0 <- 0
tau2 <- 4
sigma2 <- 4
nn <- length(observations)
xbar <- mean(observations)
vn <- 1/(1/tau2 + nn/sigma2)
mn <- vn * (m0/tau2 + nn*xbar/sigma2)
check(c(mn, vn), c(12/5, 4/5))
check(sum((observations - 2.1)^2),
      sum((observations - xbar)^2) + nn * (xbar - 2.1)^2)
normal_kernel <- function(theta) {
  dnorm(xbar, mean = theta, sd = sqrt(sigma2/nn)) *
    dnorm(theta, mean = m0, sd = sqrt(tau2))
}
normal_evidence <- quad(normal_kernel, -Inf, Inf)
test_theta <- seq(-1, 6, length.out = 20)
check(normal_kernel(test_theta)/normal_evidence,
      dnorm(test_theta, mean = mn, sd = sqrt(vn)))
check(normal_evidence, dnorm(xbar, mean = m0, sd = sqrt(tau2 + sigma2/nn)))
normal_ci <- qnorm(c(0.025, 0.975), mn, sqrt(vn))
predictive_ci <- qnorm(c(0.025, 0.975), mn, sqrt(sigma2 + vn))
check(diff(pnorm(normal_ci, mn, sqrt(vn))), 0.95)
check(diff(pnorm(predictive_ci, mn, sqrt(sigma2 + vn))), 0.95)
for (z in c(0, 2, 5)) {
  pred_density <- quad(function(theta) {
    dnorm(z, theta, sqrt(sigma2)) * dnorm(theta, mn, sqrt(vn))
  }, -Inf, Inf)
  check(pred_density, dnorm(z, mn, sqrt(sigma2 + vn)))
}
print(list(beta_credible_interval = ci_beta,
           gamma_credible_interval = gamma_ci,
           normal_credible_interval = normal_ci,
           normal_predictive_interval = predictive_ci))

# Document R block 4
# Exercises 1-5.
post_two <- c(0.5 * 0.2, 0.5 * 0.8)
post_two <- post_two / sum(post_two)
check(post_two, c(0.2, 0.8))
check(sum(post_two * c(0.2, 0.8)), 0.68)
check(c(3 + 4, 2 + 5 - 4), c(7, 3))
check(pbeta(0.5, 1, 6, lower.tail = FALSE), 1/64)
check(c(2, 3) + c(3, 1) + c(2, 4), c(7, 8))

# Exercises 7-10.
check(qbeta(0.5, 2, 1), 1/sqrt(2))
check(qbeta(4/5, 2, 1), 2/sqrt(5))
check(c(4 * (1 - 0.6), 1 * 0.6), c(1.6, 0.6))
pmf2 <- bb_pmf(0:2, 2, 7, 3)
check(pmf2[3], 28/55)
check(sum((0:2) * pmf2), 7/5)
check(sum(((0:2) - 7/5)^2 * pmf2), 126/275)
prior_pred3 <- bb_pmf(0:3, 3, 1, 1)
check(sum((0:3) * prior_pred3), 3/2)
check(sum(((0:3) - 3/2)^2 * prior_pred3), 5/4)

# Exercises 11-15.
counts <- c(1, 3, 2)
exposure <- c(0.5, 1, 0.5)
check(c(3 + sum(counts), 2 + sum(exposure)), c(9, 4))
v1 <- 1/(1/9 + 1/9)
m1 <- v1 * (10/9 + 13/9)
check(c(m1, v1, 9 + v1), c(11.5, 4.5, 13.5))
bf01 <- dbinom(6, 8, 0.5) / (1/9)
check(bf01, 63/64)
check(bf01/(1 + bf01), 63/127)
check(quad(function(t) t * (1 - t)^2, 0, 1), exp(lbeta(2, 3)))
eta_post <- function(eta) {
  theta <- plogis(eta)
  # Logarithms remain stable even at large |eta|.
  exp(log(252) + 7 * plogis(eta, log.p = TRUE) +
        3 * plogis(eta, lower.tail = FALSE, log.p = TRUE))
}
check(quad(eta_post, -Inf, Inf), 1)
eta_mode <- optimize(eta_post, c(-5, 5), maximum = TRUE, tol = 1e-10)$maximum
check(eta_mode, log(7/3), tolerance = 1e-7)

# Exercise 16: pointwise increase, average decrease.
prior_v <- beta_var(1, 9)
success_v <- beta_var(2, 9)
failure_v <- beta_var(1, 10)
average_v <- 0.1 * success_v + 0.9 * failure_v
between_v <- 0.1 * 0.9 / 11^2
check(prior_v, 9/1100)
check(success_v, 3/242)
check(as.numeric(success_v > prior_v), 1)
check(average_v, 9/1210)
check(average_v + between_v, prior_v)

# Exercise 17: Uniform(0,theta) data, Uniform(0,2) prior.
support_normalizer <- quad(function(theta) 1/(2 * theta), 1.5, 2)
check(support_normalizer, 0.5 * log(4/3))
support_post <- function(theta) 1/(theta * log(4/3))
check(quad(support_post, 1.5, 2), 1)
check(quad(function(theta) theta * support_post(theta), 1.5, 2),
      0.5/log(4/3))

# Midpoint grid posterior for the lecture's Beta(7,3) example.
# Uniform cells: their common width cancels when normalizing weights.
cells <- 20000L
width <- 1/cells
grid <- ((1:cells) - 0.5) * width
log_weight <- dbinom(6, 8, grid, log = TRUE) +
  dbeta(grid, 1, 1, log = TRUE)
shift <- max(log_weight)
weight_raw <- exp(log_weight - shift)
weight <- weight_raw / sum(weight_raw)  # Approximate cell probabilities.
grid_density <- weight / width        # Approximate density heights.
grid_evidence <- exp(shift) * sum(weight_raw) * width
check(sum(weight), 1)
check(sum(grid * weight), 0.7, tolerance = 1e-7)
check(sum(weight[grid <= 0.5]), 23/256, tolerance = 1e-7)
check(grid_evidence, 1/9, tolerance = 1e-7)
check(sum(abs(grid_density - dbeta(grid, 7, 3))) * width,
      0, tolerance = 1e-7)

cat("Passed", checks_run, "deterministic checks.\n")
cat("R version:", as.character(getRversion()), "\n")

# Document R block 5
set.seed(544)
draws <- 200000L
report_mean <- function(z, target) {
  c(estimate = mean(z), target = target,
    Monte_Carlo_SE = sd(z)/sqrt(length(z)))
}

# Posterior uncertainty about theta.
theta_draw <- rbeta(draws, shape1 = 7, shape2 = 3)
print(report_mean(theta_draw, 0.7))
print(report_mean(as.numeric(theta_draw <= 0.5), 23/256))
print(report_mean(theta_draw^2, 28/55))

# One theta per future four-trial batch, shared across its trials.
future_count <- rbinom(draws, size = 4, prob = theta_draw)
print(report_mean(future_count, 2.8))
print(report_mean(as.numeric(future_count == 4), 42/143))
print(c(empirical_variance = var(future_count), target = 294/275))

# Two future trials share each posterior theta draw.
future_y1 <- rbinom(draws, size = 1, prob = theta_draw)
future_y2 <- rbinom(draws, size = 1, prob = theta_draw)
print(c(empirical_covariance = cov(future_y1, future_y2),
        target = 21/1100))

# Gamma-Poisson predictive mixture.
lambda_draw <- rgamma(draws, shape = 11, rate = 4)
future_poisson <- rpois(draws, lambda = 2 * lambda_draw)
print(report_mean(future_poisson, 5.5))
print(report_mean(as.numeric(future_poisson == 0), (2/3)^11))
print(c(empirical_variance = var(future_poisson), target = 8.25))

# Normal predictive mixture.
mean_draw <- rnorm(draws, mean = 2.4, sd = sqrt(0.8))
new_measurement <- rnorm(draws, mean = mean_draw, sd = sqrt(4))
print(report_mean(new_measurement, 2.4))
print(c(empirical_variance = var(new_measurement), target = 4.8))

# Optional plots: set TRUE to display in an interactive graphics device.
make_plots <- FALSE
if (make_plots) {
  old_par <- par(mfrow = c(1, 2))
  curve(dbeta(x, 1, 1), from = 0, to = 1, ylim = c(0, 3),
        xlab = "Success probability", ylab = "Density",
        col = "gray45", lty = 2, main = "Prior and posterior")
  curve(dbeta(x, 7, 3), add = TRUE, col = "navy", lwd = 2)
  legend("topleft", c("Prior Beta(1,1)", "Posterior Beta(7,3)"),
         col = c("gray45", "navy"), lty = c(2, 1), bty = "n")
  barplot(rbind(bb_pmf(0:4, 4, 7, 3), dbinom(0:4, 4, 0.7)),
          beside = TRUE, names.arg = 0:4, col = c("navy", "gray70"),
          xlab = "Future successes", ylab = "Probability",
          main = "Predictive versus plug-in")
  legend("topright", c("Beta-binomial", "Plug-in binomial"),
         fill = c("navy", "gray70"), bty = "n")
  par(old_par)
}

