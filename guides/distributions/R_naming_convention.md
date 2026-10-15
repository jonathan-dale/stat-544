# R's standard naming convention

### For reference run:
```R
?Distributions
```

R’s standard “dpqr” naming convention for most probability distributions in base R’s stats package:

| Prefix | Meaning                          | Continuous distribution | Discrete distribution                                     |
| ------ | -------------------------------- | ----------------------- | --------------------------------------------------------- |
| d      | Density / mass function          | PDF f(x)f(x)f(x)        | PMF P(X=x)P(X=x)P(X=x)                                    |
| p      | Cumulative distribution function | P(X≤x)P(X \\le x)P(X≤x) | P(X≤x)P(X \\le x)P(X≤x)                                   |
| q      | Quantile function                | Inverse CDF             | Smallest integer xxx with P(X≤x)≥pP(X \\le x)\\ge pP(X≤x)≥p |
| r      | Random generation                | Simulated draws         | Simulated draws                                           |


R's distribution documentation explicitly describes the family as `dxxx`, `pxxx`, `qxxx`, and `rxxx` for density/mass, CDF, quantile, and random-variate generation, respectively.  

## Examples
### Normal Distribution: 
```R
dnorm(0)                  # density f(0) for N(0, 1)
pnorm(1.96)               # P(Z <= 1.96)
qnorm(0.975)              # 97.5th percentile of N(0, 1)
rnorm(10, mean = 5, sd = 2)  # 10 simulated N(5, 2^2) values
```

### Binomial distribution
```R
dbinom(3, size = 10, prob = 0.4)  # P(X = 3)
pbinom(3, size = 10, prob = 0.4)  # P(X <= 3)
qbinom(0.95, size = 10, prob = 0.4) # 95th percentile
rbinom(10, size = 10, prob = 0.4) # 10 simulated binomial counts
```

### Poisson distribution
```R
dpois(4, lambda = 3)      # P(X = 4)
ppois(4, lambda = 3)      # P(X <= 4)
qpois(0.90, lambda = 3)   # 90th percentile
rpois(10, lambda = 3)     # 10 simulated Poisson values
```

One Important distinction: density vs. probability:
The terminology is slightly different for continuous versus discrete random variables:

For a continuous distribution such as normal, gamma, beta, or exponential, `d...` returns a density, not the probability of an exact value. For example, `dnorm(0)` is about `0.3989`, but $$P(X=0)=0$$ for a continuous normal random variable.

For a discrete distribution such as binomial, Poisson, geometric, or negative binomial, `d...` returns a probability mass, so `dbinom(3, ...)` really is $$P(X=3)$$

For a continuous variable, probabilities over intervals come from `p...`:
```R
# P(4 < X <= 7) where X ~ N(5, 2^2)
pnorm(7, mean = 5, sd = 2) - pnorm(4, mean = 5, sd = 2)
```
For a discrete variable:
```R
# P(2 <= X <= 5), where X ~ Poisson(3)
ppois(5, lambda = 3) - ppois(1, lambda = 3)
```

### Tail probabilities
Most `p...` functions also support two particularly useful arguments:
```R
pnorm(1.96, lower.tail = FALSE)  # P(Z > 1.96)
pnorm(1.96, log.p = TRUE)        # log P(Z <= 1.96)
```
> `lower.tail = FALSE` requests an upper-tail probability, often more numerically stable than computing `1 - pnorm(...)`.
> `log.p = TRUE` returns a log-probability, useful when probabilities are extremely small.
> The corresponding `q...` functions accept `lower.tail` and `log.p` too:
```R
qnorm(0.025, lower.tail = FALSE)  # same as qnorm(0.975)
```

