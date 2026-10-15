# Why Distributions Are Everywhere

A probability distribution is a model for uncertainty: it describes which values a random quantity can take and how likely each value is.

In practical terms, distributions solve a recurring problem:

> If a process generated data like this, what would typical variation—and unusual variation—look like?

Different data-generating mechanisms lead to different distributions.

| Situation | Random variable | Typical distribution | Why |
| --- | --- | --- | --- |
| One coin flip or pass/fail outcome | 0 or 1 | Bernoulli | One binary trial |
| Number of successes in \(n\) trials | Count from 0 to \(n\) | Binomial | A fixed number of independent binary trials |
| Number of arrivals in a time window | 0, 1, 2, … | Poisson | Counts of events occurring at an average rate |
| Waiting time until an arrival | Positive continuous time | Exponential | A memoryless arrival process |
| Measurement with many small effects | Continuous value | Normal | Aggregate of many small influences |
| Waiting time for \(r\) events | Positive continuous time | Gamma | Sum of exponential waiting times |
| Uncertain probability itself | Value in \([0,1]\) | Beta | Flexible model for proportions and probabilities |

These are not unrelated formulas. They are a compact vocabulary for recurring patterns.

The goal is not to memorize 25 distributions. The goal is to recognize the process being assumed and understand what that assumption lets you calculate.

## Operational Examples

- A request succeeds or fails: **Bernoulli**.
- The number of failures in 1,000 requests: **Binomial**.
- Incidents arriving during an on-call shift: often modeled with a **Poisson**-type count process.
- The time between independent requests or failures: an **Exponential**-type waiting-time model.
- Aggregated latency or a sample mean latency: often approximately **Normal** after an appropriate transformation or sufficient averaging.
- An unknown error rate for a new service: **Beta** distributions arise naturally in Bayesian modeling.


| Layer | Question to answer | Example: **Binomial** |
| --- | --- | --- |
| Story/process | What real mechanism could generate this? | n independent requests, each succeeds/fails                                                                                           |
| Math/model | What are the assumptions and key properties? | $$X∼Binomial(n,p), \\ E[X]=np, \\ Var(X)=np(1−p)$$ |
| R/experiment | How do I calculate and simulate it? | dbinom, pbinom, qbinom, rbinom |

## A Concrete R Example

Suppose a service has a 2% error probability per request and you inspect 200 requests. Let \(X\) be the number of failed requests:

$$X ~ Binomial(n = 200, p = 0.02)$$

The model says:

$$E[X] = np = 4$$

So, four failures is the long-run expected count.

Its variance is:


$$Var(X) = np(1-p) = 3.92$$


Use R to calculate probabilities and simulate monitoring windows:

```r
n <- 200
p <- 0.02

# Theoretical expectation and variance
mean_x <- n * p
var_x <- n * p * (1 - p)

# Exact probability of exactly four failures
dbinom(4, size = n, prob = p)

# Probability of four or fewer failures
pbinom(4, size = n, prob = p)

# Simulate 10,000 monitoring windows
set.seed(42)
x <- rbinom(10000, size = n, prob = p)

mean(x) # Should be near 4
var(x)  # Should be near 3.92

hist(
  x,
  breaks = seq(-0.5, max(x) + 0.5, by = 1),
  main = "Simulated failure counts",
  xlab = "Failures per 200 requests"
)
```

The formulas are not merely dense algebra:

- `E[X]` predicts the average failure count.
- `Var(X)` describes how much ordinary monitoring windows fluctuate.
- `pbinom()` quantifies how surprising an observed count is under the assumed baseline.
- `rbinom()` lets you simulate the system and see the distribution emerge.

## What to Learn First

Do not try to learn every named distribution at equal depth. Prioritize the relationships and concepts that recur throughout statistics:

1. Random variables, PMFs/PDFs, CDFs, expectation, variance, and covariance.
2. Bernoulli and Binomial distributions.
3. Poisson and Exponential distributions, including the count/waiting-time connection.
4. Uniform and Normal distributions.
5. Gamma and Beta distributions.
6. Sampling distributions and the Central Limit Theorem.
7. Likelihood: how observed data updates which parameter values look plausible.
8. Conditional probability and conditional distributions.
9. Joint distributions, covariance, and correlation.

## A Distribution Family Tree

A useful set of relationships to internalize:

- Bernoulli trials sum to a Binomial count.
- Rare-event Binomial counts can often be approximated by Poisson counts.
- Interarrival times in a Poisson process are Exponential.
- The time until several Poisson events is Gamma.
- A Beta prior paired with Binomial data produces a Beta posterior in Bayesian analysis.
- Under broad conditions, sums and averages become approximately Normal through the Central Limit Theorem.