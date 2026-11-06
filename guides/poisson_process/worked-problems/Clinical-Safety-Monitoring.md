# Clinical safety monitoring: averages, bounds, and convergence

A clinical research team monitors the proportion of patients experiencing a mild adverse event during followup. Assume that each patient independently experiences the event with probability $p = 0.20$. Let $Xᵢ = 1$ if patient $i$ experiences an event and $Xᵢ = 0$ otherwise; let $p̂ ₙ = (X₁ + ··· + Xₙ)/n$.

(a) Find $E(p̂ ₙ)$ and $Var(p̂ ₙ)$ as functions of $n$.  
(b) For $n$ = $100$, $400$, and $1000$, use Chebyshev’s inequality to give an upper bound for $P(|p̂ ₙ − 0.20| > 0.05)$.  
(c) For each of the same sample sizes, calculate the exact probability $P(|p̂ ₙ − 0.20| > 0.05)$ using the binomial distribution (R is allowed). Compare with the bounds in (b).  
(d) State what the Weak Law of Large Numbers says about $p̂ ₙ$ in this setting. What does the Strong Law of Large Numbers say?  
(e) Explain why increasing sample size improves the stability of the observed event proportion but does not guarantee that a particular finite-sample estimate will be within 0.05 of 0.20.  


## Solution
This problem asks how accurately an observed event proportion estimates the true probability $p=0.20$. The key distinction is between an **exact finite-sample probability**, a general upper bound on that probability, and a statement about what happens as the sample size tends to infinity.

## Model and shared foundation

Each $X_i$ is a Bernoulli random variable: it equals $1$ with probability $p=0.20$ and $0$ otherwise. Because the patients’ outcomes are independent and have the same event probability, the total number of patients experiencing an event satisfies

$$
S_n=\sum_{i=1}^n X_i\sim\mathrm{Binomial}(n,0.20).
$$

The observed proportion is simply that count divided by the sample size:

$$
\widehat p_n=\frac{S_n}{n}.
$$

Before collecting data, $\widehat p_n$ is random because different samples can produce different counts. After collecting data, its realized value is the proportion you actually observe.

## Mean and variance: part (a)
Find $E(p̂ ₙ)$ and $Var(p̂ ₙ)$ as functions of $n$.  

For a Bernoulli random variable,

$$
E(X_i)=1\cdot p+0\cdot(1-p)=p.
$$

Also, $X_i^2=X_i$, because both $0^2=0$ and $1^2=1$. Therefore,

$$
\mathrm{Var}(X_i)
=E(X_i^2)-[E(X_i)]^2
=p-p^2
=p(1-p).
$$

Now use two rules: expectation is linear, and the variance of a sum of independent random variables is the sum of their variances. Dividing a random variable by $n$ divides its variance by $n^2$, so

$$
E(\widehat p_n)
=\frac1n\sum_{i=1}^n E(X_i)
=\frac{np}{n}
=p,
$$

and

$$
\mathrm{Var}(\widehat p_n)
=\frac1{n^2}\sum_{i=1}^n\mathrm{Var}(X_i)
=\frac{np(1-p)}{n^2}
=\frac{p(1-p)}n.
$$

Substituting $p=0.20$,

$$
\boxed{E(\widehat p_n)=0.20,\qquad
\mathrm{Var}(\widehat p_n)=\frac{0.16}{n}.}
$$

Thus, the estimator is **unbiased**: its average over repeated samples equals the true proportion. Increasing $n$ does not change that center; it reduces the spread around it.

## Chebyshev bounds: part (b)
For $n$ = $100$, $400$, and $1000$, use Chebyshev’s inequality to give an upper bound for $P(|p̂ ₙ − 0.20| > 0.05)$.  

Chebyshev’s inequality says that, for a random variable $Y$ with finite variance and any $\varepsilon>0$,

$$
P\bigl(|Y-E(Y)|\geq\varepsilon\bigr)
\leq\frac{\mathrm{Var}(Y)}{\varepsilon^2}.
$$

In words, variance limits the probability of being far from the mean. The same upper bound applies to a strict “$>$” event because that event is contained in the corresponding “$\geq$” event.

Here, $\varepsilon=0.05$, so

$$
P\bigl(|\widehat p_n-0.20|>0.05\bigr)
\leq
\frac{0.16/n}{(0.05)^2}
=\frac{0.16}{0.0025n}
=\boxed{\frac{64}{n}}.
$$

The requested bounds are therefore:

| Sample size $n$ | Chebyshev upper bound | As a percentage |
|---:|---:|---:|
| 100 | $0.64$ | $64\%$ |
| 400 | $0.16$ | $16\%$ |
| 1000 | $0.064$ | $6.4\%$ |

These are upper bounds, not estimates of the actual probabilities. Chebyshev uses only the mean and variance, so it does not take advantage of the full binomial model.

## Exact binomial probabilities: part (c)
For each of the same sample sizes, calculate the exact probability $P(|p̂ ₙ − 0.20| > 0.05)$ using the binomial distribution (R is allowed). Compare with the bounds in (b).  

First translate the event into a statement about the count $S_n$:

$$
|\widehat p_n-0.20|>0.05
\iff
\widehat p_n<0.15
\quad\text{or}\quad
\widehat p_n>0.25.
$$

Since $\widehat p_n=S_n/n$,

$$
|\widehat p_n-0.20|>0.05
\iff
S_n<0.15n
\quad\text{or}\quad
S_n>0.25n.
$$

The inequality is strict, so proportions exactly equal to $0.15$ or $0.25$ are **not included**. Because $S_n$ must be an integer, the events become:

| $n$ | Equivalent count event |
|---:|:---|
| 100 | $S_{100}\leq14$ or $S_{100}\geq26$ |
| 400 | $S_{400}\leq59$ or $S_{400}\geq101$ |
| 1000 | $S_{1000}\leq149$ or $S_{1000}\geq251$ |

For example, when $n=100$,

$$
\begin{aligned}
P\bigl(|\widehat p_{100}-0.20|>0.05\bigr)
&=P(S_{100}\leq14)+P(S_{100}\geq26)\\
&=\sum_{k=0}^{14}\binom{100}{k}(0.20)^k(0.80)^{100-k}\\
&\quad+\sum_{k=26}^{100}\binom{100}{k}(0.20)^k(0.80)^{100-k}.
\end{aligned}
$$

We add the two tails because a count cannot simultaneously be too low and too high. Evaluating the corresponding binomial sums gives:

| $n$ | Exact probability, rounded | Exact percentage | Chebyshev bound |
|---:|---:|---:|---:|
| 100 | $0.167919106$ | $16.7919\%$ | $64\%$ |
| 400 | $0.010304550$ | $1.03045\%$ | $16\%$ |
| 1000 | $0.0000690361$ | $0.00690361\%$ | $6.4\%$ |

Every exact probability is below its Chebyshev bound, as it must be. The bounds are conservative here, especially at $n=1000$: the full binomial distribution gives much more information than the mean and variance alone.

### Reproducing the calculation in R

In R, `pbinom(k, size = n, prob = p)` calculates $P(S_n\leq k)$. With `lower.tail = FALSE`, it calculates $P(S_n>k)$, so an upper tail beginning at $251$ uses a cutoff of $250$.

```r
n <- c(100L, 400L, 1000L)
p <- 0.20

# Strict event: S_n < 0.15*n or S_n > 0.25*n
# These boundaries are integers for all three sample sizes.
lower_last <- c(14L, 59L, 149L)
upper_boundary <- c(25L, 100L, 250L)

exact <- pbinom(lower_last, size = n, prob = p) +
         pbinom(upper_boundary, size = n, prob = p,
                lower.tail = FALSE)

chebyshev <- p * (1 - p) / (n * 0.05^2)

data.frame(
  n = n,
  chebyshev_bound = chebyshev,
  exact_probability = exact,
  exact_percent = 100 * exact
)
```

These are evaluations of the exact binomial distribution, not simulations or normal approximations. Only the displayed decimal values are rounded.

## Weak and Strong Laws: part (d)
State what the Weak Law of Large Numbers says about $p̂ ₙ$ in this setting. What does the Strong Law of Large Numbers say?  

Both laws apply because the $X_i$ are independent and identically distributed Bernoulli random variables with finite mean. They describe convergence in different senses.

### Weak Law of Large Numbers

The Weak Law says

$$
\boxed{\widehat p_n\xrightarrow{P}0.20.}
$$

By definition, this means that for every fixed $\varepsilon>0$,

$$
P\bigl(|\widehat p_n-0.20|>\varepsilon\bigr)\longrightarrow0
\qquad\text{as }n\to\infty.
$$

In plain English, the probability of an error larger than any fixed tolerance approaches zero. We can establish this directly using Chebyshev:

$$
0\leq
P\bigl(|\widehat p_n-0.20|>\varepsilon\bigr)
\leq\frac{0.16}{n\varepsilon^2}
\longrightarrow0.
$$

### Strong Law of Large Numbers

The Strong Law says

$$
\boxed{P\left(\lim_{n\to\infty}\widehat p_n=0.20\right)=1.}
$$

Imagine following one indefinitely growing sequence of patients and updating the proportion after every patient. With probability one, that sequence of proportions converges to $0.20$.

The distinction is that the Weak Law concerns the probability of a large error at sample size $n$, whereas the Strong Law concerns convergence of the entire sequence. Under the Strong Law, with probability one, any fixed tolerance is eventually satisfied and remains satisfied thereafter, but the point at which this happens can depend on the realized sequence.

## Stability is not a finite-sample guarantee: part (e)
Explain why increasing sample size improves the stability of the observed event proportion but does not guarantee that a particular finite-sample estimate will be within 0.05 of 0.20.  

The variance and standard deviation of the observed proportion are

$$
\mathrm{Var}(\widehat p_n)=\frac{0.16}{n},
\qquad
\mathrm{SD}(\widehat p_n)=\frac{0.40}{\sqrt n}.
$$

Increasing sample size therefore makes the sampling distribution more concentrated around $0.20$. For example, multiplying the sample size by four halves the standard deviation.

However, concentration does not eliminate every possible unusual sample. For any finite $n$, there is a positive probability that no patient experiences an event:

$$
P(S_n=0)=(0.80)^n>0.
$$

On that outcome, $\widehat p_n=0$, so its distance from $0.20$ is $0.20$, well above $0.05$. Thus, **a larger sample improves reliability but does not guarantee accuracy for a particular finite sample**.

Nor must a running estimate move closer to $0.20$ after every additional patient; individual updates can move it farther away. The laws of large numbers describe long-run convergence, not monotonic improvement at every step.

A reusable strategy is to identify the Bernoulli indicators, turn their sum into a binomial count, derive the sample proportion’s moments, and distinguish a probability bound from an exact probability.
