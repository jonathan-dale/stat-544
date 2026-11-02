

A factory manufactures electronic sensors used in medical devices. Each sensor independently has a probability
of 0.002 of being defective. A quality-control engineer examines 1,500 sensors produced during one shift.
(a) Let X denote the number of defective sensors. Identify the exact probability distribution of X.
(b) Use a Poisson distribution to approximate P(X = 4).
(c) Approximate P(X ≥ 5).
(d) Calculate the exact binomial probabilities for (b) and (c), and compare them with the Poisson approximations. 


The exact model is **binomial**, and the Poisson approximation uses $\lambda=np=3$. Here are the answers first; then we’ll unpack why these distributions apply.

| Part | Requested result | Answer |
|---|---|---:|
| (a) | Exact distribution | $X\sim\mathrm{Binomial}(1500,0.002)$ |
| (b) | Poisson approximation to $P(X=4)$ | $0.168031$ |
| (c) | Poisson approximation to $P(X\geq5)$ | $0.184737$ |
| (d) | Exact binomial $P(X=4)$ | $0.168200$ |
| (d) | Exact binomial $P(X\geq5)$ | $0.184568$ |

The central idea is that a Poisson distribution approximates a binomial count when there are many independent opportunities for an event, but the event is rare on each opportunity. The approximation matches the expected count, $\lambda=np$. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

## (a) Why the exact distribution is binomial

For each sensor, define an indicator random variable:

$$
I_i=
\begin{cases}
1,&\text{sensor }i\text{ is defective},\\
0,&\text{sensor }i\text{ is not defective}.
\end{cases}
$$

Each $I_i$ is a Bernoulli random variable with parameter $p=0.002$. The total number of defective sensors is

$$
X=I_1+I_2+\cdots+I_{1500}.
$$

A binomial random variable counts successes in a fixed number of independent Bernoulli trials with the same success probability. Here, “success” is simply the outcome we are counting—being defective—not a desirable outcome. Independence and a constant probability are essential to this model. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

Thus,

$$
\boxed{X\sim\mathrm{Binomial}(n=1500,p=0.002).}
$$

Its probability mass function is

$$
P(X=k)=\binom{1500}{k}(0.002)^k(0.998)^{1500-k},
\qquad k=0,1,\ldots,1500.
$$

To understand that formula, break it into two pieces:

- $(0.002)^k(0.998)^{1500-k}$ is the probability of one particular arrangement of $k$ defective and $1500-k$ nondefective sensors, by independence.
- $\binom{1500}{k}$ counts how many different arrangements have exactly $k$ defective sensors.

The binomial mean and variance are $np$ and $np(1-p)$, respectively.  In this problem, [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson24)

$$
E[X]=1500(0.002)=3,
$$

$$
\mathrm{Var}(X)=1500(0.002)(0.998)=2.994.
$$

So we expect three defective sensors per shift on average—not necessarily exactly three in any particular shift.

## Why a Poisson approximation makes sense

The relevant limiting result is

$$
n\to\infty,\qquad p\to0,\qquad np\to\lambda
\quad\Longrightarrow\quad
\mathrm{Binomial}(n,p)\text{ approaches }\mathrm{Poisson}(\lambda).
$$

In words: as the number of independent opportunities grows and the event probability shrinks, with the expected count staying finite, the binomial probabilities approach Poisson probabilities. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

Here, $n=1500$ is large and $p=0.002$ is small, so we use a separate approximating variable

$$
Y\sim\mathrm{Poisson}(3).
$$

Writing $Y$ separately helps keep the distinction clear: $X$ is exactly binomial; $Y$ provides approximate probabilities for $X$.

The Poisson probability mass function is

$$
P(Y=k)=\frac{e^{-\lambda}\lambda^k}{k!},
\qquad k=0,1,2,\ldots.
$$

Unlike the binomial formula, it depends only on the expected count $\lambda$, rather than separately on $n$ and $p$. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

### Where the Poisson formula comes from

For a fixed count $k$, substitute $p=\lambda/n$ into the binomial formula:

$$
\begin{aligned}
P(X=k)
&=\binom{n}{k}
\left(\frac{\lambda}{n}\right)^k
\left(1-\frac{\lambda}{n}\right)^{n-k}\\
&=
\frac{\lambda^k}{k!}
\left[\prod_{j=0}^{k-1}\left(1-\frac{j}{n}\right)\right]
\left(1-\frac{\lambda}{n}\right)^n
\left(1-\frac{\lambda}{n}\right)^{-k}.
\end{aligned}
$$

As $n\to\infty$, the three factors after $\lambda^k/k!$ approach $1$, $e^{-\lambda}$, and $1$, respectively. That leaves

$$
P(X=k)\longrightarrow \frac{e^{-\lambda}\lambda^k}{k!}.
$$

This is the mathematical basis of the approximation—not just a rule about “large $n$, small $p$.” [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

## (b) Approximate $P(X=4)$

Using $Y\sim\mathrm{Poisson}(3)$,

$$
\begin{aligned}
P(X=4)
&\approx P(Y=4)\\
&=\frac{e^{-3}3^4}{4!}\\
&=\frac{81e^{-3}}{24}\\
&\approx \boxed{0.168031}.
\end{aligned}
$$

Thus, the approximate probability of exactly four defective sensors is $16.8031\%$.

## (c) Approximate $P(X\geq5)$

“At least five” includes $5,6,7,\ldots$. Rather than adding all those probabilities, use the complement:

$$
P(Y\geq5)=1-P(Y\leq4).
$$

The endpoint matters: the complement of “at least five” is “at most four,” not “at most five.”

Therefore,

$$
\begin{aligned}
P(X\geq5)
&\approx 1-\sum_{k=0}^{4}\frac{e^{-3}3^k}{k!}\\
&=1-e^{-3}
\left(
1+3+\frac{9}{2}+\frac{27}{6}+\frac{81}{24}
\right)\\
&=1-16.375e^{-3}\\
&\approx \boxed{0.184737}.
\end{aligned}
$$

So the approximate probability of at least five defective sensors is $18.4737\%$.

## (d) Exact binomial probabilities and comparison

“Exact binomial” means using the original binomial model rather than approximating it with another distribution. The displayed decimals are still rounded.

For exactly four defects,

$$
\begin{aligned}
P(X=4)
&=\binom{1500}{4}(0.002)^4(0.998)^{1496}\\
&\approx \boxed{0.168200}.
\end{aligned}
$$

For at least five defects, use the same complement strategy:

$$
\begin{aligned}
P(X\geq5)
&=1-P(X\leq4)\\
&=1-\sum_{k=0}^{4}
\binom{1500}{k}(0.002)^k(0.998)^{1500-k}\\
&\approx \boxed{0.184568}.
\end{aligned}
$$

The numerical comparison is:

| Event | Exact binomial | Poisson approximation | Absolute error | Relative error |
|---|---:|---:|---:|---:|
| $X=4$ | 0.168200 | 0.168031 | 0.000168 | 0.1000% |
| $X\geq5$ | 0.184568 | 0.184737 | 0.000168 | 0.0912% |

Here,

$$
\text{relative error}
=
\frac{|\text{approximation}-\text{exact}|}{\text{exact}}\times100\%.
$$

The Poisson approximation slightly underestimates $P(X=4)$ and slightly overestimates $P(X\geq5)$. Both errors are only about $0.0168$ percentage points. This close agreement is consistent with the large-$n$, small-$p$ setting in which the approximation works well. [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson12)

### Checking the results in R

These commands reproduce all four probabilities:

```r
n <- 1500
p <- 0.002
lambda <- n * p

# Poisson approximations
dpois(4, lambda = lambda)
ppois(4, lambda = lambda, lower.tail = FALSE)

# Exact binomial probabilities
dbinom(4, size = n, prob = p)
pbinom(4, size = n, prob = p, lower.tail = FALSE)
```

For both tail calculations, `lower.tail = FALSE` with a cutoff of $4$ gives $P(X>4)=P(X\geq5)$.

The reusable pattern is: identify the exact count model first, set $\lambda=np$ for the Poisson approximation, and translate the event carefully before calculating—especially the cutoff in a tail probability.