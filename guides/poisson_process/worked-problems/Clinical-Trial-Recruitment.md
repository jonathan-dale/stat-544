# Clinical trial recruitment: Poisson arrivals and consent

A multicenter clinical trial recruits eligible patients arriving independently at three hospitals. Every eligible patient independently consents with probability 0.60. 

| Hospital | Eligible patients per week |
|---|---:|
| Hospital A | 5 |
| Hospital B | 3 |
| Hospital C | 2 |

(a) Determine the rate of the combined eligible-arrival process.  
(b) Determine the rate of the enrolled-patient process.  
(c) Calculate P(exactly 8 enrollments in one week).  
(d) Calculate P(at least 15 enrollments over two weeks).  
(e) Find expected enrollment over 12 weeks.  
(f) An investigator says 12-week recruitment must be binomial because each eligible person has a 60% probability of enrolling. Explain whether you agree.  


## Solution
The key is to separate two sources of randomness: **how many eligible patients arrive**, and **which of those patients consent**. Each patient’s consent is a Bernoulli trial, but that does not automatically make the total enrollment binomial, because the number of patients arriving is itself random.

## Set up the model

The intended assumption appears to be that the three hospitals have **independent Poisson arrival processes**, with constant rates of 5, 3, and 2 eligible patients per week. Strictly speaking, the wording you provided does not explicitly say “Poisson”; rates and independence alone are not enough to determine the count distribution, so we need that modeling assumption for the probability calculations.

A rate of 5 patients per week means 5 per week on average, not exactly 5 every week. Under a Poisson-process model with rate $\lambda$, the count in an interval of length $t$ has a Poisson distribution with parameter $\lambda t$. [probabilitycourse](https://www.probabilitycourse.com/chapter11/11_1_3_merging_and_splitting_poisson_processes.php)

It helps to distinguish these quantities:

- **Rate $\lambda$:** expected arrivals per unit of time.
- **Interval length $t$:** how long you observe the process.
- **Poisson parameter $\mu=\lambda t$:** expected count during that interval.

For a Poisson random variable $X$ with parameter $\mu$, the probability mass function is

$$
P(X=k)=\frac{e^{-\mu}\mu^k}{k!},
\qquad k=0,1,2,\ldots
$$

This is the formula we will use for exact counts. The same Poisson mass function also appears in the derivation of the thinning result below. [probabilitycourse](https://www.probabilitycourse.com/chapter11/11_1_3_merging_and_splitting_poisson_processes.php)

## Combine the hospitals and account for consent

### Part (a): Combined eligible-arrival rate

Combining independent Poisson processes gives another Poisson process whose rate is the sum of the individual rates. This operation is called **superposition**, or merging. [probabilitycourse](https://www.probabilitycourse.com/chapter11/11_1_3_merging_and_splitting_poisson_processes.php)

Therefore,

$$
\lambda_{\text{eligible}}
=\lambda_A+\lambda_B+\lambda_C
=5+3+2
=\boxed{10\text{ eligible patients per week}}.
$$

Let $N(t)$ be the total number of eligible patients arriving over $t$ weeks. Our combined model is

$$
N(t)\sim\mathrm{Poisson}(10t).
$$

### Part (b): Enrolled-patient rate

Now independently “keep” each eligible arrival with probability $p=0.60$, where keeping means enrolling. Independently retaining each arrival of a Poisson process produces another Poisson process with rate $\lambda p$; this is called **thinning**. [probabilitycourse](https://www.probabilitycourse.com/chapter11/11_1_3_merging_and_splitting_poisson_processes.php)

Thus,

$$
\lambda_{\text{enrolled}}
=10(0.60)
=\boxed{6\text{ enrolled patients per week}}.
$$

Let $X(t)$ be the number enrolled over $t$ weeks. The model we will use for all remaining parts is

$$
\boxed{X(t)\sim\mathrm{Poisson}(6t)}.
$$

As a check, we could thin each hospital separately: their enrollment rates would be $5(0.60)=3$, $3(0.60)=1.8$, and $2(0.60)=1.2$. Adding those rates gives the same total, $6$ enrollments per week.

## Calculate the probabilities and expectation

### Part (c): Exactly 8 enrollments in one week

One week means $t=1$, so the Poisson parameter is

$$
\mu=6(1)=6.
$$

Substitute $\mu=6$ and $k=8$ into the Poisson formula:

$$
\begin{aligned}
P(X(1)=8)
&=\frac{e^{-6}6^8}{8!}\\
&\approx \boxed{0.10326}.
\end{aligned}
$$

There is about a **10.33%** chance of exactly 8 enrollments in one week. Notice that 6 is the expected count, while 8 is the particular count whose probability we want.

### Part (d): At least 15 enrollments over two weeks

Two weeks means $t=2$, so

$$
X(2)\sim\mathrm{Poisson}(12).
$$

“At least 15” includes $15,16,17,\ldots$, so the complementary event is “14 or fewer.” Using that complement avoids adding infinitely many probabilities:

$$
\begin{aligned}
P(X(2)\ge 15)
&=1-P(X(2)\le 14)\\
&=1-\sum_{k=0}^{14}\frac{e^{-12}12^k}{k!}\\
&\approx \boxed{0.22798}.
\end{aligned}
$$

Thus, the probability is about **22.80%**. The cutoff is 14, not 15, because subtracting $P(X(2)\le15)$ would also remove the outcome of exactly 15.

### Part (e): Expected enrollment over 12 weeks

Expected enrollment is the enrollment rate multiplied by the observation period:

$$
E[X(12)]
=6(12)
=\boxed{72\text{ enrolled patients}}.
$$

This is an average over repeated hypothetical 12-week recruitment periods, not a guarantee. The full model for that interval is $X(12)\sim\mathrm{Poisson}(72)$.

## Part (f): Why isn’t total enrollment binomial?

The investigator has identified the consent mechanism correctly, but has overlooked the random number of eligible arrivals. **Enrollment is binomial conditional on the number of eligible patients, but Poisson before that number is known**; this distinction is the basis of the thinning theorem. [probabilitycourse](https://www.probabilitycourse.com/chapter11/11_1_3_merging_and_splitting_poisson_processes.php)

For example, if we are told that exactly $n$ eligible patients arrived during the 12 weeks, then we have $n$ independent consent trials:

$$
\boxed{X(12)\mid N(12)=n
\sim\mathrm{Binomial}(n,0.60)}.
$$

But the problem does not give a fixed $n$. Instead, under our arrival model,

$$
N(12)\sim\mathrm{Poisson}(10\cdot12)
=\mathrm{Poisson}(120).
$$

The number 120 is the **expected number of eligible arrivals**, not a fixed number of trials. Replacing the random $N(12)$ with its mean would incorrectly turn the model into $\mathrm{Binomial}(120,0.60)$.

### How the conditional binomial becomes Poisson

We can see the connection directly by averaging over every possible eligible-patient count. For $k$ enrollments, there must have been at least $k$ eligible arrivals:

$$
\begin{aligned}
P(X(12)=k)
&=\sum_{n=k}^{\infty}
P(X(12)=k\mid N(12)=n)\,P(N(12)=n)\\
&=\sum_{n=k}^{\infty}
\binom nk(0.60)^k(0.40)^{n-k}
\frac{e^{-120}120^n}{n!}.
\end{aligned}
$$

Using $\binom nk=n!/[k!(n-k)!]$ and setting $m=n-k$, this simplifies to

$$
\begin{aligned}
P(X(12)=k)
&=\frac{e^{-120}72^k}{k!}
\sum_{m=0}^{\infty}\frac{48^m}{m!}\\
&=\frac{e^{-120}72^k}{k!}e^{48}\\
&=\frac{e^{-72}72^k}{k!}.
\end{aligned}
$$

That is exactly the mass function of $\mathrm{Poisson}(72)$, illustrating the general thinning derivation. [probabilitycourse](https://www.probabilitycourse.com/chapter11/11_1_3_merging_and_splitting_poisson_processes.php)

So a concise answer for part (f) is:

> I disagree with “must be binomial.” Conditional on exactly $n$ eligible arrivals, enrollment is $\mathrm{Binomial}(n,0.60)$. However, the eligible-arrival count is random and Poisson, so the unconditional 12-week enrollment count is $\mathrm{Poisson}(72)$.

This is an exact result under the Poisson-process assumptions, not a small-probability approximation to a binomial distribution. The consent probability being $0.60$ does not cause a problem.

## Check the calculations in R

After setting up the distributions and events, the numerical checks are short. Here, `dpois` evaluates an exact-count probability, while the upper-tail call computes $P(X>14)=P(X\ge15)$.

```r
# (c) Exactly 8 enrollments in one week
dpois(8, lambda = 6)
# 0.1032577

# (d) At least 15 enrollments over two weeks
ppois(14, lambda = 12, lower.tail = FALSE)
# 0.2279755

# (e) Expected enrollment over 12 weeks
6 * 12
# 72
```

Would you like to walk through the algebra in the thinning derivation more slowly, especially why the infinite sum becomes $e^{48}$?
