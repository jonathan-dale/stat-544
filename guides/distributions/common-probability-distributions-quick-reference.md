# Common Probability Distributions — Graduate Quick Reference

A single-page reference for the distributions that appear most often in graduate probability and mathematical statistics (Casella & Berger / Blitzstein & Hwang level). Each section gives the setup in plain English, the PMF/PDF, CDF where it has a closed form, the moment generating function (MGF), mean, variance, standard deviation, notable properties, a short derivation, and the matching R functions.

**Conventions.** $X \sim \text{Dist}(\theta)$ means $X$ has that distribution. $q = 1-p$. $\Gamma(\alpha)=\int_0^\infty t^{\alpha-1}e^{-t}\,dt$ and $B(\alpha,\beta)=\Gamma(\alpha)\Gamma(\beta)/\Gamma(\alpha+\beta)$. The exponential and gamma use the **rate** parameter $\lambda$ (or $\beta$) unless noted, matching R's default `rate=`. Always check which parameterization a textbook uses — it is the most common source of off-by-a-reciprocal errors.

R naming pattern: `d*` = PMF/PDF, `p*` = CDF, `q*` = quantile, `r*` = random draws.

---

## Table of Contents

**Discrete**
1. [Bernoulli](#bernoulli)
2. [Binomial](#binomial)
3. [Geometric](#geometric)
4. [Negative Binomial](#negative-binomial)
5. [Hypergeometric](#hypergeometric)
6. [Poisson](#poisson)
7. [Discrete Uniform](#discrete-uniform)
8. [Multinomial](#multinomial)

**Continuous**

9. [Continuous Uniform](#continuous-uniform)
10. [Normal (Gaussian)](#normal-gaussian)
11. [Exponential](#exponential)
12. [Gamma](#gamma)
13. [Beta](#beta)
14. [Chi-Square](#chi-square)
15. [Student's t](#students-t)
16. [F (Snedecor)](#f-snedecor)
17. [Lognormal](#lognormal)
18. [Weibull](#weibull)
19. [Cauchy](#cauchy)
20. [Laplace (Double Exponential)](#laplace-double-exponential)
21. [Multivariate Normal](#multivariate-normal)

**Reference**

22. [Relationships Between Distributions](#relationships-between-distributions)
23. [Summary Table](#summary-table)

---

# Discrete Distributions

## Bernoulli

**Setup.** A single trial with success probability $p$. $X=1$ on success, $0$ on failure. The building block for binomial, geometric, and negative binomial.

| Quantity | Formula |
|---|---|
| Support | $x \in \{0,1\}$, $0 \le p \le 1$ |
| PMF | $P(X=x) = p^x q^{1-x}$ |
| MGF | $M(t) = q + pe^t$ |
| Mean | $E[X] = p$ |
| Variance | $\text{Var}(X) = pq$ |
| SD | $\sqrt{pq}$ |

**Derivation.** $E[X] = 0\cdot q + 1\cdot p = p$. Since $X^2 = X$ for $X\in\{0,1\}$, $E[X^2]=p$, so $\text{Var}(X)=p-p^2=pq$.

**Notable facts.**
- Indicator variables are Bernoulli: $I_A \sim \text{Bern}(P(A))$, and $E[I_A]=P(A)$ — the basis of the "indicator method" for expectations.
- Variance is maximized at $p=1/2$.
- Exponential family with natural parameter $\log\frac{p}{1-p}$ (the logit), which is why logistic regression uses it.

**R:** `dbinom(x, 1, p)`, `rbinom(n, 1, p)`

[Back to top](#table-of-contents)

---

## Binomial

**Setup.** Number of successes in $n$ independent Bernoulli($p$) trials.

| Quantity | Formula |
|---|---|
| Support | $x \in \{0,1,\dots,n\}$ |
| PMF | $P(X=x)=\binom{n}{x}p^x q^{n-x}$ |
| MGF | $M(t) = (q+pe^t)^n$ |
| Mean | $np$ |
| Variance | $npq$ |
| SD | $\sqrt{npq}$ |

**Derivation of the PMF.** Any specific sequence with $x$ successes and $n-x$ failures has probability $p^x q^{n-x}$ by independence. There are $\binom{n}{x}$ ways to choose which positions are successes, and the sequences are disjoint, so add them up. The PMF sums to 1 by the binomial theorem: $\sum_x \binom{n}{x}p^xq^{n-x}=(p+q)^n=1$.

**Moments.** Write $X=\sum_{i=1}^n X_i$ with $X_i \overset{iid}{\sim}\text{Bern}(p)$. Linearity gives $E[X]=np$; independence gives $\text{Var}(X)=\sum \text{Var}(X_i)=npq$. The MGF of a sum of independent variables is the product of MGFs: $(q+pe^t)^n$.

**Notable facts.**
- Additivity: $X\sim\text{Bin}(n,p)$, $Y\sim\text{Bin}(m,p)$ independent $\Rightarrow X+Y\sim\text{Bin}(n+m,p)$.
- Normal approximation (CLT): $X \approx N(np, npq)$ when $np$ and $nq$ are both large (rule of thumb $\ge 10$); use a continuity correction of $\pm 0.5$.
- Poisson approximation: $n\to\infty$, $p\to 0$, $np\to\lambda$ gives $\text{Pois}(\lambda)$.
- Mode is $\lfloor (n+1)p \rfloor$.
- Sampling **with** replacement; compare with the [Hypergeometric](#hypergeometric).

**R:** `dbinom(x, size=n, prob=p)`, `pbinom`, `qbinom`, `rbinom`

[Back to top](#table-of-contents)

---

## Geometric

**Setup.** Repeat independent Bernoulli($p$) trials until the first success. Two conventions:
- **Failures before the first success** (R's convention): $X\in\{0,1,2,\dots\}$.
- **Trials up to and including the first success**: $Y = X+1 \in\{1,2,\dots\}$.

| Quantity | Failures $X$ (R) | Trials $Y=X+1$ |
|---|---|---|
| PMF | $q^x p$ | $q^{y-1}p$ |
| CDF | $1-q^{x+1}$ | $1-q^{y}$ |
| MGF | $\dfrac{p}{1-qe^t}$, $t<-\log q$ | $\dfrac{pe^t}{1-qe^t}$ |
| Mean | $q/p$ | $1/p$ |
| Variance | $q/p^2$ | $q/p^2$ |
| SD | $\sqrt{q}/p$ | $\sqrt{q}/p$ |

**Derivation of the PMF.** $X=x$ means $x$ failures followed by one success: $q^x p$. The tail is $P(X\ge x)=q^x$ (the first $x$ trials all fail), which gives the CDF.

**Mean via the derivative trick.** Starting from $\sum_{k=0}^\infty q^k = \frac{1}{1-q}$ and differentiating with respect to $q$ gives $\sum_k kq^{k-1} = \frac{1}{(1-q)^2}$. So
$$E[X]=\sum_{k=0}^\infty k q^k p = pq\sum_k kq^{k-1} = \frac{pq}{p^2}=\frac{q}{p}.$$
Differentiating a second time gives $E[X(X-1)]=2q^2/p^2$, so $\text{Var}(X)=2q^2/p^2+q/p-q^2/p^2 = q/p^2$.

**Notable facts.**
- **Memoryless:** $P(X\ge m+n \mid X\ge m)=P(X\ge n)$. It is the only memoryless discrete distribution, the discrete analog of the [Exponential](#exponential).
- Special case of the negative binomial with $r=1$.

**R:** `dgeom(x, prob=p)` (counts failures)

[Back to top](#table-of-contents)

---

## Negative Binomial

**Setup.** Number of failures $X$ before the $r$-th success in independent Bernoulli($p$) trials (R's convention). The alternate convention counts total trials $Y=X+r$.

| Quantity | Formula (failures $X$) |
|---|---|
| Support | $x\in\{0,1,2,\dots\}$, $r>0$ |
| PMF | $P(X=x)=\binom{x+r-1}{x}p^r q^x$ |
| MGF | $\left(\dfrac{p}{1-qe^t}\right)^r$, $t<-\log q$ |
| Mean | $rq/p$ |
| Variance | $rq/p^2$ |
| SD | $\sqrt{rq}/p$ |

**Derivation of the PMF.** The final trial must be the $r$-th success. The first $x+r-1$ trials contain exactly $r-1$ successes and $x$ failures, which can be arranged in $\binom{x+r-1}{x}$ ways. Each arrangement has probability $p^{r}q^x$ (counting the final success).

**Moments.** $X=\sum_{i=1}^r G_i$ with $G_i\overset{iid}{\sim}\text{Geom}(p)$ (failures between consecutive successes), so mean and variance are $r$ times the geometric values.

**Notable facts.**
- **Overdispersion:** $\text{Var}(X)>E[X]$, unlike the Poisson. It is the standard model for overdispersed counts (e.g., `MASS::glm.nb`).
- **Poisson–Gamma mixture:** if $X\mid\Lambda\sim\text{Pois}(\Lambda)$ and $\Lambda\sim\text{Gamma}(r,\text{rate}=p/q)$, then $X\sim\text{NB}(r,p)$. This lets $r$ be any positive real.
- Mean parameterization: with $\mu=rq/p$, $\text{Var}(X)=\mu+\mu^2/r$.

**R:** `dnbinom(x, size=r, prob=p)` or `dnbinom(x, size=r, mu=mu)`

[Back to top](#table-of-contents)

---

## Hypergeometric

**Setup.** An urn holds $N$ items, $K$ of them "successes." Draw $n$ **without replacement**. $X$ is the number of successes drawn.

| Quantity | Formula |
|---|---|
| Support | $\max(0,n-N+K)\le x\le\min(n,K)$ |
| PMF | $P(X=x)=\dfrac{\binom{K}{x}\binom{N-K}{n-x}}{\binom{N}{n}}$ |
| Mean | $n\dfrac{K}{N}$ |
| Variance | $n\dfrac{K}{N}\dfrac{N-K}{N}\dfrac{N-n}{N-1}$ |
| SD | square root of the variance |

**Derivation of the PMF.** All $\binom{N}{n}$ subsets are equally likely. Favorable subsets choose $x$ from the $K$ successes and $n-x$ from the $N-K$ failures.

**Moments.** Write $X=\sum_{i=1}^n I_i$ where $I_i$ indicates that draw $i$ is a success. By symmetry $E[I_i]=K/N$, so $E[X]=nK/N$. The draws are dependent: for $i\ne j$, $P(I_i=I_j=1)=\frac{K(K-1)}{N(N-1)}$, so $\text{Cov}(I_i,I_j)<0$. Summing the variances and the $n(n-1)$ covariances yields the formula above.

**Notable facts.**
- $\frac{N-n}{N-1}$ is the **finite population correction**. It goes to 1 as $N\to\infty$, recovering the binomial with $p=K/N$.
- The PMF sums to 1 by Vandermonde's identity.
- Fisher's exact test for $2\times2$ tables uses this distribution.

**R:** `dhyper(x, m=K, n=N-K, k=n)` — R's argument names differ from the notation here.

[Back to top](#table-of-contents)

---

## Poisson

**Setup.** Counts of events in a fixed interval when events occur independently at a constant average rate $\lambda$.

| Quantity | Formula |
|---|---|
| Support | $x\in\{0,1,2,\dots\}$, $\lambda>0$ |
| PMF | $P(X=x)=\dfrac{e^{-\lambda}\lambda^x}{x!}$ |
| MGF | $M(t)=\exp\{\lambda(e^t-1)\}$ |
| Mean | $\lambda$ |
| Variance | $\lambda$ |
| SD | $\sqrt{\lambda}$ |

**Derivation as a limit of the binomial.** Let $X_n\sim\text{Bin}(n,\lambda/n)$. Then
$$P(X_n=x)=\frac{n(n-1)\cdots(n-x+1)}{n^x}\cdot\frac{\lambda^x}{x!}\cdot\left(1-\frac{\lambda}{n}\right)^{n}\left(1-\frac{\lambda}{n}\right)^{-x}.$$
As $n\to\infty$, the first factor goes to 1, $(1-\lambda/n)^n\to e^{-\lambda}$, and the last factor goes to 1, leaving $e^{-\lambda}\lambda^x/x!$.

**MGF and moments.**
$$M(t)=\sum_{x}e^{tx}\frac{e^{-\lambda}\lambda^x}{x!}=e^{-\lambda}\sum_x\frac{(\lambda e^t)^x}{x!}=e^{\lambda(e^t-1)}.$$
$M'(0)=\lambda$ and $M''(0)=\lambda^2+\lambda$, so the variance is $\lambda$.

**Notable facts.**
- Equidispersion: mean equals variance. A sample variance well above the sample mean suggests a [Negative Binomial](#negative-binomial).
- Additivity: independent $\text{Pois}(\lambda_1)+\text{Pois}(\lambda_2)\sim\text{Pois}(\lambda_1+\lambda_2)$.
- Conditioning: $X_1\mid X_1+X_2=n \sim \text{Bin}\!\left(n,\frac{\lambda_1}{\lambda_1+\lambda_2}\right)$.
- Thinning: if each of $N\sim\text{Pois}(\lambda)$ events is kept independently with probability $p$, the kept count is $\text{Pois}(\lambda p)$, independent of the discarded count.
- In a Poisson process with rate $\lambda$, interarrival times are [Exponential](#exponential)($\lambda$) and the waiting time to the $k$-th event is [Gamma](#gamma)($k,\lambda$).
- Normal approximation for large $\lambda$: $X\approx N(\lambda,\lambda)$. The square root $\sqrt{X}$ is variance-stabilizing.

**R:** `dpois(x, lambda)`, `ppois`, `qpois`, `rpois`

[Back to top](#table-of-contents)

---

## Discrete Uniform

**Setup.** Each integer from $a$ to $b$ is equally likely. Let $n=b-a+1$.

| Quantity | Formula |
|---|---|
| PMF | $1/n$ for $x\in\{a,\dots,b\}$ |
| CDF | $(\lfloor x\rfloor-a+1)/n$ for $a \le x \le b$ |
| MGF | $\dfrac{e^{at}-e^{(b+1)t}}{n(1-e^t)}$ |
| Mean | $(a+b)/2$ |
| Variance | $(n^2-1)/12$ |
| SD | $\sqrt{(n^2-1)/12}$ |

**Derivation (for $\{1,\dots,n\}$).** $E[X]=\frac1n\sum k=\frac{n+1}{2}$ and $E[X^2]=\frac1n\sum k^2=\frac{(n+1)(2n+1)}{6}$, so $\text{Var}(X)=\frac{(n+1)(2n+1)}{6}-\frac{(n+1)^2}{4}=\frac{n^2-1}{12}$. Shifting by a constant does not change the variance.

**Notable facts.** The German tank problem estimates $n$ from a sample; the UMVUE is $\frac{k+1}{k}M-1$, where $M$ is the sample maximum of $k$ draws without replacement.

**R:** `sample(a:b, size, replace=TRUE)`

[Back to top](#table-of-contents)

---

## Multinomial

**Setup.** $n$ independent trials, each landing in one of $k$ categories with probabilities $p_1,\dots,p_k$ ($\sum p_j=1$). $\mathbf X=(X_1,\dots,X_k)$ counts the outcomes in each category, so $\sum X_j=n$.

| Quantity | Formula |
|---|---|
| PMF | $P(\mathbf X=\mathbf x)=\dfrac{n!}{x_1!\cdots x_k!}p_1^{x_1}\cdots p_k^{x_k}$ |
| MGF | $\left(\sum_{j}p_je^{t_j}\right)^n$ |
| Mean | $E[X_j]=np_j$ |
| Variance | $\text{Var}(X_j)=np_j(1-p_j)$ |
| Covariance | $\text{Cov}(X_i,X_j)=-np_ip_j$ for $i\ne j$ |

**Derivation.** A specific ordered sequence has probability $\prod p_j^{x_j}$. The multinomial coefficient counts the arrangements of $n$ trials into groups of sizes $x_1,\dots,x_k$.

**Covariance derivation.** $X_i+X_j\sim\text{Bin}(n,p_i+p_j)$, so
$$\text{Var}(X_i+X_j)=n(p_i+p_j)(1-p_i-p_j)=\text{Var}(X_i)+\text{Var}(X_j)+2\text{Cov}(X_i,X_j).$$
Solving gives $\text{Cov}(X_i,X_j)=-np_ip_j$. The covariance is negative because the counts compete for a fixed total $n$.

**Notable facts.**
- Each marginal is binomial, and merging categories gives another multinomial.
- Covariance matrix: $n(\text{diag}(\mathbf p)-\mathbf p\mathbf p^\top)$. It is singular because $\sum X_j=n$.
- Pearson's $\chi^2$ goodness-of-fit statistic $\sum (X_j-np_j)^2/(np_j)\xrightarrow{d}\chi^2_{k-1}$.
- The Dirichlet distribution is its conjugate prior.

**R:** `dmultinom(x, prob=p)`, `rmultinom(n, size, prob)`

[Back to top](#table-of-contents)

---

# Continuous Distributions

## Continuous Uniform

**Setup.** Equally likely over the interval $[a,b]$.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{1}{b-a}$, $a\le x\le b$ |
| CDF | $F(x)=\dfrac{x-a}{b-a}$ |
| MGF | $\dfrac{e^{tb}-e^{ta}}{t(b-a)}$ ($t\ne0$), $M(0)=1$ |
| Mean | $(a+b)/2$ |
| Variance | $(b-a)^2/12$ |
| SD | $(b-a)/\sqrt{12}$ |

**Derivation.** $E[X^2]=\int_a^b\frac{x^2}{b-a}dx=\frac{a^2+ab+b^2}{3}$. Subtracting $\left(\frac{a+b}{2}\right)^2$ gives $\frac{(b-a)^2}{12}$.

**Notable facts.**
- **Probability integral transform:** if $X$ has a continuous CDF $F$, then $F(X)\sim U(0,1)$. Conversely, $F^{-1}(U)$ has CDF $F$, which is the basis of inverse-transform sampling.
- Order statistics of $U(0,1)$: $U_{(k)}\sim\text{Beta}(k,n-k+1)$.
- $-\log U\sim\text{Expo}(1)$.
- For $U(0,\theta)$, the MLE of $\theta$ is $X_{(n)}$, a classic example where regularity conditions fail because the support depends on $\theta$.

**R:** `dunif(x, min=a, max=b)`, `punif`, `qunif`, `runif`

[Back to top](#table-of-contents)

---

## Normal (Gaussian)

**Setup.** The limiting distribution of standardized sums (CLT). It is parameterized by mean $\mu$ and variance $\sigma^2$. $Z\sim N(0,1)$ is the standard normal, with PDF $\phi$ and CDF $\Phi$.

| Quantity | Formula |
|---|---|
| Support | $x\in\mathbb R$ |
| PDF | $f(x)=\dfrac{1}{\sigma\sqrt{2\pi}}\exp\left\{-\dfrac{(x-\mu)^2}{2\sigma^2}\right\}$ |
| CDF | $\Phi\left(\frac{x-\mu}{\sigma}\right)$ (no closed form) |
| MGF | $\exp\{\mu t+\tfrac12\sigma^2t^2\}$ |
| Mean / Median / Mode | $\mu$ |
| Variance | $\sigma^2$ |
| SD | $\sigma$ |
| Skewness / excess kurtosis | $0$ / $0$ |

**Derivation of the normalizing constant (polar trick).** Let $I=\int_{-\infty}^\infty e^{-x^2/2}dx$. Then
$$I^2=\iint e^{-(x^2+y^2)/2}\,dx\,dy=\int_0^{2\pi}\!\!\int_0^\infty e^{-r^2/2}\,r\,dr\,d\theta=2\pi,$$
so $I=\sqrt{2\pi}$.

**MGF derivation (complete the square).**
$$M_Z(t)=\int\frac{1}{\sqrt{2\pi}}e^{tz-z^2/2}dz=e^{t^2/2}\int\frac{1}{\sqrt{2\pi}}e^{-(z-t)^2/2}dz=e^{t^2/2}.$$
For $X=\mu+\sigma Z$: $M_X(t)=e^{\mu t}M_Z(\sigma t)$.

**Notable facts.**
- Linear combinations of independent normals are normal: $\sum a_iX_i\sim N(\sum a_i\mu_i,\sum a_i^2\sigma_i^2)$.
- Empirical rule: about 68%, 95%, and 99.7% of the mass lies within 1, 2, and 3 SD of the mean.
- For iid normal samples, $\bar X$ and $S^2$ are **independent**, and $(n-1)S^2/\sigma^2\sim\chi^2_{n-1}$. This characterizes the normal.
- Moments: $E[Z^{2k}]=(2k-1)!!$; e.g., $E[Z^4]=3$.
- Maximum entropy among distributions with a given mean and variance.
- Mills ratio tail bound: $P(Z>z)\le\phi(z)/z$ for $z>0$.

**R:** `dnorm(x, mean, sd)` — note R takes the SD, not the variance.

[Back to top](#table-of-contents)

---

## Exponential

**Setup.** Waiting time until the first event in a Poisson process with rate $\lambda$.

| Quantity | Rate $\lambda$ | Scale $\beta=1/\lambda$ |
|---|---|---|
| PDF | $\lambda e^{-\lambda x}$, $x\ge0$ | $\frac1\beta e^{-x/\beta}$ |
| CDF | $1-e^{-\lambda x}$ | $1-e^{-x/\beta}$ |
| MGF | $\dfrac{\lambda}{\lambda-t}$, $t<\lambda$ | $\dfrac{1}{1-\beta t}$ |
| Mean | $1/\lambda$ | $\beta$ |
| Variance | $1/\lambda^2$ | $\beta^2$ |
| SD | $1/\lambda$ | $\beta$ |
| Median | $\log 2/\lambda$ | $\beta\log2$ |

**Derivation from the Poisson process.** The first arrival $T$ exceeds $t$ exactly when there are zero events in $[0,t]$. With $N(t)\sim\text{Pois}(\lambda t)$, $P(T>t)=P(N(t)=0)=e^{-\lambda t}$. Differentiating $1-e^{-\lambda t}$ gives the PDF.

**Moments.** $E[X^k]=\int_0^\infty x^k\lambda e^{-\lambda x}dx=\frac{k!}{\lambda^k}$, using the gamma integral.

**Notable facts.**
- **Memoryless:** $P(X>s+t\mid X>s)=P(X>t)$. It is the only continuous memoryless distribution, and it has a constant hazard rate $\lambda$.
- The minimum of independent exponentials is exponential: $\min_i X_i\sim\text{Expo}(\sum\lambda_i)$, and $P(X_j\text{ is the min})=\lambda_j/\sum\lambda_i$.
- The sum of $n$ iid $\text{Expo}(\lambda)$ is $\text{Gamma}(n,\lambda)$.
- $\text{Expo}(1/2)=\chi^2_2$.

**R:** `dexp(x, rate=lambda)`

[Back to top](#table-of-contents)

---

## Gamma

**Setup.** Waiting time until the $\alpha$-th event in a Poisson process (for integer $\alpha$). More generally it is a flexible family on $(0,\infty)$ with shape $\alpha$ and rate $\beta$.

| Quantity | Formula (shape $\alpha$, rate $\beta$) |
|---|---|
| PDF | $f(x)=\dfrac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}$, $x>0$ |
| MGF | $\left(\dfrac{\beta}{\beta-t}\right)^\alpha$, $t<\beta$ |
| Mean | $\alpha/\beta$ |
| Variance | $\alpha/\beta^2$ |
| SD | $\sqrt\alpha/\beta$ |
| Mode | $(\alpha-1)/\beta$ for $\alpha\ge1$ |

With scale $\theta=1/\beta$: mean $\alpha\theta$ and variance $\alpha\theta^2$.

**Derivation of the PDF.** The kernel $x^{\alpha-1}e^{-\beta x}$ integrates to $\Gamma(\alpha)/\beta^\alpha$ (substitute $u=\beta x$ in the gamma function), so dividing by that constant normalizes it.

**Moments via the kernel trick.**
$$E[X^k]=\frac{\beta^\alpha}{\Gamma(\alpha)}\int x^{\alpha+k-1}e^{-\beta x}dx=\frac{\beta^\alpha}{\Gamma(\alpha)}\cdot\frac{\Gamma(\alpha+k)}{\beta^{\alpha+k}}=\frac{\Gamma(\alpha+k)}{\Gamma(\alpha)\beta^k}.$$
So $E[X]=\alpha/\beta$ and $E[X^2]=\alpha(\alpha+1)/\beta^2$.

**Notable facts.**
- Additivity in shape: independent $\text{Gamma}(\alpha_1,\beta)+\text{Gamma}(\alpha_2,\beta)\sim\text{Gamma}(\alpha_1+\alpha_2,\beta)$.
- Special cases: $\text{Gamma}(1,\lambda)=\text{Expo}(\lambda)$; $\text{Gamma}(k/2,1/2)=\chi^2_k$.
- **Beta–Gamma connection:** if $X\sim\text{Gamma}(a,\beta)$ and $Y\sim\text{Gamma}(b,\beta)$ are independent, then $\frac{X}{X+Y}\sim\text{Beta}(a,b)$, independent of $X+Y$.
- Conjugate prior for the Poisson rate and the exponential rate.
- Useful facts: $\Gamma(1/2)=\sqrt\pi$ and $\Gamma(n)=(n-1)!$.

**R:** `dgamma(x, shape=alpha, rate=beta)` — or `scale=`; always name the argument.

[Back to top](#table-of-contents)

---

## Beta

**Setup.** A flexible family on $(0,1)$, commonly used to model probabilities and proportions.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha,\beta)}$, $0<x<1$ |
| Mean | $\dfrac{\alpha}{\alpha+\beta}$ |
| Variance | $\dfrac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$ |
| SD | square root of the variance |
| Mode | $\dfrac{\alpha-1}{\alpha+\beta-2}$ for $\alpha,\beta>1$ |
| MGF | No simple closed form (confluent hypergeometric) |

**Derivation of $B(\alpha,\beta)=\Gamma(\alpha)\Gamma(\beta)/\Gamma(\alpha+\beta)$.** Write $\Gamma(\alpha)\Gamma(\beta)=\iint x^{\alpha-1}y^{\beta-1}e^{-(x+y)}dx\,dy$ and substitute $u=x+y$, $v=x/(x+y)$ (Jacobian $u$). The integral factors into $\Gamma(\alpha+\beta)\cdot\int_0^1 v^{\alpha-1}(1-v)^{\beta-1}dv$.

**Moments via the kernel trick.** $E[X^k]=\frac{B(\alpha+k,\beta)}{B(\alpha,\beta)}$, so $E[X]=\frac{\alpha}{\alpha+\beta}$ and $E[X^2]=\frac{\alpha(\alpha+1)}{(\alpha+\beta)(\alpha+\beta+1)}$.

**Notable facts.**
- $\text{Beta}(1,1)=U(0,1)$.
- **Conjugate prior** for the binomial $p$: a $\text{Beta}(\alpha,\beta)$ prior with $x$ successes in $n$ trials gives a $\text{Beta}(\alpha+x,\beta+n-x)$ posterior.
- The $k$-th order statistic of $n$ iid $U(0,1)$ variables is $\text{Beta}(k,n-k+1)$.
- Symmetric when $\alpha=\beta$. It is U-shaped when both parameters are below 1.
- Its multivariate generalization is the Dirichlet.

**R:** `dbeta(x, shape1=alpha, shape2=beta)`

[Back to top](#table-of-contents)

---

## Chi-Square

**Setup.** The sum of squares of $k$ independent standard normals: $X=\sum_{i=1}^k Z_i^2\sim\chi^2_k$, where $k$ is the degrees of freedom.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{1}{2^{k/2}\Gamma(k/2)}x^{k/2-1}e^{-x/2}$, $x>0$ |
| MGF | $(1-2t)^{-k/2}$, $t<1/2$ |
| Mean | $k$ |
| Variance | $2k$ |
| SD | $\sqrt{2k}$ |

**Derivation for $k=1$.** For $Y=Z^2$ and $y>0$,
$$F_Y(y)=P(-\sqrt y\le Z\le\sqrt y)=2\Phi(\sqrt y)-1\implies f_Y(y)=\frac{1}{\sqrt y}\phi(\sqrt y)=\frac{1}{\sqrt{2\pi}}y^{-1/2}e^{-y/2},$$
which is $\text{Gamma}(1/2,1/2)$. Gamma additivity then gives $\chi^2_k=\text{Gamma}(k/2,1/2)$.

**Notable facts.**
- For a normal sample, $(n-1)S^2/\sigma^2\sim\chi^2_{n-1}$ (one degree of freedom is lost by estimating $\mu$ with $\bar X$, via Cochran's theorem).
- Additive in degrees of freedom.
- Asymptotic distribution of the likelihood ratio statistic (Wilks): $-2\log\Lambda\xrightarrow{d}\chi^2_{r}$.
- The noncentral $\chi^2_k(\lambda)$ arises when the $Z_i$ have nonzero means. It is used in power calculations.

**R:** `dchisq(x, df=k)`

[Back to top](#table-of-contents)

---

## Student's t

**Setup.** $T=\dfrac{Z}{\sqrt{V/\nu}}$ where $Z\sim N(0,1)$ and $V\sim\chi^2_\nu$ are independent. This is the distribution of the studentized mean $\frac{\bar X-\mu}{S/\sqrt n}\sim t_{n-1}$.

| Quantity | Formula |
|---|---|
| PDF | $f(t)=\dfrac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\nu\pi}\,\Gamma\left(\frac\nu2\right)}\left(1+\dfrac{t^2}{\nu}\right)^{-(\nu+1)/2}$ |
| Mean | $0$ for $\nu>1$ (undefined for $\nu\le1$) |
| Variance | $\dfrac{\nu}{\nu-2}$ for $\nu>2$; $\infty$ for $1<\nu\le2$ |
| SD | $\sqrt{\nu/(\nu-2)}$ |
| MGF | Does not exist |

**Derivation sketch.** Condition on $V=v$: then $T\mid V=v\sim N(0,\nu/v)$. Integrate the joint density against the $\chi^2_\nu$ density of $V$:
$$f_T(t)=\int_0^\infty \sqrt{\tfrac{v}{2\pi\nu}}\,e^{-t^2v/(2\nu)}\,f_V(v)\,dv.$$
The integrand is a gamma kernel in $v$ with shape $\frac{\nu+1}{2}$ and rate $\frac12\left(1+\frac{t^2}{\nu}\right)$. Evaluating it with the gamma integral produces the formula above.

**Notable facts.**
- Heavier tails than the normal. Only moments of order $<\nu$ exist.
- $t_\nu\xrightarrow{d}N(0,1)$ as $\nu\to\infty$, and $t_1$ is the standard [Cauchy](#cauchy).
- $T^2\sim F_{1,\nu}$.
- It is a scale mixture of normals, which is why it is used for robust regression.

**R:** `dt(x, df=nu)`, `qt(0.975, df)` for confidence intervals

[Back to top](#table-of-contents)

---

## F (Snedecor)

**Setup.** $F=\dfrac{U/d_1}{V/d_2}$ where $U\sim\chi^2_{d_1}$ and $V\sim\chi^2_{d_2}$ are independent. It is the distribution of variance ratios in ANOVA and nested-model regression tests.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{1}{B\left(\frac{d_1}{2},\frac{d_2}{2}\right)}\left(\dfrac{d_1}{d_2}\right)^{d_1/2}x^{d_1/2-1}\left(1+\dfrac{d_1}{d_2}x\right)^{-(d_1+d_2)/2}$, $x>0$ |
| Mean | $\dfrac{d_2}{d_2-2}$ for $d_2>2$ |
| Variance | $\dfrac{2d_2^2(d_1+d_2-2)}{d_1(d_2-2)^2(d_2-4)}$ for $d_2>4$ |
| MGF | Does not exist |

**Derivation sketch.** $\frac{U}{U+V}\sim\text{Beta}(d_1/2,d_2/2)$ by the Beta–Gamma connection. Since $F=\frac{d_2}{d_1}\cdot\frac{B}{1-B}$ with $B=\frac{U}{U+V}$, a one-to-one transformation of the beta density gives the F density. The mean follows from independence: $E[F]=\frac{d_2}{d_1}E[U]\,E[1/V]=\frac{d_2}{d_1}\cdot d_1\cdot\frac{1}{d_2-2}$.

**Notable facts.**
- $1/F_{d_1,d_2}\sim F_{d_2,d_1}$, so lower quantiles can be computed from upper ones.
- $t_\nu^2\sim F_{1,\nu}$.
- In regression, $F=\frac{(\text{RSS}_0-\text{RSS}_1)/q}{\text{RSS}_1/(n-p)}\sim F_{q,n-p}$ under $H_0$.

**R:** `df(x, df1, df2)`, `pf`, `qf`

[Back to top](#table-of-contents)

---

## Lognormal

**Setup.** $X=e^Y$ where $Y\sim N(\mu,\sigma^2)$. This is a multiplicative analog of the normal, used for incomes, asset prices (as in geometric Brownian motion), and file sizes.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{1}{x\sigma\sqrt{2\pi}}\exp\left\{-\dfrac{(\log x-\mu)^2}{2\sigma^2}\right\}$, $x>0$ |
| CDF | $\Phi\left(\frac{\log x-\mu}{\sigma}\right)$ |
| Mean | $e^{\mu+\sigma^2/2}$ |
| Variance | $(e^{\sigma^2}-1)e^{2\mu+\sigma^2}$ |
| SD | $e^{\mu+\sigma^2/2}\sqrt{e^{\sigma^2}-1}$ |
| Median / Mode | $e^\mu$ / $e^{\mu-\sigma^2}$ |
| MGF | Does not exist (all moments are finite, but $E[e^{tX}]=\infty$ for $t>0$) |

**Derivation.** The CDF is $P(e^Y\le x)=P(Y\le\log x)$. Differentiating brings in the Jacobian $1/x$. Moments come directly from the normal MGF: $E[X^k]=E[e^{kY}]=M_Y(k)=e^{k\mu+k^2\sigma^2/2}$.

**Notable facts.**
- Products of independent lognormals are lognormal.
- Right-skewed with mean > median > mode.
- The parameters $\mu$ and $\sigma$ describe $\log X$, not $X$ — a common mistake.
- It is a classic example where the moments do not uniquely determine the distribution.

**R:** `dlnorm(x, meanlog=mu, sdlog=sigma)`

[Back to top](#table-of-contents)

---

## Weibull

**Setup.** Lifetimes with a hazard that increases, decreases, or stays constant. Shape $k>0$, scale $\lambda>0$.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{k}{\lambda}\left(\dfrac{x}{\lambda}\right)^{k-1}e^{-(x/\lambda)^k}$, $x\ge0$ |
| CDF | $1-e^{-(x/\lambda)^k}$ |
| Hazard | $h(x)=\dfrac{k}{\lambda}\left(\dfrac{x}{\lambda}\right)^{k-1}$ |
| Mean | $\lambda\,\Gamma(1+1/k)$ |
| Variance | $\lambda^2\left[\Gamma(1+2/k)-\Gamma(1+1/k)^2\right]$ |
| Median | $\lambda(\log2)^{1/k}$ |

**Derivation.** If $E\sim\text{Expo}(1)$, then $X=\lambda E^{1/k}$ satisfies $P(X>x)=P(E>(x/\lambda)^k)=e^{-(x/\lambda)^k}$. Moments: $E[X^r]=\lambda^rE[E^{r/k}]=\lambda^r\Gamma(1+r/k)$.

**Notable facts.**
- $k<1$: decreasing hazard (infant mortality). $k=1$: exponential. $k>1$: increasing hazard (wear-out).
- One of the three extreme-value limit families, describing minima of variables bounded below.

**R:** `dweibull(x, shape=k, scale=lambda)`

[Back to top](#table-of-contents)

---

## Cauchy

**Setup.** Location $x_0$, scale $\gamma$. It is the ratio of two independent standard normals, and $t_1$.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{1}{\pi\gamma\left[1+\left(\frac{x-x_0}{\gamma}\right)^2\right]}$ |
| CDF | $\frac12+\frac1\pi\arctan\left(\frac{x-x_0}{\gamma}\right)$ |
| Mean / Variance | **Undefined** / undefined |
| Median / Mode | $x_0$ |
| MGF | Does not exist. The characteristic function is $e^{ix_0t-\gamma|t|}$ |

**Derivation (standard case).** If $\Theta\sim U(-\pi/2,\pi/2)$, then $X=\tan\Theta$ has $P(X\le x)=P(\Theta\le\arctan x)=\frac12+\frac{\arctan x}{\pi}$. Differentiating gives $\frac{1}{\pi(1+x^2)}$.

**Notable facts.**
- $E|X|=\int\frac{|x|}{\pi(1+x^2)}dx=\infty$, so the law of large numbers and CLT do not apply.
- The sample mean of $n$ iid Cauchy variables is again Cauchy with the same parameters — averaging does not help. Use the sample median instead.
- A stable distribution, and a standard counterexample in asymptotics.

**R:** `dcauchy(x, location, scale)`

[Back to top](#table-of-contents)

---

## Laplace (Double Exponential)

**Setup.** Location $\mu$, scale $b$. It looks like two exponentials placed back to back.

| Quantity | Formula |
|---|---|
| PDF | $f(x)=\dfrac{1}{2b}e^{-|x-\mu|/b}$ |
| MGF | $\dfrac{e^{\mu t}}{1-b^2t^2}$, $|t|<1/b$ |
| Mean / Median / Mode | $\mu$ |
| Variance | $2b^2$ |
| SD | $\sqrt2\,b$ |

**Derivation.** The difference of two independent $\text{Expo}(1/b)$ variables is Laplace$(0,b)$. Its variance is $b^2+b^2=2b^2$.

**Notable facts.**
- The MLE of $\mu$ is the sample median, since the likelihood minimizes $\sum|x_i-\mu|$ (least absolute deviations).
- A Laplace prior on regression coefficients corresponds to the LASSO ($L_1$) penalty.
- Excess kurtosis is 3, so it has heavier tails than the normal.

**R:** no base function. Use `extraDistr::dlaplace` or `mu + b*(rexp(n) - rexp(n))`.

[Back to top](#table-of-contents)

---

## Multivariate Normal

**Setup.** $\mathbf X\sim N_p(\boldsymbol\mu,\boldsymbol\Sigma)$ with mean vector $\boldsymbol\mu$ and a positive-definite covariance matrix $\boldsymbol\Sigma$.

| Quantity | Formula |
|---|---|
| PDF | $f(\mathbf x)=(2\pi)^{-p/2}|\boldsymbol\Sigma|^{-1/2}\exp\left\{-\tfrac12(\mathbf x-\boldsymbol\mu)^\top\boldsymbol\Sigma^{-1}(\mathbf x-\boldsymbol\mu)\right\}$ |
| MGF | $\exp\{\mathbf t^\top\boldsymbol\mu+\tfrac12\mathbf t^\top\boldsymbol\Sigma\mathbf t\}$ |
| Mean | $\boldsymbol\mu$ |
| Covariance | $\boldsymbol\Sigma$ |

**Derivation.** Let $\mathbf Z$ be a vector of $p$ iid $N(0,1)$ variables, with joint density $(2\pi)^{-p/2}e^{-\mathbf z^\top\mathbf z/2}$. Set $\mathbf X=\boldsymbol\mu+\mathbf A\mathbf Z$ with $\mathbf A\mathbf A^\top=\boldsymbol\Sigma$ (e.g., the Cholesky factor). The change of variables $\mathbf z=\mathbf A^{-1}(\mathbf x-\boldsymbol\mu)$ has Jacobian $|\mathbf A|^{-1}=|\boldsymbol\Sigma|^{-1/2}$, and $\mathbf z^\top\mathbf z=(\mathbf x-\boldsymbol\mu)^\top\boldsymbol\Sigma^{-1}(\mathbf x-\boldsymbol\mu)$.

**Notable facts.**
- $\mathbf X$ is MVN if and only if every linear combination $\mathbf a^\top\mathbf X$ is univariate normal.
- Affine transformations stay MVN: $\mathbf B\mathbf X+\mathbf c\sim N(\mathbf B\boldsymbol\mu+\mathbf c,\mathbf B\boldsymbol\Sigma\mathbf B^\top)$.
- For a **jointly** normal vector, uncorrelated components are independent. Marginally normal variables alone are not enough.
- **Conditional distribution:** partition $\mathbf X=(\mathbf X_1,\mathbf X_2)$. Then
$$\mathbf X_1\mid\mathbf X_2=\mathbf x_2\sim N\left(\boldsymbol\mu_1+\boldsymbol\Sigma_{12}\boldsymbol\Sigma_{22}^{-1}(\mathbf x_2-\boldsymbol\mu_2),\ \boldsymbol\Sigma_{11}-\boldsymbol\Sigma_{12}\boldsymbol\Sigma_{22}^{-1}\boldsymbol\Sigma_{21}\right).$$
- In the bivariate case, $E[Y\mid X=x]=\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X)$ and $\text{Var}(Y\mid X)=\sigma_Y^2(1-\rho^2)$.
- The squared Mahalanobis distance satisfies $(\mathbf X-\boldsymbol\mu)^\top\boldsymbol\Sigma^{-1}(\mathbf X-\boldsymbol\mu)\sim\chi^2_p$.

**R:** `mvtnorm::dmvnorm(x, mean, sigma)`, `MASS::mvrnorm(n, mu, Sigma)`

[Back to top](#table-of-contents)

---

# Reference

## Relationships Between Distributions

| From | Operation / limit | Result |
|---|---|---|
| $n$ iid Bernoulli($p$) | sum | Binomial($n,p$) |
| Binomial($n,p$) | $n\to\infty$, $np\to\lambda$ | Poisson($\lambda$) |
| Binomial($n,p$) | $n$ large (CLT) | $\approx N(np,npq)$ |
| Hypergeometric | $N\to\infty$, $K/N\to p$ | Binomial($n,p$) |
| $r$ iid Geometric($p$) | sum | Negative Binomial($r,p$) |
| Poisson($\Lambda$), $\Lambda\sim$ Gamma | mixture | Negative Binomial |
| Poisson process | interarrival / $k$-th arrival | Exponential / Gamma($k,\lambda$) |
| $n$ iid Exponential($\lambda$) | sum | Gamma($n,\lambda$) |
| Exponential | $\min$ of independent | Exponential($\sum\lambda_i$) |
| Gamma($k/2$, rate $1/2$) | special case | $\chi^2_k$ |
| $Z_1,\dots,Z_k$ iid $N(0,1)$ | $\sum Z_i^2$ | $\chi^2_k$ |
| $Z$, $V\sim\chi^2_\nu$ independent | $Z/\sqrt{V/\nu}$ | $t_\nu$ |
| $\chi^2_{d_1}$, $\chi^2_{d_2}$ independent | ratio of $\chi^2/\text{df}$ | $F_{d_1,d_2}$ |
| $t_\nu$ | square | $F_{1,\nu}$ |
| $t_\nu$ | $\nu=1$ / $\nu\to\infty$ | Cauchy / $N(0,1)$ |
| Gamma($a$), Gamma($b$), same rate | $X/(X+Y)$ | Beta($a,b$) |
| $U(0,1)$ | $k$-th order statistic of $n$ | Beta($k,n-k+1$) |
| $U(0,1)$ | $-\log U/\lambda$ | Exponential($\lambda$) |
| Any continuous $X$ | $F_X(X)$ | $U(0,1)$ |
| $N(\mu,\sigma^2)$ | $e^X$ | Lognormal($\mu,\sigma^2$) |
| Exponential(1) | $\lambda E^{1/k}$ | Weibull($k,\lambda$) |
| Two iid Exponentials | difference | Laplace |
| $Z_1/Z_2$, iid $N(0,1)$ | ratio | Cauchy(0,1) |

**Conjugate pairs (Bayesian).** Beta–Binomial/Bernoulli/Geometric, Gamma–Poisson/Exponential, Normal–Normal (known variance), Dirichlet–Multinomial, Inverse-Gamma–Normal variance.

[Back to top](#table-of-contents)

---

## Summary Table

| Distribution | Parameters | Mean | Variance | MGF |
|---|---|---|---|---|
| [Bernoulli](#bernoulli) | $p$ | $p$ | $pq$ | $q+pe^t$ |
| [Binomial](#binomial) | $n,p$ | $np$ | $npq$ | $(q+pe^t)^n$ |
| [Geometric](#geometric) (failures) | $p$ | $q/p$ | $q/p^2$ | $\frac{p}{1-qe^t}$ |
| [Negative Binomial](#negative-binomial) | $r,p$ | $rq/p$ | $rq/p^2$ | $\left(\frac{p}{1-qe^t}\right)^r$ |
| [Hypergeometric](#hypergeometric) | $N,K,n$ | $nK/N$ | $n\frac KN\frac{N-K}{N}\frac{N-n}{N-1}$ | — |
| [Poisson](#poisson) | $\lambda$ | $\lambda$ | $\lambda$ | $e^{\lambda(e^t-1)}$ |
| [Discrete Uniform](#discrete-uniform) | $a,b$ | $\frac{a+b}{2}$ | $\frac{n^2-1}{12}$ | see section |
| [Uniform](#continuous-uniform) | $a,b$ | $\frac{a+b}{2}$ | $\frac{(b-a)^2}{12}$ | $\frac{e^{tb}-e^{ta}}{t(b-a)}$ |
| [Normal](#normal-gaussian) | $\mu,\sigma^2$ | $\mu$ | $\sigma^2$ | $e^{\mu t+\sigma^2t^2/2}$ |
| [Exponential](#exponential) | $\lambda$ (rate) | $1/\lambda$ | $1/\lambda^2$ | $\frac{\lambda}{\lambda-t}$ |
| [Gamma](#gamma) | $\alpha,\beta$ (rate) | $\alpha/\beta$ | $\alpha/\beta^2$ | $\left(\frac{\beta}{\beta-t}\right)^\alpha$ |
| [Beta](#beta) | $\alpha,\beta$ | $\frac{\alpha}{\alpha+\beta}$ | $\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$ | — |
| [Chi-Square](#chi-square) | $k$ | $k$ | $2k$ | $(1-2t)^{-k/2}$ |
| [Student's t](#students-t) | $\nu$ | $0$ ($\nu>1$) | $\frac{\nu}{\nu-2}$ ($\nu>2$) | none |
| [F](#f-snedecor) | $d_1,d_2$ | $\frac{d_2}{d_2-2}$ | see section | none |
| [Lognormal](#lognormal) | $\mu,\sigma^2$ | $e^{\mu+\sigma^2/2}$ | $(e^{\sigma^2}-1)e^{2\mu+\sigma^2}$ | none |
| [Weibull](#weibull) | $k,\lambda$ | $\lambda\Gamma(1+\frac1k)$ | $\lambda^2[\Gamma(1+\frac2k)-\Gamma^2(1+\frac1k)]$ | — |
| [Cauchy](#cauchy) | $x_0,\gamma$ | undefined | undefined | none |
| [Laplace](#laplace-double-exponential) | $\mu,b$ | $\mu$ | $2b^2$ | $\frac{e^{\mu t}}{1-b^2t^2}$ |
| [Multivariate Normal](#multivariate-normal) | $\boldsymbol\mu,\boldsymbol\Sigma$ | $\boldsymbol\mu$ | $\boldsymbol\Sigma$ | $e^{\mathbf t^\top\boldsymbol\mu+\frac12\mathbf t^\top\boldsymbol\Sigma\mathbf t}$ |

The SD is always $\sqrt{\text{Variance}}$. "—" means no convenient closed form; "none" means the MGF does not exist.

[Back to top](#table-of-contents)
