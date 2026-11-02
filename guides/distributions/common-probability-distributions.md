# Common Probability Distributions

A graduate-level quick reference with formulas, parameterization warnings, constructions, and expandable derivations.

## Table of contents

1. [Notation and calculation toolkit](#notation-and-calculation-toolkit)
2. [Bernoulli](#bernoulli)
3. [Binomial](#binomial)
4. [Categorical](#categorical)
5. [Discrete uniform](#discrete-uniform)
6. [Geometric](#geometric)
7. [Negative binomial](#negative-binomial)
8. [Hypergeometric](#hypergeometric)
9. [Poisson](#poisson)
10. [Beta-binomial](#beta-binomial)
11. [Continuous uniform](#continuous-uniform)
12. [Normal Gaussian](#normal-gaussian)
13. [Exponential](#exponential)
14. [Gamma and Erlang](#gamma-and-erlang)
15. [Beta](#beta)
16. [Chi-squared](#chi-squared)
17. [Student t](#student-t)
18. [F](#f)
19. [Lognormal](#lognormal)
20. [Weibull](#weibull)
21. [Rayleigh](#rayleigh)
22. [Laplace](#laplace)
23. [Cauchy](#cauchy)
24. [Pareto Type I](#pareto-type-i)
25. [Inverse gamma](#inverse-gamma)
26. [Multinomial](#multinomial)
27. [Dirichlet](#dirichlet)
28. [Multivariate normal](#multivariate-normal)
29. [Relationships and approximation map](#relationships-and-approximation-map)
30. [Software parameterization and numerical checks](#software-parameterization-and-numerical-checks)
31. [Source notes and further reading](#source-notes-and-further-reading)

---

## Notation and calculation toolkit

This guide emphasizes distributions commonly encountered in graduate probability, mathematical statistics, and introductory Bayesian inference. It is not an exhaustive catalog. All logarithms are natural. Unless stated otherwise, variables appearing in sums, ratios, and constructions are independent. Densities and mass functions are zero outside their stated supports.

Parameter conventions:
- Gamma and exponential use a rate, not a scale. Scale = 1 / rate.
- Normal uses variance as its second mathematical parameter; software often expects standard deviation.
- Geometric counts trials through the first success; negative binomial counts failures before a specified number of successes.
- Pareto means Type I; Weibull means the two-parameter minimum/lifetime family.
- Standard deviation (SD) is the nonnegative square root of variance, when variance exists. A scale parameter is not automatically an SD.

### What is being calculated?
A PMF assigns a probability to a discrete value. A PDF assigns density, not point probability. For an absolutely continuous variable, any single point has probability zero. The CDF works for either type:

```math
p_X(k)=P(X=k),\qquad F_X(x)=P(X\le x).
```

```math
P(a<X\le b)=F_X(b)-F_X(a),\qquad
P(a\le X\le b)=F_X(b)-F_X(a^-).
```

For integer-valued X and integer a, b:

```math
P(a\le X\le b)=F_X(b)-F_X(a-1).
```

For a PDF, calculate interval probabilities by integration. For a PMF, sum over the appropriate support points. CDF formulas below are given on the support; extend them to 0 below the support and 1 above a finite upper endpoint.

### Moments and transforms

```math
E[g(X)]=\sum_x g(x)p_X(x)\quad\text{or}\quad
E[g(X)]=\int g(x)f_X(x)\,dx.
```

```math
\operatorname{Var}(X)=E[X^2]-E[X]^2,\qquad
\operatorname{SD}(X)=\sqrt{\operatorname{Var}(X)}.
```

```math
M_X(t)=E[e^{tX}],\quad K_X(t)=\log M_X(t),\quad
E[X]=M_X'(0),\quad \operatorname{Var}(X)=K_X''(0).
```

The MGF method requires finiteness on an open interval containing zero. A characteristic function always exists; an MGF need not. Expressions below identify the MGF's finite domain when it is restricted. Multivariate analogues use a vector argument.

```math
\varphi_X(t)=E[e^{itX}],\qquad
G_X(z)=E[z^X]\quad\text{for nonnegative integer-valued }X.
```

### Transformations, sums, and mixtures
For a differentiable one-to-one transformation Y = g(X):

```math
f_Y(y)=f_X(g^{-1}(y))\left|\frac{d}{dy}g^{-1}(y)\right|.
```

For a many-to-one transformation, sum these contributions over all inverse branches. In multiple dimensions use the absolute determinant of the inverse Jacobian. For independent sums:

```math
f_{X+Y}(s)=\int f_X(x)f_Y(s-x)\,dx,\qquad
M_{X+Y}(t)=M_X(t)M_Y(t).
```

For mixtures, average the conditional distribution, not just its parameters:

```math
f_X(x)=\int f_{X\mid\Theta}(x\mid\theta)f_\Theta(\theta)\,d\theta.
```

```math
E[X]=E[E(X\mid\Theta)],\quad
\operatorname{Var}(X)=E[\operatorname{Var}(X\mid\Theta)]+\operatorname{Var}(E[X\mid\Theta]).
```

### Special functions

```math
\Gamma(a)=\int_0^\infty u^{a-1}e^{-u}\,du,\quad
\Gamma(a+1)=a\Gamma(a),\quad \Gamma(n)=(n-1)!.
```

```math
B(a,b)=\int_0^1u^{a-1}(1-u)^{b-1}\,du
=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}.
```

```math
\gamma(a,x)=\int_0^x u^{a-1}e^{-u}\,du,\qquad
I_x(a,b)=\frac{1}{B(a,b)}\int_0^x u^{a-1}(1-u)^{b-1}\,du.
```

Phi denotes the standard normal CDF. These special functions are exact definitions, not approximations. In practice use library CDFs rather than numerical integration written from scratch.

The derivations in each section explain a stochastic construction or normalization and, where useful, a moment calculation. They do not claim that a distribution has a unique derivation.

[Back to table of contents](#table-of-contents)

---

## Bernoulli

Use: one binary outcome; an indicator of an event. Parameters: success probability p in [0, 1]. Support: {0, 1}.

```math
X\sim\operatorname{Bernoulli}(p),\qquad p_X(x)=p^x(1-p)^{1-x}.
```

At p = 0 or 1 interpret the law by its point masses, avoiding ambiguous powers.

```math
E[X]=p,\quad \operatorname{Var}(X)=p(1-p),\quad
\operatorname{SD}(X)=\sqrt{p(1-p)},\quad M_X(t)=1-p+pe^t.
```

CDF: 0 for x < 0; 1 − p for 0 ≤ x < 1; 1 for x ≥ 1. A categorical distribution with two categories is Bernoulli after coding one category as 1.

<details>
<summary>Derivation: indicators and moments</summary>

Let X be 1 when event A occurs and 0 otherwise, with P(A) = p. Its two masses sum to 1. Since X squared equals X:

```math
E[X]=0(1-p)+1p=p,\qquad E[X^2]-E[X]^2=p-p^2.
```

This indicator argument is often the quickest way to derive moments of counting distributions.

</details>

[Back to table of contents](#table-of-contents)

---

## Binomial

Use: successes in n independent trials with the same success probability p. Parameters: nonnegative integer n; p in [0, 1]. Support: integers 0 through n.

```math
X\sim\operatorname{Bin}(n,p),\qquad
p_X(k)=\binom nk p^k(1-p)^{n-k}.
```

```math
F_X(x)=\sum_{k=0}^{\min(n,\lfloor x\rfloor)}\binom nk p^k(1-p)^{n-k}\quad(x\ge0).
```

```math
E[X]=np,\quad \operatorname{Var}(X)=np(1-p),\quad
\operatorname{SD}(X)=\sqrt{np(1-p)},\quad
M_X(t)=(1-p+pe^t)^n.
```

Independent binomials with the same p add by adding their trial counts. Unequal success probabilities give a Poisson-binomial distribution, not generally a binomial. Sampling without replacement calls for a hypergeometric model instead.

<details>
<summary>Derivation: count arrangements</summary>

A specified sequence containing k successes has probability equal to p to the k times (1 − p) to the n − k. There are n choose k such sequences. Summing the PMF gives 1 by the binomial theorem.

Alternatively:

```math
X=\sum_{i=1}^n B_i,\qquad B_i\overset{\mathrm{iid}}\sim\operatorname{Bernoulli}(p).
```

Linearity gives the mean; independence gives the variance and the product MGF.

</details>

[Back to table of contents](#table-of-contents)

---

## Categorical

Use: one outcome in K mutually exclusive categories. Parameters: probability vector with nonnegative entries summing to 1. Support: category labels 1 through K.

```math
P(X=j)=p_j,\qquad \sum_{j=1}^K p_j=1.
```

If numeric labels are genuinely meaningful:

```math
E[X]=\sum_j jp_j,\quad \operatorname{Var}(X)=\sum_j j^2p_j-\left(\sum_j jp_j\right)^2,\quad
\operatorname{SD}(X)=\sqrt{\operatorname{Var}(X)}.
```

For nominal categories these moments depend on arbitrary coding and are usually uninformative. Use the one-hot indicator vector Y instead:

```math
E[Y_j]=p_j,\quad \operatorname{Var}(Y_j)=p_j(1-p_j),\quad
\operatorname{Cov}(Y_i,Y_j)=-p_ip_j\quad(i\ne j).
```

<details>
<summary>Derivation: a partition of the sample space</summary>

The K possible outcomes form a disjoint partition, so their masses sum to 1. The indicator Yj equals 1 only for category j. For different categories, their indicators cannot both equal 1, making their product zero. Subtracting the product of their means gives the negative covariance. This is a multinomial experiment with n = 1.

</details>

[Back to table of contents](#table-of-contents)

---

## Discrete uniform

Use: equally likely integers. Parameters: integers a ≤ b; let m = b − a + 1. Support: integers a through b.

```math
p_X(k)=\frac1m,\qquad
F_X(x)=\frac{\lfloor x\rfloor-a+1}{m}\quad(a\le x<b).
```

```math
E[X]=\frac{a+b}{2},\quad \operatorname{Var}(X)=\frac{m^2-1}{12},\quad
\operatorname{SD}(X)=\sqrt{\frac{m^2-1}{12}}.
```

```math
M_X(t)=\frac1m\sum_{k=a}^b e^{tk}\qquad(t\in\mathbb R).
```

Do not use the continuous-uniform variance: integer endpoints include both a and b, hence the extra 1 in m.

<details>
<summary>Derivation: normalization and finite sums</summary>

Equal mass c at each of m points requires mc = 1. Shift to J = X − a, uniform on 0 through m − 1. Apply:

```math
\sum_{j=0}^{m-1}j=\frac{m(m-1)}2,\qquad
\sum_{j=0}^{m-1}j^2=\frac{m(m-1)(2m-1)}6.
```

Divide by m and subtract the squared mean to obtain the variance.

</details>

[Back to table of contents](#table-of-contents)

---

## Geometric

Use: number of trials through the first success. Parameter: 0 < p ≤ 1. Support: positive integers.

```math
p_X(k)=(1-p)^{k-1}p,\quad
F_X(x)=1-(1-p)^{\lfloor x\rfloor}\quad(x\ge1).
```

```math
E[X]=\frac1p,\quad \operatorname{Var}(X)=\frac{1-p}{p^2},\quad
\operatorname{SD}(X)=\frac{\sqrt{1-p}}p.
```

```math
M_X(t)=\frac{pe^t}{1-(1-p)e^t},\qquad t<-\log(1-p)\quad(0<p<1).
```

At p = 1, X = 1 and the MGF is finite for every real t. Memorylessness holds for nonnegative integers m, n:

```math
P(X>m+n\mid X>m)=P(X>n).
```

Alternative convention: Y = X − 1 counts failures. Its mean is (1 − p)/p, with the same variance. Always check which convention a textbook or software uses.

<details>
<summary>Derivation: a string of failures</summary>

For X = k, the first k − 1 trials must fail and the kth must succeed. Multiply independent probabilities. Normalization and the MGF follow from the geometric series:

```math
\sum_{j=0}^\infty q^j=\frac1{1-q}\quad(|q|<1).
```

Differentiate the MGF or the geometric series for the moments. The survival probability is (1 − p) to the m, which gives memorylessness by cancellation.

</details>

[Back to table of contents](#table-of-contents)

---

## Negative binomial

Use: failures before the rth success; also overdispersed counts. Parameters: r > 0, 0 < p < 1. The waiting-time interpretation requires integer r. Support: nonnegative integers.

```math
p_X(k)=\frac{\Gamma(k+r)}{\Gamma(r)k!}p^r(1-p)^k.
```

```math
F_X(x)=I_p(r,\lfloor x\rfloor+1)\quad(x\ge0).
```

```math
E[X]=\frac{r(1-p)}p,\quad \operatorname{Var}(X)=\frac{r(1-p)}{p^2},\quad
\operatorname{SD}(X)=\frac{\sqrt{r(1-p)}}p.
```

```math
M_X(t)=\left(\frac{p}{1-(1-p)e^t}\right)^r,\qquad t<-\log(1-p).
```

With mean parameter m: p = r/(r + m), variance = m + m²/r. For integer r, total trials T = X + r has mean r/p and the same variance. Geometric failure counts are the r = 1 case. The p = 1 limit is a point mass at zero.

<details>
<summary>Derivation: arrangements and a Poisson–gamma mixture</summary>

For integer r, the last trial is a success. Among the first k + r − 1 trials choose positions for r − 1 successes; multiply by the sequence probability:

```math
P(X=k)=\binom{k+r-1}{r-1}p^r(1-p)^k.
```

For general positive r, mix Poisson counts over a gamma intensity:

```math
X\mid\Lambda\sim\operatorname{Pois}(\Lambda),\quad
\Lambda\sim\operatorname{Gamma}\left(r,\text{rate}=\frac p{1-p}\right).
```

Multiplying the two densities and integrating the intensity uses the gamma integral and produces the stated PMF. Total expectation and total variance give the moments and explain the extra-Poisson variation.

</details>

[Back to table of contents](#table-of-contents)

---

## Hypergeometric

Use: successes in n draws without replacement from N objects, K of which are successes. Parameters: integers N ≥ 1, 0 ≤ K ≤ N, 0 ≤ n ≤ N.

```math
\max(0,n-N+K)\le k\le\min(n,K),\qquad
p_X(k)=\frac{\binom Kk\binom{N-K}{n-k}}{\binom Nn}.
```

CDF: sum this PMF over supported integers no larger than x. Put p = K/N. For N > 1:

```math
E[X]=np,\quad \operatorname{Var}(X)=np(1-p)\frac{N-n}{N-1},\quad
\operatorname{SD}(X)=\sqrt{np(1-p)\frac{N-n}{N-1}}.
```

For N = 1 the variable is deterministic and has zero variance. The factor (N − n)/(N − 1) is the finite-population correction. Binomial approximation becomes appropriate when the sampling fraction is small, with adequacy depending on the required accuracy.

<details>
<summary>Derivation: counting subsets and dependent indicators</summary>

There are N choose n equally likely unordered samples. A sample with k successes chooses k objects from the K successes and n − k from the remaining objects. Divide favorable samples by total samples.

For draw indicators Bi:

```math
E[B_i]=p,\quad \operatorname{Cov}(B_i,B_j)=-\frac{p(1-p)}{N-1}\quad(i\ne j).
```

Summing variances plus all pairwise covariances yields the finite-population correction. Unlike binomial trials, the draws are not independent.

</details>

[Back to table of contents](#table-of-contents)

---

## Poisson

Use: event counts in a fixed exposure under a homogeneous Poisson-process model. Parameter: intensity over that exposure λ > 0. Support: nonnegative integers. The λ = 0 limit is deterministic zero.

```math
p_X(k)=e^{-\lambda}\frac{\lambda^k}{k!},\quad
F_X(x)=e^{-\lambda}\sum_{k=0}^{\lfloor x\rfloor}\frac{\lambda^k}{k!}\quad(x\ge0).
```

```math
E[X]=\lambda,\quad \operatorname{Var}(X)=\lambda,\quad
\operatorname{SD}(X)=\sqrt\lambda,\quad
M_X(t)=\exp\{\lambda(e^t-1)\}.
```

Independent Poisson counts add by adding their intensities. Independent retention of each event with probability p gives a Poisson count with intensity pλ. For event rate ρ and exposure t, use λ = ρt. Independence of increments and constant rate are model assumptions, not consequences of observing count data.

<details>
<summary>Derivation: rare-event limit</summary>

Let Xn be binomial with n trials and success probability λ/n. For fixed k:

```math
P(X_n=k)=\frac{n(n-1)\cdots(n-k+1)}{n^k}\frac{\lambda^k}{k!}
\left(1-\frac\lambda n\right)^{n-k}
\longrightarrow e^{-\lambda}\frac{\lambda^k}{k!}.
```

The exponential series normalizes the PMF. Differentiating the log MGF gives both mean and variance λ.

</details>

[Back to table of contents](#table-of-contents)

---

## Beta-binomial

Use: successes in n trials with a shared random success probability; predictive counts under a beta prior. Parameters: integer n ≥ 0; α, β > 0. Support: integers 0 through n.

```math
P\sim\operatorname{Beta}(\alpha,\beta),\quad X\mid P\sim\operatorname{Bin}(n,P),\quad
p_X(k)=\binom nk\frac{B(k+\alpha,n-k+\beta)}{B(\alpha,\beta)}.
```

```math
E[X]=\frac{n\alpha}{\alpha+\beta},\quad
\operatorname{Var}(X)=\frac{n\alpha\beta(\alpha+\beta+n)}{(\alpha+\beta)^2(\alpha+\beta+1)},\quad
\operatorname{SD}(X)=\sqrt{\operatorname{Var}(X)}.
```

CDF: finite sum of the PMF. Given P the trials are independent; marginally they are positively correlated. The variance exceeds the binomial variance with the same n and mean when n > 1. At α = β = 1 the count is discrete uniform on 0 through n.

<details>
<summary>Derivation: integrate the success probability</summary>

```math
P(X=k)=\int_0^1\binom nk p^k(1-p)^{n-k}
\frac{p^{\alpha-1}(1-p)^{\beta-1}}{B(\alpha,\beta)}\,dp.
```

Combine exponents and recognize the beta integral. For the moments:

```math
E[X]=nE[P],\qquad
\operatorname{Var}(X)=nE[P(1-P)]+n^2\operatorname{Var}(P).
```

Insert the beta moments. A single shared P is essential; drawing an independent new P for every trial gives ordinary independent Bernoulli outcomes after marginalization.

</details>

[Back to table of contents](#table-of-contents)

---

## Continuous uniform

Use: a constant density over a bounded interval. Parameters: real a < b. Support: [a, b].

```math
f_X(x)=\frac1{b-a},\quad F_X(x)=\frac{x-a}{b-a}\quad(a\le x\le b).
```

```math
E[X]=\frac{a+b}{2},\quad \operatorname{Var}(X)=\frac{(b-a)^2}{12},\quad
\operatorname{SD}(X)=\frac{b-a}{\sqrt{12}}.
```

```math
M_X(t)=\frac{e^{tb}-e^{ta}}{t(b-a)}\quad(t\ne0),\qquad M_X(0)=1.
```

Quantile: a + (b − a)u for 0 < u < 1. If a variable has a continuous CDF F, then F(X) is uniform on (0, 1). Conversely, applying the generalized inverse CDF to a uniform draw generates the target law.

<details>
<summary>Derivation: constant density and integration</summary>

A constant c over an interval of length b − a must satisfy c(b − a) = 1. Integrating from a to x gives the CDF. Raw moments are:

```math
E[X^r]=\frac{b^{r+1}-a^{r+1}}{(r+1)(b-a)}\quad(r=0,1,2,\ldots).
```

Use r = 1 and 2 for mean and variance.

</details>

[Back to table of contents](#table-of-contents)

---

## Normal Gaussian

Use: additive measurements and errors; exact normal models and central-limit approximations. Parameters: real μ, variance σ² > 0. Support: all real numbers.

```math
X\sim N(\mu,\sigma^2),\qquad
f_X(x)=\frac1{\sigma\sqrt{2\pi}}\exp\left[-\frac{(x-\mu)^2}{2\sigma^2}\right].
```

```math
F_X(x)=\Phi\left(\frac{x-\mu}{\sigma}\right),\quad
E[X]=\mu,\quad \operatorname{Var}(X)=\sigma^2,\quad \operatorname{SD}(X)=\sigma.
```

```math
M_X(t)=e^{\mu t+\sigma^2t^2/2},\qquad
\frac{X-\mu}{\sigma}\sim N(0,1).
```

Mean = median = mode. Affine transformations and sums of independent normal variables remain normal. A central limit theorem describes a limiting standardized sum under assumptions; it does not imply every observed variable is normal.

<details>
<summary>Derivation: Gaussian normalization and completing the square</summary>

Let I be the integral of exp(−x²/2) over the real line. Square the integral and change to polar coordinates:

```math
I^2=\int_0^{2\pi}\int_0^\infty e^{-r^2/2}r\,dr\,d\theta=2\pi.
```

Thus I = square root of 2π, giving the normalizing constant. Shift and scale the standard density to obtain the stated PDF. For the MGF, complete the square:

```math
tx-\frac{(x-\mu)^2}{2\sigma^2}
=-\frac{(x-\mu-\sigma^2t)^2}{2\sigma^2}+\mu t+\frac{\sigma^2t^2}{2}.
```

The remaining normal density integrates to 1. Derivatives of the log MGF give the moments.

</details>

[Back to table of contents](#table-of-contents)

---

## Exponential

Use: waiting time to the next event in a homogeneous Poisson process; constant-hazard lifetimes. Parameter: rate λ > 0. Support: x ≥ 0.

```math
f_X(x)=\lambda e^{-\lambda x},\quad F_X(x)=1-e^{-\lambda x},\quad
S_X(x)=e^{-\lambda x}.
```

```math
E[X]=\frac1\lambda,\quad \operatorname{Var}(X)=\frac1{\lambda^2},\quad
\operatorname{SD}(X)=\frac1\lambda,\quad
M_X(t)=\frac\lambda{\lambda-t}\quad(t<\lambda).
```

Hazard = λ; median = log(2)/λ; quantile = −log(1 − u)/λ. Memorylessness:

```math
P(X>s+t\mid X>s)=P(X>t)\quad(s,t\ge0).
```

Scale convention: θ = 1/λ. Multiplying a waiting time by c > 0 divides its rate by c.

<details>
<summary>Derivation: no events before time x</summary>

If N(x) is Poisson with mean λx and X is the first arrival time:

```math
P(X>x)=P(N(x)=0)=e^{-\lambda x}.
```

Differentiate 1 minus this survival function to obtain the PDF. Integrating exp(tx) times the PDF gives the MGF when t < λ. The constant hazard also follows from density divided by survival.

</details>

[Back to table of contents](#table-of-contents)

---

## Gamma and Erlang

Use: positive waiting times, sums of exponentials, and positive latent rates. Parameters: shape α > 0 and rate β > 0. Support: x > 0.

```math
X\sim\operatorname{Gamma}(\alpha,\text{rate}=\beta),\quad
f_X(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x},\quad
F_X(x)=\frac{\gamma(\alpha,\beta x)}{\Gamma(\alpha)}.
```

```math
E[X]=\frac\alpha\beta,\quad \operatorname{Var}(X)=\frac\alpha{\beta^2},\quad
\operatorname{SD}(X)=\frac{\sqrt\alpha}\beta,\quad
M_X(t)=\left(\frac\beta{\beta-t}\right)^\alpha\quad(t<\beta).
```

Mode = (α − 1)/β for α > 1. For α ≤ 1 the density decreases from its boundary (unbounded there if α < 1). Integer shape gives Erlang; shape 1 gives exponential. Independent gammas with the same rate add by adding shapes. With scale θ = 1/β, mean = αθ and variance = αθ².

<details>
<summary>Derivation: waiting times and the gamma integral</summary>

For integer α = n, the waiting time Tn to the nth Poisson-process event exceeds x exactly when fewer than n events have arrived:

```math
P(T_n>x)=e^{-\beta x}\sum_{j=0}^{n-1}\frac{(\beta x)^j}{j!}.
```

Differentiate to obtain the Erlang density. For noninteger shape, normalize the kernel using substitution u = βx and the gamma integral. Raw moments follow by the same substitution:

```math
E[X^q]=\frac{\Gamma(\alpha+q)}{\Gamma(\alpha)\beta^q}\quad(q>-\alpha).
```

The cases q = 1 and 2 yield the mean and variance.

</details>

[Back to table of contents](#table-of-contents)

---

## Beta

Use: probabilities and proportions bounded between 0 and 1. Parameters: α, β > 0. Support: 0 < x < 1.

```math
f_X(x)=\frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha,\beta)},\qquad
F_X(x)=I_x(\alpha,\beta).
```

```math
E[X]=\frac\alpha{\alpha+\beta},\quad
\operatorname{Var}(X)=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)},\quad
\operatorname{SD}(X)=\sqrt{\operatorname{Var}(X)}.
```

Interior mode = (α − 1)/(α + β − 2) when both parameters exceed 1. Both below 1 give a U-shaped density; both equal 1 give uniform. Swapping the parameters gives the law of 1 − X. The MGF is finite for all real arguments because the support is bounded.

Beta–binomial conjugacy: after s successes and f failures, a Beta(α, β) prior becomes Beta(α + s, β + f). This updates a distribution on the probability, not a binomial distribution on counts.

<details>
<summary>Derivation: normalized gamma variables and beta moments</summary>

Let U and V be independent gammas with shapes α and β and a common rate c. Transform to X = U/(U + V), T = U + V. The inverse is U = XT, V = (1 − X)T and its absolute Jacobian is T. The joint density factors into beta in X and gamma in T. Thus X and T are independent and X has the stated law.

Alternatively the beta integral directly normalizes the density. Moments are:

```math
E[X^q]=\frac{B(\alpha+q,\beta)}{B(\alpha,\beta)}\quad(q>-\alpha).
```

Insert q = 1 and 2. Conjugacy follows by multiplying the beta kernel by the binomial likelihood and collecting exponents.

</details>

[Back to table of contents](#table-of-contents)

---

## Chi-squared

Use: normal-sample variance, sums of squared standard normals, and classical test statistics. Parameter: degrees of freedom ν > 0. Support: x > 0.

```math
X\sim\chi^2_\nu=\operatorname{Gamma}(\nu/2,\text{rate}=1/2),\quad
f_X(x)=\frac{x^{\nu/2-1}e^{-x/2}}{2^{\nu/2}\Gamma(\nu/2)}.
```

```math
F_X(x)=\frac{\gamma(\nu/2,x/2)}{\Gamma(\nu/2)},\quad
E[X]=\nu,\quad \operatorname{Var}(X)=2\nu,\quad \operatorname{SD}(X)=\sqrt{2\nu}.
```

```math
M_X(t)=(1-2t)^{-\nu/2}\quad(t<1/2).
```

Independent chi-squares add their degrees of freedom. For an iid normal sample of size n ≥ 2, with unbiased sample variance S²:

```math
\frac{(n-1)S^2}{\sigma^2}\sim\chi^2_{n-1}.
```

This exact sampling law requires normal observations; it is not distribution-free.

<details>
<summary>Derivation: square a standard normal</summary>

For Z standard normal and Y = Z², the inverse branches are plus and minus square root of y:

```math
f_Y(y)=\frac{\phi(\sqrt y)+\phi(-\sqrt y)}{2\sqrt y}
=\frac{y^{-1/2}e^{-y/2}}{\sqrt{2\pi}}.
```

This is Gamma(1/2, rate 1/2). Summing ν independent copies for positive integer ν gives Gamma(ν/2, rate 1/2). The gamma definition extends degrees of freedom to every positive real value.

</details>

[Back to table of contents](#table-of-contents)

---

## Student t

Use: inference on a normal mean with unknown variance; symmetric heavy-tailed models. Parameter: degrees of freedom ν > 0. Support: all real numbers. This section gives the standard, centered, unit-scale law.

```math
f_T(t)=\frac{\Gamma((\nu+1)/2)}{\sqrt{\nu\pi}\,\Gamma(\nu/2)}
\left(1+\frac{t^2}\nu\right)^{-(\nu+1)/2}.
```

```math
E[T]=0\quad(\nu>1),\qquad
\operatorname{Var}(T)=\frac\nu{\nu-2},\quad
\operatorname{SD}(T)=\sqrt{\frac\nu{\nu-2}}\quad(\nu>2).
```

For ν ≤ 1 the mean is undefined; for 1 < ν ≤ 2 the variance is infinite. No finite variance or SD exists for ν ≤ 2. The MGF is infinite at every nonzero real argument. For t ≥ 0:

```math
F_T(t)=1-\frac12 I_{\nu/(\nu+t^2)}(\nu/2,1/2),\qquad F_T(-t)=1-F_T(t).
```

Median and mode = 0. At ν = 1 this is standard Cauchy; as ν grows it converges to standard normal. For X = μ + sT, s > 0, s is a scale, not its SD. For ν > 2 the SD is s times the standard-law SD.

<details>
<summary>Derivation: normal divided by an estimated scale</summary>

```math
T=\frac{Z}{\sqrt{V/\nu}},\qquad Z\sim N(0,1),\quad V\sim\chi^2_\nu,\quad Z\perp V.
```

Conditional on V = v, T is normal with variance ν/v. Multiply that conditional density by the chi-square density and integrate v. The gamma integral gives the stated PDF. For ν > 2:

```math
E[T^2]=E[Z^2]E[\nu/V]=\frac\nu{\nu-2}.
```

More generally, the absolute moment of order q ≥ 0 exists exactly when q < ν. For an iid normal sample:

```math
\frac{\bar X-\mu}{S/\sqrt n}\sim t_{n-1}.
```

Independence of the normal sample mean and sample variance is essential to this exact ratio construction.

</details>

[Back to table of contents](#table-of-contents)

---

## F

Use: ratios of independent variance estimates; normal-theory ANOVA and regression tests. Parameters: numerator degrees of freedom d1 > 0 and denominator degrees of freedom d2 > 0. Support: x > 0.

```math
f_X(x)=\frac{(d_1/d_2)^{d_1/2}x^{d_1/2-1}}
{B(d_1/2,d_2/2)(1+d_1x/d_2)^{(d_1+d_2)/2}},\quad
F_X(x)=I_{d_1x/(d_1x+d_2)}(d_1/2,d_2/2).
```

```math
E[X]=\frac{d_2}{d_2-2}\quad(d_2>2).
```

```math
\operatorname{Var}(X)=\frac{2d_2^2(d_1+d_2-2)}{d_1(d_2-2)^2(d_2-4)},\quad
\operatorname{SD}(X)=\sqrt{\operatorname{Var}(X)}\quad(d_2>4).
```

Mean is infinite when d2 ≤ 2; second moment is infinite when d2 ≤ 4. There is no finite SD for d2 ≤ 4. The MGF is finite for t ≤ 0, not for positive t, so no MGF exists on a neighborhood of zero.

Reciprocal: 1/X has F(d2, d1) law. Squared standard Student t with ν degrees of freedom has F(1, ν) law. Reversing numerator and denominator changes the distribution.

<details>
<summary>Derivation: ratio of independent chi-squares</summary>

```math
X=\frac{U/d_1}{V/d_2},\qquad U\sim\chi^2_{d_1},\quad V\sim\chi^2_{d_2},\quad U\perp V.
```

Transform from (U, V) to (X, V), with U = (d1/d2)XV and absolute Jacobian (d1/d2)V. Integrate V using the gamma integral. Moments follow from independence and gamma moments:

```math
E[X^q]=\left(\frac{d_2}{d_1}\right)^q
\frac{\Gamma(d_1/2+q)\Gamma(d_2/2-q)}{\Gamma(d_1/2)\Gamma(d_2/2)},\qquad
-d_1/2<q<d_2/2.
```

Use q = 1 and 2.

</details>

[Back to table of contents](#table-of-contents)

---

## Lognormal

Use: multiplicative effects and positive right-skewed quantities. Parameters: real μ and σ² > 0, which describe log X, not X. Support: x > 0.

```math
\log X\sim N(\mu,\sigma^2),\quad
f_X(x)=\frac1{x\sigma\sqrt{2\pi}}\exp\left[-\frac{(\log x-\mu)^2}{2\sigma^2}\right],\quad
F_X(x)=\Phi\left(\frac{\log x-\mu}{\sigma}\right).
```

```math
E[X]=e^{\mu+\sigma^2/2},\quad
\operatorname{Var}(X)=(e^{\sigma^2}-1)e^{2\mu+\sigma^2},\quad
\operatorname{SD}(X)=e^{\mu+\sigma^2/2}\sqrt{e^{\sigma^2}-1}.
```

Median = exp(μ); mode = exp(μ − σ²). All positive integer moments exist, but the MGF is infinite for t > 0; therefore it does not exist on a neighborhood of zero. Products of independent lognormals are lognormal; sums are not generally lognormal.

<details>
<summary>Derivation: exponentiate a normal variable</summary>

Let Y be normal and X = exp(Y). The inverse y = log x has derivative 1/x, which supplies the PDF's Jacobian factor. Raw moments use the normal MGF:

```math
E[X^q]=E[e^{qY}]=\exp(q\mu+q^2\sigma^2/2)\quad(q\in\mathbb R).
```

Apply q = 1 and 2 for the moments. For positive t, exp(tx) grows faster in the upper tail than the lognormal density decays, explaining why having every polynomial moment does not imply MGF existence.

</details>

[Back to table of contents](#table-of-contents)

---

## Weibull

Use: positive lifetimes with increasing or decreasing hazard. Parameters: shape k > 0 and scale η > 0. Support: x > 0.

```math
f_X(x)=\frac{k}{\eta}\left(\frac{x}{\eta}\right)^{k-1}e^{-(x/\eta)^k},\quad
F_X(x)=1-e^{-(x/\eta)^k}.
```

```math
E[X]=\eta\Gamma(1+1/k),\quad
\operatorname{Var}(X)=\eta^2\left[\Gamma(1+2/k)-\Gamma(1+1/k)^2\right],\quad
\operatorname{SD}(X)=\sqrt{\operatorname{Var}(X)}.
```

```math
h_X(x)=\frac{k}{\eta}(x/\eta)^{k-1},\qquad
Q(u)=\eta[-\log(1-u)]^{1/k}.
```

Hazard decreases for k < 1, is constant for k = 1, and increases for k > 1. The k = 1 case is exponential with rate 1/η. MGF: finite for all real t if k > 1; finite for t < 1/η if k = 1; finite only for t ≤ 0 if 0 < k < 1.

<details>
<summary>Derivation: a power-transformed exponential</summary>

Take Y exponential with rate 1 and set X = ηY to the power 1/k. Then:

```math
P(X>x)=P(Y>(x/\eta)^k)=e^{-(x/\eta)^k}.
```

Differentiate the CDF for the PDF. Since X to the q is η to the q times Y to the q/k:

```math
E[X^q]=\eta^q\Gamma(1+q/k)\quad(q>-k).
```

Insert q = 1 and 2. The hazard follows by dividing the PDF by survival.

</details>

[Back to table of contents](#table-of-contents)

---

## Rayleigh

Use: magnitude of a two-dimensional isotropic zero-mean Gaussian vector. Parameter: scale σ > 0. Support: x ≥ 0.

```math
f_X(x)=\frac{x}{\sigma^2}e^{-x^2/(2\sigma^2)},\qquad
F_X(x)=1-e^{-x^2/(2\sigma^2)}.
```

```math
E[X]=\sigma\sqrt{\pi/2},\quad
\operatorname{Var}(X)=\frac{4-\pi}{2}\sigma^2,\quad
\operatorname{SD}(X)=\sigma\sqrt{\frac{4-\pi}{2}}.
```

Mode = σ. Equivalent to Weibull shape 2 and scale σ√2. The MGF is finite for every real argument. The Rayleigh scale is the SD of either underlying normal component, not the SD of their magnitude.

<details>
<summary>Derivation: polar coordinates</summary>

For independent Z1, Z2 distributed N(0, σ²), let R be their Euclidean norm. The planar joint density is proportional to exp(−r²/(2σ²)); polar coordinates introduce Jacobian r. Integrating the angle from 0 to 2π gives the stated PDF. Equivalently:

```math
R^2/\sigma^2\sim\chi^2_2,\qquad R^2/(2\sigma^2)\sim\operatorname{Exp}(1).
```

Use the Weibull moment formula to obtain mean and variance.

</details>

[Back to table of contents](#table-of-contents)

---

## Laplace

Use: symmetric, sharper-peaked and heavier-tailed errors than normal; absolute-deviation loss. Parameters: real location μ, scale b > 0. Support: all real numbers.

```math
f_X(x)=\frac1{2b}e^{-|x-\mu|/b},\quad
F_X(x)=\begin{cases}\tfrac12e^{(x-\mu)/b},&x<\mu,\\1-\tfrac12e^{-(x-\mu)/b},&x\ge\mu.\end{cases}
```

```math
E[X]=\mu,\quad \operatorname{Var}(X)=2b^2,\quad
\operatorname{SD}(X)=\sqrt2b,\quad
M_X(t)=\frac{e^{\mu t}}{1-b^2t^2}\quad(|t|<1/b).
```

Mean = median = mode = μ. For fixed b, maximizing a Laplace likelihood over μ is equivalent to minimizing the sum of absolute residuals, giving a sample median (possibly a nonunique interval).

<details>
<summary>Derivation: signed exponential or difference of exponentials</summary>

Let E be exponential with rate 1/b, and choose an independent sign with probability 1/2 each. Then X = μ ± E has the stated two-sided density. Alternatively:

```math
X=\mu+E_1-E_2,\qquad E_i\overset{\mathrm{iid}}\sim\operatorname{Exp}(1/b).
```

The means cancel, variances add, and multiplying their MGFs gives the formula.

</details>

[Back to table of contents](#table-of-contents)

---

## Cauchy

Use: a canonical heavy-tailed law and ratio of centered normal variables. Parameters: real location x0, scale γ > 0. Support: all real numbers.

```math
f_X(x)=\frac1{\pi\gamma[1+((x-x_0)/\gamma)^2]},\quad
F_X(x)=\frac12+\frac1\pi\arctan\left(\frac{x-x_0}{\gamma}\right).
```

Mean, variance, and SD are undefined. Median and mode = x0; interquartile range = 2γ. The scale is not an SD. The MGF is infinite for every nonzero real t. Its characteristic function is:

```math
\varphi_X(t)=e^{ix_0t-\gamma|t|}.
```

The average of iid Cauchy(x0, γ) variables has exactly the same Cauchy(x0, γ) distribution: averaging does not reduce its scale. A symmetric principal-value integral of x times the standard density is zero, but that is not an expectation.

<details>
<summary>Derivation: ratio of standard normals</summary>

For independent standard normals Z1 and Z2, transform Y = Z1/Z2, V = Z2. The inverse has Z1 = YV and absolute Jacobian |V|:

```math
f_Y(y)=\int_{-\infty}^{\infty}\frac{|v|}{2\pi}e^{-(1+y^2)v^2/2}\,dv
=\frac1{\pi(1+y^2)}.
```

Shift and scale for the general law. The absolute first-moment integrand behaves as a constant divided by |x| in the tails and is not integrable. This is why symmetry cannot establish a mean.

</details>

[Back to table of contents](#table-of-contents)

---

## Pareto Type I

Use: power-law upper tails above a minimum threshold. Parameters: minimum xm > 0 and shape α > 0. Support: x ≥ xm.

```math
f_X(x)=\frac{\alpha x_m^\alpha}{x^{\alpha+1}},\quad
F_X(x)=1-\left(\frac{x_m}{x}\right)^\alpha,\quad
S_X(x)=\left(\frac{x_m}{x}\right)^\alpha.
```

```math
E[X]=\frac{\alpha x_m}{\alpha-1}\quad(\alpha>1),\quad
\operatorname{Var}(X)=\frac{\alpha x_m^2}{(\alpha-1)^2(\alpha-2)},\quad
\operatorname{SD}(X)=\frac{x_m\sqrt\alpha}{(\alpha-1)\sqrt{\alpha-2}}\quad(\alpha>2).
```

Mean is infinite for α ≤ 1; second moment is infinite for α ≤ 2. No finite SD exists for α ≤ 2. Median = xm times 2 to the power 1/α; mode = xm. The MGF is finite for t ≤ 0, not positive t. Smaller α means a heavier upper tail. A Pareto tail model does not imply the entire data distribution is Pareto.

<details>
<summary>Derivation: exponentiate an exponential</summary>

Let Y be exponential with rate α and set X = xm exp(Y). Then:

```math
P(X>x)=P(Y>\log(x/x_m))=(x_m/x)^\alpha.
```

Differentiate for the density. Raw moments are:

```math
E[X^q]=\frac{\alpha x_m^q}{\alpha-q}\quad(q<\alpha).
```

This establishes both the moment formulas and their existence thresholds.

</details>

[Back to table of contents](#table-of-contents)

---

## Inverse gamma

Use: positive scale uncertainty, including priors on variance. Parameters: shape α > 0, inverse-scale coefficient β > 0. Support: x > 0. Define by the density rather than relying on inconsistent naming of β across sources.

```math
f_X(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{-\alpha-1}e^{-\beta/x},\quad
F_X(x)=\frac{\Gamma(\alpha,\beta/x)}{\Gamma(\alpha)}.
```

Here the two-argument capital Gamma is the upper incomplete gamma integral from its second argument to infinity.

```math
E[X]=\frac\beta{\alpha-1}\quad(\alpha>1),\quad
\operatorname{Var}(X)=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)},\quad
\operatorname{SD}(X)=\frac\beta{(\alpha-1)\sqrt{\alpha-2}}\quad(\alpha>2).
```

Mean is infinite for α ≤ 1; second moment is infinite for α ≤ 2. No finite SD exists for α ≤ 2. Mode = β/(α + 1). The MGF is finite for t ≤ 0, not positive t.

<details>
<summary>Derivation: reciprocal of a gamma variable</summary>

If Y is Gamma(α, rate β) and X = 1/Y, the inverse y = 1/x has absolute derivative 1/x². Substitution produces the stated PDF. Thus:

```math
E[X^q]=\beta^q\frac{\Gamma(\alpha-q)}{\Gamma(\alpha)}\quad(q<\alpha).
```

For normal observations with known mean μ, an inverse-gamma prior on variance updates as:

```math
\sigma^2\sim\operatorname{InvGamma}(\alpha,\beta)
\Longrightarrow
\sigma^2\mid x\sim\operatorname{InvGamma}\left(\alpha+\frac n2,\beta+\frac12\sum_i(x_i-\mu)^2\right).
```

This update assumes μ is known; estimating or integrating an unknown mean changes the derivation.

</details>

[Back to table of contents](#table-of-contents)

---

## Multinomial

Use: category counts from n independent categorical trials. Parameters: integer n ≥ 0; K probabilities summing to 1. Support: nonnegative integer vectors with total n.

```math
P(\mathbf X=\mathbf x)=\frac{n!}{\prod_{j=1}^Kx_j!}\prod_{j=1}^Kp_j^{x_j},\qquad \sum_j x_j=n.
```

```math
E[X_j]=np_j,\quad \operatorname{Var}(X_j)=np_j(1-p_j),\quad
\operatorname{SD}(X_j)=\sqrt{np_j(1-p_j)},\quad
\operatorname{Cov}(X_i,X_j)=-np_ip_j\ (i\ne j).
```

```math
\operatorname{Cov}(\mathbf X)=n[\operatorname{diag}(\mathbf p)-\mathbf p\mathbf p^T],\quad
M_{\mathbf X}(\mathbf t)=\left(\sum_jp_je^{t_j}\right)^n.
```

Every individual count is binomial, but counts are not independent because their sum is fixed. The covariance matrix is singular. With K = 2, specifying one count determines the other.

<details>
<summary>Derivation: count category assignments</summary>

A specified ordered sequence with counts x has probability equal to the product of category probabilities to their counts. The multinomial coefficient counts distinct sequences with those counts. The multinomial theorem normalizes the PMF.

Write each count as the sum of n one-hot indicators. Add their means and covariance matrices across independent trials; indicators within each trial are mutually exclusive, which gives negative cross-category covariance.

</details>

[Back to table of contents](#table-of-contents)

---

## Dirichlet

Use: random probability vectors; the multivariate extension of beta. Parameters: K ≥ 2 positive concentrations αj; write α0 for their sum. Support: positive vectors summing to 1 (the simplex).

```math
f(\mathbf p)=\frac{\Gamma(\alpha_0)}{\prod_j\Gamma(\alpha_j)}\prod_jp_j^{\alpha_j-1},\qquad
\alpha_0=\sum_j\alpha_j.
```

This is a density in K − 1 free coordinates, with the last coordinate determined by the sum constraint, not a density with respect to K-dimensional volume.

```math
E[P_j]=\frac{\alpha_j}{\alpha_0},\quad
\operatorname{Var}(P_j)=\frac{\alpha_j(\alpha_0-\alpha_j)}{\alpha_0^2(\alpha_0+1)},\quad
\operatorname{SD}(P_j)=\sqrt{\operatorname{Var}(P_j)}.
```

```math
\operatorname{Cov}(P_i,P_j)=-\frac{\alpha_i\alpha_j}{\alpha_0^2(\alpha_0+1)}\quad(i\ne j),\qquad
P_j\sim\operatorname{Beta}(\alpha_j,\alpha_0-\alpha_j).
```

For a fixed mean vector, larger α0 concentrates the distribution more tightly around that mean. Multinomial conjugacy updates each αj by its observed category count. K = 2 reduces to beta. The MGF is finite for all finite real vector arguments because the simplex is bounded.

<details>
<summary>Derivation: normalize independent gamma variables</summary>

Let Gj be independent Gamma(αj, rate c) variables with a common rate. Put T equal to their sum and Pj = Gj/T. The transformation from the K − 1 free proportions and T has absolute Jacobian T to the K − 1. The joint density factors into:

```math
T\sim\operatorname{Gamma}(\alpha_0,\text{rate}=c),\qquad
\mathbf P\sim\operatorname{Dirichlet}(\boldsymbol\alpha),\qquad T\perp\mathbf P.
```

Integrating out T produces the normalizing constant. Mixed raw moments are:

```math
E\left[\prod_jP_j^{r_j}\right]=\frac{\Gamma(\alpha_0)}{\Gamma(\alpha_0+\sum_jr_j)}
\prod_j\frac{\Gamma(\alpha_j+r_j)}{\Gamma(\alpha_j)},\qquad r_j>-\alpha_j.
```

First and second moments yield the stated covariance structure.

</details>

[Back to table of contents](#table-of-contents)

---

## Multivariate normal

Use: jointly Gaussian vectors, linear models, Gaussian conditioning, and asymptotic inference. Parameters: d-dimensional mean vector μ and symmetric positive-definite covariance matrix Σ. Support: all of real d-space.

```math
f_{\mathbf X}(\mathbf x)=\frac{\exp[-\tfrac12(\mathbf x-\boldsymbol\mu)^T\Sigma^{-1}(\mathbf x-\boldsymbol\mu)]}{(2\pi)^{d/2}|\Sigma|^{1/2}}.
```

```math
E[\mathbf X]=\boldsymbol\mu,\quad \operatorname{Cov}(\mathbf X)=\Sigma,\quad
\operatorname{SD}(X_j)=\sqrt{\Sigma_{jj}},\quad
M_{\mathbf X}(\mathbf t)=e^{\mathbf t^T\boldsymbol\mu+\mathbf t^T\Sigma\mathbf t/2}.
```

A vector has a covariance matrix and componentwise SDs, not a single scalar variance. Affine transformations remain jointly normal:

```math
A\mathbf X+\mathbf b\sim N(A\boldsymbol\mu+\mathbf b,A\Sigma A^T).
```

Uncorrelated jointly Gaussian components are independent; uncorrelatedness alone does not imply independence for arbitrary distributions. Having normal marginals alone does not establish joint normality. Positive-semidefinite singular Σ also defines a Gaussian law, but it has no full-dimensional PDF of the form above.

<details>
<summary>Derivation and conditional distribution</summary>

Let Z have independent standard normal components and take a nonsingular L with LL transpose = Σ (for example a Cholesky factor). Set X = μ + LZ. The inverse Jacobian determinant has magnitude 1/|L| and |Σ| = |L|²; substitution produces the density.

Partition X into blocks Xa and Xb, and partition μ and Σ accordingly. Completing the square gives:

```math
\mathbf X_a\mid\mathbf X_b=\mathbf b\sim
N\left(\boldsymbol\mu_a+\Sigma_{ab}\Sigma_{bb}^{-1}(\mathbf b-\boldsymbol\mu_b),
\Sigma_{aa}-\Sigma_{ab}\Sigma_{bb}^{-1}\Sigma_{ba}\right).
```

This formula assumes the conditioning covariance block is invertible. The conditional covariance does not depend on the observed value b.

</details>

[Back to table of contents](#table-of-contents)

---

## Relationships and approximation map

### Exact relationships
Each statement below assumes independence where more than one variable is combined.

| Starting point | Result | Conditions / parameter mapping |
|---|---|---|
| Sum of Bernoulli variables | Binomial | Common success probability p |
| Sum of geometric trial counts | Negative-binomial total trials | Integer number r of terms, common p; subtract r for failure count |
| Sum of Poisson variables | Poisson | Add intensities |
| Sum of exponentials | Erlang | Common rate; shape = number of terms |
| Sum of gamma variables | Gamma | Common rate; add shapes |
| Squared standard normal | Chi-square | 1 degree of freedom |
| Sum of squared standard normals | Chi-square | Degrees of freedom = number of terms |
| Standard normal / square root of scaled chi-square | Student t | Denominator chi-square divided by its degrees of freedom |
| Ratio of scaled chi-squares | F | Each chi-square divided by its own degrees of freedom |
| Squared Student t | F | Numerator degrees of freedom 1 |
| Ratio of two standard normals | Cauchy | Both centered, unit variance |
| Exponential of normal | Lognormal | Normal mean and variance become log parameters |
| Reciprocal of gamma | Inverse gamma | Gamma rate becomes inverse-gamma coefficient β |
| Gamma / sum of two gammas | Beta | Common rate |
| Normalized vector of gammas | Dirichlet | Common rate |
| Poisson with gamma-random intensity | Negative binomial | Gamma shape r, rate p/(1 − p) |
| Binomial with beta-random probability | Beta-binomial | A shared random probability |

### Asymptotic approximations

- Binomial to Poisson: n grows, p shrinks, and np converges to a finite intensity. A finite-sample approximation should be checked at the relevant tail or interval.
- Binomial to normal: center at np, variance np(1 − p), when both expected successes and failures are sufficiently large. There is no universal cutoff that guarantees all tail accuracies.
- Poisson to normal: center at λ, variance λ, as λ grows. Relative accuracy can still be poor far into a tail.
- Hypergeometric to binomial: when sampling changes population proportions negligibly; check the sampling fraction and accuracy required.
- Student t to standard normal: degrees of freedom grow.
- Standardized gamma to standard normal: shape grows.

For normal approximation to integer-valued X, a continuity correction for inclusive integer endpoints is:

```math
P(a\le X\le b)\approx
\Phi\left(\frac{b+1/2-\mu}{\sigma}\right)-
\Phi\left(\frac{a-1/2-\mu}{\sigma}\right).
```

Approximations are not identities. Prefer exact library PMFs/CDFs when readily available, especially for small samples or tail probabilities.

[Back to table of contents](#table-of-contents)

---

## Software parameterization and numerical checks

The table maps this guide's mathematical parameters to familiar R and SciPy APIs. In R, replace the density prefix `d` with `p`, `q`, or `r` for CDF, quantile, or sampling. In SciPy, distribution objects expose `pmf`/`pdf`, `cdf`, `ppf`, `rvs`, `mean`, `var`, and `std` where supported.

| Distribution | R density / mass | SciPy distribution object |
|---|---|---|
| Bernoulli(p) | `dbinom(x, size=1, prob=p)` | `stats.bernoulli(p)` |
| Binomial(n, p) | `dbinom(x, size=n, prob=p)` | `stats.binom(n, p)` |
| Categorical(p) | Probability-vector lookup | `stats.rv_discrete(values=(labels, p))` |
| Discrete uniform on a through b | Probability-vector lookup | `stats.randint(a, b+1)` |
| Geometric trials | `dgeom(x-1, prob=p)` | `stats.geom(p)` |
| Negative-binomial failures | `dnbinom(x, size=r, prob=p)` | `stats.nbinom(r, p)` |
| Hypergeometric(N, K, n) | `dhyper(x, m=K, n=N-K, k=n)` | `stats.hypergeom(N, K, n)` |
| Poisson(λ) | `dpois(x, lambda=lam)` | `stats.poisson(mu=lam)` |
| Beta-binomial(n, α, β) | No base-R density | `stats.betabinom(n, alpha, beta)` |
| Uniform(a, b) | `dunif(x, min=a, max=b)` | `stats.uniform(loc=a, scale=b-a)` |
| Normal(μ, σ²) | `dnorm(x, mean=mu, sd=sigma)` | `stats.norm(loc=mu, scale=sigma)` |
| Exponential(rate λ) | `dexp(x, rate=lam)` | `stats.expon(scale=1/lam)` |
| Gamma(shape α, rate β) | `dgamma(x, shape=alpha, rate=beta)` | `stats.gamma(a=alpha, scale=1/beta)` |
| Beta(α, β) | `dbeta(x, alpha, beta)` | `stats.beta(alpha, beta)` |
| Chi-square(ν) | `dchisq(x, df=nu)` | `stats.chi2(df=nu)` |
| Standard Student t(ν) | `dt(x, df=nu)` | `stats.t(df=nu)` |
| F(d1, d2) | `df(x, df1=d1, df2=d2)` | `stats.f(dfn=d1, dfd=d2)` |
| Lognormal(μ, σ²) | `dlnorm(x, meanlog=mu, sdlog=sigma)` | `stats.lognorm(s=sigma, scale=exp(mu))` |
| Weibull(k, η) | `dweibull(x, shape=k, scale=eta)` | `stats.weibull_min(c=k, scale=eta)` |
| Rayleigh(σ) | No base-R density | `stats.rayleigh(scale=sigma)` |
| Laplace(μ, b) | No base-R density | `stats.laplace(loc=mu, scale=b)` |
| Cauchy(x0, γ) | `dcauchy(x, location=x0, scale=gam)` | `stats.cauchy(loc=x0, scale=gam)` |
| Pareto I(xm, α) | No base-R density | `stats.pareto(b=alpha, scale=xm)` |
| Inverse gamma(α, β) | No base-R density | `stats.invgamma(a=alpha, scale=beta)` |
| Multinomial(n, p) | `dmultinom(x, size=n, prob=p)` | `stats.multinomial(n=n, p=p)` |
| Dirichlet(α vector) | No base-R density | `stats.dirichlet(alpha)` |
| Multivariate normal(μ, Σ) | No base-R density | `stats.multivariate_normal(mean=mu, cov=Sigma)` |

SciPy entries assume `from scipy import stats`; `exp(mu)` denotes the exponential function. R's geometric distribution counts failures, unlike the convention in this guide. SciPy's integer-uniform upper endpoint is exclusive. Inverse gamma takes `scale=beta`, even though gamma takes `scale=1/beta` in this guide's rate convention.

### Probability calculation examples

```python
import numpy as np
from scipy import stats

# P(3 <= X <= 7) for X ~ Binomial(20, 0.25)
x = stats.binom(n=20, p=0.25)
p_interval = x.cdf(7) - x.cdf(2)

# P(X > 10) for Gamma(shape=3, rate=0.5)
y = stats.gamma(a=3, scale=1 / 0.5)
p_tail = y.sf(10)

# Geometric: 4 trials means 3 failures before the first success.
p_four_trials = stats.geom(p=0.2).pmf(4)
p_three_failures = stats.nbinom(n=1, p=0.2).pmf(3)
assert np.isclose(p_four_trials, p_three_failures)

# Check that a rate/scale mapping agrees with the analytic moments.
assert np.isclose(y.mean(), 3 / 0.5)
assert np.isclose(y.var(), 3 / 0.5**2)
```

```r
# Same binomial interval.
pbinom(7, size = 20, prob = 0.25) -
  pbinom(2, size = 20, prob = 0.25)

# Upper gamma tail; avoid subtracting a CDF close to 1.
pgamma(10, shape = 3, rate = 0.5, lower.tail = FALSE)

# R geometric counts failures: trial 4 corresponds to failure count 3.
dgeom(3, prob = 0.2)
```

### Numerical habits
- Use survival functions or upper-tail options rather than subtracting a nearly-one CDF from 1.
- Use `logpmf`/`logpdf` in SciPy or `log=TRUE` for R densities when multiplying many small probabilities; add log probabilities instead of multiplying probabilities.
- For large combinatorial coefficients use log-gamma functions or library PMFs rather than factorials.
- Distinguish existence of moments from numerical software output: undefined means, infinite moments, and NaN are not interchangeable concepts.
- Verify support, endpoint inclusivity, units, and rate-versus-scale conventions before substituting numbers.

[Back to table of contents](#table-of-contents)

---

## Source notes and further reading

The guide synthesizes standard distribution identities, with short derivations included so the formulas can be checked independently. Parameterizations are fixed explicitly above rather than copied implicitly from software. Useful verification entry points:

- SciPy continuous-distribution tutorial: https://docs.scipy.org/doc/scipy-1.13.0/tutorial/stats/continuous.html
- SciPy discrete-distribution tutorial: https://docs.scipy.org/doc/scipy-1.13.0/tutorial/stats/discrete.html
- SciPy statistical-distribution API, including multivariate families: https://docs.scipy.org/doc/scipy-1.13.0/reference/stats.html
- R distribution documentation is available locally with `?Distributions` and individual help pages such as `?Gamma`, `?Geometric`, and `?NegBinomial`.

These are documentation entry points, not claims that every identity uses the sources' parameter names. For any implementation, check its specific distribution page and installed version. The software mapping is a convenience reference, not a substitute for those definitions.

[Back to table of contents](#table-of-contents)
