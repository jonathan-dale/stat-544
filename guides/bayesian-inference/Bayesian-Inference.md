# STAT 544 Study Guide and Lecture Prep: Bayesian Inference

Bayesian inference uses a probability model to learn about an unknown parameter from observed data. The central progression is **sampling model → likelihood × prior → normalized posterior → probabilities, estimates, and predictions**. The posterior is a conditional distribution, so the ideas from conditional distributions are not being replaced; they are being applied to an unknown parameter.

This single-part guide follows the foundation-first, explicit-derivation approach of the attached conditional-distribution and Poisson-process guides. It is designed for graduate-level study without assuming the usual undergraduate statistics sequence. Theory and hand calculations come before R, and optional extensions are separated from the lecture's main material.

## How to use this guide

Read the prerequisite review before trying to memorize conjugate updates. On a second pass, cover the solutions and reconstruct the likelihood, support, normalization, and interpretation yourself.

- **Lecture core:** Classical versus Bayesian inference; Bayes' theorem; prior, likelihood, and posterior; uniform and beta priors for binomial data; posterior means; loss functions; estimation of functions of a parameter; and the posterior mean as a weighted average.
- **Foundational completion:** Credible intervals, sequential updating, prior sensitivity, and posterior predictive distributions. These explain how to use and interpret a posterior rather than stopping at its formula.
- **Clearly labeled extensions:** Poisson-gamma and normal-normal models, Bayes factors for a point null, improper priors, transformations, and a short formal probability reference. These are useful graduate-level connections, not claims about additional assigned lecture topics.
- **Practice and verification:** Seven worked examples, seventeen practice problems with separate solutions, a compact formula sheet, and executable base-R checks with simulations reported separately.

### Lecture alignment and source conventions

The supplied `4L544F2611BAY.pdf` contains 27 lecture slides arranged four per PDF page, except for the last page. References below use the slide numbers printed on the slides, not the seven PDF-page numbers.

| Lecture slides | Main topic | Guide location |
|---|---|---|
| 3–4 | Classical versus Bayesian inference | Prerequisite review and the Bayesian model |
| 5–9, 13 | Bayes reversal and manuscript authorship | Bayes' theorem and worked example A |
| 10–12 | Joint, marginal, and conditional densities | The Bayesian model and posterior derivation |
| 14–20 | Uniform/beta prior and binomial likelihood | Beta-binomial derivations and worked examples B–C |
| 21–25 | Posterior mean and loss-based estimation | Bayes estimation and worked example D |
| 26–27 | Posterior mean as a weighted average | Shrinkage, sequential updating, and prior sensitivity |

The manuscript and eight-trial examples retain the lecture's assumptions and data. Other numerical examples and exercises are instructional additions. The lecture writes the prior and posterior as $g(\theta)$ and $g(\theta\mid x)$; this guide uses $\pi(\theta)$ and $\pi(\theta\mid x)$, with exactly the same meaning.

### Contents

- [Prerequisite review](#prerequisite-review)
- [The Bayesian model and the meaning of inference](#the-bayesian-model-and-the-meaning-of-inference)
- [Bayes theorem from joint to marginal to posterior](#bayes-theorem-from-joint-to-marginal-to-posterior)
- [Beta and gamma integrals](#beta-and-gamma-integrals)
- [The beta-binomial model](#the-beta-binomial-model)
- [Sequential updating and prior sensitivity](#sequential-updating-and-prior-sensitivity)
- [Bayes estimation and loss functions](#bayes-estimation-and-loss-functions)
- [Credible intervals and posterior prediction](#credible-intervals-and-posterior-prediction)
- [Extension: Poisson-gamma inference](#extension-poisson-gamma-inference)
- [Extension: normal-normal inference](#extension-normal-normal-inference)
- [A reliable problem-solving procedure](#a-reliable-problem-solving-procedure)
- [Seven fully worked examples](#seven-fully-worked-examples)
- [Graduate reference and modeling cautions](#graduate-reference-and-modeling-cautions)
- [Practice problems](#practice-problems)
- [Solutions and checks](#solutions-and-checks)
- [Executable R checks](#executable-r-checks)
- [Common errors and lecture preparation](#common-errors-and-lecture-preparation)
- [Compact reference sheet](#compact-reference-sheet)

## Prerequisite review

### Random quantities, observations, and parameters

A random variable is a numerical description of an uncertain outcome. Its distribution describes how probability is assigned to its possible values. Uppercase letters usually denote random quantities, while lowercase letters denote particular values.

For example, $X$ may be the number of successful requests among eight trials. Before observing the trials, $X$ is uncertain; after observing six successes, its realized value is $x=6$.

A parameter indexes a family of probability models. In

$$
X\mid\Theta=\theta\sim\mathrm{Binomial}(8,\theta),
$$

$\theta$ is a possible success probability, not the observed success count. Different choices of $\theta$ describe different distributions for $X$.

This guide writes $\Theta$ when emphasizing the random parameter under a prior or posterior, and $\theta$ for its possible values. Many statistics books and the lecture use lowercase $\theta$ for both roles. Read $X\mid\theta$ as shorthand for $X\mid\Theta=\theta$.

| Object | Example | What it represents |
|---|---|---|
| Random data | $X$ | What could be observed |
| Observed data | $x=6$ | What actually was observed |
| Parameter value | $\theta=0.7$ | A possible setting of the sampling model |
| Random parameter | $\Theta$ | Uncertainty about the parameter in the Bayesian model |
| Hyperparameters | $\alpha,\beta$ | Known settings that specify a prior distribution |
| Statistic | $T(X)=X/8$ | A function of the data that does not use the unknown parameter |
| Estimator | $\delta(X)$ | A data-dependent rule for estimating a target |
| Realized estimate | $\delta(6)$ | The numerical output of the rule after observing data |

A parameter need not physically change between trials for a prior to be useful. The Bayesian distribution can represent uncertainty about a single persistent value. Randomizing uncertainty is different from claiming that nature independently redraws the parameter for each observation.

### Probability mass versus probability density

A probability mass function, or PMF, gives probabilities at individual points of a discrete distribution:

$$
p_X(x)=P(X=x),\qquad
p_X(x)\ge0,\qquad
\sum_xp_X(x)=1.
$$

A probability density function, or PDF, gives probabilities through integration:

$$
P(\Theta\in A)=\int_A\pi(\theta)\,d\theta,\qquad
\pi(\theta)\ge0,\qquad
\int\pi(\theta)\,d\theta=1.
$$

For a parameter with an ordinary continuous density, $P(\Theta=\theta)=0$ at every fixed point. A density height can exceed 1 because a height is not itself a probability. Whenever a PDF appears, identify the variable over which it integrates to 1.

The support is the set of allowed values for the calculation. A success probability lies in $[0,1]$, a positive event rate lies in $(0,\infty)$, and an unrestricted normal mean lies in $\mathbb R$. For continuous densities, including or excluding an endpoint does not change interval probabilities, but it can matter when discussing a mode or an allowed decision.

### Conditional probability, partitions, and independence

An event is a set of outcomes, such as “the manuscript was written by $H$” or “the count equals six.” For events $A,B$ with $P(B)>0$, conditional probability is

$$
P(A\mid B)=\frac{P(A\cap B)}{P(B)}.
$$

The intersection $A\cap B$ means both events occur. Conditioning restricts attention to outcomes in $B$ and rescales their probabilities to sum to 1; it does not by itself assert that $B$ caused $A$.

Rearranging gives the product rule $P(A\cap B)=P(A\mid B)P(B)$. Applying the same rule in the other order is the basis for Bayes' theorem.

A partition $H_1,\ldots,H_k$ divides the outcome space into disjoint, exhaustive events. Every outcome belongs to exactly one member, so

$$
P(D)=\sum_{j=1}^kP(D\cap H_j)
=\sum_{j=1}^kP(D\mid H_j)P(H_j),
$$

when the conditional probabilities are defined. This is the law of total probability: calculate the observation probability within each possible source, then average using the source probabilities.

Independence instead means $P(A\cap B)=P(A)P(B)$. Disjoint events cannot both occur, whereas independent events do not change one another's probabilities; two disjoint events with positive probabilities cannot be independent. Conditional independence requires this factorization within the distribution given the conditioning information, not necessarily in the unconditional distribution.

### Joint, marginal, and conditional distributions

A joint distribution describes two quantities together. A marginal distribution removes one quantity by summing or integrating over it. A conditional distribution fixes information about one quantity and describes uncertainty about the other.

For a continuous parameter and discrete data, the Bayesian model has a mixed joint mass-density:

$$
q(x,\theta)=p(x\mid\theta)\pi(\theta).
$$

Its precise meaning is

$$
P(X=x,\Theta\in A)=\int_Aq(x,\theta)\,d\theta.
$$

It is not a point probability for the continuous coordinate. Marginalizing the parameter gives the data PMF:

$$
m(x)=\int p(x\mid\theta)\pi(\theta)\,d\theta.
$$

Dividing this joint mass-density by $m(x)>0$ gives a density in $\theta$. This is the same joint-to-marginal-to-conditional procedure used in the earlier guides.

### Expectation, variance, and quantiles

Expectation is a probability-weighted average. For an integrable function $h$,

$$
E[h(\Theta)]
=\sum_jh(\theta_j)P(\Theta=\theta_j)
\quad\text{or}\quad
\int h(\theta)\pi(\theta)\,d\theta.
$$

The mean is $E[\Theta]$. With a finite second moment, variance and standard deviation are

$$
\mathrm{Var}(\Theta)
=E[(\Theta-E[\Theta])^2]
=E[\Theta^2]-E[\Theta]^2,
\qquad
\mathrm{SD}(\Theta)=\sqrt{\mathrm{Var}(\Theta)}.
$$

The second moment $E[\Theta^2]$ is not the variance. Also, $E[h(\Theta)]$ generally differs from $h(E[\Theta])$ when $h$ is nonlinear. This distinction will determine the correct Bayes estimate of a transformed parameter.

A cumulative distribution function, or CDF, is $F(t)=P(\Theta\le t)$. A $p$-quantile can be defined as

$$
Q(p)=\inf\{t:F(t)\ge p\},\qquad 0<p<1.
$$

The median is a $0.5$-quantile, while the mode is a value where a chosen PMF or density is maximized, if such a value exists. Mean, median, and mode summarize different features and need not coincide.

### Bernoulli trials and binomial counts

A Bernoulli variable $Y$ takes values 1 and 0, usually labeled success and failure. Conditional on a success probability $\theta$,

$$
P(Y=y\mid\theta)=\theta^y(1-\theta)^{1-y},
\qquad y\in\{0,1\}.
$$

If $Y_1,\ldots,Y_n$ are conditionally independent with this same probability, their sum $X=\sum_{i=1}^nY_i$ has PMF

$$
p(x\mid\theta)
=\binom nx\theta^x(1-\theta)^{n-x},
\qquad x=0,\ldots,n.
$$

Here

$$
n!=n(n-1)\cdots1,\qquad 0!=1,\qquad
\binom nx=\frac{n!}{x!(n-x)!}.
$$

Why does the binomial coefficient appear? Each particular ordered sequence with $x$ successes has probability $\theta^x(1-\theta)^{n-x}$, and there are $\binom nx$ sequences with that count. A sequence and its success count are different data descriptions, even though their likelihoods are proportional as functions of $\theta$.

The moments are

$$
E[X\mid\theta]=n\theta,\qquad
\mathrm{Var}(X\mid\theta)=n\theta(1-\theta).
$$

For the sample proportion $\widehat\theta=X/n$, with $n>0$,

$$
E[\widehat\theta\mid\theta]=\theta,\qquad
\mathrm{Var}(\widehat\theta\mid\theta)
=\frac{\theta(1-\theta)}n.
$$

The first identity follows by linearity. The variance follows by adding the variances of conditionally independent Bernoulli trials and then dividing by $n^2$. These describe repeated-sample behavior at a fixed parameter, not posterior uncertainty about the parameter.

### Conditional independence is the important assumption

In a Bayesian Bernoulli model, all trials share one $\Theta$. They are independent given $\Theta$, so their conditional joint PMF factors:

$$
p(y_1,\ldots,y_n\mid\theta)
=\prod_{i=1}^np(y_i\mid\theta).
$$

After averaging over an uncertain common $\Theta$, the trials are generally dependent. For two distinct trials,

$$
\begin{aligned}
\mathrm{Cov}(Y_i,Y_j)
&=E[\mathrm{Cov}(Y_i,Y_j\mid\Theta)]
 +\mathrm{Cov}(E[Y_i\mid\Theta],E[Y_j\mid\Theta])\\
&=0+\mathrm{Var}(\Theta).
\end{aligned}
$$

The first term vanishes by conditional independence, but the second need not vanish. Observing one success provides information about the common parameter and therefore about other trials. This is the same shared-latent-rate mechanism discussed in the Poisson-process guide.

### The algebra and calculus used below

The main operations are multiplication of powers, normalization by an integral, differentiation of a log likelihood, and completing the square. Keep fixed quantities outside the relevant sum or integral, but do not discard a factor merely because its name looks like a constant.

- **Powers:** $\theta^r\theta^s=\theta^{r+s}$.
- **Logarithms:** $\log(ab)=\log a+\log b$, and $\log(a^r)=r\log a$ for positive arguments.
- **Power integration:** $\int t^r\,dt=t^{r+1}/(r+1)$ for $r\ne-1$, while $\int t^{-1}\,dt=\log t$ on $t>0$.
- **Square completion:** $A\theta^2-2B\theta=A(\theta-B/A)^2-B^2/A$, for $A>0$.
- **Dummy variables:** In $\int f(x\mid u)\pi(u)\,du$, $x$ is fixed and $u$ disappears after integration.

When we write $a(\theta)\propto b(\theta)$, we mean that $a(\theta)=C b(\theta)$ for a positive factor $C$ that does not depend on $\theta$. That factor may depend on the already observed data. Proportionality identifies shape; it does not by itself give a normalized distribution.

## The Bayesian model and the meaning of inference

### A likelihood is a function of the unknown parameter

The sampling model $f(x\mid\theta)$ describes possible data at a fixed parameter. After observing $x$, the likelihood is the same numerical expression viewed as a function of the parameter:

$$
L(\theta;x)=f(x\mid\theta),\qquad x\text{ fixed}.
$$

For example, if $x=6$ of $n=8$ trials succeeded,

$$
L(\theta;6)=\binom86\theta^6(1-\theta)^2,
\qquad 0\le\theta\le1.
$$

For fixed $\theta$, summing $p(x\mid\theta)$ over $x=0,\ldots,8$ gives 1. For fixed $x=6$, integrating this likelihood over $\theta$ need not give 1. A likelihood is not automatically a probability density for the parameter.

A maximum likelihood estimate, or MLE, maximizes $L(\theta;x)$ over allowed values of $\theta$. For $0<x<n$, differentiate the binomial log likelihood:

$$
\ell(\theta)=C+x\log\theta+(n-x)\log(1-\theta),
$$

$$
\ell'(\theta)=\frac{x}{\theta}
-\frac{n-x}{1-\theta}=0
\quad\Longrightarrow\quad
\widehat\theta_{\mathrm{MLE}}=\frac xn.
$$

The second derivative is negative in the interior, establishing a maximum. If $x=0$ or $x=n$, the maximum lies at 0 or 1 when those endpoints are allowed. On an open parameter space, it is instead a boundary supremum.

### Prior, posterior, evidence, and predictive distribution

A Bayesian model specifies a prior distribution for the parameter and a conditional distribution for data given the parameter. Together they determine a joint law; the posterior is the conditional law in the reverse direction ([MIT comparison of Bayesian and frequentist inference](https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_class20-prep.pdf)).

| Quantity | Notation | Question answered |
|---|---|---|
| Prior | $\pi(\theta)$ | What parameter values are plausible before these data? |
| Sampling model | $f(x\mid\theta)$ | What data are plausible if the parameter equals $\theta$? |
| Likelihood | $L(\theta;x)$ | How do the observed data compare different parameter values? |
| Evidence or marginal likelihood | $m(x)$ | How plausible are these data after averaging over the prior? |
| Posterior | $\pi(\theta\mid x)$ | What parameter values are plausible after these data? |
| Posterior predictive | $p(\widetilde x\mid x)$ | What future observations are plausible after these data? |

The evidence is also called the prior predictive PMF or density for the data. The posterior predictive distribution uses the posterior instead of the prior and describes a new observation. Do not confuse these two different averages.

### Classical versus Bayesian inference without caricatures

In the classical, or frequentist, setup, an unknown parameter is held fixed while the sampling distribution describes randomness in the data and in estimators. In a Bayesian model, probability also represents uncertainty about that parameter, and conditioning on the observed data gives the posterior; this is the central distinction emphasized in both the lecture and [MIT's comparison notes](https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_class20-prep.pdf).

Both approaches use probability models and require assumptions. A prior does not remove the need to justify a likelihood, and not using a prior does not make modeling choices disappear. For this guide, focus on understanding which distribution is being used for which question rather than declaring one school universally superior.

A frequentist confidence procedure has a repeated-sampling coverage property. A Bayesian credible set has a specified posterior probability conditional on the observed data. The intervals can sometimes coincide numerically, but those interpretations are not interchangeable.

## Bayes theorem from joint to marginal to posterior

### Finite hypotheses: normalize a weighted likelihood table

Suppose $H_1,\ldots,H_k$ are mutually exclusive and exhaustive hypotheses, with prior probabilities $w_j=P(H_j)$. Exhaustive means that exactly one is true within the model, and $\sum_jw_j=1$.

For an observation $D$ having positive total probability,

$$
P(D)=\sum_{j=1}^kP(D\mid H_j)w_j.
$$

The product rule gives $P(H_j,D)=P(D\mid H_j)w_j$. Dividing by the total probability yields

$$
\boxed{
P(H_j\mid D)=
\frac{P(D\mid H_j)w_j}
{\sum_{\ell=1}^kP(D\mid H_\ell)w_\ell}.
}
$$

Multiply each likelihood by its prior, then divide every resulting weight by their sum. Using likelihoods alone is equivalent to assuming equal prior weights, not to doing an assumption-free calculation.

For two hypotheses with positive prior weights and positive likelihoods, division of their posterior probabilities gives

$$
\frac{P(H_1\mid D)}{P(H_2\mid D)}
=
\frac{P(H_1)}{P(H_2)}
\frac{P(D\mid H_1)}{P(D\mid H_2)}.
$$

The likelihood ratio here is the Bayes factor for these two simple hypotheses. It multiplies prior odds to produce posterior odds. Evidence can favor $H_2$ while the posterior still favors $H_1$, if the starting odds for $H_1$ were sufficiently large.

### A continuous parameter: build, integrate, divide

Let $\Theta$ have prior density $\pi(\theta)$, and let $f(x\mid\theta)$ be a conditional PMF or density for the data. The joint mass-density or density is

$$
q(x,\theta)=f(x\mid\theta)\pi(\theta).
$$

Integrating out the parameter gives the marginal PMF or density of the observed data:

$$
m(x)=\int_{\mathcal T}f(x\mid u)\pi(u)\,du,
$$

where $\mathcal T$ is the parameter support. For $0<m(x)<\infty$, divide the joint by this marginal:

$$
\boxed{
\pi(\theta\mid x)
=\frac{f(x\mid\theta)\pi(\theta)}
{\int_{\mathcal T}f(x\mid u)\pi(u)\,du}.
}
$$

The normalization check is immediate:

$$
\int_{\mathcal T}\pi(\theta\mid x)\,d\theta
=\frac{\int_{\mathcal T}f(x\mid\theta)\pi(\theta)\,d\theta}
{m(x)}=1.
$$

If the parameter is discrete, replace the integral over parameter values with a sum. Whether the data are discrete or continuous does not determine that choice; the type of the parameter being marginalized does.

For continuous data, $m(x)$ is a density value, not $P(X=x)$. The formula is a conditional-density construction and does not divide by a probability-zero point event. As in the conditional-distribution guides, continuous conditional laws are determined only marginal-almost everywhere unless a particular version is selected.

### Kernel recognition and what may be omitted

With data $x$ fixed, the evidence does not depend on $\theta$. Thus

$$
\pi(\theta\mid x)\propto
f(x\mid\theta)\pi(\theta).
$$

For posterior-shape identification, omit only factors independent of $\theta$. For binomial data, $\binom nx$ may be omitted, but $(1-\theta)^{n-x}$ may not. Support restrictions depending on $\theta$ are also part of the likelihood and must remain.

When computing $m(x)$, probabilities of datasets, or a Bayes factor between models, restore any factors needed for the full sampling PMF or density. A factor common to the compared quantities may cancel, but that cancellation must be shown rather than assumed.

## Beta and gamma integrals

### Why normalizing constants deserve attention

A formula proportional to $\theta^6(1-\theta)^2$ specifies relative heights on $(0,1)$, not the total probability. The beta and gamma integrals provide the constants needed to normalize common posterior shapes and calculate their moments.

These functions extend factorial calculations to positive noninteger arguments. You do not need a special-functions course to use them here; the integral definitions and recurrence relations explain the steps we need.

### The gamma function and its recurrence

For $a>0$, define

$$
\Gamma(a)=\int_0^\infty t^{a-1}e^{-t}\,dt.
$$

Integration by parts, with the boundary terms vanishing for $a>0$, gives

$$
\begin{aligned}
\Gamma(a+1)
&=\int_0^\infty t^ae^{-t}\,dt\\
&=\left[-t^ae^{-t}\right]_0^\infty
 +a\int_0^\infty t^{a-1}e^{-t}\,dt\\
&=a\Gamma(a).
\end{aligned}
$$

Because $\Gamma(1)=1$, this implies $\Gamma(k)=(k-1)!$ for positive integers $k$. In particular, $\Gamma(3)=2!$, not $3!$.

### The beta integral and its relation to gamma

For $a,b>0$, define

$$
B(a,b)=\int_0^1t^{a-1}(1-t)^{b-1}\,dt.
$$

The identity needed in the lecture is

$$
\boxed{
B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}.
}
$$

To see why, write the product of two gamma integrals as a nonnegative double integral. Set $s=u+v$, $t=u/(u+v)$, so $u=st$, $v=s(1-t)$, and the absolute Jacobian is $s$:

$$
\begin{aligned}
\Gamma(a)\Gamma(b)
&=\int_0^\infty\int_0^\infty
u^{a-1}v^{b-1}e^{-(u+v)}\,du\,dv\\
&=\int_0^\infty s^{a+b-1}e^{-s}\,ds
\int_0^1t^{a-1}(1-t)^{b-1}\,dt\\
&=\Gamma(a+b)B(a,b).
\end{aligned}
$$

Nonnegativity justifies separating the integrals. For positive integers,

$$
B(a,b)=\frac{(a-1)!(b-1)!}{(a+b-1)!}.
$$

The gamma recurrence also gives the useful ratios

$$
\frac{B(a+1,b)}{B(a,b)}=\frac a{a+b},
\qquad
\frac{B(a+2,b)}{B(a,b)}
=\frac{a(a+1)}{(a+b)(a+b+1)}.
$$

These ratios will calculate beta means and second moments directly from their defining integrals. Recognizing a distribution family is convenient, but these equations explain where its moment formulas come from.

## The beta-binomial model

### The beta prior

For $\alpha,\beta>0$, a beta prior has density

$$
\pi(\theta)
=\frac{\theta^{\alpha-1}(1-\theta)^{\beta-1}}
{B(\alpha,\beta)}
\mathbf1_{\{0<\theta<1\}}.
$$

Its integral is 1 by the definition of $B$. The beta family is flexible on the probability scale, and its conjugacy with binomial data is also developed in [MIT's beta and normal conjugate-prior notes](https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_class15-prep.pdf).

The shape parameters control where prior mass is placed. A $\mathrm{Beta}(1,1)$ density is uniform; equal parameters produce symmetry around $1/2$; parameters both greater than 1 give an interior peak; and parameters both less than 1 can give an integrable density that diverges near the endpoints. A divergent height is not a point mass at an endpoint.

The prior moments follow from the ratios above:

$$
E[\Theta]=\frac{\alpha}{\alpha+\beta},
\qquad
\mathrm{Var}(\Theta)
=\frac{\alpha\beta}
{(\alpha+\beta)^2(\alpha+\beta+1)}.
$$

Writing $\kappa=\alpha+\beta$ and $p_0=\alpha/\kappa$ gives

$$
\alpha=\kappa p_0,\qquad
\beta=\kappa(1-p_0),\qquad
\mathrm{Var}(\Theta)=\frac{p_0(1-p_0)}{\kappa+1}.
$$

For a fixed mean, a larger $\kappa$ means a more concentrated prior. This parameterization is often easier to interpret than selecting $\alpha,\beta$ separately.

### Derive the posterior by matching exponents

Assume

$$
\Theta\sim\mathrm{Beta}(\alpha,\beta),
\qquad
X\mid\Theta=\theta\sim\mathrm{Binomial}(n,\theta).
$$

For an observed count $x\in\{0,\ldots,n\}$, multiplication gives

$$
\begin{aligned}
p(x\mid\theta)\pi(\theta)
&=\frac{\binom nx}{B(\alpha,\beta)}
\theta^x(1-\theta)^{n-x}
\theta^{\alpha-1}(1-\theta)^{\beta-1}\\
&=\frac{\binom nx}{B(\alpha,\beta)}
\theta^{\alpha+x-1}(1-\theta)^{\beta+n-x-1}.
\end{aligned}
$$

Define

$$
A=\alpha+x,\qquad B=\beta+n-x,\qquad S=A+B=\alpha+\beta+n.
$$

Here $B$ without arguments is a posterior shape parameter, while $B(a,b)$ with arguments is the beta function. The notation is conventional but worth checking when reading a formula.

The kernel is $\theta^{A-1}(1-\theta)^{B-1}$ on $(0,1)$, whose integral is $B(A,B)$. Consequently,

$$
\boxed{
\Theta\mid X=x\sim\mathrm{Beta}(A,B)
=\mathrm{Beta}(\alpha+x,\beta+n-x).
}
$$

The normalized posterior density is

$$
\boxed{
\pi(\theta\mid x)
=\frac{\theta^{A-1}(1-\theta)^{B-1}}{B(A,B)}
\mathbf1_{\{0<\theta<1\}}.
}
$$

The prior family is conjugate because the posterior remains in the beta family, not because prior and posterior parameters are unchanged. This derivation works for positive real $\alpha,\beta$, not only integer shapes.

### Derive the evidence instead of forgetting the denominator

Integrating the full prior-times-likelihood expression gives

$$
\boxed{
m(x)
=\binom nx\frac{B(\alpha+x,\beta+n-x)}
{B(\alpha,\beta)},
\qquad x=0,\ldots,n.
}
$$

This is a beta-binomial PMF. It describes the count before observing data, after averaging over the prior. It is generally not a binomial PMF with success probability replaced by the prior mean.

Its normalization also follows without a special summation identity:

$$
\begin{aligned}
\sum_{x=0}^nm(x)
&=\int_0^1
\left[\sum_{x=0}^np(x\mid\theta)\right]\pi(\theta)\,d\theta\\
&=\int_0^1 1\cdot\pi(\theta)\,d\theta=1.
\end{aligned}
$$

The finite sum and integral can be interchanged. The inner sum is 1 because it is a binomial PMF at fixed $\theta$.

### Uniform prior: the lecture's special case

With $\alpha=\beta=1$, the prior density equals 1 on $(0,1)$. The evidence reduces to

$$
\begin{aligned}
m(x)
&=\binom nx B(x+1,n-x+1)\\
&=\frac{n!}{x!(n-x)!}
\frac{x!(n-x)!}{(n+1)!}\\
&=\boxed{\frac1{n+1}}.
\end{aligned}
$$

Thus the prior predictive count is discrete uniform on $0,\ldots,n$. This does not say that $X\mid\theta$ is uniform, and it does not say that the posterior of $\Theta$ is uniform.

Dividing the likelihood by $1/(n+1)$ gives

$$
\pi(\theta\mid x)
=(n+1)\binom nx\theta^x(1-\theta)^{n-x},
\qquad 0<\theta<1.
$$

Equivalently,

$$
\Theta\mid x\sim\mathrm{Beta}(x+1,n-x+1).
$$

### Derive posterior moments

For a beta posterior with positive shapes $A,B$,

$$
\begin{aligned}
E[\Theta\mid x]
&=\int_0^1\theta
\frac{\theta^{A-1}(1-\theta)^{B-1}}{B(A,B)}\,d\theta\\
&=\frac{B(A+1,B)}{B(A,B)}
=\frac A{A+B}.
\end{aligned}
$$

Likewise,

$$
E[\Theta^2\mid x]
=\frac{B(A+2,B)}{B(A,B)}
=\frac{A(A+1)}{S(S+1)}.
$$

Subtract the square of the mean to obtain

$$
\begin{aligned}
\mathrm{Var}(\Theta\mid x)
&=\frac{A(A+1)}{S(S+1)}-\frac{A^2}{S^2}\\
&=\frac{AS(A+1)-A^2(S+1)}{S^2(S+1)}\\
&=\boxed{\frac{AB}{S^2(S+1)}}.
\end{aligned}
$$

More generally, whenever the moment exists,

$$
E[\Theta^r\mid x]
=\frac{B(A+r,B)}{B(A,B)},
\qquad A+r>0.
$$

For nonnegative integer $r$, this ratio is a product of $r$ successive factors. Negative powers require a condition near zero, so do not apply this formula without checking existence.

### Posterior mean, MLE, and MAP are different objects

The maximum a posteriori estimate, or MAP, maximizes a posterior density in a specified parameterization. For $A>1,B>1$, differentiate its log density:

$$
\frac{d}{d\theta}\log\pi(\theta\mid x)
=\frac{A-1}{\theta}-\frac{B-1}{1-\theta}=0.
$$

This gives

$$
\widehat\theta_{\mathrm{MAP}}
=\frac{A-1}{A+B-2}.
$$

Under a uniform prior and $0<x<n$, MAP equals the MLE $x/n$, but the posterior mean is $(x+1)/(n+2)$. With nonuniform priors, even MAP need not equal MLE.

The displayed MAP formula assumes both shapes exceed 1. Boundary maxima, boundary suprema, a flat density, or two endpoint singularities need separate treatment. For example, a $\mathrm{Beta}(1,6)$ posterior decreases toward the right and has its density supremum at the left endpoint, not at a value obtained by blindly using the interior formula.

### Posterior mean as a weighted average

For $n>0$, rearrange the beta-binomial mean:

$$
\begin{aligned}
E[\Theta\mid x]
&=\frac{\alpha+x}{\alpha+\beta+n}\\
&=\frac{\alpha+\beta}{\alpha+\beta+n}
\frac{\alpha}{\alpha+\beta}
+\frac{n}{\alpha+\beta+n}\frac xn.
\end{aligned}
$$

Let

$$
w=\frac{\kappa}{\kappa+n},\qquad
p_0=\frac{\alpha}{\kappa},\qquad
\widehat p=\frac xn.
$$

Then

$$
\boxed{E[\Theta\mid x]=wp_0+(1-w)\widehat p.}
$$

The mean is pulled toward the prior mean, an effect called shrinkage. Because the weights are nonnegative and sum to 1, the posterior mean lies between $p_0$ and $\widehat p$. If they are equal, the mean does not move, although the distribution can become more concentrated.

For fixed positive prior shapes, $w\to0$ as $n\to\infty$. More explicitly,

$$
E[\Theta\mid x]-\frac xn
=\frac{\kappa}{\kappa+n}\left(p_0-\frac xn\right).
$$

This is an algebraic statement about this model with a fixed prior. It is not a universal theorem that any prior becomes irrelevant under any model.

The quantity $\kappa=\alpha+\beta$ acts as a prior sample size in this mean-weighting identity. Pseudo-count language has multiple conventions: the density exponents are $\alpha-1,\beta-1$, while the posterior-mean weights use $\alpha,\beta$. These are interpretations of prior information, not literal extra observed trials.

## Sequential updating and prior sensitivity

### Why yesterday's posterior can be today's prior

Let $D_1,D_2$ be data batches conditionally independent given a common $\Theta$. First update using $D_1$:

$$
\pi(\theta\mid D_1)\propto
p(D_1\mid\theta)\pi(\theta).
$$

Then update that posterior using only the new batch:

$$
\begin{aligned}
\pi(\theta\mid D_1,D_2)
&\propto p(D_2\mid\theta)\pi(\theta\mid D_1)\\
&\propto p(D_2\mid\theta)p(D_1\mid\theta)\pi(\theta).
\end{aligned}
$$

This equals a single update using both batches. If conditional independence is absent, the second likelihood must be $p(D_2\mid D_1,\theta)$ instead. Reusing $D_1$'s likelihood after already including it in the prior would double-count those data.

For binomial batches $(x_1,n_1)$ and $(x_2,n_2)$, the posterior shapes are

$$
A=\alpha+x_1+x_2,\qquad
B=\beta+(n_1-x_1)+(n_2-x_2).
$$

The shared parameter and stated sampling assumptions are essential. If the success probability changes between batches, this pooled update answers the wrong model.

### Choosing a prior transparently

A practical prior specification should identify its support, central location, concentration, and rationale. Previous measurements can inform it, and so can deliberately weak information, but “uniform” is not a synonym for “no assumptions.”

One way to choose a beta prior is to specify a plausible prior mean $p_0$ and concentration $\kappa$. Set $\alpha=\kappa p_0$ and $\beta=\kappa(1-p_0)$, then inspect its intervals and prior predictive counts. A prior that produces implausible datasets deserves reconsideration before it is used for an important conclusion.

Sensitivity analysis means repeating the same likelihood calculation under several defensible priors and comparing the posterior quantities that matter. It is especially important when data are sparse or a conclusion lies near a decision threshold; [MIT's inference-comparison notes](https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_class20-prep.pdf) emphasize assessing robustness to alternate priors.

Do not choose a strong prior after inspecting the same data and then report an ordinary update as though that prior were independent information. Data-dependent prior construction is a separate modeling issue and needs disclosure.

## Bayes estimation and loss functions

### A posterior is not yet a decision

An inference describes uncertainty; a decision selects an action. The same posterior can support different sensible decisions when the costs of errors differ.

Let $d$ be a proposed estimate or action, and let $L(\theta,d)\ge0$ be its loss if the true parameter value is $\theta$. After observing data $x$, the posterior expected loss is

$$
\rho(d\mid x)
=E[L(\Theta,d)\mid x]
=\int L(\theta,d)\pi(\theta\mid x)\,d\theta,
$$

with a sum for a discrete parameter. A Bayes action is a minimizer

$$
\boxed{
\delta_{\mathrm B}(x)\in
\mathop{\mathrm{arg\,min}}_{d\in\mathcal A}
\rho(d\mid x),
}
$$

where $\mathcal A$ is the allowed action set. Existence and uniqueness need not hold automatically; the examples below use cases where the minimizer exists.

Notice that $x$ and $d$ are fixed during this posterior average. We average over possible parameter values, not over possible datasets. The action becomes a data-dependent rule when the calculation is performed for every possible $x$.

### Squared error: derive the posterior mean

For squared error $L(\theta,d)=(\theta-d)^2$, assume $E[\Theta^2\mid x]<\infty$. Let $\mu_x=E[\Theta\mid x]$, and write

$$
\Theta-d=(\Theta-\mu_x)+(\mu_x-d).
$$

Squaring and taking the posterior expectation gives

$$
\begin{aligned}
E[(\Theta-d)^2\mid x]
&=E[(\Theta-\mu_x)^2\mid x]
+2(\mu_x-d)E[\Theta-\mu_x\mid x]
+(\mu_x-d)^2\\
&=\mathrm{Var}(\Theta\mid x)+(\mu_x-d)^2.
\end{aligned}
$$

The cross term is zero because $\mu_x$ is the posterior mean. The variance term does not depend on $d$, and the remaining square is uniquely minimized at $d=\mu_x$ if this action is allowed.

Therefore,

$$
\boxed{\delta_{\mathrm B}(x)=E[\Theta\mid x]}
$$

under squared error. The minimized posterior expected loss equals the posterior variance. The mean is not best under every possible loss; its optimality is tied to this criterion.

### Estimating a function is not plugging into the mean

Suppose the target is $\psi=h(\Theta)$, and loss is squared error in that target:

$$
L(\theta,d)=(h(\theta)-d)^2.
$$

Applying the same decomposition to $h(\Theta)$, assuming a finite posterior second moment of that target, gives

$$
\boxed{
\widehat\psi_{\mathrm B}
=E[h(\Theta)\mid x].
}
$$

In general, this is not $h(E[\Theta\mid x])$. For $h(\theta)=\theta^2$,

$$
E[\Theta^2\mid x]
=E[\Theta\mid x]^2+\mathrm{Var}(\Theta\mid x).
$$

The difference is exactly the posterior variance. This identity is a useful sanity check and a direct application of the lecture's function-of-a-parameter result.

### Absolute error: derive the posterior median

For $L(\theta,d)=|\theta-d|$, assume a finite posterior first absolute moment. With a continuous posterior density and an interior action $d$, split the expected loss at $d$:

$$
\rho(d\mid x)
=\int_{-\infty}^d(d-\theta)\pi(\theta\mid x)\,d\theta
+\int_d^\infty(\theta-d)\pi(\theta\mid x)\,d\theta.
$$

Differentiating the two pieces gives

$$
\rho'(d\mid x)
=F_{\Theta\mid x}(d)-[1-F_{\Theta\mid x}(d)]
=2F_{\Theta\mid x}(d)-1.
$$

The expected loss is convex. A value with CDF $1/2$ minimizes it, so the Bayes action is a posterior median.

For a posterior with atoms, the derivative argument is replaced by the median condition

$$
P(\Theta\le d\mid x)\ge\frac12,
\qquad
P(\Theta\ge d\mid x)\ge\frac12.
$$

Any posterior median satisfying these inequalities minimizes absolute error. There can be more than one.

### Asymmetric error: derive a posterior quantile

Suppose underestimating costs $c_U>0$ per unit and overestimating costs $c_O>0$ per unit:

$$
L(\theta,d)
=c_U(\theta-d)_+ +c_O(d-\theta)_+,
\qquad z_+=\max(z,0).
$$

For a continuous posterior, differentiation gives

$$
\rho'(d\mid x)
=-c_U[1-F_{\Theta\mid x}(d)]
+c_OF_{\Theta\mid x}(d).
$$

Setting this equal to zero yields

$$
\boxed{
F_{\Theta\mid x}(d)=\frac{c_U}{c_U+c_O}.
}
$$

If underestimation is more costly, the selected quantile is above $1/2$. This direction check helps prevent reversing the cost labels. Equal costs recover the median.

### Discrete hypotheses and binary decisions

If a parameter takes finitely many possible values and loss is 0 for a correct choice and 1 for an incorrect one, the posterior expected loss of selecting $\theta_j$ is

$$
1-P(\Theta=\theta_j\mid x).
$$

Selecting a highest-posterior-probability value minimizes this loss. For a continuous parameter, literal exact-match 0–1 loss gives loss 1 almost surely for every fixed $d$, so it does not single out the density mode. Continuous MAP is a density maximizer, not a universally loss-optimal estimate.

For a yes/no decision, let $H$ be the event of interest and $q=P(H\mid x)$. Suppose a false positive costs $c_{\mathrm{FP}}>0$, a false negative costs $c_{\mathrm{FN}}>0$, and correct decisions cost zero.

The posterior losses are

$$
\rho(\text{yes}\mid x)=c_{\mathrm{FP}}(1-q),
\qquad
\rho(\text{no}\mid x)=c_{\mathrm{FN}}q.
$$

Choose yes when

$$
\boxed{
q>\frac{c_{\mathrm{FP}}}
{c_{\mathrm{FP}}+c_{\mathrm{FN}}}.
}
$$

At equality, either action has the same expected loss. A threshold of $1/2$ is justified by equal error costs, not by probability theory alone.

## Credible intervals and posterior prediction

### Credible intervals are posterior probability statements

A $100(1-\gamma)\%$ credible set $C(x)$ satisfies

$$
P(\Theta\in C(x)\mid x)=1-\gamma,
$$

or at least this probability when discreteness prevents exact equality. Its probability refers to the posterior conditional on the observed data and the specified model.

For a continuous posterior with CDF $F_x$, the equal-tailed interval is

$$
\boxed{
\left[
F_x^{-1}(\gamma/2),
F_x^{-1}(1-\gamma/2)
\right].
}
$$

The two excluded tails each have probability $\gamma/2$. A highest-posterior-density region instead selects points above an appropriate density threshold; it can be disconnected for a multimodal posterior and need not equal the equal-tailed interval. Neither construction should be automatically replaced by mean plus or minus two posterior standard deviations for a skewed distribution.

A frequentist $95\%$ confidence procedure instead aims to contain a fixed parameter in $95\%$ of repeated samples under its assumptions. A particular Bayesian $95\%$ credible interval does not automatically have $95\%$ frequentist coverage for every parameter value.

### Posterior predictive distribution: average over parameter uncertainty

Let $\widetilde X$ denote a future observation. If future and observed data are conditionally independent given $\Theta$, then

$$
\boxed{
p(\widetilde x\mid x)
=\int p(\widetilde x\mid\theta)
\pi(\theta\mid x)\,d\theta.
}
$$

This calculation averages the future sampling distribution over plausible parameter values. Replacing $\theta$ by a single estimate is a plug-in approximation and generally gives a different distribution.

For a single future Bernoulli trial,

$$
P(\widetilde Y=1\mid x)
=E[\Theta\mid x].
$$

This equality holds because the conditional success probability is linear in $\theta$. It does not imply that every nonlinear predictive probability is obtained by plugging in the posterior mean.

### Future binomial counts have a beta-binomial law

If

$$
\Theta\mid x\sim\mathrm{Beta}(A,B),
\qquad
\widetilde X\mid\Theta=\theta\sim\mathrm{Binomial}(m,\theta),
$$

with the same shared success probability and conditional independence from past data, then

$$
\begin{aligned}
P(\widetilde X=k\mid x)
&=\int_0^1\binom mk\theta^k(1-\theta)^{m-k}
\frac{\theta^{A-1}(1-\theta)^{B-1}}{B(A,B)}\,d\theta\\
&=\boxed{\binom mk\frac{B(A+k,B+m-k)}{B(A,B)}},
\end{aligned}
$$

for $k=0,\ldots,m$. This is the posterior predictive beta-binomial distribution. The same formula with prior shapes describes prior prediction, but the information being conditioned on differs.

Let $\mu=A/S$, $v=AB/[S^2(S+1)]$, and $S=A+B$. Total expectation and total variance give

$$
E[\widetilde X\mid x]=m\mu,
$$

$$
\begin{aligned}
\mathrm{Var}(\widetilde X\mid x)
&=E[m\Theta(1-\Theta)\mid x]
+\mathrm{Var}(m\Theta\mid x)\\
&=m\{\mu(1-\mu)-v\}+m^2v\\
&=m\mu(1-\mu)+m(m-1)v\\
&=\boxed{\frac{mAB(S+m)}{S^2(S+1)}}.
\end{aligned}
$$

The first line separates remaining trial randomness from parameter uncertainty. For $m>1$, this exceeds the variance $m\mu(1-\mu)$ of a plug-in binomial unless $v=0$. For one future Bernoulli observation, the predictive variance is exactly $\mu(1-\mu)$; do not add an extra $v$ to that value.

For distinct future trials,

$$
\mathrm{Cov}(\widetilde Y_i,\widetilde Y_j\mid x)
=\mathrm{Var}(\Theta\mid x).
$$

Their shared uncertain parameter creates posterior predictive dependence, even though they are independent conditional on the parameter. Drawing a new unrelated parameter independently for each future trial would describe a different model.

## Extension: Poisson-gamma inference

This section connects Bayesian inference to the earlier Poisson-process guides. It is an extension beyond the supplied lecture's beta-binomial and loss-function core, and all gamma distributions below use shape and rate, not shape and scale.

### Sampling model, exposure, and prior

Suppose disjoint interval counts are conditionally independent given a common rate $\Lambda=\lambda$:

$$
X_i\mid\lambda\sim\mathrm{Poisson}(\lambda t_i),
\qquad t_i>0.
$$

The known $t_i$ are observation durations or exposures in a common unit. Define

$$
s=\sum_{i=1}^nx_i,\qquad T=\sum_{i=1}^nt_i.
$$

The full likelihood and its parameter-dependent part are

$$
\begin{aligned}
L(\lambda;\boldsymbol x)
&=\prod_{i=1}^n
\frac{e^{-\lambda t_i}(\lambda t_i)^{x_i}}{x_i!}\\
&=\left(\prod_{i=1}^n\frac{t_i^{x_i}}{x_i!}\right)
\lambda^se^{-T\lambda}
\propto\lambda^se^{-T\lambda}.
\end{aligned}
$$

Take a gamma prior

$$
\pi(\lambda)
=\frac{b^a}{\Gamma(a)}
\lambda^{a-1}e^{-b\lambda}\mathbf1_{\{\lambda>0\}},
\qquad a,b>0.
$$

Its mean is $a/b$ and variance $a/b^2$. Substituting $z=b\lambda$ in the defining integrals gives these moments using the gamma recurrence.

### Posterior and units

Multiplication gives

$$
\pi(\lambda\mid\boldsymbol x)
\propto\lambda^{a+s-1}e^{-(b+T)\lambda}.
$$

This is a gamma kernel with shape $a+s$ and rate $b+T$. Hence

$$
\boxed{
\Lambda\mid\boldsymbol x
\sim\mathrm{Gamma}(a+s,b+T).
}
$$

Writing $a_*=a+s$, $b_*=b+T$,

$$
E[\Lambda\mid\boldsymbol x]=\frac{a_*}{b_*},
\qquad
\mathrm{Var}(\Lambda\mid\boldsymbol x)
=\frac{a_*}{b_*^2}.
$$

The posterior mean is a weighted average of the prior mean and observed event rate:

$$
\frac{a+s}{b+T}
=\frac b{b+T}\frac ab
+\frac T{b+T}\frac sT.
$$

If $\lambda$ is measured per hour, $b,T$ have hour units, and $a,s$ are dimensionless. The parameter called the gamma rate $b$ is not itself the Poisson event rate $\lambda$; they play different roles.

### Predict a count in a future interval

Let a new interval have duration $u>0$. Given the rate, its count satisfies

$$
\widetilde X\mid\lambda\sim\mathrm{Poisson}(\lambda u).
$$

Integrating over the posterior gives

$$
\begin{aligned}
P(\widetilde X=k\mid\boldsymbol x)
&=\frac{u^kb_*^{a_*}}{k!\Gamma(a_*)}
\int_0^\infty
\lambda^{a_*+k-1}e^{-(b_*+u)\lambda}\,d\lambda\\
&=\boxed{
\frac{\Gamma(a_*+k)}{\Gamma(a_*)k!}
\left(\frac{b_*}{b_*+u}\right)^{a_*}
\left(\frac{u}{b_*+u}\right)^k}.
\end{aligned}
$$

This is a negative-binomial law on $k=0,1,\ldots$, using R's convention `size = a_star`, `prob = b_star / (b_star + u)`. The formula remains valid for noninteger positive $a_*$.

Total expectation and total variance give

$$
E[\widetilde X\mid\boldsymbol x]=\frac{ua_*}{b_*},
\qquad
\mathrm{Var}(\widetilde X\mid\boldsymbol x)
=\frac{ua_*}{b_*}+\frac{u^2a_*}{b_*^2}.
$$

The extra term comes from uncertainty about the rate. A plug-in Poisson distribution omits it.

For the next waiting time $W$, conditional on $\lambda$ the survival function is $e^{-\lambda w}$. Averaging gives

$$
P(W>w\mid\boldsymbol x)
=E[e^{-\Lambda w}\mid\boldsymbol x]
=\left(\frac{b_*}{b_*+w}\right)^{a_*},
\qquad w\ge0.
$$

This predictive waiting law is not exponential. When $a_*>1$, its mean is $b_* /(a_*-1)$, obtained by integrating the survival function; it is not $1/E[\Lambda\mid\boldsymbol x]$. Thus exponential memorylessness conditional on a known rate need not survive mixing over an uncertain rate.

## Extension: normal-normal inference

This section provides a second conjugate model with a different algebraic shape. The unknown quantity is a mean, and the observation variance is assumed known; the same known-variance conjugate pair is covered in [MIT's beta and normal notes](https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_class15-prep.pdf).

### Normal distribution and parameter conventions

The notation $Z\sim N(m,v)$ means mean $m$ and variance $v>0$. Its density is

$$
f_Z(z)=\frac1{\sqrt{2\pi v}}
\exp\left[-\frac{(z-m)^2}{2v}\right],
\qquad z\in\mathbb R.
$$

The standard deviation is $\sqrt v$, which is the quantity R's `dnorm`, `pnorm`, and `qnorm` require as `sd`. Precision is the reciprocal variance, $1/v$, not the reciprocal standard deviation.

Assume

$$
X_1,\ldots,X_n\mid\Theta=\theta
\ \text{are iid }N(\theta,\sigma^2),
\qquad \sigma^2>0\text{ known},
$$

and take

$$
\Theta\sim N(m_0,\tau_0^2),
\qquad \tau_0^2>0.
$$

Here $\overline x=n^{-1}\sum_i x_i$. Because the data share the same conditional mean,

$$
\sum_{i=1}^n(x_i-\theta)^2
=\sum_{i=1}^n(x_i-\overline x)^2
+n(\overline x-\theta)^2.
$$

To verify this identity, expand $x_i-\theta=(x_i-\overline x)+(\overline x-\theta)$. The cross terms sum to zero because $\sum_i(x_i-\overline x)=0$.

### Complete the square in the parameter

Ignoring factors independent of $\theta$, the posterior is

$$
\pi(\theta\mid\boldsymbol x)
\propto
\exp\left[-\frac12\left\{
\frac{n(\theta-\overline x)^2}{\sigma^2}
+\frac{(\theta-m_0)^2}{\tau_0^2}
\right\}\right].
$$

The quadratic expression inside the braces has the form

$$
A\theta^2-2B\theta+C,
$$

where

$$
A=\frac1{\tau_0^2}+\frac n{\sigma^2},
\qquad
B=\frac{m_0}{\tau_0^2}
+\frac{n\overline x}{\sigma^2}.
$$

Complete the square:

$$
A\theta^2-2B\theta+C
=A\left(\theta-\frac BA\right)^2
+C-\frac{B^2}{A}.
$$

The final two terms do not depend on $\theta$, so the normalized posterior is

$$
\boxed{
\Theta\mid\boldsymbol x\sim N(m_n,v_n),
\quad
v_n=\left(\frac1{\tau_0^2}+\frac n{\sigma^2}\right)^{-1},
\quad
m_n=v_n\left(\frac{m_0}{\tau_0^2}
+\frac{n\overline x}{\sigma^2}\right).
}
$$

Prior precision and data precision add. The posterior mean is a precision-weighted average of $m_0$ and $\overline x$. Within this model, $v_n$ is smaller than both $\tau_0^2$ and $\sigma^2/n$; do not generalize that pointwise variance-reduction statement to all Bayesian models.

### Predict a new observation, not just its mean

For a future observation from the same process, conditionally independent of past data given $\Theta$,

$$
\boxed{
\widetilde X\mid\boldsymbol x
\sim N(m_n,\sigma^2+v_n).
}
$$

One way to see the distribution is to write $\widetilde X=\Theta+\varepsilon$, where $\varepsilon\sim N(0,\sigma^2)$ is independent of $\Theta$ and the observed data. Conditional on the data, this is a sum of independent normal quantities.

The variance has two parts:

$$
\mathrm{Var}(\widetilde X\mid\boldsymbol x)
=E[\sigma^2\mid\boldsymbol x]
+\mathrm{Var}(\Theta\mid\boldsymbol x)
=\sigma^2+v_n.
$$

A credible interval for the mean uses $v_n$, while a predictive interval for a new observation uses $\sigma^2+v_n$. If $\sigma^2$ is unknown, these formulas are not a full Bayesian treatment of that unknown variance; an additional prior and a different derivation are required.

## A reliable problem-solving procedure

Start by saying what the unknown parameter means and what has actually been observed. A correct conjugate formula used with the wrong experiment or support still gives the wrong answer.

1. **Identify the target:** A posterior distribution, probability, estimate, interval, future count, or decision.
2. **State the sampling model:** Include conditional independence, fixed sample size or exposure, parameter support, and parameterization.
3. **Write the likelihood:** Fix the observed data and retain all parameter-dependent factors, including support indicators.
4. **Write the prior:** Include its support and hyperparameters. Distinguish it from the sampling distribution.
5. **Multiply and simplify:** Form the posterior kernel as a function of the parameter.
6. **Normalize:** Recognize a proper distribution or compute a finite positive normalizing integral or sum.
7. **Answer the actual question:** Integrate the posterior, use its moments or quantiles, minimize the specified loss, or integrate a future sampling model.
8. **Interpret and check:** Verify support, total mass 1, probability bounds, parameter units, and an independent moment or numerical calculation.

An exam-ready template is

$$
\begin{gathered}
L(\theta;x)=f(x\mid\theta),\qquad
\pi(\theta)=\text{prior on }\mathcal T,\\
\pi(\theta\mid x)\propto L(\theta;x)\pi(\theta),\\
\pi(\theta\mid x)=
\frac{L(\theta;x)\pi(\theta)}
{\int_{\mathcal T}L(u;x)\pi(u)\,du},\\
\text{requested summary}
=\text{an operation under this posterior}.
\end{gathered}
$$

For a discrete parameter, use a sum. For an estimate, name the loss before calling it optimal. For a future observation, average its sampling model over the posterior instead of treating an estimate as certain.

## Seven fully worked examples

### Worked example A: manuscript authorship from the lecture

An old manuscript has 1,000 words. The two possible authors are $H$ and $M$, with prior probabilities $P(H)=3/4$ and $P(M)=1/4$. The marker-word count satisfies

$$
X\mid H\sim\mathrm{Poisson}(3),
\qquad
X\mid M\sim\mathrm{Poisson}(6).
$$

The observed manuscript contains six uses of the marker word. Find $P(H\mid X=6)$, and explain how the data changed the prior judgment.

The Poisson model is an explicit assumption from the lecture. Its means refer to counts in this 1,000-word manuscript; comparing texts of substantially different lengths would require an exposure adjustment rather than blindly reusing these means.

The two likelihoods are

$$
\begin{aligned}
P(X=6\mid H)
&=e^{-3}\frac{3^6}{6!}
\approx0.0504094067,\\
P(X=6\mid M)
&=e^{-6}\frac{6^6}{6!}
\approx0.1606231410.
\end{aligned}
$$

The observation is more likely under $M$. However, likelihoods are not posterior probabilities, so include the starting prior weights.

| Author | Prior | Likelihood at $x=6$ | Prior × likelihood | Posterior |
|---|---:|---:|---:|---:|
| $H$ | 0.75 | 0.0504094067 | 0.0378070550 | 0.4849368609 |
| $M$ | 0.25 | 0.1606231410 | 0.0401557853 | 0.5150631391 |

The denominator is

$$
m(6)=\frac34P(X=6\mid H)+\frac14P(X=6\mid M)
\approx0.0779628403.
$$

Thus

$$
\boxed{
P(H\mid X=6)
=\frac{0.0378070550}{0.0779628403}
\approx0.4849368609.
}
$$

As an independent check, the likelihood ratio is

$$
\mathrm{BF}_{H:M}
=\frac{e^{-3}3^6/6!}{e^{-6}6^6/6!}
=\frac{e^3}{2^6}
\approx0.3138365144.
$$

Multiplying by prior odds 3 gives posterior odds about $0.9415095433$. Converting odds $o$ to probability $o/(1+o)$ recovers $0.4849368609$.

The posterior favors $M$ only slightly, despite starting with a 3:1 preference for $H$. This is a model-conditional conclusion, not proof of authorship. Treating the unknown Poisson mean as a parameter gives an equivalent two-point posterior: mass $0.4849368609$ at $\lambda=3$ and mass $0.5150631391$ at $\lambda=6$.

### Worked example B: the lecture's uniform prior and eight trials

Assume $\Theta\sim\mathrm{Uniform}(0,1)$, and $X\mid\theta\sim\mathrm{Binomial}(8,\theta)$. We observe $x=6$. Find the posterior, $P(\Theta\le1/2\mid X=6)$, its mean and variance, and a $95\%$ equal-tailed credible interval.

The binomial model requires eight conditionally independent trials sharing a fixed success probability. The uniform prior has density 1 on $(0,1)$, so the posterior kernel is

$$
\pi(\theta\mid6)\propto
\theta^6(1-\theta)^2.
$$

Because a beta kernel has exponents one less than its shapes, the shapes are 7 and 3, not 6 and 2. The evidence and beta integral are

$$
m(6)=\frac19,\qquad
B(7,3)=\frac{6!2!}{9!}=\frac1{252}.
$$

Therefore,

$$
\boxed{
\Theta\mid6\sim\mathrm{Beta}(7,3),\qquad
\pi(\theta\mid6)=252\theta^6(1-\theta)^2
\ \text{on }(0,1).
}
$$

To compute the requested probability by hand, expand the polynomial:

$$
\begin{aligned}
P(\Theta\le1/2\mid6)
&=252\int_0^{1/2}
(\theta^6-2\theta^7+\theta^8)\,d\theta\\
&=252\left[
\frac{\theta^7}{7}
-\frac{\theta^8}{4}
+\frac{\theta^9}{9}
\right]_0^{1/2}\\
&=252\left(\frac1{896}-\frac1{1024}+\frac1{4608}\right)\\
&=\boxed{\frac{23}{256}=0.08984375}.
\end{aligned}
$$

The prior probability was $1/2$; the posterior probability is about $0.08984$. The observed successes moved substantial mass to values above $1/2$. The slide reports this probability approximately; the fraction above is the exact value.

The posterior moments are

$$
E[\Theta\mid6]=\frac7{10}=0.7,
\qquad
\mathrm{Var}(\Theta\mid6)
=\frac{7\cdot3}{10^2\cdot11}
=\frac{21}{1100}\approx0.01909091.
$$

The posterior standard deviation is approximately $0.13817$. The sample proportion and interior MAP are $6/8=0.75$, so they are not the posterior mean.

The equal-tailed interval is

$$
\left[
Q_{\mathrm{Beta}(7,3)}(0.025),
Q_{\mathrm{Beta}(7,3)}(0.975)
\right]
\approx\boxed{[0.399906,\ 0.925145]}.
$$

Its posterior probability is $0.95$. The R check evaluates the two quantiles and verifies the included CDF mass. This is an interval for the success probability, not a predicted count of future successes.

### Worked example C: informative prior, shrinkage, and sequential batches

Suppose the prior is $\mathrm{Beta}(2,8)$, and ten conditionally independent trials yield six successes. Find the posterior mean and explain how it reconciles the prior mean $0.2$ with the observed proportion $0.6$.

The likelihood has six success factors and four failure factors. Multiplying by the prior gives

$$
\pi(\theta\mid x)
\propto
\theta^{6+2-1}(1-\theta)^{4+8-1},
$$

so

$$
\Theta\mid x\sim\mathrm{Beta}(8,12).
$$

The mean and variance are

$$
E[\Theta\mid x]=\frac8{20}=0.4,\qquad
\mathrm{Var}(\Theta\mid x)
=\frac{8\cdot12}{20^2\cdot21}
=\frac2{175}\approx0.01142857.
$$

Here the prior concentration and new sample size both equal 10:

$$
E[\Theta\mid x]
=\frac{10}{20}(0.2)+\frac{10}{20}(0.6)=0.4.
$$

The mean lies between the prior mean and sample proportion, as it should. This prior is not weak relative to a sample of ten; its influence is visible.

Suppose the same data arrive as three successes in four trials, followed by three successes in six new trials. The first posterior is $\mathrm{Beta}(5,9)$, and the second update gives

$$
\mathrm{Beta}(5+3,\ 9+3)=\mathrm{Beta}(8,12).
$$

The batch and pooled answers agree because both use six successes and four failures exactly once. This agreement does not justify pooling batches that actually have different success probabilities.

For the same six-of-ten likelihood, consider three priors:

| Prior | Prior mean | Posterior | Posterior mean | $P(\Theta>1/2\mid x)$ |
|---|---:|---|---:|---:|
| $\mathrm{Beta}(1,1)$ | 0.5 | $\mathrm{Beta}(7,5)$ | 0.583333 | 0.725586 |
| $\mathrm{Beta}(2,8)$ | 0.2 | $\mathrm{Beta}(8,12)$ | 0.400000 | 0.179642 |
| $\mathrm{Beta}(20,80)$ | 0.2 | $\mathrm{Beta}(26,84)$ | 0.236364 | approximately $6.08\times10^{-9}$ |

The likelihood is identical in all rows. The strong last prior is a deliberate stress test of sensitivity, not a recommended default. If different defensible priors change the conclusion materially, report that dependence instead of hiding it behind a single posterior mean.

### Worked example D: one posterior, several losses and targets

Reuse the $\mathrm{Beta}(7,3)$ posterior from example B. Compare squared-error and absolute-error estimates of $\Theta$, estimate $\Theta^2$ under squared error, and decide whether to assert $\Theta>1/2$ under asymmetric decision costs.

Under squared error, the optimal estimate of $\Theta$ is $d=0.7$. For any candidate $d$,

$$
\rho(d\mid x)=\frac{21}{1100}+(d-0.7)^2.
$$

For example, using the MAP $d=0.75$ instead adds $0.05^2=0.0025$ to posterior expected squared loss. This does not make MAP a wrong summary; it shows that it is not the squared-loss optimum.

Under absolute error, the Bayes estimate is the posterior median:

$$
d=Q_{\mathrm{Beta}(7,3)}(0.5)\approx0.713763.
$$

The median differs from both the mean $0.7$ and mode $0.75$. A request for “the Bayes estimate” is incomplete unless a loss or decision criterion is specified.

For the target $\psi=\Theta^2$, squared-loss optimality gives

$$
\widehat\psi_{\mathrm B}
=E[\Theta^2\mid x]
=\frac{7\cdot8}{10\cdot11}
=\boxed{\frac{28}{55}\approx0.509091}.
$$

Squaring the posterior mean would give $0.49$, which is smaller by $21/1100$, the posterior variance. This is a concrete reason not to replace posterior expectations of nonlinear functions by plug-in transformations.

For the event $H=\{\Theta>1/2\}$,

$$
q=P(H\mid x)=1-\frac{23}{256}=\frac{233}{256}
\approx0.91015625.
$$

If a false positive costs 9 and a false negative costs 1, the yes threshold is $9/10$. The posterior losses are

$$
\rho(\text{yes}\mid x)=9\frac{23}{256}
=\frac{207}{256},
\qquad
\rho(\text{no}\mid x)=\frac{233}{256}.
$$

Choose yes because $207/256<233/256$. If the false-positive cost rises to 19 while the false-negative cost remains 1, the threshold is $19/20=0.95$, and the same posterior now favors no. The inference has not changed; the decision costs have.

### Worked example E: predict four future trials

Again take $\Theta\mid x\sim\mathrm{Beta}(7,3)$. Let $\widetilde X$ count successes in four new conditionally independent trials sharing this same probability. Find its distribution, probability of four successes, mean, and variance.

Integrating the binomial likelihood over the posterior gives

$$
P(\widetilde X=k\mid x)
=\binom4k\frac{B(7+k,7-k)}{B(7,3)},
\qquad k=0,\ldots,4.
$$

The second beta argument is $3+4-k=7-k$. The resulting PMF is

| $k$ | 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|
| $P(\widetilde X=k\mid x)$ | $3/143$ | $14/143$ | $168/715$ | $252/715$ | $42/143$ |

Converting to denominator 715 gives numerators $15,70,168,252,210$, which sum to 715. This verifies normalization without rounding.

For four successes,

$$
\begin{aligned}
P(\widetilde X=4\mid x)
&=E[\Theta^4\mid x]\\
&=\frac{7\cdot8\cdot9\cdot10}{10\cdot11\cdot12\cdot13}\\
&=\boxed{\frac{42}{143}\approx0.293706}.
\end{aligned}
$$

The plug-in calculation $0.7^4=0.2401$ is different. Although a single future success probability equals the posterior mean, the probability of four successes is a fourth moment.

The predictive mean is

$$
E[\widetilde X\mid x]=4(0.7)=2.8.
$$

Its variance is

$$
\begin{aligned}
\mathrm{Var}(\widetilde X\mid x)
&=4(0.7)(0.3)+4(3)\frac{21}{1100}\\
&=\boxed{\frac{294}{275}\approx1.069091}.
\end{aligned}
$$

The plug-in binomial variance is only $0.84$. Both distributions have the same mean, but the posterior predictive law accounts for the shared uncertainty and dependence among the four future outcomes.

### Worked example F: update an uncertain Poisson rate

This is an extension connecting to Poisson processes. Suppose $\Lambda\sim\mathrm{Gamma}(2,1)$ in shape-rate form, with the event rate measured per hour, and we observe nine events over a fixed three-hour window.

The count likelihood is

$$
L(\lambda;9)
=e^{-3\lambda}\frac{(3\lambda)^9}{9!}
\propto \lambda^9e^{-3\lambda}.
$$

The prior kernel is $\lambda^{2-1}e^{-\lambda}$. Their product is $\lambda^{11-1}e^{-4\lambda}$, so

$$
\boxed{\Lambda\mid x\sim\mathrm{Gamma}(11,4)}.
$$

The posterior mean, variance, and interior mode are

$$
E[\Lambda\mid x]=\frac{11}{4}=2.75,
\qquad
\mathrm{Var}(\Lambda\mid x)=\frac{11}{16}=0.6875,
\qquad
\lambda_{\mathrm{MAP}}=\frac{10}{4}=2.5.
$$

The MLE is $9/3=3$ events per hour. The posterior mean can be checked through weighting:

$$
\frac14(2)+\frac34(3)=2.75.
$$

Posterior uncertainty about exceeding three events per hour is

$$
P(\Lambda>3\mid x)\approx0.347229,
$$

and the $95\%$ equal-tailed credible interval is approximately $[1.372790,4.597589]$ events per hour. The R checks use gamma tail and quantile functions with `rate = 4`.

For the number of events in the next two hours,

$$
\widetilde X\mid x
\sim\mathrm{NegBin}\left(
\mathrm{size}=11,\ \mathrm{prob}=\frac46
\right).
$$

Its mean and variance are

$$
E[\widetilde X\mid x]=2\frac{11}{4}=5.5,
\qquad
\mathrm{Var}(\widetilde X\mid x)
=5.5+4\frac{11}{16}=8.25.
$$

The probability of no future events is

$$
P(\widetilde X=0\mid x)
=\left(\frac46\right)^{11}
\approx0.0115610.
$$

Using a Poisson model with fixed mean 5.5 would miss the extra variance from rate uncertainty. For the next event's waiting time, the posterior predictive mean is $4/(11-1)=0.4$ hour, not $1/2.75\approx0.363636$ hour.

### Worked example G: normal mean and a future measurement

This known-variance extension uses prior $\Theta\sim N(0,4)$ and conditionally iid observations $X_i\mid\theta\sim N(\theta,4)$. Four observed values are $2,3,3,4$, so $\overline x=3$.

The prior precision is $1/4$, and the data precision is $n/\sigma^2=4/4=1$. Therefore,

$$
v_n=\left(\frac14+1\right)^{-1}=\frac45,
$$

$$
m_n=\frac45\left(\frac04+\frac{4\cdot3}{4}\right)
=\frac{12}{5}=2.4.
$$

Thus

$$
\boxed{\Theta\mid\boldsymbol x\sim N(2.4,0.8)}.
$$

The posterior mean is $0.2(0)+0.8(3)=2.4$. Data have greater weight because their precision is larger, and the result lies between the prior mean and sample mean.

The $95\%$ equal-tailed credible interval for $\Theta$ is

$$
2.4\pm z_{0.975}\sqrt{0.8}
\approx[0.646955,\ 4.153045],
$$

where $z_{0.975}\approx1.959964$ is the $0.975$-quantile of a standard normal distribution. Here symmetry makes this interval calculation exact up to numerical quantiles.

A new individual measurement has posterior predictive distribution

$$
\widetilde X\mid\boldsymbol x\sim N(2.4,4.8),
$$

and a $95\%$ predictive interval approximately

$$
[-1.894066,\ 6.694066].
$$

The predictive interval is wider because an individual measurement has observation noise in addition to uncertainty about the common mean. The known variance is 4, while its standard deviation is 2; interchanging them in R would change the model.

## Graduate reference and modeling cautions

### Bayes' theorem as an update of a prior measure

This subsection is a second-pass connection to the earlier measure-theoretic material. It is not needed to carry out the lecture's beta-binomial computations.

Let $\Pi$ be a proper prior probability measure on the parameter space. Suppose a sampling family is dominated by a common reference measure for the data and has a jointly measurable nonnegative density $f(x\mid\theta)$.

For observed $x$ with

$$
0<m(x)=\int f(x\mid\theta)\,\Pi(d\theta)<\infty,
$$

the posterior measure is

$$
\Pi(A\mid x)
=\frac{\int_Af(x\mid\theta)\,\Pi(d\theta)}
{m(x)}.
$$

Its Radon-Nikodym derivative relative to the prior is

$$
\frac{d\Pi(\cdot\mid x)}{d\Pi}(\theta)
=\frac{f(x\mid\theta)}{m(x)}.
$$

This formulation covers discrete priors, continuous priors, and mixtures with atoms. It also makes prior support explicit: under this finite-normalizer construction, a set with zero prior probability has zero posterior probability.

The posterior is a regular conditional law in the usual real or Euclidean settings, and statements involving versions retain the almost-everywhere qualification from the conditional-distribution guides. An observed $x$ with zero or infinite displayed normalizer cannot be handled by blindly dividing.

### Sufficient statistics and discarded data

A statistic $T(X)$ is sufficient for a parameter when it preserves the parameter-related information in the model. In the dominated models used here, a factorization

$$
f(\boldsymbol x\mid\theta)
=h(\boldsymbol x)g(T(\boldsymbol x),\theta)
$$

shows why a posterior depends on $\boldsymbol x$ only through $T(\boldsymbol x)$, whenever the normalizing calculation is valid. The factor $h$, independent of the parameter, cancels.

For iid Bernoulli data, the success count is sufficient. For the known-variance normal mean model with fixed $n$, the sample mean is sufficient. For a homogeneous Poisson process on a fixed observation window, the count contains the rate information; conditional event locations do not add information about the constant rate.

Sufficiency does not mean the remaining data are useless for every purpose. Event times can help reveal a time-varying rate that violates homogeneity, and a sequence can reveal dependence that a success count hides. A statistic sufficient within a model does not validate that model.

### Posterior loss, frequentist risk, and Bayes risk

These are three different averages. Keeping their conditioning explicit prevents the word “risk” from hiding which quantity is random.

For a rule $\delta$,

$$
\begin{aligned}
\rho(\delta(x)\mid x)
&=E[L(\Theta,\delta(x))\mid X=x],
&&\text{average over the posterior},\\
R(\theta,\delta)
&=E_{X\mid\theta}[L(\theta,\delta(X))],
&&\text{average over repeated data},\\
r(\Pi,\delta)
&=\int R(\theta,\delta)\,\Pi(d\theta),
&&\text{average repeated-data risk over the prior}.
\end{aligned}
$$

For nonnegative loss, Tonelli's theorem allows changing integration order. With a proper prior,

$$
\begin{aligned}
r(\Pi,\delta)
&=\int\!\int
L(\theta,\delta(x))f(x\mid\theta)\,\Pi(d\theta)\,dx\\
&=\int \rho(\delta(x)\mid x)m(x)\,dx.
\end{aligned}
$$

Use a sum instead of $dx$ for discrete data. A measurable rule that minimizes posterior expected loss at almost every data value minimizes this average Bayes risk, provided the minimization and integrals are well defined. This explains the relationship between the pointwise posterior calculation in the lecture and the broader decision-theoretic term “Bayes estimator.”

### Posterior variance need not decrease at every observed value

For a proper joint model with finite parameter variance, total variance gives

$$
\mathrm{Var}(\Theta)
=E_X[\mathrm{Var}(\Theta\mid X)]
+\mathrm{Var}_X(E[\Theta\mid X]).
$$

Therefore,

$$
E_X[\mathrm{Var}(\Theta\mid X)]
\le\mathrm{Var}(\Theta).
$$

This is an average over possible data under the prior predictive distribution. It does not imply $\mathrm{Var}(\Theta\mid X=x)\le\mathrm{Var}(\Theta)$ for every surprising $x$. Practice problem 16 gives an explicit beta-Bernoulli counterexample.

### Improper priors require a propriety check

A function called an improper prior does not integrate to 1, so it is not a prior probability distribution. It can sometimes be used as a formal weight if multiplying by the likelihood produces a finite positive posterior normalizing integral, but this must be checked.

For example, the weight

$$
\pi(\theta)\propto\theta^{-1}(1-\theta)^{-1},
\qquad 0<\theta<1,
$$

is improper. With binomial data, the posterior kernel is

$$
\theta^{x-1}(1-\theta)^{n-x-1}.
$$

This produces a proper beta posterior only when $0<x<n$. At $x=0$ or $x=n$, an endpoint integral diverges. Writing “$\mathrm{Beta}(0,0)$” as though it were an ordinary distribution conceals the problem.

An arbitrary multiplicative constant in an improper prior prevents an ordinary marginal likelihood or Bayes factor from being defined without additional construction. A proper posterior alone does not repair that missing prior normalization.

### Continuous priors and exact point hypotheses

Under an absolutely continuous prior and likelihood update with finite normalizer,

$$
P(\Theta=1/2\mid x)=0.
$$

This does not mean the data disproved a fair-coin hypothesis. The model assigned that individual point zero probability in advance. To compare an exact point null with a continuous alternative using posterior model probabilities, give the point hypothesis positive prior mass.

For instance, a prior with probability $w$ at $\theta=1/2$ and probability $1-w$ distributed by a proper alternative density yields

$$
P(H_0\mid x)
=\frac{w\,p(x\mid\theta=1/2)}
{w\,p(x\mid\theta=1/2)+(1-w)m_1(x)},
$$

where $m_1(x)$ is the alternative's integrated likelihood. This is finite-hypothesis Bayes' theorem applied to models; the alternative likelihood is averaged over its prior rather than maximized.

### Transformations and MAP

For a differentiable one-to-one change $\eta=h(\theta)$ with inverse $\theta=g(\eta)$, the posterior density transforms as

$$
\pi_\eta(\eta\mid x)
=\pi_\Theta(g(\eta)\mid x)|g'(\eta)|.
$$

The Jacobian changes density heights, so a density maximizer need not transform into a maximizer in the new coordinate. Posterior probabilities of corresponding sets agree, because the change-of-variables formula preserves probability.

For example, the logit is

$$
\eta=\log\frac{\theta}{1-\theta},\qquad
\theta=\frac1{1+e^{-\eta}},\qquad
\frac{d\theta}{d\eta}=\theta(1-\theta).
$$

A uniform prior in $\theta$ is therefore not uniform in $\eta$. The Jacobian is not optional when constructing a density for the transformed parameter.

### Model checks and computation beyond conjugacy

A posterior is conditional on both likelihood and prior. Check whether trial independence, a shared success probability, exposure handling, or known normal variance are reasonable before treating a narrow posterior as strong evidence.

Prior predictive checks ask whether the prior-plus-sampling-model combination generates plausible data before seeing the study's results. Posterior predictive checks simulate replicated observations from the fitted model and compare relevant features with the observed data. Reusing the data to fit the model means these checks are not an independent validation experiment.

When no conjugate form exists, the same Bayes formula still applies. One-dimensional problems may be handled by numerical integration or a carefully normalized grid; higher-dimensional problems often need more advanced approximation or sampling methods. Independent direct draws from the beta, gamma, and normal posteriors in this guide are not Markov chain Monte Carlo, so they do not require chain convergence diagnostics.

The R section includes an interior midpoint grid and log weights. Grid weights approximate probability masses at cells, while density heights are masses divided by cell width. For nonuniform grids or endpoint-singular densities, quadrature weights and integration behavior require additional care.

## Practice problems

Try these before opening the solutions. For each calculation, state what is fixed, what remains uncertain, the support, and any required moment or independence assumptions.

1. **Finite hypotheses:** A success probability is $0.2$ or $0.8$, each with prior probability $1/2$. Observe one success. Find the posterior probabilities and the probability the next trial succeeds, assuming the same common parameter and conditional independence.
2. **Likelihood versus posterior:** For $X\mid\theta\sim\mathrm{Binomial}(8,\theta)$ and $x=6$, calculate $\int_0^1L(\theta;6)\,d\theta$. Explain why the likelihood is not a density for $\theta$, and give the posterior under a uniform prior.
3. **Conjugate update:** Use a $\mathrm{Beta}(3,2)$ prior and four successes in five trials. Derive the posterior, its mean, variance, and interior MAP.
4. **No observed successes:** Under a uniform prior, observe zero successes in five trials. Find the posterior mean and $P(\Theta>1/2\mid x)$. Explain why a posterior probability need not become exactly zero despite observing no successes.
5. **Sequential learning:** Start from $\mathrm{Beta}(2,3)$. Observe three successes in four trials, then two successes in six new trials. Compute both sequential posteriors and verify the final answer by pooling.
6. **A transformed target:** For posterior $\mathrm{Beta}(7,3)$, estimate $\Theta^2$ under squared error. Compare it with the square of the posterior mean and explain the difference.
7. **Asymmetric estimation:** Suppose the posterior is $\mathrm{Beta}(2,1)$, with density $2\theta$ on $(0,1)$. Find the mean, median, and Bayes estimate when underestimation costs four times as much per unit as overestimation.
8. **Binary decision:** The posterior probability of $H$ is 0.6. A false positive costs 4 and a false negative costs 1, with zero loss for correct decisions. Calculate both expected losses and choose an action.
9. **Two future trials:** For posterior $\mathrm{Beta}(7,3)$, find the probability both of two new trials succeed, and the mean and variance of the future success count. Compare the variance with a plug-in binomial.
10. **Prior prediction:** Under a uniform prior, find the prior predictive distribution, mean, and variance of a count in three conditionally independent Bernoulli trials. Why is it not $\mathrm{Binomial}(3,1/2)$?
11. **Poisson exposure extension:** The prior for a common rate is $\mathrm{Gamma}(3,2)$. Observe counts $1,3,2$ in disjoint windows of lengths $0.5,1,0.5$. Derive the posterior mean and variance using shape-rate notation.
12. **Normal extension:** Take prior $N(10,9)$ and one observation $x=13$ with known sampling variance 9. Find the parameter posterior and the posterior predictive distribution for one new observation.
13. **Point-null extension:** Give $H_0:\theta=1/2$ prior probability $1/2$, with a uniform prior on $\theta$ under the alternative $H_1$. For six successes in eight trials, compute the Bayes factor $H_0:H_1$ and $P(H_0\mid x)$. Contrast this with a purely continuous uniform prior.
14. **Improper-prior extension:** Use the formal prior weight $\theta^{-1}(1-\theta)^{-1}$. For $n=5$, determine whether a proper posterior exists when $x=0$ and when $x=2$.
15. **Transformation extension:** For posterior $\mathrm{Beta}(7,3)$, derive the density of the logit $\eta=\log(\theta/(1-\theta))$. Find its mode on the $\eta$-scale, and compare with the logit of the original posterior's mode.
16. **Does uncertainty always decrease?:** Start from $\mathrm{Beta}(1,9)$ and observe one Bernoulli trial. Compute the prior variance and the posterior variance after a success. Then average the success and failure posterior variances using their prior predictive probabilities, and verify total variance.
17. **Parameter-dependent support extension:** Let $\Theta\sim\mathrm{Uniform}(0,2)$, and let $X\mid\theta\sim\mathrm{Uniform}(0,\theta)$. Observe $x=1.5$. Derive the posterior support, normalized density, and mean; explain why the likelihood's support indicator cannot be discarded.

## Solutions and checks

### Finite hypotheses and basic conjugacy

- **Problem 1:** The prior-times-likelihood weights are $0.5(0.2)=0.1$ and $0.5(0.8)=0.4$, which sum to 0.5. The posterior weights are therefore $0.2$ and $0.8$. The next-success probability is the posterior average $0.2(0.2)+0.8(0.8)=0.68$; it is not just the most probable parameter value 0.8.
- **Problem 2:** The integral is $\binom86B(7,3)=28/252=1/9$, not 1. A uniform-prior posterior rescales the likelihood by 9, giving density $252\theta^6(1-\theta)^2$ and distribution $\mathrm{Beta}(7,3)$. Summing the sampling PMF over possible counts is a different operation from integrating the likelihood over parameter values.
- **Problem 3:** Multiplication gives kernel $\theta^{3+4-1}(1-\theta)^{2+1-1}$, so the posterior is $\mathrm{Beta}(7,3)$. Its mean is $7/10$, variance $21/1100$, and interior MAP $(7-1)/(10-2)=3/4$. The observed proportion is $4/5$, illustrating that a nonuniform prior can make both posterior mean and MAP differ from MLE.
- **Problem 4:** The posterior is $\mathrm{Beta}(1,6)$, with density $6(1-\theta)^5$. Its mean is $1/7$, and $\int_{1/2}^1 6(1-\theta)^5\,d\theta=(1/2)^6=1/64$. No successes provide evidence against large success probabilities, but finite data do not eliminate every positive probability value that had positive prior density.
- **Problem 5:** The first update gives $\mathrm{Beta}(2+3,3+1)=\mathrm{Beta}(5,4)$. The second gives $\mathrm{Beta}(5+2,4+4)=\mathrm{Beta}(7,8)$, with mean $7/15$. Pooling gives five successes and five failures, producing the same shapes $2+5=7$, $3+5=8$.

### Losses, targets, and prediction

For problem 6, apply squared-loss optimality to the target, not to an intermediate estimate:

$$
E[\Theta^2\mid x]=\frac{7\cdot8}{10\cdot11}
=\frac{28}{55}.
$$

The squared mean is $49/100$, so the difference is

$$
\frac{28}{55}-\frac{49}{100}
=\frac{21}{1100}
=\mathrm{Var}(\Theta\mid x).
$$

This is positive because the posterior is nondegenerate. It is an exact algebraic check, not a numerical approximation.

For problem 7, integrate $2\theta$ to get posterior CDF $F(t)=t^2$ for $0<t<1$. Thus

$$
E[\Theta\mid x]=2\int_0^1\theta^2\,d\theta=\frac23,
\qquad
\mathrm{median}=\frac1{\sqrt2}.
$$

With $c_U=4,c_O=1$, the relevant quantile is $4/5$, giving

$$
d=F^{-1}(4/5)=\sqrt{4/5}=\frac2{\sqrt5}.
$$

It exceeds the median because underestimation is more costly. Using the $1/5$-quantile would reverse the loss interpretation.

For problem 8, choosing yes has posterior expected loss $4(1-0.6)=1.6$, while choosing no has loss $1(0.6)=0.6$. Choose no; the yes threshold is $4/(4+1)=0.8$, not $0.5$.

For problem 9, the probability of both successes is a posterior second moment:

$$
P(\widetilde X=2\mid x)=E[\Theta^2\mid x]
=\frac{28}{55}.
$$

The mean and variance of the future count are

$$
E[\widetilde X\mid x]=2\frac7{10}=\frac75,
$$

$$
\mathrm{Var}(\widetilde X\mid x)
=2\frac7{10}\frac3{10}
+2(1)\frac{21}{1100}
=\frac{126}{275}\approx0.458182.
$$

The plug-in binomial variance is $0.42$, smaller by $42/1100$. Both future outcomes share the same uncertain probability.

For problem 10, the uniform-prior evidence calculation gives

$$
P(X=k)=\frac14,\qquad k=0,1,2,3.
$$

Therefore,

$$
E[X]=\frac{0+1+2+3}{4}=\frac32,\qquad
E[X^2]=\frac{0+1+4+9}{4}=\frac72,
$$

$$
\mathrm{Var}(X)=\frac72-\left(\frac32\right)^2=\frac54.
$$

A $\mathrm{Binomial}(3,1/2)$ has the same mean but variance $3/4$. Averaging over one shared uncertain success probability introduces positive dependence among the trials.

### Conjugate-model extensions

For problem 11, total count is $s=6$ and total exposure is $T=2$. The likelihood kernel is $\lambda^6e^{-2\lambda}$, so multiplication by $\lambda^{3-1}e^{-2\lambda}$ gives

$$
\Lambda\mid\boldsymbol x\sim\mathrm{Gamma}(9,4),
\quad
E[\Lambda\mid\boldsymbol x]=\frac94,
\quad
\mathrm{Var}(\Lambda\mid\boldsymbol x)=\frac9{16}.
$$

Adding the counts but forgetting exposures would give an incorrect posterior rate. Shape and rate receive different kinds of information.

For problem 12, prior and observation precisions are both $1/9$. Hence

$$
v_1=\left(\frac19+\frac19\right)^{-1}=\frac92,
\qquad
m_1=\frac92\left(\frac{10}{9}+\frac{13}{9}\right)
=\frac{23}{2}=11.5.
$$

The posterior is $N(11.5,4.5)$, and a new observation has distribution $N(11.5,9+4.5)=N(11.5,13.5)$. The extra 9 in predictive variance is observation noise, not additional uncertainty about the mean.

### Point nulls, propriety, and transformations

For problem 13, the null likelihood and alternative evidence are

$$
m_0(x)=\binom86(1/2)^8=\frac7{64},
\qquad
m_1(x)=\frac19.
$$

Thus

$$
\mathrm{BF}_{0:1}=\frac{7/64}{1/9}=\frac{63}{64}.
$$

Equal prior model odds yield

$$
P(H_0\mid x)=\frac{63/64}{1+63/64}
=\boxed{\frac{63}{127}\approx0.496063}.
$$

Under a purely continuous uniform prior, the exact point $\theta=1/2$ has posterior probability zero regardless of this comparison. The two analyses use different prior models; that difference is not a contradiction.

For problem 14 with $x=0$, the posterior kernel is $\theta^{-1}(1-\theta)^4$. Near zero, it behaves as $\theta^{-1}$, whose integral diverges, so there is no proper posterior from this formal weight.

With $x=2$, the kernel is $\theta^{1}(1-\theta)^2$, which normalizes to a $\mathrm{Beta}(2,3)$ posterior. The distinction comes from integrability at the endpoints, not from whether the algebra resembles a beta formula.

For problem 15, let $\theta(\eta)=1/(1+e^{-\eta})$. Multiplying by the Jacobian gives

$$
\pi_\eta(\eta\mid x)
=252\,\theta(\eta)^7[1-\theta(\eta)]^3,
\qquad \eta\in\mathbb R.
$$

The log derivative is

$$
\frac{d}{d\eta}\log\pi_\eta(\eta\mid x)
=7[1-\theta(\eta)]-3\theta(\eta)
=7-10\theta(\eta).
$$

It vanishes at $\theta(\eta)=0.7$, giving mode

$$
\eta_{\mathrm{mode}}=\log(7/3).
$$

The original theta-density mode is $0.75$, whose logit is $\log3$, a different value. The transformed mode is found after including the Jacobian; posterior set probabilities remain consistent.

### Variance reduction and parameter-dependent support

For problem 16, the prior mean is $1/10$ and prior variance is

$$
\mathrm{Var}(\Theta)=\frac{1\cdot9}{10^2\cdot11}
=\frac9{1100}\approx0.008182.
$$

A success gives posterior $\mathrm{Beta}(2,9)$, with variance

$$
v_1=\frac{2\cdot9}{11^2\cdot12}
=\frac3{242}\approx0.012397.
$$

This exceeds the prior variance. A failure gives posterior $\mathrm{Beta}(1,10)$, with variance $v_0=5/726$, and prior predictive success probability is $1/10$.

The average posterior variance is

$$
E_X[\mathrm{Var}(\Theta\mid X)]
=\frac1{10}\frac3{242}
+\frac9{10}\frac5{726}
=\frac9{1210}.
$$

The posterior mean is $(1+X)/11$, so

$$
\mathrm{Var}_X(E[\Theta\mid X])
=\frac{\mathrm{Var}(X)}{121}
=\frac{(1/10)(9/10)}{121}
=\frac9{12100}.
$$

Adding gives $9/1210+9/12100=9/1100$, verifying total variance. Conditioning reduced variance on average, not after every possible observation.

For problem 17, the full likelihood is

$$
f(x\mid\theta)=\frac1\theta
\mathbf1_{\{0<x<\theta\}}.
$$

At $x=1.5$, the data eliminate prior values below 1.5. Combining with the prior support gives $1.5<\theta<2$, and the unnormalized posterior is $1/(2\theta)$ there.

Its integral is

$$
m(1.5)=\frac12\int_{1.5}^2\frac1\theta\,d\theta
=\frac12\log(4/3).
$$

Therefore,

$$
\boxed{
\pi(\theta\mid1.5)=
\frac1{\theta\log(4/3)}
\mathbf1_{\{1.5<\theta<2\}}.
}
$$

The posterior mean is

$$
E[\Theta\mid1.5]
=\frac1{\log(4/3)}\int_{1.5}^2 1\,d\theta
=\frac{1/2}{\log(4/3)}
\approx1.738030.
$$

It lies between 1.5 and 2, as required. Omitting the parameter-dependent support indicator would retain impossible parameter values and would even make the integral near zero diverge.

## Executable R checks

### Run instructions and mathematical mapping

All code below uses base R and is intended to run in order in a clean session. The companion `Bayesian_Inference_checks.R` contains the same R blocks in the same order, so you can run `Rscript --vanilla Bayesian_Inference_checks.R` from its directory or use `source("Bayesian_Inference_checks.R")` in R.

R's distribution functions have consistent prefixes. The `d` functions evaluate a PMF or density, `p` functions evaluate CDF probabilities, `q` functions find quantiles, and `r` functions generate random draws.

| Mathematical task | R call |
|---|---|
| Beta posterior density | `dbeta(theta, shape1 = A, shape2 = B)` |
| $P(\Theta\le t\mid x)$, beta posterior | `pbeta(t, A, B)` |
| $P(\Theta>t\mid x)$, beta posterior | `pbeta(t, A, B, lower.tail = FALSE)` |
| Beta credible interval | `qbeta(c(0.025, 0.975), A, B)` |
| Gamma posterior tail | `pgamma(t, shape = a_star, rate = b_star, lower.tail = FALSE)` |
| Normal posterior interval | `qnorm(c(0.025, 0.975), mean = m_n, sd = sqrt(v_n))` |
| Negative-binomial predictive PMF | `dnbinom(k, size = a_star, prob = b_star / (b_star + u))` |
| Log beta function for stable ratios | `lbeta(A, B)` |

An upper-tail function avoids unnecessary subtraction of nearly equal numbers. Gamma calls explicitly specify `rate`, and normal calls explicitly convert variance to standard deviation.

### Setup and lecture examples

The helper checks numerical agreement between two calculations and counts successful assertions. Integrations provide numerical checks of formulas, not mathematical proofs.

```r
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
```

### Shrinkage, losses, and beta-binomial prediction

This block compares closed forms with direct integration and numerical minimization. It also prints the sensitivity table without treating any one prior as automatically appropriate.

```r
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
```

### Poisson-gamma and normal-normal checks

The predictive checks integrate over parameter uncertainty rather than substituting a posterior mean. Numerical integration is sufficient for these small examples but should not be assumed reliable for every parameter scale or tail.

```r
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
```

The normal evidence in this block is for the sample mean's sampling model, not for the full four-dimensional dataset. Those likelihoods are proportional in the parameter for fixed data, but their marginal densities refer to different data objects.

### Practice-problem checks and a stable grid approximation

The improper-prior exercise is resolved analytically by integrability, not by asking a numerical integrator to certify divergence. The finite integrals below verify the proper case and the parameter-dependent-support example.

```r
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
```

The grid approximation is intentionally compared with a known exact posterior. Subtracting the maximum log weight prevents numerical underflow without changing normalized weights; using `exp(log_weight)` directly can fail for much larger datasets.

### Independent posterior and predictive simulation

These simulations use independent direct draws and report Monte Carlo standard errors for means and indicator probabilities. The variance and covariance comparisons are descriptive, not deterministic equality assertions.

```r
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
```

There is no requirement for a simulation estimate to equal its target exactly, even with a fixed seed. The Monte Carlo standard error measures numerical sampling uncertainty in the estimate, not posterior uncertainty about a parameter; increasing simulation size reduces the former without changing the latter.

### Verification record

The code blocks above are the source for the companion R script. A clean-session run in R 4.5.2 passed all 113 deterministic assertions, and the separate seeded simulations completed. The companion README records the verification scope, including local Markdown navigation, code-fence, and mathematical-syntax checks.

These computational checks support reproducibility, but do not replace the derivations or prove that a statistical model is appropriate. No GitHub publication or verification inside GitHub's own renderer is implied.

## Common errors and lecture preparation

### Errors to catch before finalizing a solution

| Error | Repair |
|---|---|
| Treat $p(x\mid\theta)$ as $p(\theta\mid x)$ | Multiply by the prior and normalize over parameter values. |
| Average likelihoods without prior weights | Use the actual prior weights or state an equal-prior assumption. |
| Call a density height a probability | Integrate over a parameter set, or sum posterior masses for a discrete parameter. |
| Confuse a likelihood's variable | Hold data fixed and identify what varies in the calculation. |
| Drop a factor involving the unknown parameter | Omit only parameter-independent factors when finding posterior shape. |
| Drop parameter-dependent support | Retain indicators such as $x<\theta$. |
| Match beta shapes to exponents directly | Add 1 to each exponent. |
| Use the beta MAP formula at all shapes | Check interior-mode conditions and endpoints separately. |
| Call the mean universally optimal | State squared-error loss and the relevant moment assumption. |
| Estimate $h(\Theta)$ by $h(E[\Theta\mid x])$ | Under squared loss, compute $E[h(\Theta)\mid x]$. |
| Claim continuous MAP minimizes literal 0–1 loss | Exact-match 0–1 loss is degenerate for continuous posteriors. |
| Treat future trials as independent after mixing | They share an uncertain parameter and can be predictively dependent. |
| Use a plug-in predictive distribution as exact | Integrate the future likelihood over the posterior. |
| Interpret a credible interval as automatic coverage | Distinguish conditional posterior probability from repeated-sample coverage. |
| Claim every observation reduces variance | Total variance gives reduction on average, not pointwise. |
| Use gamma scale where rate is intended | With shape-rate $(a,b)$, mean is $a/b$; specify `rate=` in R. |
| Pass a variance as normal `sd` | Take its square root. |
| Reuse earlier data in a new likelihood | Update with new data only, unless explicitly refitting from the original prior. |
| Treat an improper prior as an ordinary distribution | Check posterior integrability and avoid undefined marginal likelihoods. |
| Interpret zero point probability as evidence against a null | Put positive prior mass on an exact point null if comparing it as a model. |

### A practical study sequence

Use three passes through this one-part guide. The goal is to make the workflow reusable rather than to memorize a different formula for every example.

- **Foundations pass:** Explain PMF versus PDF, likelihood versus posterior, and conditional independence. Rework examples A and B, then problems 1–5.
- **Lecture mastery pass:** Derive beta moments, shrinkage, squared-loss optimality, and the absolute-loss median. Rework examples C and D, then problems 6–8 and 16.
- **Synthesis pass:** Derive the predictive beta-binomial and explore the Poisson and normal connections. Work the extension problems, then use R to check calculations after doing them by hand.

### Questions worth bringing to lecture

- **Assessment scope:** Are credible intervals and prediction included, or primarily posterior distributions and loss-based point estimation?
- **Derivation expectations:** Should the beta-gamma identity be proved, or may it be used as a reference formula?
- **Notation:** Will the same lowercase $\theta$ denote both a random parameter and its possible values?
- **Loss conventions:** Are asymmetric losses or decisions included, or only squared and absolute error?
- **Prior assumptions:** Are all priors proper, and are hyperparameters treated as known?
- **Extensions:** Are Poisson-gamma, normal-normal, or point-null comparisons assigned, or enrichment only?

Before class, practice explaining why the posterior mean differs from the sample proportion and why it minimizes one specific loss. If you can derive those answers, the last third of the lecture becomes a consequence of conditional expectation rather than a new list of rules.

## Compact reference sheet

### General update and decision formulas

Keep parameter support and a finite positive normalizer with every update. Continuous densities require integration, while discrete parameter laws require summation.

| Task | Formula | Conditions or interpretation |
|---|---|---|
| Posterior | $\pi(\theta\mid x)=f(x\mid\theta)\pi(\theta)/m(x)$ | $0<m(x)<\infty$ |
| Evidence | $m(x)=\int f(x\mid u)\pi(u)\,du$ | Full prior predictive PMF or density |
| Posterior probability | $P(\Theta\in A\mid x)=\int_A\pi(\theta\mid x)\,d\theta$ | Use posterior masses for discrete parameters |
| Posterior odds | Prior odds × Bayes factor | Composite-model likelihoods are integrated over their priors |
| Squared-loss action | $E[\Theta\mid x]$ | Finite posterior second moment; action allowed |
| Squared-loss action for $h(\Theta)$ | $E[h(\Theta)\mid x]$ | Finite posterior second moment of the target |
| Absolute-loss action | Posterior median | Finite posterior first absolute moment |
| Asymmetric absolute-loss action | Quantile $c_U/(c_U+c_O)$ | $c_U$: underestimation cost; $c_O$: overestimation cost |
| Binary yes threshold | $q>c_{\mathrm{FP}}/(c_{\mathrm{FP}}+c_{\mathrm{FN}})$ | $q=P(H\mid x)$, zero loss for correct decisions |
| Equal-tailed credible interval | $[Q_x(\gamma/2),Q_x(1-\gamma/2)]$ | Continuous posterior, credibility $1-\gamma$ |
| Posterior predictive | $p(\widetilde x\mid x)=\int p(\widetilde x\mid\theta)\pi(\theta\mid x)\,d\theta$ | Future/past conditionally independent given parameter |

### Beta-binomial formulas

Assume $\Theta\sim\mathrm{Beta}(\alpha,\beta)$, $X\mid\theta\sim\mathrm{Binomial}(n,\theta)$, and observed count $x$. Define $A=\alpha+x$, $B=\beta+n-x$, $S=A+B$, with positive prior shapes.

$$
\Theta\mid x\sim\mathrm{Beta}(A,B),
\qquad
m(x)=\binom nx\frac{B(A,B)}{B(\alpha,\beta)}.
$$

$$
E[\Theta\mid x]=\frac AS,
\qquad
\mathrm{Var}(\Theta\mid x)=\frac{AB}{S^2(S+1)},
\qquad
E[\Theta^r\mid x]=\frac{B(A+r,B)}{B(A,B)}
\quad[A+r>0].
$$

$$
\theta_{\mathrm{MAP}}=\frac{A-1}{S-2}
\quad[A>1,\ B>1].
$$

$$
E[\Theta\mid x]
=\frac{\alpha+\beta}{\alpha+\beta+n}
\frac{\alpha}{\alpha+\beta}
+\frac{n}{\alpha+\beta+n}\frac xn
\quad[n>0].
$$

For $m$ new trials sharing the same probability,

$$
P(\widetilde X=k\mid x)
=\binom mk\frac{B(A+k,B+m-k)}{B(A,B)},
\qquad k=0,\ldots,m,
$$

$$
E[\widetilde X\mid x]=\frac{mA}{S},
\qquad
\mathrm{Var}(\widetilde X\mid x)
=\frac{mAB(S+m)}{S^2(S+1)}.
$$

### Poisson-gamma extension

For conditionally independent $X_i\mid\lambda\sim\mathrm{Poisson}(\lambda t_i)$, let $s=\sum_i x_i$, $T=\sum_i t_i$, and use shape-rate prior $\mathrm{Gamma}(a,b)$.

$$
\Lambda\mid\boldsymbol x\sim\mathrm{Gamma}(a_*,b_*),
\qquad a_*=a+s,\quad b_*=b+T,
$$

$$
E[\Lambda\mid\boldsymbol x]=\frac{a_*}{b_*},
\qquad
\mathrm{Var}(\Lambda\mid\boldsymbol x)=\frac{a_*}{b_*^2}.
$$

For future exposure $u$,

$$
\widetilde X\mid\boldsymbol x
\sim\mathrm{NegBin}\left(
\mathrm{size}=a_*,
\mathrm{prob}=\frac{b_*}{b_*+u}
\right),
$$

$$
E[\widetilde X\mid\boldsymbol x]=\frac{ua_*}{b_*},
\qquad
\mathrm{Var}(\widetilde X\mid\boldsymbol x)
=\frac{ua_*}{b_*}+\frac{u^2a_*}{b_*^2}.
$$

### Normal-normal extension

For conditionally iid $X_i\mid\theta\sim N(\theta,\sigma^2)$, known $\sigma^2$, and prior $\Theta\sim N(m_0,\tau_0^2)$,

$$
\Theta\mid\boldsymbol x\sim N(m_n,v_n),
$$

$$
v_n=\left(\frac1{\tau_0^2}+\frac n{\sigma^2}\right)^{-1},
\qquad
m_n=v_n\left(\frac{m_0}{\tau_0^2}
+\frac{n\overline x}{\sigma^2}\right),
$$

$$
\widetilde X\mid\boldsymbol x\sim N(m_n,\sigma^2+v_n).
$$

The final diagnostic is **What is fixed, what remains uncertain, and which distribution answers the question?** If those are correct, posterior calculation is ordinary conditional probability, and loss-based estimation is ordinary conditional expectation with a clearly stated objective.
