## Conditional Distribution I

**Core intuition:** a conditional distribution answers, “Once $X$ is known, how should the uncertainty about $Y$ change?” Think of starting with a joint distribution over $(X,Y)$, taking the slice at a fixed $X=x$ and renormalizing to get the distribution of $Y$ in that slice.

### What to know first

You will usually begin with a **joint distribution**, which describes two variables together.

- **Joint:** $p_{X,Y}(x,y)$ or $f_{X,Y}(x,y)$: how likely the pair $(X,Y)$ is.
- **Marginal:** $p_X(x)$ or $f_X(x)$: the distribution of $X$ alone, after summing or integrating out $Y$.
- **Conditional:** $p_{Y\mid X}(y\mid x)$ or $f_{Y\mid X}(y\mid x)$: the distribution of $Y$ when $X=x$ is treated as known.

For discrete variables, marginalize by summing:

$$
p_X(x)=\sum_y p_{X,Y}(x,y)
$$

For continuous variables, marginalize by integrating:

$$
f_X(x)=\int_{-\infty}^{\infty} f_{X,Y}(x,y)\,dy
$$

## Prereq. Definitions

### Standard (undergraduate-level) definitions

**Discrete case.** For discrete $X, Y$ with joint PMF $p_{X,Y}(x,y) = P(X=x, Y=y)$:

- Marginal: $p_X(x) = \sum_y p_{X,Y}(x,y)$.
- Conditional PMF of $Y$ given $X=x$, defined only when $p_X(x) > 0$:
$$
p_{Y|X}(y|x) = P(Y=y \mid X=x) = \frac{p_{X,Y}(x,y)}{p_X(x)}.
$$
- Validity: $p_{Y|X}(y|x) \ge 0$ and $\sum_y p_{Y|X}(y|x) = 1$ for each fixed $x$.

**Continuous case.** For jointly continuous $X, Y$ with joint PDF $f_{X,Y}$:

- Marginal: $f_X(x) = \int_{-\infty}^{\infty} f_{X,Y}(x,y)\,dy$.
- Conditional PDF of $Y$ given $X=x$, defined when $f_X(x) > 0$:
$$
f_{Y|X}(y|x) = \frac{f_{X,Y}(x,y)}{f_X(x)}.
$$
- Validity: $f_{Y|X}(y|x)\ge 0$ and $\int f_{Y|X}(y|x)\,dy = 1$.
- Probabilities: $P(Y\in B \mid X=x) = \int_B f_{Y|X}(y|x)\,dy$.

**Product rule and independence.**
$$
f_{X,Y}(x,y) = f_{Y|X}(y|x)\,f_X(x) = f_{X|Y}(x|y)\,f_Y(y).
$$
$X$ and $Y$ are independent iff $f_{Y|X}(y|x) = f_Y(y)$ for all $x$ with $f_X(x)>0$. The support must be a product set $A\times B$, so a triangular support always means dependence.

**Bayes' rule for densities.**
$$
f_{X|Y}(x|y) = \frac{f_{Y|X}(y|x)\,f_X(x)}{\int f_{Y|X}(y|x')\,f_X(x')\,dx'}.
$$
This is the engine of Bayesian inference later in the course. The denominator is just $f_Y(y)$.

**Why the continuous definition is not "just conditional probability."** Here $P(X=x)=0$, so the elementary formula $P(A|B)=P(A\cap B)/P(B)$ is $0/0$. The density ratio is a definition chosen so that the product rule holds. The measure-theoretic definitions below show it is the right one, and also where it can go wrong.

### Measure-theoretic definitions

Setting: a probability space $(\Omega,\mathcal F,P)$; random variables $X:\Omega\to\mathbb R^d$ and $Y:\Omega\to\mathbb R$ are measurable maps. The law of $X$ is the pushforward $P_X(A) = P(X\in A)$.

**Joint, marginal, and densities as Radon-Nikodym derivatives.** The joint law $P_{X,Y}$ is a measure on $\mathbb R^{d+1}$. If $P_{X,Y}\ll \mu\otimes\nu$ for $\sigma$-finite $\mu,\nu$, then the density is
$$
f_{X,Y} = \frac{dP_{X,Y}}{d(\mu\otimes\nu)}.
$$
- Taking $\mu,\nu$ as counting measure gives the PMF.
- Taking $\mu,\nu$ as Lebesgue measure gives the PDF.
- Taking a mix gives mixed discrete and continuous pairs under one framework.

Fubini's theorem gives the marginal density $f_X(x)=\int f_{X,Y}(x,y)\,\nu(dy)$.

**Conditional expectation given a sigma-algebra (Kolmogorov).** Let $Y\in L^1(P)$ and $\mathcal G\subseteq\mathcal F$ a sub-sigma-algebra. $E[Y|\mathcal G]$ is the (a.s. unique) random variable $Z$ with
1. $Z$ is $\mathcal G$-measurable, and
2. $\int_G Z\,dP = \int_G Y\,dP$ for every $G\in\mathcal G$.

Existence comes from Radon-Nikodym applied to the signed measure $G\mapsto\int_G Y\,dP$, which is absolutely continuous with respect to $P|_{\mathcal G}$. Uniqueness is only almost sure.

**Conditioning on a random variable.** $E[Y|X] := E[Y|\sigma(X)]$. By the Doob-Dynkin lemma there is a Borel function $h$ with $E[Y|X]=h(X)$. We then define
$$
E[Y|X=x] := h(x),
$$
determined only $P_X$-almost everywhere. Different versions of $h$ can differ on a $P_X$-null set of $x$ values. For continuous $X$, the value at a single point $x$ is therefore not canonically defined. The continuous formulas below give a good version.

**Conditional probability** is the special case $P(A|\mathcal G)=E[\mathbf 1_A|\mathcal G]$.

**Regular conditional distribution (probability kernel).** A function $\kappa:\mathbb R^d\times\mathcal B(\mathbb R)\to[0,1]$ such that
1. for each $x$, $\kappa(x,\cdot)$ is a probability measure;
2. for each Borel $B$, $x\mapsto\kappa(x,B)$ is measurable;
3. disintegration holds:
$$
P(X\in A,\,Y\in B)=\int_A \kappa(x,B)\,P_X(dx)\quad\text{for all Borel }A,B.
$$
We write $\kappa(x,B)=P(Y\in B|X=x)$. A regular version exists whenever $Y$ takes values in a Polish (Borel) space, so $\mathbb R$ and $\mathbb R^k$ are covered. For general $Y$, conditional probabilities exist for each $B$ but may fail to assemble into a measure.

**Link to densities.** If $P_{X,Y}\ll\mu\otimes\nu$ with density $f_{X,Y}$, a valid kernel is
$$
\kappa(x,B)=\int_B \frac{f_{X,Y}(x,y)}{f_X(x)}\,\nu(dy)\qquad\text{for }f_X(x)>0.
$$
Check: for any $g\ge 0$, $E[g(X)\mathbf 1_B(Y)]=\iint g(x)\mathbf 1_B(y)f_{X,Y}\,d\nu\,d\mu=\int g(x)\,\kappa(x,B)\,f_X(x)\,\mu(dx)$. This is exactly the disintegration identity, so the elementary definition is a theorem here.

**Borel-Kolmogorov caution.** Conditioning on a null event such as $\{X=Y\}$ depends on how the event is approximated. For example, conditioning on $Y-X=0$ and on $Y/X=1$ give different answers. The statement "condition on $X=x$" is meaningful only relative to the specific random variable $X$ and the sigma-algebra $\sigma(X)$.


## Part 1: likely foundations

Expect this first portion to focus on constructing and interpreting conditional pmfs/pdf’s.

### Discrete case

When $p_X(x)>0$, use the **conditional-pmf formula**:

$$
p_{Y\mid X}(y\mid x) = \frac{p_{X,Y}(x,y)}{p_X(x)}
$$

This is the usual conditional-probability rule applied to random variables. The denominator is the total probability in the selected $X=x$ slice, so dividing by it makes the conditional probabilities sum to 1 within that slice.

A reliable workflow:

1. Identify the support: which $(x,y)$ pairs have nonzero joint probability?
2. Find the marginal $p_X(x)$ by summing over $y$.
3. Divide the joint probability by that marginal.
4. Check that, for each fixed $x$, $\sum_y p_{Y\mid X}(y\mid x)=1$.

### Continuous case

For a joint pdf, the conditional density is

$$
f_{Y\mid X}(y\mid x) = \frac{f_{X,Y}(x,y)}{f_X(x)}
$$

provided $f_X(x)>0$.

The notation “given $X=x$” can feel strange because a continuous random variable has probability zero of equaling one exact number. Do not interpret the formula as ordinary division by $P(X=x)$; instead, it is a density defined on the slice $X=x$, which is constructed by conditioning in the limit.

### The reconstruction identity

A central identity is:

$$
p_{X,Y}(x,y) = p_{Y\mid X}(y\mid x)p_X(x)
$$

or, for densities,

$$
f_{X,Y}(x,y) = f_{Y\mid X}(y\mid x)f_X(x)
$$

This matters because many statistical models are built in exactly this order: model $X$, then model the distribution of $Y$ conditional on $X$.

## Part 2: likely extensions

The second section commonly builds from “find the conditional distribution” to “use it to compute quantities and reason about dependence.”

- **Conditional CDF:** $F_{Y\mid X}(y\mid x)=P(Y\le y\mid X=x)$.
- **Conditional expectation:** $E[Y\mid X=x]$, the conditional mean of $Y$ at a fixed value of $X$.
- **Conditional variance:** $\mathrm{Var}(Y\mid X=x)$, the remaining spread in $Y$ after $X$ is known.
- **Independence:** if $X$ and $Y$ are independent, conditioning does not change the distribution:

$$
f_{Y\mid X}(y\mid x) = f_Y(y)
$$

Equivalently, the joint density or pmf factorizes into marginals.

- **Bayes’ rule:** reverse the conditioning direction:

$$
p_{X\mid Y}(x\mid y) = \frac{p_{Y\mid X}(y\mid x)p_X(x)}{p_Y(y)}
$$

Bayes’ rule becomes especially important later in statistical inference: prior distribution $\times$ likelihood, normalized into a posterior distribution.

## A useful mental model

For your ML background, frame this as a **generative factorization**:

$$
p(x,y)=p(y\mid x)p(x)
$$

For supervised learning, $p(y\mid x)$ is often the object of greatest interest: it describes the target’s uncertainty after observing features. Regression estimates a conditional mean such as $E[Y\mid X=x]$.

Before class, try explaining in your own words why the denominator in $p_{Y\mid X}(y\mid x)$ must be $p_X(x)$, rather than $p_Y(y)$.
