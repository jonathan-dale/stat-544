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
