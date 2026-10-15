# STAT 544 Study Guide and Lecture Prep: Conditional Distributions I

## How to use this guide

Read Part 1 for the vocabulary, Part 2 for the derivation recipe, and Part 3 for conditional expectation. Then try each worked example in Part 4 yourself before reading the solution. Part 5 lists pitfalls and self-test questions. Part 6 gives R checks.

The one idea to hold onto: a conditional distribution is the joint distribution with one variable frozen, then rescaled so it integrates (or sums) to 1. Everything below is bookkeeping around that sentence.

---

## Part 1. Definitions

### 1.1 Standard (undergraduate-level) definitions

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

### 1.2 Measure-theoretic definitions

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

### 1.3 Translation table

| Standard | Measure-theoretic |
|---|---|
| $p_{Y\|X}(y\|x)$ or $f_{Y\|X}(y\|x)$ | density of $\kappa(x,\cdot)$ with respect to $\nu$ |
| $E[Y\|X=x]$ | $h(x)$ where $E[Y\|\sigma(X)]=h(X)$, unique $P_X$-a.e. |
| "$f_X(x)>0$" | $x$ outside a $P_X$-null set |
| "integrating out" | Fubini on $\mu\otimes\nu$ |

---

## Part 2. From joint to marginal to conditional: the recipe

Given a joint distribution, follow these steps every time.

1. **Draw the support.** Write the region where $f_{X,Y}>0$ as inequalities. Note whether bounds of one variable depend on the other.
2. **Check normalization** (sum or integral equals 1). A cheap sanity check on the problem statement.
3. **Marginalize the conditioning variable.** For $f_X(x)$, integrate (or sum) over $y$, using $y$-limits that depend on $x$ if the support is not a rectangle.
4. **State where $f_X(x)>0$.** The conditional is defined only there.
5. **Divide:** $f_{Y|X}(y|x)=f_{X,Y}(x,y)/f_X(x)$. Keep the $y$-support, which may depend on $x$.
6. **Verify** that $\int f_{Y|X}(y|x)\,dy=1$. If it is not 1, step 3 has an error.
7. **Compute what is asked:** a conditional probability (integrate over $B$), a conditional mean or variance (Part 3), or the other direction $f_{X|Y}$ by swapping roles.
8. **Cross-check** with the tower property $E[E[Y|X]]=E[Y]$, computed independently.

Shortcut: in step 5, the conditional density is proportional to the joint in $y$, with $x$ held fixed. Anything depending only on $x$ is part of the normalizing constant, so you can often identify a named distribution by inspecting the $y$-dependent kernel and skip the explicit marginal.

---

## Part 3. Conditional expectation $E[Y|X=x]$

### 3.1 Formulas

$$
E[Y|X=x]=\sum_y y\,p_{Y|X}(y|x)\qquad\text{or}\qquad \int y\,f_{Y|X}(y|x)\,dy.
$$
For a function $g$: $E[g(Y)|X=x]=\int g(y)f_{Y|X}(y|x)\,dy$.

$E[Y|X=x]$ is a number for each $x$. $E[Y|X]=h(X)$ is a random variable. Confusing these is the most common mistake.

### 3.2 Properties (all hold almost surely)

1. **Linearity:** $E[aY_1+bY_2|X]=aE[Y_1|X]+bE[Y_2|X]$.
2. **Taking out what is known:** $E[g(X)Y|X]=g(X)E[Y|X]$.
3. **Tower / law of total expectation:** $E[E[Y|X]]=E[Y]$. More generally, for $\mathcal H\subseteq\mathcal G$, $E[E[Y|\mathcal G]|\mathcal H]=E[Y|\mathcal H]$.
4. **Independence:** if $Y\perp X$, then $E[Y|X]=E[Y]$. The converse is false: zero dependence in the mean does not imply independence.
5. **Law of total variance:**
$$
\mathrm{Var}(Y)=E[\mathrm{Var}(Y|X)]+\mathrm{Var}(E[Y|X]).
$$
6. **Conditional variance:** $\mathrm{Var}(Y|X=x)=E[Y^2|X=x]-(E[Y|X=x])^2$.
7. **Jensen:** $\varphi(E[Y|X])\le E[\varphi(Y)|X]$ for convex $\varphi$.

### 3.3 Why it matters: the best predictor

If $E[Y^2]<\infty$, then $E[Y|X]$ is the orthogonal projection of $Y$ onto the closed subspace of square-integrable functions of $X$. It uniquely minimizes
$$
E\big[(Y-g(X))^2\big]\quad\text{over all measurable }g.
$$
This is the theoretical target of regression: the regression function is $E[Y|X=x]$. Linear regression is the special case where this function happens to be linear, as with the bivariate normal below.

### 3.4 Reference result: bivariate normal

If $(X,Y)\sim N(\mu_X,\mu_Y,\sigma_X^2,\sigma_Y^2,\rho)$, then
$$
Y|X=x\sim N\!\Big(\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X),\ \sigma_Y^2(1-\rho^2)\Big).
$$
The conditional mean is linear in $x$, and the conditional variance is constant. This is a special feature of the normal.

---

## Part 4. Four worked examples

### Example 1 (discrete): a joint PMF table

Let $X\in\{1,2,3\}$ and $Y\in\{0,1\}$ with joint PMF

| | $x=1$ | $x=2$ | $x=3$ |
|---|---|---|---|
| $y=0$ | 0.10 | 0.20 | 0.10 |
| $y=1$ | 0.20 | 0.15 | 0.25 |

Find the marginals, $p_{Y|X}$, $E[Y|X=x]$, verify the tower property, and find $E[X|Y=1]$.

**Step 1: normalization.** $0.10+0.20+0.10+0.20+0.15+0.25=1.00$. Valid.

**Step 2: marginals.** Sum columns for $X$:
$p_X(1)=0.10+0.20=0.30$, $p_X(2)=0.20+0.15=0.35$, $p_X(3)=0.10+0.25=0.35$.
Sum rows for $Y$: $p_Y(0)=0.40$, $p_Y(1)=0.60$.

**Step 3: conditional PMF of $Y$ given $X=x$.** Divide each column by its total:

| | $x=1$ | $x=2$ | $x=3$ |
|---|---|---|---|
| $P(Y=0\|X=x)$ | $0.10/0.30=1/3$ | $0.20/0.35=4/7$ | $0.10/0.35=2/7$ |
| $P(Y=1\|X=x)$ | $0.20/0.30=2/3$ | $0.15/0.35=3/7$ | $0.25/0.35=5/7$ |

Each column sums to 1.

**Step 4: conditional mean.** Since $Y$ is 0/1, $E[Y|X=x]=P(Y=1|X=x)$:
$h(1)=2/3,\ h(2)=3/7,\ h(3)=5/7$.

**Step 5: tower check.**
$$
E[h(X)]=0.30\cdot\tfrac23+0.35\cdot\tfrac37+0.35\cdot\tfrac57=0.20+0.15+0.25=0.60,
$$
which equals $E[Y]=p_Y(1)=0.60$. Verified.

**Step 6: the other direction.** $p_{X|Y}(x|1)=p(x,1)/0.60$: $1/3,\ 1/4,\ 5/12$ (sum = 1). Then
$$
E[X|Y=1]=1\cdot\tfrac13+2\cdot\tfrac14+3\cdot\tfrac5{12}=\tfrac{4+6+15}{12}=\tfrac{25}{12}\approx 2.083.
$$

**Takeaway.** Conditioning on a column means renormalizing that column. $X$ and $Y$ are dependent because the conditional probabilities differ across $x$, and they differ from $P(Y=1)=0.6$.

---

### Example 2 (continuous): $f(x,y)=x+y$ on the unit square

Let $f_{X,Y}(x,y)=x+y$ for $0<x<1,\ 0<y<1$. Find $f_X$, $f_{Y|X}$, $E[Y|X=x]$, and $P(Y<\tfrac12\mid X=\tfrac12)$.

**Step 1: normalization.** $\int_0^1\int_0^1 (x+y)\,dy\,dx=\int_0^1 (x+\tfrac12)\,dx=\tfrac12+\tfrac12=1$. Valid. The support is a rectangle, but the density does not factor, so $X,Y$ are dependent.

**Step 2: marginal of $X$.**
$$
f_X(x)=\int_0^1 (x+y)\,dy=x+\tfrac12,\qquad 0<x<1.
$$
It is positive on the whole interval, so the conditional is defined for every $x\in(0,1)$.

**Step 3: conditional density.**
$$
f_{Y|X}(y|x)=\frac{x+y}{x+\tfrac12},\qquad 0<y<1.
$$
Check: $\int_0^1\frac{x+y}{x+1/2}dy=\frac{x+1/2}{x+1/2}=1$.

**Step 4: conditional mean.**
$$
E[Y|X=x]=\frac{1}{x+\tfrac12}\int_0^1 y(x+y)\,dy=\frac{x/2+1/3}{x+1/2}=\frac{3x+2}{6x+3}.
$$
Sanity values: $x\to0$ gives $2/3$, and $x\to1$ gives $5/9$. The mean of $Y$ falls as $x$ rises, which seems odd, since the joint density $x+y$ increases in both. It is correct: a larger $x$ makes the $x$ term dominate and flattens the $y$-shape toward uniform, whose mean is $1/2$.

**Step 5: tower check.**
$$
E[Y]=\int_0^1 (x+\tfrac12)\cdot\frac{x/2+1/3}{x+1/2}\,dx=\int_0^1\left(\tfrac x2+\tfrac13\right)dx=\tfrac14+\tfrac13=\tfrac7{12}.
$$
Direct: $f_Y(y)=y+\tfrac12$ by symmetry, so $E[Y]=\int_0^1 y(y+\tfrac12)dy=\tfrac13+\tfrac14=\tfrac7{12}$. Matches.

**Step 6: a conditional probability at $x=\tfrac12$.** Here $f_X(\tfrac12)=1$, so $f_{Y|X}(y|\tfrac12)=\tfrac12+y$:
$$
P(Y<\tfrac12|X=\tfrac12)=\int_0^{1/2}(\tfrac12+y)\,dy=\tfrac14+\tfrac18=\tfrac38.
$$

---

### Example 3 (continuous, non-rectangular support): $f(x,y)=8xy$ on $0<x<y<1$

Find $f_X$, $f_Y$, both conditional densities, both conditional means, and verify with $E[X]$.

**Step 1: support and normalization.** The region is the triangle above the diagonal. The bounds interact, so $X,Y$ are dependent regardless of the density.
$$
\int_0^1\int_0^y 8xy\,dx\,dy=\int_0^1 8y\cdot\tfrac{y^2}{2}\,dy=\int_0^1 4y^3\,dy=1.
$$
Valid.

**Step 2: marginals.** For fixed $x$, $y$ runs from $x$ to 1:
$$
f_X(x)=\int_x^1 8xy\,dy=4x(1-x^2),\qquad 0<x<1.
$$
For fixed $y$, $x$ runs from 0 to $y$:
$$
f_Y(y)=\int_0^y 8xy\,dx=4y^3,\qquad 0<y<1.
$$

**Step 3: conditional densities.** Keep the support limits.
$$
f_{X|Y}(x|y)=\frac{8xy}{4y^3}=\frac{2x}{y^2},\quad 0<x<y,
$$
$$
f_{Y|X}(y|x)=\frac{8xy}{4x(1-x^2)}=\frac{2y}{1-x^2},\quad x<y<1.
$$
Checks: $\int_0^y 2x/y^2\,dx=1$. $\int_x^1 2y/(1-x^2)\,dy=(1-x^2)/(1-x^2)=1$.

**Step 4: conditional means.**
$$
E[X|Y=y]=\int_0^y x\cdot\frac{2x}{y^2}\,dx=\frac{2}{y^2}\cdot\frac{y^3}{3}=\frac{2y}{3}.
$$
$$
E[Y|X=x]=\int_x^1 y\cdot\frac{2y}{1-x^2}\,dy=\frac{2(1-x^3)}{3(1-x^2)}=\frac{2(1+x+x^2)}{3(1+x)}.
$$
The last step uses $1-x^3=(1-x)(1+x+x^2)$ and $1-x^2=(1-x)(1+x)$. Sanity: at $x\to0$ it is $2/3$; at $x\to1$ it is $2\cdot3/(3\cdot2)=1$, which is correct, since $y$ is squeezed to 1.

**Step 5: tower check.**
$$
E[X]=E[E[X|Y]]=\tfrac23E[Y]=\tfrac23\int_0^1 y\cdot4y^3\,dy=\tfrac23\cdot\tfrac45=\tfrac8{15}.
$$
Direct: $E[X]=\int_0^1 x\cdot4x(1-x^2)\,dx=4(\tfrac13-\tfrac15)=\tfrac8{15}$. Matches.

**Takeaway.** The conditional support depends on the conditioning value. Forgetting this and integrating over $(0,1)$ is the usual error.

---

### Example 4 (hierarchical, discrete): Poisson thinning

Let $N\sim\text{Poisson}(\lambda)$, and given $N=n$, $X|N=n\sim\text{Binomial}(n,p)$. (Think of $N$ arrivals, each kept independently with probability $p$.) Let $q=1-p$. Find the marginal of $X$, the conditional distribution of $N$ given $X=x$, and $E[N|X=x]$. Evaluate for $\lambda=10$, $p=0.3$, $x=4$.

**Step 1: joint PMF from the product rule.** For $0\le x\le n$:
$$
p(n,x)=P(N=n)\,P(X=x|N=n)=\frac{e^{-\lambda}\lambda^n}{n!}\binom nx p^xq^{n-x}=\frac{e^{-\lambda}\lambda^n\,p^xq^{n-x}}{x!\,(n-x)!}.
$$

**Step 2: marginal of $X$.** Sum over $n\ge x$ and let $m=n-x$:
$$
p_X(x)=\frac{e^{-\lambda}(\lambda p)^x}{x!}\sum_{m=0}^\infty\frac{(\lambda q)^m}{m!}=\frac{e^{-\lambda}(\lambda p)^x}{x!}e^{\lambda q}=\frac{e^{-\lambda p}(\lambda p)^x}{x!}.
$$
So $X\sim\text{Poisson}(\lambda p)$, using $-\lambda+\lambda q=-\lambda p$.

**Step 3: conditional of $N$ given $X=x$.** Divide:
$$
p(n|x)=\frac{e^{-\lambda}\lambda^np^xq^{n-x}/(x!(n-x)!)}{e^{-\lambda p}(\lambda p)^x/x!}=\frac{e^{-\lambda q}(\lambda q)^{n-x}}{(n-x)!},\qquad n=x,x+1,\dots
$$
(using $\lambda^n p^x = \lambda^{n-x}(\lambda p)^x$ and $e^{-\lambda+\lambda p}=e^{-\lambda q}$). This says $N-X\,|\,X=x\sim\text{Poisson}(\lambda q)$, regardless of $x$. The number of discarded arrivals is independent of the number kept.

**Step 4: conditional expectation.**
$$
E[N|X=x]=x+E[N-X|X=x]=x+\lambda(1-p).
$$
This is linear in $x$.

**Step 5: numbers.** With $\lambda=10,p=0.3,x=4$: $E[N|X=4]=4+10(0.7)=11$. Tower check: $E[E[N|X]]=E[X]+\lambda q=\lambda p+\lambda q=\lambda=10=E[N]$. Verified.

**Takeaway.** Recognizing the kernel in $n$ lets you skip normalization: after dividing, the only terms left in $n$ are $(\lambda q)^{n-x}/(n-x)!$, which is a shifted Poisson kernel.

---

## Part 5. Pitfalls and self-test

### Common mistakes

1. Treating $E[Y|X=x]$ (a number) and $E[Y|X]$ (a random variable) as the same object.
2. Dropping the dependence of the conditional support on the conditioning value (Example 3).
3. Dividing by the wrong marginal: $f_{Y|X}$ divides by $f_X$, not $f_Y$.
4. Treating $f_{Y|X}(y|x)$ as a probability. It is a density in $y$, and can exceed 1.
5. Claiming independence from $E[Y|X]=E[Y]$ alone. That shows only mean independence.
6. Conditioning on null events such as $\{X=Y\}$ without specifying the limiting random variable (Borel-Kolmogorov).
7. Using the factorization test with a non-rectangular support.

### Self-test (answers are short; derive them yourself first)

1. For $f(x,y)=2$ on $0<x<y<1$, find $f_{Y|X}$ and $E[Y|X=x]$. Answer: $Y|X=x\sim\text{Unif}(x,1)$, mean $(1+x)/2$.
2. Show from the definition in 1.2 that $E[E[Y|\mathcal G]]=E[Y]$. Hint: take $G=\Omega$.
3. For the bivariate normal, derive $E[Y|X=x]$ from $f_{X,Y}/f_X$ by completing the square.
4. Use total variance to compute $\mathrm{Var}(N)$ in Example 4 from $N\sim\text{Poisson}(\lambda)$ and verify against the conditional structure of $X$.
5. Explain in one paragraph why a regular conditional distribution exists for real-valued $Y$, and what could fail otherwise.

### Lecture prep questions to bring

- Which version of $E[Y|X=x]$ does the course adopt when $x$ is a single point of a continuous distribution?
- How does the Bayes' rule form here connect to posterior densities in the next unit?
- How is conditioning handled when $X$ is mixed discrete and continuous?

---

## Part 6. R reference checks

```r
# Example 1: discrete table
P <- matrix(c(.10,.20,.10,
              .20,.15,.25), nrow = 2, byrow = TRUE,
            dimnames = list(y = 0:1, x = 1:3))
px <- colSums(P)                 # marginal of X
cond <- sweep(P, 2, px, "/")     # P(Y=y | X=x), columns sum to 1
cond["1", ]                      # E[Y|X=x] = P(Y=1|X=x)
sum(px * cond["1", ])            # tower: should be 0.60

# Example 2: continuous, numerical conditional mean
f  <- function(x, y) x + y
fx <- function(x) integrate(function(y) f(x, y), 0, 1)$value
h  <- function(x) integrate(function(y) y * f(x, y), 0, 1)$value / fx(x)
c(h(0.5), (3*0.5 + 2) / (6*0.5 + 3))       # both 7/12
integrate(Vectorize(function(x) fx(x) * h(x)), 0, 1)$value   # 7/12

# Example 3: triangular support
hx <- function(x) integrate(function(y) y * 2*y/(1 - x^2), x, 1)$value
c(hx(0.4), 2*(1 + .4 + .4^2) / (3*(1 + .4)))   # should agree

# Example 4: simulation check of thinning
set.seed(544)
lambda <- 10; p <- 0.3
N <- rpois(2e6, lambda); X <- rbinom(2e6, N, p)
mean(N[X == 4])                  # approx 11 = 4 + lambda*(1-p)
mean(X); var(X)                  # both approx lambda*p = 3
```

---

## One-page summary

- Conditional = joint with one variable fixed, rescaled by the marginal: $f_{Y|X}=f_{X,Y}/f_X$.
- Measure-theoretically, it is a probability kernel $\kappa(x,\cdot)$ satisfying $P(X\in A,Y\in B)=\int_A\kappa(x,B)P_X(dx)$, unique $P_X$-a.e.
- $E[Y|X]$ is the $\sigma(X)$-measurable variable matching integrals of $Y$ over every set in $\sigma(X)$, and it is the $L^2$-best predictor of $Y$ from $X$.
- Recipe: support, normalize, marginalize, divide, verify the conditional integrates to 1, compute, cross-check with the tower property.
- Always carry the conditional support, and always run the tower check.
