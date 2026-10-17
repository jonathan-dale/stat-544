# STAT 544 Study Guide and Lecture Prep: Conditional Distributions I

A conditional distribution describes the uncertainty about one random variable after information about another becomes available. For a joint PMF or PDF, the working idea is **fix the conditioning variable, retain the correct slice of the support, and renormalize**.

This guide develops that idea through prerequisite review, definitions, a repeatable calculation procedure, and the original four worked examples. The formal section explains why conditioning still makes sense when an event such as $X=x$ has probability zero.

## How to use this guide

On a first pass, read the prerequisite review and standard definitions, then work through the examples. On a second pass, study conditional expectation as a random variable and connect the calculations to the formal definitions.

- **Core Part I skills:** Find marginals, determine conditional support, construct conditional PMFs/PDFs, and compute conditional probabilities and means.
- **Bridge to Part II:** Conditional CDFs, total expectation, conditional variance, Bayes reversal, and regression interpretation.
- **Advanced reference:** Sigma-algebras, Radon–Nikodym derivatives, probability kernels, and almost-sure uniqueness. Check with the instructor about the expected proof depth.

### Contents

- [Prerequisite review](#prerequisite-review)
- [Standard definitions](#standard-definitions)
- [Joint-to-marginal-to-conditional procedure](#joint-to-marginal-to-conditional-procedure)
- [Conditional expectation](#conditional-expectation)
- [Four worked examples](#four-worked-examples)
- [Formal measure-theoretic reference](#formal-measure-theoretic-reference)
- [Pitfalls and lecture preparation](#pitfalls-and-lecture-preparation)
- [R verification](#r-verification)
- [Compact reference sheet](#compact-reference-sheet)

## Prerequisite review

### Random variables, events, and observed values

A probability space $(\Omega,\mathcal F,P)$ consists of possible outcomes, measurable events, and their probabilities. A random variable is a measurable function on this space; the sigma-algebra framework becomes useful when specifying exactly what information conditioning reveals ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)).

Use uppercase letters for random variables and lowercase letters for their possible values. For example, $X$ is random before observation, while $x$ is held fixed in a calculation conditional on $X=x$.

- **Joint distribution:** Describes the pair $(X,Y)$.
- **Marginal distribution:** Describes one variable alone, ignoring the other.
- **Conditional distribution:** Describes one variable after fixing information about the other.
- **Support for calculation:** The values or region on which the mass or density is nonzero.
- **Almost surely (a.s.):** A property holds except on an event with probability zero.

### Mass is not density

A discrete PMF assigns probability to a point:

$$
p_X(x)=P(X=x),\qquad p_X(x)\ge0,\qquad \sum_xp_X(x)=1.
$$

A PDF assigns probability through integration:

$$
P(X\in A)=\int_Af_X(x)\,dx,\qquad f_X(x)\ge0,\qquad \int_{\mathbb R}f_X(x)\,dx=1.
$$

A density height need not lie between 0 and 1; its integral over an event must do so ([Penn State probability distributions](https://online.stat.psu.edu/stat500/Lesson03)). For an absolutely continuous variable, $P(X=x)=0$ at every fixed point, even where $f_X(x)>0$.

Throughout the density calculations here, “continuous” means **absolutely continuous with respect to Lebesgue measure**. A continuous CDF alone does not guarantee that a PDF exists.

For a pair, the analogous definitions are

$$
p_{X,Y}(x,y)=P(X=x,Y=y),\qquad \sum_x\sum_y p_{X,Y}(x,y)=1,
$$

or

$$
P((X,Y)\in A)=\iint_Af_{X,Y}(x,y)\,dx\,dy,\qquad
\iint_{\mathbb R^2}f_{X,Y}(x,y)\,dx\,dy=1.
$$

All PMFs and PDFs must be nonnegative. A joint PDF exists only under an assumption about the joint law, not merely because each variable separately has a PDF.

### Support determines the integration limits

Write the whole model, including where it is zero:

$$
f_{X,Y}(x,y)=g(x,y)\mathbf1_S(x,y).
$$

For fixed $x$, the allowed values of $y$ form the slice $S_x=\{y:(x,y)\in S\}$. For fixed $y$, the allowed values of $x$ form $S^y=\{x:(x,y)\in S\}$.

For example, if $S=\{(x,y):0<x<y<1\}$, then:

| Quantity wanted | Conditioning or retained value | Treatment of the other variable | Bounds |
|---|---|---|---|
| $f_X(x)$ | $x$ | $y$ | $x<y<1$ |
| $f_Y(y)$ | $y$ | $x$ | $0<x<y$ |
| Law of $Y\mid X=x$ | Fixed $x$ | $y$ remains variable | $x<y<1$ |
| Law of $X\mid Y=y$ | Fixed $y$ | $x$ remains variable | $0<x<y$ |

The distinction is between removing a variable by averaging and fixing a variable by conditioning. Both operations start from the same joint law, but they answer different questions.

### Expectation, variance, and elementary conditioning

Expectation is a probability-weighted average, calculated using the appropriate PMF or PDF ([Penn State probability distributions](https://online.stat.psu.edu/stat500/Lesson03)):

$$
E[g(Y)]=\sum_y g(y)p_Y(y)
\quad\text{or}\quad
E[g(Y)]=\int_{\mathbb R}g(y)f_Y(y)\,dy.
$$

Assume $E|g(Y)|<\infty$ when a finite signed expectation is needed. With $E[Y^2]<\infty$, variance is

$$
\operatorname{Var}(Y)=E[Y^2]-(E[Y])^2.
$$

Ordinary conditioning on an event with positive probability is

$$
P(A\mid B)=\frac{P(A\cap B)}{P(B)},\qquad P(B)>0.
$$

The conditional PMF follows by setting $A=\{Y=y\}$ and $B=\{X=x\}$ ([MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf)). The continuous case needs a different justification because the event $\{X=x\}$ has probability zero.

## Standard definitions

### Discrete conditional distributions

Start with the joint PMF and obtain the marginal of the conditioning variable:

$$
p_X(x)=\sum_{y\in S_x}p_{X,Y}(x,y).
$$

For $p_X(x)>0$, the conditional PMF is ([MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf))

$$
\boxed{p_{Y\mid X}(y\mid x)=P(Y=y\mid X=x)=\frac{p_{X,Y}(x,y)}{p_X(x)}.}
$$

It is nonnegative, zero outside the conditional support, and sums to 1:

$$
\sum_y p_{Y\mid X}(y\mid x)
=\frac{\sum_y p_{X,Y}(x,y)}{p_X(x)}
=1.
$$

The denominator is constant while summing over $y$, because $x$ is fixed. If $p_X(x)=0$, this ratio does not define a conditional law at $x$; the formal framework allows an arbitrary probability law there without changing the joint distribution.

### Continuous conditional distributions

Assume the pair $(X,Y)$ has a joint PDF. First marginalize, then divide, using a positive finite denominator ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)):

$$
f_X(x)=\int_{S_x}f_{X,Y}(x,y)\,dy,
$$

$$
\boxed{f_{Y\mid X}(y\mid x)=\frac{f_{X,Y}(x,y)}{f_X(x)},\qquad 0<f_X(x)<\infty.}
$$

For each such $x$, the conditional density integrates to 1:

$$
\int_{S_x}f_{Y\mid X}(y\mid x)\,dy
=\frac{\int_{S_x}f_{X,Y}(x,y)\,dy}{f_X(x)}=1.
$$

Consequently,

$$
P(Y\in B\mid X=x)=\int_{B\cap S_x}f_{Y\mid X}(y\mid x)\,dy.
$$

The formula specifies a convenient version of the conditional law; uniqueness is only $P_X$-almost everywhere. It is not obtained by dividing by $P(X=x)$.

A useful intuition comes from first conditioning on a narrow interval:

$$
P(Y\in B\mid x<X<x+\varepsilon)
=\frac{\int_x^{x+\varepsilon}\int_B f_{X,Y}(u,y)\,dy\,du}
{\int_x^{x+\varepsilon}f_X(u)\,du}.
$$

If the numerator's inner integral and $f_X$ are continuous at $x$, with $f_X(x)>0$, dividing by $\varepsilon$ and taking a limit gives the density-ratio formula. The formal kernel construction below supplies the general justification without requiring this pointwise continuity argument.

### Conditional CDFs and reconstruction

For a chosen conditional law, define

$$
F_{Y\mid X}(t\mid x)=P(Y\le t\mid X=x).
$$

Calculate this by summing conditional masses with $y\le t$ or integrating the conditional PDF up to $t$. Include the values 0 below the support and 1 above it, just as for an ordinary CDF.

The product rule reconstructs the joint from the conditional and marginal laws ([MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf)):

$$
p_{X,Y}(x,y)=p_{Y\mid X}(y\mid x)p_X(x),
$$

$$
f_{X,Y}(x,y)=f_{Y\mid X}(y\mid x)f_X(x)
=f_{X\mid Y}(x\mid y)f_Y(y).
$$

For densities, identities are understood almost everywhere. This is also the basis of a generative model: draw $X$ from its marginal, then draw $Y$ from its conditional law.

### Independence and Bayes reversal

Independence means $P(X\in A,Y\in B)=P_X(A)P_Y(B)$ for every pair of measurable sets. For PMFs it is equivalent to pointwise factorization; for PDFs the equivalent factorization holds almost everywhere ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)):

$$
p_{X,Y}(x,y)=p_X(x)p_Y(y),
\qquad
f_{X,Y}(x,y)=f_X(x)f_Y(y).
$$

Equivalently, the conditional law of $Y$ equals its marginal law for $P_X$-almost every conditioning value. Conditioning then provides no information about $Y$.

Check the complete density, including its support indicator. A rectangular support does not prove independence; positive density throughout a nondegenerate triangle is incompatible with independence, but merely describing a larger containing region as triangular is not enough.

Bayes' rule reverses the direction of conditioning:

$$
f_{X\mid Y}(x\mid y)
=\frac{f_{Y\mid X}(y\mid x)f_X(x)}
{\int f_{Y\mid X}(y\mid u)f_X(u)\,du},
$$

where the denominator is $f_Y(y)$ and must be positive and finite. The analogous discrete formula uses PMFs and a sum; both are instances of reconstructing the joint and dividing by the new conditioning marginal ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)).

## Joint-to-marginal-to-conditional procedure

Use the same sequence each time. State the target before doing algebra so that the direction of conditioning stays clear.

1. **Identify the model:** Is it a joint PMF, joint PDF, mixed law, or something without an ordinary joint density?
2. **Write the support:** For $Y\mid X=x$, determine the allowed values $S_x$ of $y$.
3. **Check normalization:** Verify nonnegativity and total mass 1; solve for any unknown normalizing constant.
4. **Obtain the conditioning marginal:** To find $p_X$ or $f_X$, sum or integrate out $Y$ while retaining $x$.
5. **State admissible conditioning values:** Require $p_X(x)>0$ or $0<f_X(x)<\infty$ for the displayed ratio.
6. **Divide and keep the support:** The conditional law is joint divided by the conditioning marginal.
7. **Check conditional normalization:** A failed check may indicate incorrect bounds, a wrong marginal, or an algebraic error.
8. **Compute and cross-check:** Find probabilities or moments, then compare with the joint law, marginal law, or total expectation.

The mnemonic is **retain, remove, divide**: retain the conditioning value, remove the other variable to obtain its marginal normalizer, then divide the joint by that normalizer. During the final conditional calculation, the unconditioned variable still ranges over its allowed values.

A shortcut is to recognize the conditional density as a function of $y$ proportional to the joint density, with $x$ fixed. This can avoid an explicit marginal calculation, but only if the support and the normalized distribution family have been identified correctly.

## Conditional expectation

### Fixed-value mean versus a random conditional mean

For $E|Y|<\infty$, calculate a conditional mean using ordinary expectation under the conditional law ([MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf)):

$$
h(x)=E[Y\mid X=x]
=\sum_y y\,p_{Y\mid X}(y\mid x)
\quad\text{or}\quad
\int y\,f_{Y\mid X}(y\mid x)\,dy.
$$

For $E[g(Y)\mid X=x]$, replace the multiplier $y$ with $g(y)$, assuming the required integrability. Do not change the conditional distribution merely because the quantity being averaged has changed.

| Expression | Object | Meaning |
|---|---|---|
| $E[Y]$ | Number | Unconditional mean |
| $h(x)=E[Y\mid X=x]$ | Number at fixed $x$ | Mean after observing that value |
| $h$ | Function | Maps an observed value to a mean |
| $h(X)=E[Y\mid X]$ | Random variable | Conditional mean before $X$ is observed |

If $Y$ takes only values 0 and 1, then $E[Y\mid X=x]=P(Y=1\mid X=x)$. This simplification is specific to indicator-valued variables, not arbitrary numerical outcomes.

### Total expectation and its derivation

The law of total expectation is ([MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf))

$$
\boxed{E[E[Y\mid X]]=E[Y].}
$$

For the discrete case, substitution and cancellation show why:

$$
\begin{aligned}
E[h(X)]
&=\sum_{x:p_X(x)>0}\left[\sum_y y\frac{p_{X,Y}(x,y)}{p_X(x)}\right]p_X(x)\\
&=\sum_x\sum_y y\,p_{X,Y}(x,y)
=\sum_y y\,p_Y(y)=E[Y].
\end{aligned}
$$

The continuous proof replaces sums with integrals. The outer average must weight $h(x)$ by the marginal distribution of $X$; it is not generally an equally weighted average of conditional means.

### Properties and a Part II preview

The following expectation identities hold almost surely under integrability assumptions; the variance identities require $E[Y^2]<\infty$ ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf); [MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf)).

- **Linearity:** $E[aY_1+bY_2\mid X]=aE[Y_1\mid X]+bE[Y_2\mid X]$ for constants $a,b$.
- **Taking out what is known:** $E[g(X)Y\mid X]=g(X)E[Y\mid X]$, for example when $g$ is bounded and $Y$ integrable.
- **Nested conditioning:** If $\mathcal H\subseteq\mathcal G$, then $E[E[Y\mid\mathcal G]\mid\mathcal H]=E[Y\mid\mathcal H]$.
- **Independence:** If $Y$ is independent of $X$, then $E[Y\mid X]=E[Y]$. A constant conditional mean alone does not imply independence.
- **Conditional Jensen:** $\varphi(E[Y\mid X])\le E[\varphi(Y)\mid X]$ for convex $\varphi$, when the relevant variables are integrable.

Conditional variance and total variance are

$$
v(x)=\operatorname{Var}(Y\mid X=x)
=E[Y^2\mid X=x]-h(x)^2,
$$

$$
\boxed{\operatorname{Var}(Y)=E[v(X)]+\operatorname{Var}(h(X)).}
$$

To see the decomposition, take expectations of $E[Y^2\mid X]=v(X)+h(X)^2$, then subtract $E[Y]^2$. This separates average remaining variation from variation among the conditional means.

Total variance implies $E[v(X)]\le\operatorname{Var}(Y)$, not necessarily $v(x)\le\operatorname{Var}(Y)$ at every individual $x$. Likewise, the tower property requires nested information, not two arbitrary conditioning variables.

### Why this is the target of regression

For $Y\in L^2$, $h(X)=E[Y\mid X]$ minimizes mean squared error among square-integrable predictors based on $X$, uniquely up to almost-sure equality ([MIT conditional expectation lecture](https://ocw.mit.edu/courses/18-440-probability-and-random-variables-spring-2014/45b2b86ffc4234ef33d78b9cbda159e9_MIT18_440S14_Lecture26.pdf)). The key identity is

$$
E[(Y-g(X))^2]
=E[(Y-h(X))^2]+E[(h(X)-g(X))^2].
$$

The cross term vanishes because $E[Y-h(X)\mid X]=0$. The conditional mean need not be linear; linear regression may approximate it rather than equal it.

For a jointly normal pair with $\sigma_X,\sigma_Y>0$ and $|\rho|<1$, the familiar reference result is ([Penn State conditional normal distributions](https://online.stat.psu.edu/stat505/Lesson06))

$$
Y\mid X=x\sim
N\left(\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X),\
\sigma_Y^2(1-\rho^2)\right).
$$

Here $N(m,v)$ means mean $m$ and variance $v$, and $\rho$ is the correlation. Joint normality produces a linear conditional mean and constant conditional variance, but those features are not exclusive to normal models; normal marginals alone also do not establish joint normality.

## Four worked examples

### Discrete joint PMF table

Example 1. Let $X\in\{1,2,3\}$ and $Y\in\{0,1\}$ have the joint PMF below. Find the marginals, conditional law of $Y$, conditional mean, and $E[X\mid Y=1]$.

|  | $x=1$ | $x=2$ | $x=3$ |
|---|---:|---:|---:|
| $y=0$ | 0.10 | 0.20 | 0.10 |
| $y=1$ | 0.20 | 0.15 | 0.25 |

1. **Normalize.** All entries are nonnegative, and their sum is $0.10+0.20+0.10+0.20+0.15+0.25=1$.
2. **Find marginals.** Column totals give $p_X(1)=0.30$, $p_X(2)=0.35$, and $p_X(3)=0.35$. Row totals give $p_Y(0)=0.40$ and $p_Y(1)=0.60$.
3. **Normalize each column.** Divide the entries in an $X=x$ column by that column's total.

|  | $x=1$ | $x=2$ | $x=3$ |
|---|---:|---:|---:|
| $P(Y=0\mid X=x)$ | $0.10/0.30=1/3$ | $0.20/0.35=4/7$ | $0.10/0.35=2/7$ |
| $P(Y=1\mid X=x)$ | $0.20/0.30=2/3$ | $0.15/0.35=3/7$ | $0.25/0.35=5/7$ |

Each column sums to 1. The $y$ values remain 0 and 1; only their probabilities have changed.

4. **Calculate the mean.** Since $Y$ is 0 or 1,

$$
h(x)=0\cdot P(Y=0\mid X=x)+1\cdot P(Y=1\mid X=x),
$$

so $h(1)=2/3$, $h(2)=3/7$, and $h(3)=5/7$. The random variable $h(X)$ takes these values with probabilities $0.30,0.35,0.35$, respectively.

5. **Check total expectation.**

$$
E[h(X)]=0.30\left(\frac23\right)+0.35\left(\frac37\right)
+0.35\left(\frac57\right)=0.20+0.15+0.25=0.60.
$$

Directly, $E[Y]=0(0.40)+1(0.60)=0.60$. Both calculations agree.

6. **Reverse the conditioning.** For $Y=1$, normalize the second row:

$$
\bigl(p_{X\mid Y}(1\mid1),p_{X\mid Y}(2\mid1),p_{X\mid Y}(3\mid1)\bigr)
=\left(\frac13,\frac14,\frac5{12}\right).
$$

$$
E[X\mid Y=1]
=1\left(\frac13\right)+2\left(\frac14\right)+3\left(\frac5{12}\right)
=\frac{25}{12}\approx2.0833.
$$

The variables are dependent because the conditional probability of $Y=1$ changes with $x$. As a variance preview, $Y^2=Y$, so $\operatorname{Var}(Y\mid X=x)=h(x)(1-h(x))$, giving $2/9,12/49,10/49$.

### Continuous density on the unit square

Example 2. Let $f_{X,Y}(x,y)=x+y$ for $0<x<1$ and $0<y<1$, and zero otherwise. Find the marginals, conditional law of $Y$, conditional mean, and $P(Y<1/2\mid X=1/2)$.

1. **Normalize over the whole square.**

$$
\int_0^1\int_0^1(x+y)\,dy\,dx
=\int_0^1\left(x+\frac12\right)\,dx
=\frac12+\frac12=1.
$$

2. **Find marginals.** Integrate out $y$ to retain $x$:

$$
f_X(x)=\int_0^1(x+y)\,dy=x+\frac12,\qquad 0<x<1.
$$

By symmetry, $f_Y(y)=y+1/2$ for $0<y<1$. Both marginals are zero outside their stated intervals.

3. **Divide and verify.**

$$
f_{Y\mid X}(y\mid x)=\frac{x+y}{x+1/2},\qquad 0<y<1,\quad0<x<1,
$$

$$
\int_0^1\frac{x+y}{x+1/2}\,dy
=\frac{x+1/2}{x+1/2}=1.
$$

4. **Compute the conditional mean.**

$$
\begin{aligned}
h(x)
&=\int_0^1y\,\frac{x+y}{x+1/2}\,dy\\
&=\frac1{x+1/2}\left[x\frac{y^2}{2}+\frac{y^3}{3}\right]_0^1\\
&=\frac{x/2+1/3}{x+1/2}
=\frac{3x+2}{6x+3}.
\end{aligned}
$$

The mean decreases from the limit $2/3$ as $x\downarrow0$ to the limit $5/9$ as $x\uparrow1$. This is consistent with the conditional density: increasing $x$ increases its uniform component relative to the component proportional to $y$.

5. **Compute the conditional CDF and requested probability.**

$$
F_{Y\mid X}(t\mid x)=
\begin{cases}
0,&t\le0,\\
\dfrac{xt+t^2/2}{x+1/2},&0<t<1,\\
1,&t\ge1.
\end{cases}
$$

At $x=1/2$, the density is $1/2+y$. Therefore

$$
P(Y<1/2\mid X=1/2)
=\int_0^{1/2}\left(\frac12+y\right)\,dy
=\frac14+\frac18=\frac38.
$$

6. **Check total expectation.**

$$
\begin{aligned}
E[h(X)]
&=\int_0^1\frac{x/2+1/3}{x+1/2}\left(x+\frac12\right)\,dx\\
&=\int_0^1\left(\frac x2+\frac13\right)\,dx
=\frac14+\frac13=\frac7{12}.
\end{aligned}
$$

Directly, $E[Y]=\int_0^1y(y+1/2)\,dy=1/3+1/4=7/12$. Notice that the outer expectation uses the marginal density $f_X$, not equal weighting over $x$.

The support is rectangular, but $f_{Y\mid X}$ depends on $x$, so the variables are dependent. For instance, at $x=y=1/2$, the natural continuous versions give $f_{Y\mid X}(1/2\mid1/2)=1$; equality at one point is not an independence test.

### Continuous density on a triangle

Example 3. Let $f_{X,Y}(x,y)=8xy$ on $0<x<y<1$, and zero otherwise. Find both marginals, both conditional densities and means, and check the marginal means in both directions.

1. **Interpret the support and normalize.** For fixed $y$, $x$ ranges from 0 to $y$:

$$
\int_0^1\int_0^y8xy\,dx\,dy
=\int_0^1 4y^3\,dy=1.
$$

The factorization $8xy=(8x)y$ does not establish independence because it omits the triangular support indicator. The full density is zero below the diagonal.

2. **Calculate both marginals using different slices.**

$$
f_X(x)=\int_x^1 8xy\,dy
=4x(1-x^2),\qquad0<x<1,
$$

$$
f_Y(y)=\int_0^y8xy\,dx
=4y^3,\qquad0<y<1.
$$

Both marginals are positive throughout $(0,1)$, but the joint density is zero for $0<y<x<1$. Their product therefore cannot equal the joint density almost everywhere, establishing dependence.

3. **Calculate and normalize the conditional densities.**

$$
f_{X\mid Y}(x\mid y)=\frac{8xy}{4y^3}
=\frac{2x}{y^2},\qquad0<x<y,\quad0<y<1,
$$

$$
f_{Y\mid X}(y\mid x)=\frac{8xy}{4x(1-x^2)}
=\frac{2y}{1-x^2},\qquad x<y<1,\quad0<x<1.
$$

Outside each conditional support the corresponding density is zero. The checks are

$$
\int_0^y\frac{2x}{y^2}\,dx=1,\qquad
\int_x^1\frac{2y}{1-x^2}\,dy=\frac{1-x^2}{1-x^2}=1.
$$

4. **Calculate both conditional means.**

$$
E[X\mid Y=y]
=\int_0^yx\frac{2x}{y^2}\,dx
=\frac2{y^2}\frac{y^3}{3}
=\frac{2y}{3}.
$$

$$
\begin{aligned}
E[Y\mid X=x]
&=\int_x^1y\frac{2y}{1-x^2}\,dy\\
&=\frac{2(1-x^3)}{3(1-x^2)}
=\frac{2(1+x+x^2)}{3(1+x)}.
\end{aligned}
$$

The last step cancels $1-x$ using the difference-of-powers factorizations. For $0<x<1$, this mean lies strictly between $x$ and 1, as required by the conditional support.

5. **Practice intersecting an event with the conditional support.** At $x=1/2$, the conditional support is $(1/2,1)$, so

$$
P(Y>3/4\mid X=1/2)
=\int_{3/4}^1\frac{2y}{1-1/4}\,dy
=\frac{1-9/16}{3/4}=\frac7{12}.
$$

In contrast, $P(Y<1/4\mid X=1/2)=0$ because the event is outside the conditional support. Always intersect the requested event with the allowed slice before integrating.

6. **Check both total-expectation identities.**

$$
E[Y]=\int_0^1y(4y^3)\,dy=\frac45,
$$

$$
E[E[X\mid Y]]
=\frac23E[Y]=\frac8{15}
=\int_0^1x\,4x(1-x^2)\,dx.
$$

For the reverse direction,

$$
\begin{aligned}
E[E[Y\mid X]]
&=\int_0^1\frac{2(1-x^3)}{3(1-x^2)}4x(1-x^2)\,dx\\
&=\frac83\int_0^1(x-x^4)\,dx
=\frac83\left(\frac12-\frac15\right)=\frac45.
\end{aligned}
$$

The limits of the natural conditional mean are $2/3$ as $x\downarrow0$ and 1 as $x\uparrow1$. These limits do not mean the original ratio is defined at the endpoints, where the displayed marginal density is zero.

### Hierarchical discrete model: Poisson thinning

Example 4. Let $N\sim\operatorname{Poisson}(\lambda)$ with $\lambda>0$, and let $X\mid N=n\sim\operatorname{Binomial}(n,p)$ with $0<p<1$. Interpret $N$ as arrivals and $X$ as the number retained independently with probability $p$; write $q=1-p$.

Find the marginal of $X$, the law of $N\mid X=x$, and its conditional moments. The numerical case is $\lambda=10$, $p=0.3$, and $x=4$.

1. **Build the joint PMF.** For nonnegative integers $0\le x\le n$,

$$
\begin{aligned}
p_{N,X}(n,x)
&=P(N=n)P(X=x\mid N=n)\\
&=\frac{e^{-\lambda}\lambda^n}{n!}\binom nx p^xq^{n-x}\\
&=\frac{e^{-\lambda}\lambda^np^xq^{n-x}}{x!(n-x)!}.
\end{aligned}
$$

For each $n$, summing over $x=0,\ldots,n$ gives $P(N=n)$. Summing over $n$ then gives 1, verifying joint normalization.

2. **Sum over the possible totals.** Fixing $X=x$ forces $n\ge x$; use $m=n-x$:

$$
\begin{aligned}
p_X(x)
&=\sum_{n=x}^{\infty}p_{N,X}(n,x)\\
&=\frac{e^{-\lambda}(\lambda p)^x}{x!}
\sum_{m=0}^{\infty}\frac{(\lambda q)^m}{m!}\\
&=\frac{e^{-\lambda}(\lambda p)^x}{x!}e^{\lambda q}
=\frac{e^{-\lambda p}(\lambda p)^x}{x!}.
\end{aligned}
$$

The series is the exponential series $\sum_{m=0}^{\infty}a^m/m!=e^a$. Therefore $X\sim\operatorname{Poisson}(\lambda p)$.

3. **Reverse the conditioning.**

$$
\begin{aligned}
p_{N\mid X}(n\mid x)
&=\frac{p_{N,X}(n,x)}{p_X(x)}\\
&=\frac{e^{-\lambda q}(\lambda q)^{n-x}}{(n-x)!},
\qquad n=x,x+1,\ldots.
\end{aligned}
$$

The factors $(\lambda p)^x/x!$ cancel, and $e^{-\lambda}/e^{-\lambda p}=e^{-\lambda q}$. Substituting $m=n-x$ shows that these probabilities sum to 1.

4. **Interpret the shifted Poisson law.** The result says

$$
N-X\mid X=x\sim\operatorname{Poisson}(\lambda q).
$$

It is the number discarded that has a conditional Poisson law, not the total $N$ without its shift. Indeed, if $D=N-X$, the transformed joint PMF factors as

$$
p_{X,D}(x,d)
=p_{N,X}(x+d,x)
=\frac{e^{-\lambda p}(\lambda p)^x}{x!}
\frac{e^{-\lambda q}(\lambda q)^d}{d!},
$$

for $x,d\ge0$. This directly proves that retained and discarded counts are independent.

5. **Calculate conditional moments.**

$$
E[N\mid X=x]=x+\lambda q,\qquad
\operatorname{Var}(N\mid X=x)=\lambda q.
$$

With $\lambda=10$, $p=0.3$, and $x=4$, the conditional law is a shift of a Poisson distribution with mean 7. Its mean is 11 and its variance is 7.

6. **Check total expectation and total variance.**

$$
E[E[N\mid X]]=E[X]+\lambda q=\lambda p+\lambda q=\lambda.
$$

$$
\begin{aligned}
E[\operatorname{Var}(N\mid X)]
+\operatorname{Var}(E[N\mid X])
&=\lambda q+\operatorname{Var}(X+\lambda q)\\
&=\lambda q+\lambda p=\lambda.
\end{aligned}
$$

Both match the original Poisson distribution of $N$. Recognizing a distribution family helps with moments, but the support and normalizing constant must still be correct.

## Formal measure-theoretic reference

This section places the preceding formulas in a common framework. The essential distinction is between the existence of a conditional law and the existence of an ordinary density for that law.

### Probability spaces, laws, and reference measures

Let $X:\Omega\to\mathbb R^d$ and $Y:\Omega\to\mathbb R$ be measurable maps on $(\Omega,\mathcal F,P)$. Their laws are pushforward measures, for example $P_X(A)=P(X\in A)$ ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)).

Here $\mathcal F$ is a sigma-algebra: it contains $\Omega$ and is closed under complements and countable unions. The Borel sigma-algebra $\mathcal B(\mathbb R)$ is generated by open subsets of $\mathbb R$.

If $P_{X,Y}$ is absolutely continuous with respect to a product of sigma-finite measures, write

$$
P_{X,Y}\ll\mu\otimes\nu,\qquad
r(x,y)=\frac{dP_{X,Y}}{d(\mu\otimes\nu)}(x,y).
$$

The symbol $\ll$ means every reference-measure-null set is also null for the joint law. The usual density derivations extend to this setting ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)).

- **PMF:** Use counting measure on finite or countable value sets.
- **PDF:** Use Lebesgue measure on the real coordinates.
- **Mixed mass-density:** Use counting measure for a discrete coordinate and Lebesgue measure for a continuous coordinate.

Since $r\ge0$, Tonelli's theorem gives the marginal density with respect to $\mu$:

$$
r_X(x)=\int r(x,y)\,\nu(dy).
$$

Fubini's theorem applies to signed integrands when absolute integrability holds. Distinguishing these conditions is useful when moving from densities to expectations.

### Conditional expectation given a sigma-algebra

Let $\mathcal G\subseteq\mathcal F$ represent the available information, and assume $Y\in L^1(P)$, meaning $E|Y|<\infty$. A version of $E[Y\mid\mathcal G]$ is an integrable random variable $Z$ satisfying two conditions ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)):

1. **Measurability:** $Z$ is $\mathcal G$-measurable.
2. **Integral matching:** For every $G\in\mathcal G$,

$$
\int_G Z\,dP=\int_G Y\,dP.
$$

Measurability says $Z$ uses only the available information. Integral matching says it reproduces the expected contribution of $Y$ on every event distinguishable using that information.

Existence follows from Radon–Nikodym applied to the finite signed measure $G\mapsto\int_GY\,dP$, which is absolutely continuous with respect to $P|_{\mathcal G}$; uniqueness is only almost sure ([MIT measure-theoretic probability lecture](https://ocw.mit.edu/courses/18-175-theory-of-probability-spring-2014/6dc0eeda7f7a0cbb95f9b7a16cd1f276_MIT18_175S14_Lecture26.pdf)). Taking $G=\Omega$ immediately gives total expectation.

Conditional probability is the special case

$$
P(A\mid\mathcal G)=E[\mathbf1_A\mid\mathcal G].
$$

### Conditioning on a random variable

The information revealed by $X$ is $\sigma(X)$, the sigma-algebra generated by $X$. Define

$$
E[Y\mid X]:=E[Y\mid\sigma(X)].
$$

The Doob–Dynkin lemma gives a Borel function $h$ such that $E[Y\mid X]=h(X)$ a.s.; a selected version lets us write $E[Y\mid X=x]=h(x)$ ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)). Only the $P_X$-almost-everywhere values of $h$ are determined.

For an absolutely continuous $X$, changing $h$ at a single point does not change $h(X)$ a.s. A smooth formula may provide a natural version, but the abstract definition does not uniquely determine each isolated point value.

### Regular conditional distributions

A probability kernel $\kappa:\mathbb R^d\times\mathcal B(\mathbb R)\to[0,1]$ represents the conditional law of $Y$ given $X$ when it satisfies the following conditions ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)):

1. **Probability measure:** For each $x$, $\kappa(x,\cdot)$ is a probability measure.
2. **Measurability:** For each Borel $B$, $x\mapsto\kappa(x,B)$ is measurable.
3. **Disintegration:** For all Borel $A,B$,

$$
P(X\in A,Y\in B)=\int_A\kappa(x,B)\,P_X(dx).
$$

We write $\kappa(x,B)=P(Y\in B\mid X=x)$. Averaging the conditional laws over $X$ reconstructs the joint law.

For the real-valued and Euclidean-valued variables considered here, regular conditional distributions exist; more generally, standard Borel state spaces support this framework ([CMU regular conditional distributions lecture](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Scribed_Lectures/Mar6.pdf)). Versions agree as probability measures for $P_X$-almost every $x$; arbitrary measurable spaces require additional care.

When $0<r_X(x)<\infty$, the density construction is

$$
\kappa(x,B)=\int_B\frac{r(x,y)}{r_X(x)}\,\nu(dy).
$$

On the exceptional $P_X$-null set, assign a fixed probability law to complete the kernel. The construction satisfies disintegration because

$$
\begin{aligned}
\int_A\kappa(x,B)\,P_X(dx)
&=\int_A\left[\int_B\frac{r(x,y)}{r_X(x)}\,\nu(dy)\right]r_X(x)\,\mu(dx)\\
&=\int_A\int_Br(x,y)\,\nu(dy)\,\mu(dx)\\
&=P(X\in A,Y\in B),
\end{aligned}
$$

with the exceptional set omitted in the ratio calculation. Thus “joint divided by marginal” follows from the general construction rather than division by a probability-zero event.

When $E|Y|<\infty$, the conditional expectation can be recovered as

$$
h(x)=\int y\,\kappa(x,dy)
$$

for $P_X$-almost every $x$ ([CMU probability notes](https://www.stat.cmu.edu/~arinaldo/Teaching/36752/S18/Notes/lec_notes_6.pdf)). A sum, a density integral, and a point-mass expectation are all instances of this formula.

### Null events and the limits of density formulas

For a simple illustration, let $X$ be uniform on $(0,1)$ and set $Y=X$. Both marginals have PDFs, but the joint law lies on a diagonal of two-dimensional Lebesgue measure zero, so there is no ordinary two-dimensional joint PDF.

The natural conditional law nevertheless exists: $\kappa(x,\cdot)=\delta_x$, where $\delta_x$ denotes unit mass at $x$. Consequently $E[Y\mid X=x]=x$ in this version.

More generally, a probability-zero event by itself does not select a unique conditional law. Different conditioning variables or limiting schemes can yield different versions or limits, but it is too strong to say that conditioning on $Y-X=0$ and on $Y/X=1$ always gives different answers.

The practical rule is to specify the conditioning variable or sigma-algebra. When a particular null-event value matters, also state the version, continuity convention, or limiting construction being used.

### Translation table

| Computational language | Formal interpretation |
|---|---|
| Conditional PMF or PDF | Density of $\kappa(x,\cdot)$ relative to counting or Lebesgue measure, when such a density exists |
| $E[Y\mid X=x]$ | $h(x)$ for a chosen version with $E[Y\mid X]=h(X)$ a.s. |
| $0<f_X(x)<\infty$ | The density ratio is valid at that value for the chosen density construction |
| Almost-everywhere uniqueness | Conditional laws or means may differ on a $P_X$-null set |
| Integrating out a variable | Marginalization using Tonelli, or Fubini when absolute integrability is needed |

## Pitfalls and lecture preparation

### Common mistakes

- **Wrong marginal:** Divide $f_{Y\mid X}$ by $f_X$, not $f_Y$.
- **Wrong direction:** To obtain $f_X$, integrate out $y$; do not integrate out the variable being retained.
- **Wrong support:** Recompute the allowed slice after conditioning, especially on triangles or piecewise regions.
- **Density as probability:** Integrate the conditional density over an event; do not interpret its height as a point probability.
- **Mean as a distribution:** Equal conditional and marginal means do not establish independence.
- **Variance as a second moment:** Subtract the squared conditional mean from the conditional second moment.
- **Unweighted outer average:** Total expectation uses the marginal distribution of the conditioning variable.
- **Incomplete factorization:** Include support indicators when testing independence.
- **Unjustified pointwise claim:** Continuous conditional laws are unique only almost everywhere.
- **Unnormalized kernel recognition:** Identifying a familiar algebraic shape does not excuse incorrect support or a missing constant.

### Self-test

Try the questions before reading the short answers. The normal-model and kernel questions can serve as a second-pass or Part II preview.

1. For $f(x,y)=2$ on $0<x<y<1$, find $f_X(x)$, $f_{Y\mid X}(y\mid x)$, and $E[Y\mid X=x]$.
2. Why is the joint density not enough to calculate a conditional probability without normalization?
3. Show from the formal definition that $E[E[Y\mid\mathcal G]]=E[Y]$.
4. Is $F_{X,Y}(x,t)/F_X(x)$ generally the CDF of $Y\mid X=x$?
5. For the jointly normal model, identify the algebraic identity needed to derive the conditional density.
6. Use $N\mid X$ in the thinning example to recover $\operatorname{Var}(N)$.
7. Does having marginal PDFs for $X$ and $Y$ guarantee a joint PDF?
8. What does a regular conditional distribution provide beyond separate conditional probabilities for individual events?

### Short answers

1. **Uniform slice:** $f_X(x)=2(1-x)$ for $0<x<1$; $f_{Y\mid X}(y\mid x)=1/(1-x)$ for $x<y<1$; mean $(1+x)/2$.
2. **Total conditional mass:** A joint slice generally has total mass or integrated height equal to its marginal, not 1.
3. **Whole-space event:** Choose $G=\Omega$ in the defining integral identity.
4. **Different conditioning event:** Where the denominator is positive, that ratio gives $P(Y\le t\mid X\le x)$, not generally $P(Y\le t\mid X=x)$.
5. **Complete the square:** With standardized $u=(x-\mu_X)/\sigma_X$ and $v=(y-\mu_Y)/\sigma_Y$, use $u^2-2\rho uv+v^2=(v-\rho u)^2+(1-\rho^2)u^2$.
6. **Total variance:** $\operatorname{Var}(N)=\lambda q+\operatorname{Var}(X+\lambda q)=\lambda q+\lambda p=\lambda$.
7. **No:** The diagonal example $Y=X$ has continuous marginals but a singular joint law.
8. **A coherent law:** At each conditioning value, the event probabilities assemble into a countably additive probability measure, with measurable dependence on the conditioning value.

### Questions to bring to lecture

- **Scope:** Which results are required for Part I, and which are being deferred to Part II?
- **Formal definitions:** Are Radon–Nikodym and regular conditional distribution proofs examinable, or only their statements?
- **Pointwise conventions:** Does the course implicitly use the natural continuous version of a conditional density?
- **Bounds:** What level of support sketching or piecewise setup is expected in solutions?
- **Model extensions:** Will mixed discrete–continuous laws or singular laws be included?

## R verification

The following blocks are designed to run in order in a clean base-R session and were tested with R 4.5.2; all 30 deterministic checks passed. Deterministic checks use `stopifnot`, while the final simulation reports sampling uncertainty rather than pretending that a random estimate must equal its population target.

### Discrete table checks

```r
checks_run <- 0L
check <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(
    as.numeric(actual), as.numeric(expected), tolerance = tolerance
  )))
  checks_run <<- checks_run + 1L
}

P <- matrix(
  c(.10, .20, .10,
    .20, .15, .25),
  nrow = 2, byrow = TRUE,
  dimnames = list(y = 0:1, x = 1:3)
)
px <- colSums(P)
py <- rowSums(P)
cond <- sweep(P, 2, px, "/")
h <- cond["1", ]
v <- h * (1 - h)

check(sum(P), 1)
check(px, c(.30, .35, .35))
check(py, c(.40, .60))
check(colSums(cond), rep(1, 3))
check(h, c(2/3, 3/7, 5/7))
check(v, c(2/9, 12/49, 10/49))
check(sum(px * h), .60)
check(sum((1:3) * P["1", ] / py["1"]), 25/12)
check(sum(px * v) + sum(px * h^2) - sum(px * h)^2, .60 * .40)
```

### Square-density checks

```r
joint2 <- function(x, y) x + y  # Used only on the unit square.
fx2 <- function(x) x + 1/2
cond2 <- function(y, x) (x + y) / fx2(x)
h2 <- function(x) {
  integrate(function(y) y * cond2(y, x), 0, 1)$value
}

check(integrate(fx2, 0, 1)$value, 1)
check(integrate(function(y) joint2(.5, y), 0, 1)$value, fx2(.5))
check(integrate(function(y) cond2(y, .5), 0, 1)$value, 1)
check(h2(.5), 7/12)
check(integrate(function(y) cond2(y, .5), 0, .5)$value, 3/8)
check(integrate(function(x) Vectorize(h2)(x) * fx2(x), 0, 1)$value, 7/12)
```

### Triangular-support checks

```r
fx3 <- function(x) 4 * x * (1 - x^2)
fy3 <- function(y) 4 * y^3
hy3 <- function(x) {
  stopifnot(length(x) == 1L, x > 0, x < 1)
  integrate(function(y) y * 2 * y / (1 - x^2), x, 1)$value
}
hx3 <- function(y) {
  stopifnot(length(y) == 1L, y > 0, y < 1)
  integrate(function(x) x * 2 * x / y^2, 0, y)$value
}

check(integrate(fx3, 0, 1)$value, 1)
check(integrate(fy3, 0, 1)$value, 1)
check(integrate(function(y) 2 * y / (1 - .4^2), .4, 1)$value, 1)
check(integrate(function(x) 2 * x / .6^2, 0, .6)$value, 1)
check(hy3(.4), 2 * (1 + .4 + .4^2) / (3 * (1 + .4)))
check(hx3(.6), 2 * .6 / 3)
check(integrate(function(y) 2 * y / (1 - .5^2), .75, 1)$value, 7/12)
check(integrate(function(x) Vectorize(hy3)(x) * fx3(x), 0, 1)$value, 4/5)
check(integrate(function(y) Vectorize(hx3)(y) * fy3(y), 0, 1)$value, 8/15)
```

### Thinning checks and simulation

```r
lambda <- 10
p <- .3
q <- 1 - p
x <- 4
n <- x + 0:100
joint <- dpois(n, lambda) * dbinom(x, size = n, prob = p)
conditional <- joint / dpois(x, lambda * p)

# Finite sums approximate infinite tails; 100 extra counts is ample here.
check(sum(joint), dpois(x, lambda * p))
check(conditional, dpois(n - x, lambda * q))
check(sum(conditional), 1)
check(sum(n * conditional), 11)
check(sum((n - 11)^2 * conditional), 7)
check(lambda * p + lambda * q, lambda)
cat("Passed", checks_run, "deterministic checks.\n")

set.seed(544)
draws <- 200000L
N <- rpois(draws, lambda)
X <- rbinom(draws, size = N, prob = p)
selected <- N[X == x]
stopifnot(length(selected) > 1L)

estimate <- mean(selected)
mc_se <- sd(selected) / sqrt(length(selected))
print(c(
  selected_draws = length(selected),
  conditional_mean_estimate = estimate,
  conditional_mean_target = x + lambda * q,
  Monte_Carlo_SE = mc_se,
  conditional_variance_estimate = var(selected),
  conditional_variance_target = lambda * q,
  retained_mean_estimate = mean(X),
  retained_variance_estimate = var(X),
  retained_mean_and_variance_target = lambda * p
))
```

The density functions above are evaluated only on the supports indicated by the integration limits. In a general-purpose distribution implementation, return zero outside the support and explicitly handle invalid conditioning values.

## Compact reference sheet

For the ratio formulas, assume a joint PMF or joint PDF exists and that the conditioning marginal is positive, and finite in the density case. Always retain the conditional support.

$$
p_X(x)=\sum_y p_{X,Y}(x,y),\qquad
f_X(x)=\int f_{X,Y}(x,y)\,dy.
$$

$$
p_{Y\mid X}(y\mid x)=\frac{p_{X,Y}(x,y)}{p_X(x)},\qquad
f_{Y\mid X}(y\mid x)=\frac{f_{X,Y}(x,y)}{f_X(x)}.
$$

$$
h(x)=E[Y\mid X=x]
=\sum_y y\,p_{Y\mid X}(y\mid x)
\quad\text{or}\quad
\int y\,f_{Y\mid X}(y\mid x)\,dy.
$$

$$
E[Y\mid X]=h(X),\qquad E[h(X)]=E[Y].
$$

$$
\operatorname{Var}(Y\mid X=x)=E[Y^2\mid X=x]-h(x)^2.
$$

$$
\operatorname{Var}(Y)
=E[\operatorname{Var}(Y\mid X)]
+\operatorname{Var}(E[Y\mid X]).
$$

Formally, $E[Y\mid X]$ is the integrable $\sigma(X)$-measurable random variable matching the integrals of $Y$ over every event in $\sigma(X)$. A regular conditional distribution is a probability kernel satisfying disintegration, determined only $P_X$-almost everywhere.

Before finalizing a solution, ask: **What is fixed, what still varies, and does the conditional law have total mass 1?** Then use an independently calculated marginal quantity to check the answer whenever practical.
