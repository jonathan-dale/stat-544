# STAT 544 Study Guide and Lecture Prep: Poisson Process I

A Poisson process connects three descriptions of the same random phenomenon: how many events occur by a fixed time, how long we wait between events, and when the successive events occur. The central progression in this guide is **process assumptions → independent increments → Poisson counts → exponential waiting times → gamma arrival times**.

This guide follows the prerequisite review, explicit derivations, worked problems, formal reference, and verification approach of your [Conditional Distributions I](https://github.com/jonathan-dale/stat-544/blob/main/guides/conditional_distributions/cd_1/Conditional_Distributions_I.md) and [Conditional Distributions II](https://github.com/jonathan-dale/stat-544/blob/main/guides/conditional_distributions/cd_2/Conditional_Distributions_II.md) guides. It is a graduate-level preparation document, not an instructor-issued syllabus; the exact lecture notes for “Poisson process I” and “Poisson process II” were not supplied, so the I/II division here is a study organization rather than a verified course boundary.

## How to use this guide

On a first pass, focus on the count–waiting-time dictionary and the worked examples. On a second pass, reproduce the derivations without looking, then use the formal reference to connect the calculations to stochastic-process theory.

- **Core foundations:** Counting processes, rate and units, stationary and independent increments, Poisson probabilities, exponential waiting times, and gamma arrival times.
- **Graduate derivations:** Small-interval characterization, forward equations, finite-dimensional distributions, covariance, and joint arrival-time densities.
- **Optional advanced reference:** Filtrations, compensated martingales, stopping times, and the generator. These are enrichment topics unless your instructor explicitly includes them.
- **Bridge to Part II:** Conditioning on a total count, uniform order statistics, thinning, superposition, and time-varying rates.

### Contents

- [Prerequisite review](#prerequisite-review)
- [Defining a homogeneous Poisson process](#defining-a-homogeneous-poisson-process)
- [Deriving the count distribution](#deriving-the-count-distribution)
- [Counts, waiting times, and arrival epochs](#counts-waiting-times-and-arrival-epochs)
- [Joint laws, dependence, and prediction](#joint-laws-dependence-and-prediction)
- [A reliable problem-solving procedure](#a-reliable-problem-solving-procedure)
- [Six fully worked examples](#six-fully-worked-examples)
- [Formal graduate reference](#formal-graduate-reference)
- [Practice problems](#practice-problems)
- [Solutions and checks](#solutions-and-checks)
- [R verification and simulation](#r-verification-and-simulation)
- [Compact reference sheet](#compact-reference-sheet)

## Prerequisite review

### A random variable versus a stochastic process

A random variable assigns a numerical value to an outcome; a stochastic process assigns a family of random variables indexed by time. For an arrival process, the count at a fixed time and the entire evolving count path are different mathematical objects, a distinction used throughout [Gallager’s development of Poisson processes](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process).

Write

$$
N(t,\omega)=\text{number of events in }(0,t]
$$

for a particular outcome $\omega$. Usually we suppress $\omega$ and write $N(t)$.

- **Fixed time $t$:** $N(t)$ is a nonnegative integer-valued random variable.
- **Fixed outcome $\omega$:** $t\mapsto N(t,\omega)$ is a sample path.
- **Whole process:** $\{N(t):t\ge 0\}$ describes counts across all times, including their joint dependence.
- **Increment:** $N(b)-N(a)$ counts events in $(a,b]$, for $0\le a<b$.

Knowing the distribution of $N(t)$ separately at every time is not enough to know the process. One must also specify how counts at different times depend on one another.

### Notation and parameter conventions

Throughout this guide, $\lambda>0$ denotes a constant event rate. Time must be measured in a unit compatible with that rate.

| Symbol | Meaning | Type or units |
|---|---|---|
| $\lambda$ | Arrival rate | Events per unit time |
| $N(t)$ | Number of events in $(0,t]$ | Count |
| $T_n$ | Time of the $n$ th event | Time |
| $T_0=0$ | Starting point, not an actual arrival | Time |
| $X_n=T_n-T_{n-1}$ | $n$ th interarrival time | Time |
| $\lambda t$ | Expected count by time $t$ | Dimensionless Poisson parameter |
| $\mathcal F_t$ | Information available through time $t$ | Sigma-algebra |

Some books use $S_n$ for the $n$ th arrival time and $T_n$ for an interarrival time. Translate notation before using a theorem; here $T_n$ always means an arrival epoch, not a gap.

The distribution $\mathrm{Gamma}(\alpha,\lambda)$ below uses **shape $\alpha$ and rate $\lambda$**, not scale $\lambda$. Its scale is $1/\lambda$.

### Three prerequisite distributions

The following distributions form the basic count–waiting-time dictionary in [Gallager’s Poisson-process notes](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process). The displayed moments and transforms can also be recovered directly from these mass and density formulas.

For $Y\sim\mathrm{Poisson}(\mu)$,

$$
P(Y=k)=e^{-\mu}\frac{\mu^k}{k!},
\qquad k=0,1,2,\ldots,
$$

$$
E[Y]=\mathrm{Var}(Y)=\mu,
\qquad
M_Y(\theta)=\exp\{\mu(e^\theta-1)\}.
$$

For $X\sim\mathrm{Exp}(\lambda)$,

$$
f_X(x)=\lambda e^{-\lambda x}\mathbf 1_{\{x>0\}},
\qquad
P(X>x)=e^{-\lambda x},\quad x\ge 0,
$$

$$
E[X]=\frac1\lambda,
\qquad
\mathrm{Var}(X)=\frac1{\lambda^2},
\qquad
M_X(\theta)=\frac{\lambda}{\lambda-\theta},\quad \theta<\lambda.
$$

For $G\sim\mathrm{Gamma}(\alpha,\lambda)$,

$$
f_G(t)=\frac{\lambda^\alpha}{\Gamma(\alpha)}
t^{\alpha-1}e^{-\lambda t}\mathbf 1_{\{t>0\}},
$$

$$
E[G]=\frac{\alpha}{\lambda},
\qquad
\mathrm{Var}(G)=\frac{\alpha}{\lambda^2}.
$$

When $\alpha=n$ is a positive integer, $\Gamma(n)=(n-1)!$, and this gamma distribution is also called an Erlang distribution. The exponential is the special case $\alpha=1$.

### Expectation and conditional information

For square-integrable variables,

$$
\mathrm{Var}(Y)=E[Y^2]-(E[Y])^2,
\qquad
\mathrm{Cov}(Y,Z)=E[YZ]-E[Y]E[Z].
$$

Linearity of expectation does not require independence. Adding variances without covariance terms does require suitable uncorrelatedness, which independence supplies.

The conditional-distribution tools from the earlier guides remain central:

$$
P(A\mid B)=\frac{P(A\cap B)}{P(B)},\qquad P(B)>0,
$$

$$
E[Y]=E[E[Y\mid Z]],
$$

$$
\mathrm{Var}(Y)
=E[\mathrm{Var}(Y\mid Z)]
+\mathrm{Var}(E[Y\mid Z]).
$$

In process problems, first identify whether the information concerns the past, a future interval, or a fixed total over a larger interval. These are not interchangeable conditioning statements.

## Defining a homogeneous Poisson process

### The count-and-increments definition

The standard homogeneous Poisson-process definition combines Poisson interval counts with independent increments; equivalent constructions use independent exponential interarrival times or small-interval probabilities ([Gallager’s definitions and equivalence discussion](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process)). Here we use the right-continuous counting-process version.

A homogeneous Poisson process of rate $\lambda>0$ satisfies:

1. **Initial condition:** $N(0)=0$ almost surely.
2. **Counting paths:** Paths are nondecreasing, integer-valued, right-continuous, and have unit jumps, with finitely many jumps on every bounded time interval almost surely.
3. **Independent increments:** Counts in disjoint time intervals are mutually independent.
4. **Poisson increments:** For $0\le s<t$,

$$
N(t)-N(s)\sim\mathrm{Poisson}(\lambda(t-s)).
$$

The last property gives stationary increments: an interval count depends in distribution only on the interval length, not on where the interval starts. We write $N\sim\mathrm{HPP}(\lambda)$ as shorthand, where HPP means homogeneous Poisson process.

### Independent increments are not independent cumulative counts

For $0<s<t$, the variables $N(s)$ and $N(t)-N(s)$ concern disjoint intervals. They are independent by definition.

But

$$
N(t)=N(s)+[N(t)-N(s)]
$$

contains the earlier count. Consequently, $N(s)$ and $N(t)$ are generally dependent; later we calculate their covariance explicitly.

Stationarity also needs precise wording. The process itself is not stationary, since the distribution of $N(t)$ changes with $t$; its increments are stationary.

### The small-interval characterization

An equivalent local characterization starts with a counting process with stationary independent increments and the following behavior as $h\downarrow0$ ([Gallager’s infinitesimal characterization](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process)):

$$
P(N(h)=0)=1-\lambda h+o(h),
$$

$$
P(N(h)=1)=\lambda h+o(h),
$$

$$
P(N(h)\ge2)=o(h).
$$

The notation $r(h)=o(h)$ means

$$
\lim_{h\downarrow0}\frac{r(h)}h=0.
$$

For example, $h^2=o(h)$, but $0.01h$ is not $o(h)$. The multiple-arrival condition rules out batches arriving with probability proportional to $h$.

For an actual Poisson increment, put $x=\lambda h$ and expand:

$$
P(N(h)\ge2)
=1-e^{-x}(1+x)
=\frac{x^2}{2}+O(x^3).
$$

Thus two or more arrivals can occur in a short interval, but their probability is smaller-order than the interval length. “Negligible to first order” does not mean impossible.

### What the model assumes, and what it does not

Use the model because its assumptions are reasonable for a particular observation scale, not merely because the data are counts. These are modeling checks rather than automatic consequences of seeing a roughly stable average.

- **Constant rate:** Equal-length intervals have the same expected count under an HPP.
- **No event-triggered dependence:** An event does not raise or lower the subsequent rate in this model.
- **No simultaneous batches:** Jumps are single events; batch counts lead to a different process.
- **No deterministic spacing:** Constant average rate does not mean events occur every $1/\lambda$ time units.
- **No “one event per interval” restriction:** Any finite interval of positive length can contain arbitrarily many events.

Equality of sample mean and sample variance is not a proof of a Poisson process. A process model also makes claims about dependence, timing, and stability of the rate.

## Deriving the count distribution

### Forward equations from small intervals

Let

$$
p_n(t)=P(N(t)=n).
$$

To derive its evolution, decompose the count at $t+h$ according to what happens in $(t,t+h]$. Independent increments allow the probabilities for the old count and the new increment to multiply.

For $n=0$,

$$
p_0(t+h)=p_0(t)(1-\lambda h+o(h)),
$$

so

$$
p_0'(t)=-\lambda p_0(t),
\qquad p_0(0)=1.
$$

Therefore

$$
p_0(t)=e^{-\lambda t}.
$$

For $n\ge1$, either the old count is $n$ and no arrival occurs, or the old count is $n-1$ and one occurs. All possibilities with at least two new arrivals contribute an $o(h)$ remainder:

$$
p_n(t+h)
=p_n(t)(1-\lambda h)
+p_{n-1}(t)\lambda h+o(h).
$$

Subtract $p_n(t)$, divide by $h$, and take the limit:

$$
\boxed{
p_n'(t)=-\lambda p_n(t)+\lambda p_{n-1}(t),
\quad n\ge1.
}
$$

The initial condition is $p_n(0)=0$ for $n\ge1$. These are the forward equations of a pure-birth process with constant birth rate $\lambda$.

### Solving the recursion

Multiplying the $n$ th equation by $e^{\lambda t}$ gives

$$
\frac{d}{dt}\{e^{\lambda t}p_n(t)\}
=\lambda e^{\lambda t}p_{n-1}(t).
$$

Assuming inductively that

$$
p_{n-1}(t)=e^{-\lambda t}
\frac{(\lambda t)^{n-1}}{(n-1)!},
$$

we obtain

$$
e^{\lambda t}p_n(t)
=\int_0^t\frac{\lambda^n u^{n-1}}{(n-1)!}\,du
=\frac{(\lambda t)^n}{n!}.
$$

Hence

$$
\boxed{
P(N(t)=n)=e^{-\lambda t}\frac{(\lambda t)^n}{n!},
\quad n=0,1,2,\ldots.
}
$$

The normalization check uses the exponential series:

$$
\sum_{n=0}^\infty e^{-\lambda t}\frac{(\lambda t)^n}{n!}
=e^{-\lambda t}e^{\lambda t}=1.
$$

This establishes the count formula from the local assumptions. If Poisson increments were taken as the starting definition, the same calculation is a consistency check rather than a new assumption.

### An alternative derivation using the probability generating function

For $0\le z\le1$, define

$$
G_t(z)=E[z^{N(t)}].
$$

Independent stationary increments imply

$$
G_{t+h}(z)=G_t(z)G_h(z).
$$

Using the local probabilities,

$$
G_h(z)=1+\lambda h(z-1)+o(h).
$$

The resulting differential equation is

$$
\frac{\partial G_t(z)}{\partial t}
=\lambda(z-1)G_t(z),
\qquad G_0(z)=1,
$$

with solution

$$
G_t(z)=\exp\{\lambda t(z-1)\}.
$$

Expanding this expression in powers of $z$ recovers the same Poisson PMF. This second method highlights the additive-in-time structure rather than solving each probability separately.

## Counts, waiting times, and arrival epochs

### Translating events before calculating

Define

$$
T_n=\inf\{t\ge0:N(t)\ge n\},
\qquad n\ge1.
$$

The fundamental event identities are

$$
\boxed{\{T_n\le t\}=\{N(t)\ge n\}},
$$

$$
\boxed{\{T_n>t\}=\{N(t)\le n-1\}},
$$

$$
\{N(t)=n\} = \{T_n\le t \le T_{n+1}\}.
$$

Each identity describes the same sample paths in two ways. These count–arrival relations underpin the Poisson/gamma connection in [Gallager’s arrival-time derivation](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process).

The safest exam habit is to write the event identity before evaluating a CDF. “The third event occurs after $t$” means at most two events by $t$, not at most three.

### The first waiting time is exponential

No event by time $t$ is equivalent to waiting longer than $t$ for the first event:

$$
P(T_1>t)=P(N(t)=0)=e^{-\lambda t},
\qquad t\ge0.
$$

Therefore

$$
T_1\sim\mathrm{Exp}(\lambda).
$$

Differentiating $F_{T_1}(t)=1-e^{-\lambda t}$ gives its density. Notice the units: $1/\lambda$ is a mean waiting time, whereas $\lambda t$ is a mean count.

### Memorylessness

For $X\sim\mathrm{Exp}(\lambda)$ and $a,b\ge0$,

$$
\begin{aligned}
P(X>a+b\mid X>a)
&=\frac{P(X>a+b)}{P(X>a)}\\
&=\frac{e^{-\lambda(a+b)}}{e^{-\lambda a}}
=e^{-\lambda b}.
\end{aligned}
$$

Thus

$$
X-a\mid\{X>a\}\sim\mathrm{Exp}(\lambda).
$$

This says that the remaining wait has the original distribution after learning that the event has not yet occurred. It does not say that the total elapsed-plus-future wait is exponential with the same origin.

For example,

$$
E[X\mid X>a]=a+\frac1\lambda,
$$

whereas

$$
E[X-a\mid X>a]=\frac1\lambda.
$$

### All interarrival times

For an HPP, the successive gaps satisfy

$$
\boxed{X_1,X_2,\ldots\ \text{are iid }\mathrm{Exp}(\lambda).}
$$

One can either construct the process from these gaps or derive the statement from the count-based definition using the restart property at arrival times ([Gallager’s renewal construction](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process); [University of Toronto’s stopping-time formulation](https://utstat.utoronto.ca/mikevans/stac63/lecture6b.pdf)). Restarting at a random arrival time is a strong Markov statement, not merely a direct substitution into deterministic-time independent increments.

Conversely, start with iid exponential gaps and set

$$
T_n=\sum_{j=1}^n X_j,
\qquad
N(t)=\max\{n\ge0:T_n\le t\}.
$$

The strong law gives $T_n/n\to1/\lambda$ almost surely, so $T_n\to\infty$ and the number of events on bounded intervals is finite. Exponential memorylessness then yields stationary independent increments, producing the HPP.

### The $n$ th arrival has a gamma distribution

Since $T_n$ is a sum of $n$ independent exponential gaps,

$$
M_{T_n}(\theta)
=\left(\frac{\lambda}{\lambda-\theta}\right)^n,
\qquad \theta<\lambda.
$$

This is the MGF of $\mathrm{Gamma}(n,\lambda)$. Thus

$$
\boxed{
f_{T_n}(t)=\frac{\lambda^n t^{n-1}e^{-\lambda t}}{(n-1)!},
\quad t>0,
}
$$

$$
E[T_n]=\frac n\lambda,
\qquad
\mathrm{Var}(T_n)=\frac n{\lambda^2}.
$$

The same law follows without MGFs:

$$
P(T_n>t)
=P(N(t)\le n-1)
=e^{-\lambda t}\sum_{k=0}^{n-1}\frac{(\lambda t)^k}{k!}.
$$

Differentiating the complement of this survival function makes all but one term cancel. The remaining term is exactly the gamma density above.

### A direct convolution proof for the second arrival

For $t>0$, the support of $X_1$ in the convolution is $0<x<t$:

$$
\begin{aligned}
f_{T_2}(t)
&=\int_0^t f_{X_1}(x)f_{X_2}(t-x)\,dx\\
&=\int_0^t
\lambda e^{-\lambda x}\lambda e^{-\lambda(t-x)}\,dx\\
&=\lambda^2e^{-\lambda t}\int_0^t1\,dx\\
&=\lambda^2te^{-\lambda t}.
\end{aligned}
$$

The exponential terms cancel their dependence on $x$, but the integration interval still depends on $t$. Forgetting its length would lose the factor $t$.

### Waiting from a fixed observation time

For deterministic $s\ge0$, let

$$
W_s=\inf\{u>0:N(s+u)>N(s)\}.
$$

The entire future increment process is independent of the history through $s$, so

$$
P(W_s>u\mid\mathcal F_s)
=P(N(s+u)-N(s)=0)
=e^{-\lambda u}.
$$

Therefore $W_s\sim\mathrm{Exp}(\lambda)$ and is independent of the past. This deterministic-time restart result is also established directly in [Gallager’s memoryless-process theorem](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process).

The time until the $r$ th future event after $s$ is $\mathrm{Gamma}(r,\lambda)$. But if the question is about an absolute event number $n$, first subtract the number already observed.

Given $N(s)=m<n$,

$$
T_n-s\mid\{N(s)=m\}\sim\mathrm{Gamma}(n-m,\lambda).
$$

If $m\ge n$, then $T_n\le s$ already. A future-waiting gamma formula would answer the wrong question.

## Joint laws, dependence, and prediction

### Finite-dimensional count distributions

For $0=t_0<t_1<\cdots<t_k$ and integers $0=n_0\le n_1\le\cdots\le n_k$, rewrite cumulative counts as disjoint increments:

$$
\begin{aligned}
&P(N(t_1)=n_1,\ldots,N(t_k)=n_k)\\
&\quad=
\prod_{j=1}^k
e^{-\lambda(t_j-t_{j-1})}
\frac{[\lambda(t_j-t_{j-1})]^{n_j-n_{j-1}}}
{(n_j-n_{j-1})!}.
\end{aligned}
$$

The probability is zero for an invalid count ordering. This product follows from the defining independence of increments and specifies the joint law rather than only its one-time marginals.

In particular, for $0<s<t$ and $0\le i\le j$,

$$
P(N(s)=i,N(t)=j)
=e^{-\lambda t}
\frac{(\lambda s)^i}{i!}
\frac{[\lambda(t-s)]^{j-i}}{(j-i)!}.
$$

### Covariance of cumulative counts

For $0\le s\le t$,

$$
\begin{aligned}
\mathrm{Cov}(N(s),N(t))
&=\mathrm{Cov}(N(s),N(s)+N(t)-N(s))\\
&=\mathrm{Var}(N(s))+0\\
&=\lambda s.
\end{aligned}
$$

Hence, for arbitrary nonnegative $s,t$,

$$
\boxed{\mathrm{Cov}(N(s),N(t))=\lambda\min(s,t).}
$$

For $s,t>0$, the correlation is

$$
\mathrm{Corr}(N(s),N(t))
=\frac{\min(s,t)}{\sqrt{st}}.
$$

For two bounded time intervals $A$ and $B$, decomposition into their shared and disjoint portions similarly gives

$$
\mathrm{Cov}(N(A),N(B))=\lambda |A\cap B|,
$$

where $N(A)$ counts arrivals in $A$ and $|A\cap B|$ is the overlap length. Disjoint intervals have zero overlap and independent counts.

### Conditional prediction of a future count

For $t>s$, write $Z=N(t)-N(s)$. Since $Z$ is independent of $\mathcal F_s$ and has mean and variance $\lambda(t-s)$,

$$
N(t)\mid\mathcal F_s
\quad\text{has the law of}\quad
N(s)+Z,
\qquad
Z\sim\mathrm{Poisson}(\lambda(t-s)).
$$

Consequently,

$$
E[N(t)\mid\mathcal F_s]=N(s)+\lambda(t-s),
$$

$$
\mathrm{Var}(N(t)\mid\mathcal F_s)=\lambda(t-s).
$$

If only $N(s)=m$ is given, replace $N(s)$ by $m$. The shifted Poisson support starts at $m$, not at zero.

### Joint arrival-time density

Using $X_1=T_1$ and $X_j=T_j-T_{j-1}$, the transformation from arrival epochs to gaps has a triangular Jacobian with determinant 1. Multiplying the exponential gap densities gives

$$
\begin{aligned}
f_{T_1,\ldots,T_n}(t_1,\ldots,t_n)
&=\lambda^n e^{-\lambda[t_1+(t_2-t_1)+\cdots+(t_n-t_{n-1})]}\\
&=\lambda^n e^{-\lambda t_n},
\end{aligned}
$$

for $0<t_1<\cdots<t_n$.


This ordered-support density is the arrival-time formula developed in [Gallager’s joint-density argument](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process). Its simple algebra does not make the arrival times independent; the support and the shared gaps matter.

For example,

$$
\mathrm{Cov}(T_i,T_j)=\frac{\min(i,j)}{\lambda^2},
$$

because $T_i$ and $T_j$ share exactly the first $\min(i,j)$ independent gaps. Integrating the joint density over $0<t_1<\cdots<t_{n-1}<t_n$ also recovers the gamma marginal of $T_n$.

## A reliable problem-solving procedure

Translate the event and state the information before selecting a distribution. Most errors in this topic occur before the arithmetic.

1. **Name the object:** Is the question about an interval count, one gap, a numbered arrival, or several times jointly?
2. **Standardize units:** Convert the rate and all times to the same time unit.
3. **Write the event:** For instance, $\{T_4>t\}=\{N(t)\le3\}$.
4. **Identify conditioning:** Distinguish observed history from a known total over a larger window.
5. **Choose the representation:** Use independent increments for counts and exponential sums for arrival times.
6. **State support:** Counts are integers; waits are nonnegative; arrival epochs are ordered.
7. **Calculate:** Give the symbolic setup before a numerical approximation or R call.
8. **Cross-check:** Check units, probability bounds, monotonicity, moments, or a second distributional representation.

The main warning is **do not condition away your independence assumption**. Part II shows that interval counts independent before fixing a total become dependent after that total is fixed.

## Six fully worked examples

### Counts in disjoint intervals

Suppose requests arrive according to an HPP with rate $\lambda=6$ per hour. Find the probability of exactly two requests in the first 20 minutes and at least one in the next 10 minutes.

Convert the durations:

$$
t_1=\frac13\text{ hour},
\qquad
t_2=\frac16\text{ hour}.
$$

Define

$$
A=N(1/3),
\qquad
B=N(1/2)-N(1/3).
$$

These count disjoint intervals, so

$$
A\sim\mathrm{Poisson}(2),
\qquad
B\sim\mathrm{Poisson}(1),
\qquad A\perp B.
$$

The requested probability is

$$
\begin{aligned}
P(A=2,B\ge1)
&=P(A=2)P(B\ge1)\\
&=\left(e^{-2}\frac{2^2}{2!}\right)(1-e^{-1})\\
&=2e^{-2}(1-e^{-1})
\approx0.171096.
\end{aligned}
$$

The parameter for the second interval is 1, not the mean count 3 for the entire first half-hour. If the total half-hour count were also specified, this unconditional multiplication would no longer be the right calculation.

### The third arrival after a threshold

Let $\lambda=2$ per hour. Find $P(T_3>2.5)$ and the mean and variance of $T_3$.

Translate the event:

$$
\{T_3>2.5\}=\{N(2.5)\le2\}.
$$

Since $N(2.5)\sim\mathrm{Poisson}(5)$,

$$
\begin{aligned}
P(T_3>2.5)
&=e^{-5}\left(1+5+\frac{5^2}{2!}\right)\\
&=18.5e^{-5}
\approx0.124652.
\end{aligned}
$$

Alternatively, $T_3\sim\mathrm{Gamma}(3,2)$ and the probability is its upper tail. The moments are

$$
E[T_3]=\frac32\text{ hours},
\qquad
\mathrm{Var}(T_3)=\frac34\text{ hours}^2.
$$

The standard deviation is $\sqrt3/2$ hours. Variance has squared-time units, not time units.

### Updating after a period with no arrivals

Let $\lambda=4$ per hour and suppose no event occurs in the first 15 minutes. Find the probability of waiting at least another 10 minutes, and find the conditional expected time of the first event measured from time zero.

With $a=1/4$ hour and $b=1/6$ hour,

$$
\begin{aligned}
P(T_1>a+b\mid T_1>a)
&=e^{-\lambda b}\\
&=e^{-2/3}\approx0.513417.
\end{aligned}
$$

The remaining mean wait is $1/\lambda=1/4$ hour. The conditional expected absolute arrival time is therefore

$$
E[T_1\mid T_1>1/4]
=\frac14+\frac14
=\frac12\text{ hour}.
$$

One expectation is “another 15 minutes”; the other is “30 minutes after time zero.” Keeping the time origin visible prevents a common memorylessness mistake.

### Joint cumulative counts and future prediction

Let $\lambda=3$ per hour. Find $P(N(1)=2,N(2)=5)$, $\mathrm{Cov}(N(1),N(2))$, and the conditional law of $N(2)$ given $N(1)=2$.

Rewrite the event as

$$
\{N(1)=2,\ N(2)-N(1)=3\}.
$$

The two increments are independent and each has mean 3:

$$
\begin{aligned}
P(N(1)=2,N(2)=5)
&=\left(e^{-3}\frac{3^2}{2!}\right)
\left(e^{-3}\frac{3^3}{3!}\right)\\
&=\frac{81}{4}e^{-6}
\approx0.050195.
\end{aligned}
$$

The covariance is

$$
\mathrm{Cov}(N(1),N(2))=3.
$$

For prediction,

$$
N(2)\mid\{N(1)=2\}\ \overset d=\ 2+Z,
\qquad Z\sim\mathrm{Poisson}(3).
$$

Thus its conditional mean is 5, its conditional variance is 3, and its support is $\{2,3,4,\ldots\}$. The cumulative counts are dependent even though the new increment is independent of the old count.

### A joint arrival-time density on a triangle

Let $\lambda=2$ per hour. Derive the joint density of $(T_1,T_2)$, then calculate $E[T_1\mid T_2=t]$.

The independent gaps have joint density

$$
f_{X_1,X_2}(x_1,x_2)=4e^{-2(x_1+x_2)},
\qquad x_1,x_2>0.
$$

Use $x_1=t_1$, $x_2=t_2-t_1$, with Jacobian determinant 1. The result is

$$
f_{T_1,T_2}(t_1,t_2)
=4e^{-2t_2},
\qquad 0<t_1<t_2,
$$

and zero elsewhere. Integrating out $t_1$ gives

$$
f_{T_2}(t_2)
=\int_0^{t_2}4e^{-2t_2}\,dt_1
=4t_2e^{-2t_2},\qquad t_2>0.
$$

Divide the joint density by this marginal:

$$
f_{T_1\mid T_2}(u\mid t)
=\frac{4e^{-2t}}{4te^{-2t}}
=\frac1t,\qquad 0<u<t.
$$

Therefore, in the natural density-ratio version,

$$
T_1\mid\{T_2=t\}\sim\mathrm{Uniform}(0,t),
\qquad
E[T_1\mid T_2=t]=\frac t2.
$$

Check total expectation:

$$
E[E[T_1\mid T_2]]
=\frac12E[T_2]
=\frac12(1)=\frac12
=E[T_1].
$$

This directly reuses the earlier joint → marginal → conditional workflow. Conditioning on $T_2=t$ is a density-based conditional law, not division by the point probability $P(T_2=t)=0$.

Do not confuse it with conditioning on $N(t)=2$. In Part II, $E[T_1\mid N(t)=2]=t/3$, because two arrivals are merely known to occur before the endpoint rather than the second being fixed at it.

### The number of occupied observation windows

Suppose false alarms follow an HPP with rate $1/5$ per day. Over 30 days, let $D$ be the number of days containing at least one alarm; find its distribution, mean, and variance.

For $j=1,\ldots,30$, define

$$
I_j=\mathbf1_{\{N(j)-N(j-1)\ge1\}}.
$$

The daily increments are independent and each is $\mathrm{Poisson}(1/5)$. Thus the indicators are iid Bernoulli with

$$
p=P(I_j=1)=1-e^{-1/5}.
$$

Since $D=\sum_{j=1}^{30}I_j$,

$$
D\sim\mathrm{Binomial}(30,1-e^{-1/5}),
$$

$$
E[D]=30(1-e^{-1/5})\approx5.438077,
$$

$$
\mathrm{Var}(D)
=30(1-e^{-1/5})e^{-1/5}
\approx4.452321.
$$

The total alarm count is instead $N(30)\sim\mathrm{Poisson}(6)$. Counting occupied days is not counting individual alarms, and $D$ cannot exceed 30.

## Formal graduate reference

### Probability space, filtration, and adaptedness

Work on a probability space $(\Omega,\mathcal F,P)$. The natural filtration

$$
\mathcal F_t^N=\sigma(N(u):0\le u\le t)
$$

represents the observed count history, and $N(t)$ is adapted because it is measurable with respect to this information.

The deterministic-time Markov and compensated-martingale properties are presented in [Evans’s University of Toronto notes](https://utstat.utoronto.ca/mikevans/stac63/lecture6b.pdf). The calculations below make their connection to conditional expectation explicit.

For $t>s$ and integers $j\ge i$,

$$
P(N(t)=j\mid\mathcal F_s^N)
=e^{-\lambda(t-s)}
\frac{[\lambda(t-s)]^{j-i}}{(j-i)!}
\quad\text{on }\{N(s)=i\}.
$$

It is zero for $j<i$. Thus the current state and elapsed time determine the conditional law of the future state; the full earlier history adds no predictive information under the model.

### The compensated Poisson martingale

Define

$$
M(t)=N(t)-\lambda t.
$$

A martingale is an adapted, integrable process satisfying $E[M(t)\mid\mathcal F_s]=M(s)$ for $s\le t$. Here

$$
\begin{aligned}
E[M(t)\mid\mathcal F_s^N]
&=E[N(t)\mid\mathcal F_s^N]-\lambda t\\
&=N(s)+\lambda(t-s)-\lambda t\\
&=N(s)-\lambda s=M(s).
\end{aligned}
$$

Integrability follows from $E[N(t)]=\lambda t<\infty$. The martingale property concerns the centered process, not the increasing count $N(t)$ itself.

There is also a square-integrable martingale

$$
M(t)^2-\lambda t.
$$

To verify it, write $M(t)=M(s)+Z$ where $Z$ is an independent centered increment with variance $\lambda(t-s)$. Conditional expansion gives

$$
E[M(t)^2\mid\mathcal F_s^N]
=M(s)^2+\lambda(t-s),
$$

and subtracting $\lambda t$ establishes the claim. This is a useful bridge from the mean-rate correction to second-moment calculations.

### Stopping times and the strong Markov property

A random time $\tau$ is a stopping time if

$$
\{\tau\le t\}\in\mathcal F_t
\quad\text{for every }t\ge0.
$$

This means that by time $t$, one can decide whether stopping has already occurred using only information available by then. In particular, $T_n$ is a stopping time because $\{T_n\le t\}=\{N(t)\ge n\}$.

Under the usual right-continuous, complete filtration, if $\tau<\infty$ almost surely is a stopping time, the process

$$
\{N(\tau+u)-N(\tau):u\ge0\}
$$

is an HPP of rate $\lambda$, independent of the information at $\tau$ ([University of Toronto’s strong Markov statement](https://utstat.utoronto.ca/mikevans/stac63/lecture6b.pdf)). This justifies restarting at an observed arrival time.

Not every random time is a stopping time. For example, selecting the last arrival before a fixed future horizon requires future information, so one cannot automatically apply the same restart argument there.

Optional stopping also has hypotheses. A martingale identity proved at deterministic times cannot simply be evaluated at an arbitrary random time without checking an applicable stopping theorem.

### The infinitesimal generator

For a bounded function $g$ on $\{0,1,2,\ldots\}$, define the generator through the short-time conditional change:

$$
Ag(n)
=\lim_{h\downarrow0}
\frac{E[g(N(t+h))\mid N(t)=n]-g(n)}h.
$$

Using the local zero- and one-jump probabilities,

$$
E[g(N(t+h))\mid N(t)=n]
=g(n)(1-\lambda h)+g(n+1)\lambda h+o(h),
$$

so

$$
\boxed{Ag(n)=\lambda[g(n+1)-g(n)].}
$$

For bounded $g$, the contribution of two or more jumps is bounded by a constant times $o(h)$, which justifies the remainder. The expression captures a pure-birth chain whose rate is $\lambda$ at every state, rather than proportional to the number already counted.

## Practice problems

Try these without the solutions first. State the event, its support, the relevant independence, and any parameterization assumptions before calculating.

1. **Units and one increment:** Events arrive at rate 12 per hour. Find the probability of at least two arrivals in five minutes.
2. **Separated intervals:** At rate 2 per hour, find the probability of no arrivals in $(0,1]$ and exactly three in $(1,3]$.
3. **Arrival threshold:** At rate 3 per hour, find $P(T_4\le1)$ and identify the distribution of $T_4$.
4. **Memorylessness:** If $X\sim\mathrm{Exp}(5)$, calculate $P(X>0.7\mid X>0.2)$ and $E[X\mid X>0.2]$.
5. **A new waiting origin:** Given $N(2)=4$ in an HPP of rate 3, find the conditional distribution and mean of $T_7-2$.
6. **Covariance:** For rate 2, find $\mathrm{Cov}(N(2),N(5))$ and $\mathrm{Corr}(N(2),N(5))$.
7. **Overlapping windows:** For rate 2, define $A=N(3)-N(1)$ and $B=N(4)-N(2)$. Find their covariance and correlation.
8. **Conditional density:** Derive the conditional density of $T_1$ given $T_3=t$, including its support and mean.
9. **Occupied windows:** Divide an hour into six disjoint 10-minute intervals for rate 6 per hour. Find the distribution of the number of nonempty intervals.
10. **A model counterexample:** Let $P(t)$ be an HPP of rate $\lambda$ and set $C(t)=2P(t)$. Which HPP properties does $C$ retain, and which does it violate?
11. **Martingale proof:** Derive $E[N(t)^2\mid\mathcal F_s]$ and use it to verify that $(N(t)-\lambda t)^2-\lambda t$ is a martingale.
12. **Local-to-global derivation:** Starting from the small-interval assumptions, derive the forward equation for $p_2(t)$ and solve it using $p_1(t)=\lambda te^{-\lambda t}$.

## Solutions and checks

### Count and waiting-time problems

- **Problem 1:** The mean five-minute count is $12(5/60)=1$, so the probability is $1-P(N=0)-P(N=1)=1-2e^{-1}\approx0.264241$.
- **Problem 2:** The interval means are 2 and 4, respectively. Independence gives $e^{-2}(e^{-4}4^3/3!)=(32/3)e^{-6}\approx0.026440$.
- **Problem 3:** $T_4\sim\mathrm{Gamma}(4,3)$ in shape–rate form. The probability is $1-e^{-3}(1+3+9/2+27/6)=1-13e^{-3}\approx0.352768$.
- **Problem 4:** The additional waiting duration is $0.7-0.2=0.5$, so the probability is $e^{-2.5}\approx0.082085$. The conditional total mean is $0.2+1/5=0.4$.
- **Problem 5:** Three more events are needed, so $T_7-2\mid\{N(2)=4\}\sim\mathrm{Gamma}(3,3)$. Its mean is 1 hour.
- **Problem 6:** The covariance is $2\min(2,5)=4$. The correlation is $4/\sqrt{4\cdot10}=\sqrt{2/5}\approx0.632456$.
- **Problem 7:** The overlap $(2,3]$ has length 1, so the covariance is 2. Each interval count has variance 4, giving correlation $2/4=1/2$.
- **Problem 9:** Each interval has Poisson mean 1, so its nonempty probability is $1-e^{-1}$. Independent increments give $\mathrm{Binomial}(6,1-e^{-1})$.

### Conditional-density and counterexample problems

For Problem 8, integrate the first-three-arrival joint density over $u<t_2<t$:

$$
f_{T_1,T_3}(u,t)
=\int_u^t\lambda^3e^{-\lambda t}\,dt_2
=\lambda^3e^{-\lambda t}(t-u),
\qquad 0<u<t.
$$

Since

$$
f_{T_3}(t)=\frac{\lambda^3t^2e^{-\lambda t}}2,
$$

the conditional density and mean are

$$
f_{T_1\mid T_3}(u\mid t)
=\frac{2(t-u)}{t^2}\mathbf1_{\{0<u<t\}},
\qquad
E[T_1\mid T_3=t]=\frac{t}{3}.
$$

The integral of the density is 1. The rate cancels, but the support still depends on the conditioning time.

For Problem 10, $C$ retains nondecreasing right-continuous paths and stationary independent increments. It has jumps of size 2, and

$$
P(C(h)\ge2)=P(P(h)\ge1)=\lambda h+o(h),
$$

which violates the HPP multiple-arrival condition. Its mean is $2\lambda t$ and its variance is $4\lambda t$, so even its one-time count law is not Poisson when $t>0$.

### Derivation problems

For Problem 11, set $a=\lambda(t-s)$ and write $N(t)=N(s)+Z$ with $Z\sim\mathrm{Poisson}(a)$ independent of $\mathcal F_s$. Since $E[Z^2]=a+a^2$,

$$
E[N(t)^2\mid\mathcal F_s]
=N(s)^2+2aN(s)+a+a^2.
$$

Expanding $(N(t)-\lambda t)^2$, substituting this expression, and simplifying gives

$$
E[(N(t)-\lambda t)^2-\lambda t\mid\mathcal F_s]
=(N(s)-\lambda s)^2-\lambda s.
$$

For Problem 12,

$$
p_2'(t)+\lambda p_2(t)=\lambda p_1(t)
=\lambda^2te^{-\lambda t}.
$$

Multiplication by $e^{\lambda t}$ followed by integration from 0 to $t$ gives

$$
p_2(t)=\frac{\lambda^2t^2}{2}e^{-\lambda t}.
$$

The initial value $p_2(0)=0$ is essential. Without it, the differential equation alone admits an extra integration constant.

## R verification and simulation

### How the R functions map to the mathematics

The examples below use base R only and explicitly set `rate=` for exponential and gamma functions. Run the blocks in order in a clean session; the deterministic checks compare independently calculated expressions rather than treating a simulation as proof.

| Mathematical operation | R expression |
|---|---|
| $P(N(t)=k)$ | `dpois(k, lambda = rate * t)` |
| $P(N(t)\le k)$ | `ppois(k, lambda = rate * t)` |
| $P(N(t)\ge k)$ for integer $k\ge1$ | `ppois(k - 1, rate * t, lower.tail = FALSE)` |
| $P(T_n>t)$ | `pgamma(t, shape = n, rate = rate, lower.tail = FALSE)` |
| $f_{T_n}(t)$ | `dgamma(t, shape = n, rate = rate)` |
| Simulate interarrival gaps | `rexp(n, rate = rate)` |

For a continuous waiting time, strict versus non-strict inequalities at a fixed finite endpoint have the same probability. For discrete counts, changing $\ge k$ to $>k$ changes the cutoff by one.

### Deterministic checks

```r
checks_run <- 0L
check <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(
    as.numeric(actual), as.numeric(expected), tolerance = tolerance
  )))
  checks_run <<- checks_run + 1L
}

# Six worked examples.
check(dpois(2, 2) * ppois(0, 1, lower.tail = FALSE),
      2 * exp(-2) * (1 - exp(-1)))
check(ppois(2, 5), 18.5 * exp(-5))
check(pgamma(2.5, shape = 3, rate = 2, lower.tail = FALSE),
      ppois(2, 5))
check(exp(-4 * (1/4 + 1/6)) / exp(-4 * 1/4), exp(-2/3))
check(1/4 + 1/4, 1/2)
check(dpois(2, 3) * dpois(3, 3), (81/4) * exp(-6))
check(3 * min(1, 2), 3)
check(integrate(function(u) 4 * exp(-2 * 1.5) + 0 * u,
                0, 1.5)$value,
      dgamma(1.5, shape = 2, rate = 2))
check(integrate(function(u) rep(1/1.5, length(u)), 0, 1.5)$value, 1)
check(integrate(function(u) u / 1.5, 0, 1.5)$value, 1.5/2)
check(30 * (1 - exp(-1/5)), 30 * ppois(0, 1/5, lower.tail = FALSE))
check(30 * (1 - exp(-1/5)) * exp(-1/5),
      30 * ppois(0, 1/5, lower.tail = FALSE) * dpois(0, 1/5))

# Poisson/gamma duality for several parameters and thresholds.
for (rate in c(0.5, 2, 7)) {
  for (n in c(1L, 3L, 6L)) {
    t <- 1.2
    check(pgamma(t, shape = n, rate = rate, lower.tail = FALSE),
          ppois(n - 1, lambda = rate * t))
  }
}

# Normalization, gamma moments, and conditional-density support.
check(sum(dpois(0:100, 6)), 1)
check(integrate(function(t) dgamma(t, shape = 3, rate = 2),
                0, Inf)$value, 1)
check(integrate(function(t) t * dgamma(t, shape = 3, rate = 2),
                0, Inf)$value, 3/2)
check(integrate(function(t) (t - 3/2)^2 *
                  dgamma(t, shape = 3, rate = 2), 0, Inf)$value, 3/4)
t_fixed <- 2.4
check(integrate(function(u) 2 * (t_fixed - u) / t_fixed^2,
                0, t_fixed)$value, 1)
check(integrate(function(u) u * 2 * (t_fixed - u) / t_fixed^2,
                0, t_fixed)$value, t_fixed/3)

# A small-interval calculation with a stable R upper-tail function.
rate <- 3
h <- 1e-4
exact_multiple <- ppois(1, lambda = rate * h, lower.tail = FALSE)
check(exact_multiple / ((rate * h)^2 / 2), 1, tolerance = 3e-4)

cat("Passed", checks_run, "deterministic checks.\n")
```

The finite PMF sum above truncates an infinite tail, but the omitted mass is below the stated tolerance for that parameter. The density integrations only evaluate points on their declared support.

### Exact path simulation on a finite horizon

Generating exponential gaps until the next arrival passes the horizon avoids guessing a fixed number of gaps. It also handles a horizon with no arrivals by returning an empty vector.

```r
simulate_hpp <- function(rate, horizon) {
  stopifnot(length(rate) == 1L, is.finite(rate), rate > 0,
            length(horizon) == 1L, is.finite(horizon), horizon >= 0)
  if (horizon == 0) return(numeric(0))
  times <- numeric(0)
  current <- 0
  repeat {
    current <- current + rexp(1, rate = rate)
    if (current > horizon) break
    times <- c(times, current)
  }
  times
}

set.seed(544)
arrivals <- simulate_hpp(rate = 3, horizon = 4)
print(arrivals)

# With nondecreasing arrival times, findInterval counts arrivals <= each time.
grid <- seq(0, 4, length.out = 401)
counts <- findInterval(grid, arrivals)
stopifnot(all(diff(counts) >= 0), counts[1] == 0,
          tail(counts, 1) == length(arrivals))

# Optional plot, not needed for the checks.
# plot(stepfun(arrivals, 0:length(arrivals)), xlim = c(0, 4),
#      do.points = FALSE, verticals = TRUE,
#      xlab = "Time (hours)", ylab = "N(t)")
```

Repeated vector concatenation is adequate for a teaching example with a small expected event count. For very large simulations, use chunking or the count-plus-uniform construction introduced in Part II.

### Monte Carlo checks with uncertainty

The next block simulates disjoint increments, not independent cumulative counts. Monte Carlo standard errors describe uncertainty in the empirical mean and event-probability estimates; the sample variances and covariance are descriptive comparisons.

```r
set.seed(1544)
B <- 200000L
rate <- 3
A <- rpois(B, lambda = rate * 1)
D <- rpois(B, lambda = rate * 1)
N1 <- A
N2 <- A + D

report_mean <- function(x, target) {
  c(estimate = mean(x), target = target,
    Monte_Carlo_SE = sd(x) / sqrt(length(x)))
}

print(report_mean(N1, 3))
print(report_mean(N2, 6))
print(c(variance_N1 = var(N1), target = 3,
        covariance_N1_N2 = cov(N1, N2), covariance_target = 3))

# T3 as the sum of three independent exponential gaps.
T3 <- rowSums(matrix(rexp(3L * B, rate = 2), nrow = B))
indicator <- as.numeric(T3 > 2.5)
print(report_mean(T3, 3/2))
print(report_mean(indicator, ppois(2, 5)))
print(c(variance_T3 = var(T3), variance_target = 3/4))
```

These simulation estimates should be close to their targets relative to Monte Carlo uncertainty, but they are not required to equal them exactly. A correct theoretical result does not depend on a particular random seed.

### Verification record

All three R blocks were extracted from this Markdown file and run together in a clean R 4.5.2 session; all 28 deterministic numerical checks passed, and the path and Monte Carlo examples completed without errors. Separately, all 114 display equations and 293 inline math expressions passed a strict KaTeX syntax check, and the Markdown parser recognized the same 407 math expressions.

The document uses `$...$` for inline mathematics and separate-line `$$...$$` display blocks, with `\mathrm` for named distributions and operators. These are syntax and implementation checks, not a claim that this file was published to GitHub or tested inside GitHub's own renderer.

## Common errors and lecture preparation

### Errors to catch before submitting

- **Rate versus mean count:** Use $\lambda t$, not $\lambda$, as the Poisson parameter for an interval of length $t$.
- **Rate versus scale:** In a shape–rate gamma law, the mean is $n/\lambda$; using $\lambda$ as a scale reverses it.
- **Cumulative-count independence:** $N(s)$ and $N(t)$ share arrivals. The independent pair is $N(s)$ and $N(t)-N(s)$.
- **Wrong threshold:** $T_n>t$ corresponds to $N(t)\le n-1$, not $N(t)\le n$.
- **Wrong clock:** Distinguish a future wait from an absolute arrival time.
- **Overextended memorylessness:** Individual exponential gaps are memoryless; $T_n$ is not exponential when $n>1$.
- **Forgotten support:** Arrival times must be ordered; a density formula without its support is incomplete.
- **Random-time substitution:** Independent increments at deterministic times do not by themselves prove restart at arbitrary random times.
- **Unverified model claim:** A fitted Poisson histogram is not evidence by itself for independent stationary increments.
- **Simulation as proof:** Numerical agreement checks a calculation or implementation, not a mathematical theorem.

### Questions to bring to lecture

- **Definition choice:** Does the course start from exponential gaps, Poisson increments, or infinitesimal probabilities?
- **Proof expectations:** Are forward equations and equivalence of definitions required, or only their use?
- **Notation:** Does the instructor use $T_n$ for gaps or epochs, and rate or scale for gamma distributions?
- **Conditioning:** Will questions condition on a fixed total, exact arrival times, or the full past history?
- **Advanced scope:** Are filtrations, martingales, and strong Markov arguments examinable or enrichment?

## Compact reference sheet

Assume $N$ is a homogeneous Poisson process with rate $\lambda>0$. All times below use the same unit as that rate.

| Object | Distribution or identity | Mean | Variance |
|---|---|---|---|
| $N(t)$ | $\mathrm{Poisson}(\lambda t)$ | $\lambda t$ | $\lambda t$ |
| $N(t)-N(s)$, $s<t$ | $\mathrm{Poisson}(\lambda(t-s))$ | $\lambda(t-s)$ | $\lambda(t-s)$ |
| $X_n$ | $\mathrm{Exp}(\lambda)$, iid | $1/\lambda$ | $1/\lambda^2$ |
| $T_n$ | $\mathrm{Gamma}(n,\lambda)$, shape–rate | $n/\lambda$ | $n/\lambda^2$ |
| $W_s$, deterministic $s$ | $\mathrm{Exp}(\lambda)$, independent of past | $1/\lambda$ | $1/\lambda^2$ |

$$
P(N(t)=k)=e^{-\lambda t}\frac{(\lambda t)^k}{k!},
\qquad k=0,1,\ldots.
$$

$$
\{T_n\le t\}=\{N(t)\ge n\},
\qquad
P(T_n>t)=e^{-\lambda t}
\sum_{k=0}^{n-1}\frac{(\lambda t)^k}{k!}.
$$

$$
f_{T_1,\ldots,T_n}(t_1,\ldots,t_n)
=\lambda^ne^{-\lambda t_n},
\qquad 0<t_1<\cdots<t_n.
$$

$$
\mathrm{Cov}(N(s),N(t))=\lambda\min(s,t),
\qquad
\mathrm{Cov}(T_i,T_j)=\frac{\min(i,j)}{\lambda^2}.
$$

$$
E[N(t)\mid\mathcal F_s]=N(s)+\lambda(t-s),
\qquad
\mathrm{Var}(N(t)\mid\mathcal F_s)=\lambda(t-s).
$$

$$
N(t)-\lambda t\text{ is a martingale},
\qquad
Ag(n)=\lambda[g(n+1)-g(n)].
$$

Before choosing a formula, ask: **Am I counting events in a fixed interval, or waiting for a fixed number of events?** Then identify exactly what has been observed or conditioned on.
