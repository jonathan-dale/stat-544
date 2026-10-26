# STAT 544 Study Guide and Lecture Prep: Poisson Process II

Part I connected Poisson counts, exponential interarrival times, and gamma arrival epochs. Part II develops the structural tools that make the process useful: **condition on a total → allocate events across time → describe ordered arrival times → split or combine streams → extend the rate model**.

This guide continues the explicit joint-to-marginal-to-conditional calculations in your [Conditional Distributions I](https://github.com/jonathan-dale/stat-544/blob/main/guides/conditional_distributions/cd_1/Conditional_Distributions_I.md) and [Conditional Distributions II](https://github.com/jonathan-dale/stat-544/blob/main/guides/conditional_distributions/cd_2/Conditional_Distributions_II.md) guides. The precise lecture notes for “Poisson process II” were not supplied, so the core/extension labels below are proposed study priorities rather than a claim about the instructor’s syllabus.

## How to use this guide

Start by deriving the conditional binomial distribution and the conditional joint density of arrival times. Those two calculations explain much of the rest of the guide, including why conditioning changes dependence and why the unknown constant rate disappears from the timing distribution once the total count is fixed.

- **Core Part II skills:** Conditional binomial and multinomial counts, uniform order statistics, beta arrival times, independent thinning, superposition, and competing streams.
- **Graduate derivations:** Joint-density normalization, conditional covariance, bounded-window spacings, conditional versus unconditional independence, and process-level independence.
- **Extensions to confirm against lecture:** Nonhomogeneous processes, time change, compound Poisson sums, fixed-window likelihood, and random-rate mixtures.
- **Practice structure:** Seven worked problems, additional exercises with solutions, reproducible R checks, and a compact formula sheet.

### Contents

- [Prerequisite bridge from Part I](#prerequisite-bridge-from-part-i)
- [Conditioning on a fixed total count](#conditioning-on-a-fixed-total-count)
- [Conditional arrival times and order statistics](#conditional-arrival-times-and-order-statistics)
- [Thinning and marking](#thinning-and-marking)
- [Superposition and competing streams](#superposition-and-competing-streams)
- [Nonhomogeneous Poisson processes](#nonhomogeneous-poisson-processes)
- [Compound Poisson sums](#compound-poisson-sums)
- [Statistical extensions](#statistical-extensions)
- [Seven fully worked problems](#seven-fully-worked-problems)
- [Practice problems](#practice-problems)
- [Solutions and checks](#solutions-and-checks)
- [R verification and simulation](#r-verification-and-simulation)
- [Compact reference sheet](#compact-reference-sheet)

## Prerequisite bridge from Part I

### Keep the three random objects distinct

Let $N$ be a homogeneous Poisson process of rate $\lambda>0$. We retain the notation $T_n$ for the $n$th arrival time and $X_n=T_n-T_{n-1}$ for the $n$th gap, with $T_0=0$.

The basic dictionary is

$$
N(t)-N(s)\sim\mathrm{Poisson}(\lambda(t-s)),
\qquad 0\le s<t,
$$

$$
X_1,X_2,\ldots\ \text{are iid }\mathrm{Exp}(\lambda),
\qquad
T_n\sim\mathrm{Gamma}(n,\lambda).
$$

Gamma parameters are shape and rate. The event equivalence

$$
\{T_n\le t\}=\{N(t)\ge n\}
$$

translates between timing and counting questions, as developed in [Gallager’s basic Poisson-process theory](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.02:_Definition_and_Properties_of_a_Poisson_Process).

### The conditional-distribution workflow still applies

For discrete $A,B$ with $P(B=b)>0$,

$$
P(A=a\mid B=b)
=\frac{P(A=a,B=b)}{P(B=b)}.
$$

The support comes first. If the total count in a window is $n$, every subwindow count is restricted to $\{0,\ldots,n\}$.

For density-based conditioning, divide the appropriate joint density by the conditioning density. For a continuous arrival time, a statement such as $T_k=t$ is handled by a conditional density version, not by dividing by $P(T_k=t)=0$.

We will repeatedly use

$$
E[Y]=E[E[Y\mid Z]],
$$

$$
\mathrm{Var}(Y)
=E[\mathrm{Var}(Y\mid Z)]
+\mathrm{Var}(E[Y\mid Z]),
$$

and the covariance counterpart

$$
\mathrm{Cov}(Y,W)
=E[\mathrm{Cov}(Y,W\mid Z)]
+\mathrm{Cov}(E[Y\mid Z],E[W\mid Z]).
$$

The covariance identity follows by conditioning $E[YW]$ and subtracting $E[Y]E[W]$. It will explain how negative conditional dependence can average to unconditional independence.

## Conditioning on a fixed total count

### Deriving the conditional binomial distribution

Fix $0<s<t$ and suppose $N(t)=n$. Write

$$
A=N(s),\qquad B=N(t)-N(s).
$$

Unconditionally, $A$ and $B$ are independent Poisson variables with means $\lambda s$ and $\lambda(t-s)$. For $k=0,\ldots,n$,

$$
\begin{aligned}
P(N(s)=k\mid N(t)=n)
&=\frac{P(A=k,B=n-k)}{P(A+B=n)}\\
&=\frac{
e^{-\lambda s}(\lambda s)^k/k!\;
e^{-\lambda(t-s)}[\lambda(t-s)]^{n-k}/(n-k)!
}{
e^{-\lambda t}(\lambda t)^n/n!
}\\
&=\binom nk
\left(\frac st\right)^k
\left(1-\frac st\right)^{n-k}.
\end{aligned}
$$

Therefore

$$
\boxed{
N(s)\mid\{N(t)=n\}
\sim\mathrm{Binomial}\left(n,\frac st\right).
}
$$

This is the fixed-total allocation principle underlying [Gallager’s conditional-arrival analysis](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.05:_Conditional_Arrival_Densities_and_Order_Statistics). Its moments are

$$
E[N(s)\mid N(t)=n]=n\frac st,
$$

$$
\mathrm{Var}(N(s)\mid N(t)=n)
=n\frac st\left(1-\frac st\right).
$$

The constant rate $\lambda$ cancels. Given the total number of events, their relative allocation over a homogeneous time window does not depend on how large that constant rate was.

### Conditioning backward is not predicting forward

Compare two questions for $s<t$:

$$
N(t)\mid\{N(s)=m\}
\ \overset d=\ m+\mathrm{Poisson}(\lambda(t-s)),
$$

$$
N(s)\mid\{N(t)=n\}
\sim\mathrm{Binomial}(n,s/t).
$$

The first adds independent future events to a known past count. The second allocates a fixed final total between an earlier interval and its complement.

The supports are different: the future cumulative count is at least $m$, while the earlier count cannot exceed the specified final total $n$. Checking the support often reveals the correct distribution immediately.

### Several subintervals give a multinomial law

Partition $(0,t]$ at $0=t_0<t_1<\cdots<t_m=t$ and define

$$
C_j=N(t_j)-N(t_{j-1}),
\qquad
p_j=\frac{t_j-t_{j-1}}t.
$$

For nonnegative integers $c_1,\ldots,c_m$ summing to $n$, division of the product Poisson PMF by the total-count PMF gives

$$
\boxed{
P(C_1=c_1,\ldots,C_m=c_m\mid N(t)=n)
=\frac{n!}{\prod_{j=1}^m c_j!}\prod_{j=1}^m p_j^{c_j}.
}
$$

Thus

$$
(C_1,\ldots,C_m)\mid\{N(t)=n\}
\sim\mathrm{Multinomial}(n;p_1,\ldots,p_m).
$$

Representing the count vector as the sum of $n$ independent one-hot allocation vectors yields

$$
E[C_j\mid N(t)=n]=np_j,
$$

$$
\mathrm{Var}(C_j\mid N(t)=n)=np_j(1-p_j),
$$

$$
\mathrm{Cov}(C_i,C_j\mid N(t)=n)
=-np_ip_j,\qquad i\ne j.
$$

The negative covariance reflects the fixed total: assigning an event to one interval removes it from every other interval. Independent increments is an unconditional property, not an assertion that survives every conditioning operation.

### Recovering unconditional independence

Let $K=N(t)\sim\mathrm{Poisson}(\lambda t)$. For two distinct partition cells,

$$
\begin{aligned}
\mathrm{Cov}(C_i,C_j)
&=E[-Kp_ip_j]+\mathrm{Cov}(Kp_i,Kp_j)\\
&=-\lambda t\,p_ip_j+p_ip_j\mathrm{Var}(K)\\
&=0.
\end{aligned}
$$

This cancellation agrees with their unconditional independence. Zero covariance alone would not prove independence for arbitrary variables; here the independent-increments definition gives the stronger result.

## Conditional arrival times and order statistics

### Deriving the conditional joint density

Suppose $N(t)=n\ge1$. To get the joint density of $(T_1,\ldots,T_n)$ together with this event, use the unconditional first-$n$-arrival density and require that the next gap extends beyond $t$.

On $0<u_1<\cdots<u_n<t$,

$$
f_{T_1,\ldots,T_n}(u_1,\ldots,u_n)
=\lambda^ne^{-\lambda u_n}.
$$

The independent next gap satisfies

$$
P(X_{n+1}>t-u_n)=e^{-\lambda(t-u_n)}.
$$

Therefore the mixed density–mass for those arrival times and $N(t)=n$ is

$$
\lambda^ne^{-\lambda u_n}e^{-\lambda(t-u_n)}
=\lambda^ne^{-\lambda t}.
$$

Divide by $P(N(t)=n)=e^{-\lambda t}(\lambda t)^n/n!$:

$$
\boxed{
f_{T_1,\ldots,T_n\mid N(t)}
(u_1,\ldots,u_n\mid n)
=\frac{n!}{t^n},
\quad 0<u_1<\cdots<u_n<t.
}
$$

This is the conditional density theorem in [Gallager’s treatment of arrival densities and order statistics](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.05:_Conditional_Arrival_Densities_and_Order_Statistics). Its value is zero outside the ordered region.

The ordered simplex has volume $t^n/n!$. Multiplying that volume by $n!/t^n$ verifies normalization.

The no-more-arrivals factor is essential. Dividing $\lambda^ne^{-\lambda u_n}$ alone by $P(N(t)=n)$ would not give the correct conditional law.

### Uniform points versus ordered uniform points

Let $U_1,\ldots,U_n$ be iid $\mathrm{Uniform}(0,t)$, and write their order statistics as $U_{(1)}<\cdots<U_{(n)}$. Sorting folds $n!$ equally likely orderings onto the same ordered region, giving density $n!/t^n$ there.

Hence

$$
(T_1,\ldots,T_n)\mid\{N(t)=n\}
\ \overset d=\
(U_{(1)},\ldots,U_{(n)}).
$$

The statement is about ordered points, not independent ordered arrival times. The chronological order forces dependence.

- **Unsorted representation:** $U_1,\ldots,U_n$ are independent uniforms before sorting.
- **Ordered representation:** $T_1,\ldots,T_n$ have the order-statistic joint law.
- **Randomly selected event:** A uniformly selected one of the $n$ events has a uniform location.
- **First event:** The minimum is not uniform unless $n=1$.

For $n=0$, there are no arrival coordinates in the window. One must not apply a formula for $T_k$ with $k>n$.

### Individual arrival times have scaled beta laws

For $1\le k\le n$ and $0<u<t$, integrate the constant joint density over the coordinates before and after $u$:

$$
\begin{aligned}
f_{T_k\mid N(t)}(u\mid n)
&=\frac{n!}{t^n}
\frac{u^{k-1}}{(k-1)!}
\frac{(t-u)^{n-k}}{(n-k)!}.
\end{aligned}
$$

Equivalently,

$$
\boxed{
\frac{T_k}{t}\mid\{N(t)=n\}
\sim\mathrm{Beta}(k,n+1-k).
}
$$

Consequently,

$$
E[T_k\mid N(t)=n]=\frac{kt}{n+1},
$$

$$
\mathrm{Var}(T_k\mid N(t)=n)
=\frac{k(n+1-k)t^2}{(n+1)^2(n+2)}.
$$

The denominator is $n+1$, not $n$. There are $n+1$ gaps when the two window endpoints are included.

Useful special cases are

$$
P(T_1>u\mid N(t)=n)=\left(1-\frac ut\right)^n,
$$

$$
P(T_n\le u\mid N(t)=n)=\left(\frac ut\right)^n,
$$

for $0\le u\le t$. More generally, count–arrival duality gives

$$
P(T_k\le u\mid N(t)=n)
=P\!\left(\mathrm{Binomial}(n,u/t)\ge k\right).
$$

### Bounded-window spacings

Given $N(t)=n\ge1$, define the $n+1$ window spacings

$$
D_1=T_1,\qquad
D_j=T_j-T_{j-1}\ (2\le j\le n),
\qquad
D_{n+1}=t-T_n.
$$

They sum exactly to $t$. The last spacing is the unfinished interval from the final observed arrival to the window endpoint, not the full next interarrival time.

The unit-Jacobian transformation from ordered arrivals gives constant density $n!/t^n$ for the first $n$ spacings on

$$
d_1,\ldots,d_n>0,\qquad
\sum_{j=1}^n d_j<t,
$$

with $D_{n+1}=t-\sum_{j=1}^nD_j$. In standard notation,

$$
\left(\frac{D_1}t,\ldots,\frac{D_{n+1}}t\right)
\mid\{N(t)=n\}
\sim\mathrm{Dirichlet}(1,\ldots,1).
$$

This is the simplex form of the conditional-gap symmetry described in [Gallager’s spacing calculation](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.05:_Conditional_Arrival_Densities_and_Order_Statistics). The full $n+1$-vector lies on a hyperplane and does not have an ordinary $(n+1)$-dimensional Lebesgue density.

Integrating the simplex density gives, for distinct $i,j$,

$$
E[D_i\mid N(t)=n]=\frac{t}{n+1},
$$

$$
\mathrm{Var}(D_i\mid N(t)=n)
=\frac{nt^2}{(n+1)^2(n+2)},
$$

$$
\mathrm{Cov}(D_i,D_j\mid N(t)=n)
=-\frac{t^2}{(n+1)^2(n+2)}.
$$

The spacings are exchangeable but not independent. Their fixed sum makes independence impossible in the nondegenerate case.

Summing the appropriate spacing covariances also gives, for $1\le i\le j\le n$,

$$
\mathrm{Cov}(T_i,T_j\mid N(t)=n)
=\frac{i(n+1-j)t^2}{(n+1)^2(n+2)}.
$$

The covariance of ordered arrival times is positive even though the distinct spacings have negative covariance. These are different random quantities.

### Conditioning on an exact arrival epoch is different

In the natural conditional-density version, the first $n-1$ arrival times given $T_n=t$ are the ordered values of $n-1$ iid uniforms on $(0,t)$. To see it, divide the joint first-$n$ density by the gamma density of $T_n$:

$$
\frac{\lambda^ne^{-\lambda t}}
{\lambda^nt^{n-1}e^{-\lambda t}/(n-1)!}
=\frac{(n-1)!}{t^{n-1}}.
$$

Thus, for $1\le k<n$,

$$
E[T_k\mid T_n=t]=\frac{kt}{n}.
$$

Compare

$$
E[T_k\mid N(t)=n]=\frac{kt}{n+1}.
$$

In the first problem the $n$th arrival is fixed at the endpoint. In the second, all $n$ arrivals occur somewhere before the endpoint, and none is fixed there.

## Thinning and marking

### Independent Bernoulli thinning

Let $N$ be an HPP of rate $\lambda$. Independently for each event, retain it with probability $p\in(0,1)$ and reject it with probability $q=1-p$, with the labels also independent of all arrival times.

Write $R(t)$ and $D(t)$ for retained and discarded counts. Given $N(t)=n$,

$$
R(t)\mid\{N(t)=n\}\sim\mathrm{Binomial}(n,p),
\qquad
D(t)=n-R(t).
$$

For $r,d\ge0$,

$$
\begin{aligned}
P(R(t)=r,D(t)=d)
&=P(N(t)=r+d)\binom{r+d}{r}p^rq^d\\
&=e^{-\lambda t}
\frac{(\lambda t)^{r+d}}{r!d!}p^rq^d\\
&=\left[e^{-\lambda pt}\frac{(\lambda pt)^r}{r!}\right]
\left[e^{-\lambda qt}\frac{(\lambda qt)^d}{d!}\right].
\end{aligned}
$$

This factors into independent Poisson PMFs. Applying the same calculation separately to the cells of any finite time partition, then using independence across disjoint cells, establishes the process-level thinning theorem ([Gallager’s splitting proof](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.03:_Combining_and_Splitting_Poisson_Processes)).

Therefore

$$
\boxed{
R\sim\mathrm{HPP}(\lambda p),\qquad
D\sim\mathrm{HPP}(\lambda q),
\qquad R\perp D\text{ as processes}.
}
$$

Independence at a single time would not, by itself, establish independence of entire processes. The finite-partition argument supplies that additional step.

### Why conditioning changes the independence conclusion

Unconditionally, $R(t)$ and $D(t)$ are independent. Given $N(t)=n$, their sum is fixed, so

$$
\mathrm{Cov}(R(t),D(t)\mid N(t)=n)
=-npq.
$$

Given the full original event path, the retained and discarded paths also determine each other. There is no contradiction: the independence statement refers to the unconditional marked-process model.

The covariance decomposition makes the cancellation explicit:

$$
\begin{aligned}
\mathrm{Cov}(R(t),D(t))
&=E[-N(t)pq]
+\mathrm{Cov}(pN(t),qN(t))\\
&=-\lambda tpq+\lambda tpq=0.
\end{aligned}
$$

Again, the PMF factorization proves independence; the zero-covariance calculation is a check. It would not be sufficient on its own.

### Reverse conditioning and the link to Conditional Distributions I

Since $D(t)$ is independent of $R(t)$,

$$
\boxed{
N(t)\mid\{R(t)=r\}
\ \overset d=\ r+\mathrm{Poisson}(\lambda qt).
}
$$

The conditional mean and variance are

$$
E[N(t)\mid R(t)=r]=r+\lambda qt,
\qquad
\mathrm{Var}(N(t)\mid R(t)=r)=\lambda qt.
$$

This is the process version of the hierarchical Poisson-thinning example in your [Conditional Distributions I guide](https://github.com/jonathan-dale/stat-544/blob/main/guides/conditional_distributions/cd_1/Conditional_Distributions_I.md). The unseen rejected count is Poisson; the total is a shifted Poisson variable.

### Several marks and non-independent selection rules

With iid labels taking values $1,\ldots,m$ with probabilities $p_1,\ldots,p_m$, independent of the event times, the marked streams are independent HPPs with rates $\lambda p_j$ ([Gallager’s splitting construction](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.03:_Combining_and_Splitting_Poisson_Processes)). Conditional on the total, the mark counts instead have a multinomial law.

Keeping every second event is not independent thinning. Starting with retained events $T_2,T_4,\ldots$, successive retained gaps are

$$
X_1+X_2,\quad X_3+X_4,\quad\ldots,
$$

which are iid $\mathrm{Gamma}(2,\lambda)$, not exponential. The resulting renewal process is not an HPP, even though its long-run rate is $\lambda/2$.

## Superposition and competing streams

### Combining independent processes

Let $N_1,\ldots,N_m$ be independent HPPs with rates $\lambda_1,\ldots,\lambda_m$. Define

$$
N(t)=\sum_{j=1}^mN_j(t),
\qquad
\lambda_\ast=\sum_{j=1}^m\lambda_j.
$$

For an interval of length $h$, the PGF of the summed increment is

$$
\prod_{j=1}^m\exp\{\lambda_jh(z-1)\}
=\exp\{\lambda_\ast h(z-1)\}.
$$

This is a Poisson increment with mean $\lambda_\ast h$. Independent increments persist across disjoint intervals, so $N$ is an HPP with rate $\lambda_\ast$, as established in [Gallager’s superposition theorem](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.03:_Combining_and_Splitting_Poisson_Processes).

The independence assumption matters. Adding a process to itself doubles its jumps, rather than producing a new HPP with twice the rate.

### Which stream produces the next event?

Let $W_j\sim\mathrm{Exp}(\lambda_j)$ be independent waits for the next events in the separate streams. Their minimum satisfies

$$
P(\min_jW_j>u)
=\prod_je^{-\lambda_ju}
=e^{-\lambda_\ast u}.
$$

Thus the next combined wait is $\mathrm{Exp}(\lambda_\ast)$. The probability that stream $j$ wins is

$$
\begin{aligned}
P(W_j=\min_iW_i)
&=\int_0^\infty
\lambda_je^{-\lambda_ju}
\prod_{i\ne j}e^{-\lambda_iu}\,du\\
&=\frac{\lambda_j}{\lambda_\ast}.
\end{aligned}
$$

The joint density of waiting time $u$ and winning label $j$ is

$$
\lambda_je^{-\lambda_\ast u}
=\frac{\lambda_j}{\lambda_\ast}
\lambda_\ast e^{-\lambda_\ast u}.
$$

Hence the winning label and the combined waiting time are independent. By repeated memoryless restarting, combined event labels are iid with probabilities $\lambda_j/\lambda_\ast$.

### More general arrival races

For two independent streams A and B, let

$$
p=\frac{\lambda_A}{\lambda_A+\lambda_B}.
$$

The event that A’s $r$th arrival precedes B’s $s$th is equivalent to at least $r$ A-labels among the first $r+s-1$ combined events. Therefore

$$
\boxed{
P(T_r^A<T_s^B)
=\sum_{k=r}^{r+s-1}
\binom{r+s-1}{k}p^k(1-p)^{r+s-1-k}.
}
$$

This event-label method avoids integrating two gamma densities, and is developed in [Gallager’s examples with independent processes](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.03:_Combining_and_Splitting_Poisson_Processes). Ties have probability zero for these continuous arrival times.

## Nonhomogeneous Poisson processes

### Deterministic rate functions and integrated intensity

This section is an extension unless your lecture explicitly includes time-varying rates. Let $\lambda(u)\ge0$ be deterministic, measurable, and locally integrable, and define

$$
\Lambda(t)=\int_0^t\lambda(u)\,du.
$$

A nonhomogeneous Poisson process, or NHPP, has independent increments and

$$
\boxed{
N(t)-N(s)
\sim\mathrm{Poisson}(\Lambda(t)-\Lambda(s)),
\qquad 0\le s<t.
}
$$

This integrated-intensity definition is given in [Siegrist’s nonhomogeneous Poisson-process treatment](https://www.randomservices.org/random/poisson/Nonhomogeneous.html). Stationary increments generally disappear because equal-length intervals can contain different integrated intensities.

The basic moments become

$$
E[N(t)]=\mathrm{Var}(N(t))=\Lambda(t),
$$

$$
\mathrm{Cov}(N(s),N(t))=\Lambda(\min(s,t)).
$$

The parameter for an interval is its integrated rate, not its endpoint rate multiplied by the duration unless that rate is constant there. The compensated process $N(t)-\Lambda(t)$ is a martingale by the same conditional-increment calculation as in Part I.

### Waiting times and arrival densities

For deterministic $s\ge0$, let $W_s$ be the wait for the next event. Independent increments yield

$$
P(W_s>u\mid\mathcal F_s)
=\exp\{-[\Lambda(s+u)-\Lambda(s)]\},
\qquad u\ge0.
$$

The residual wait generally depends on the clock time $s$, so it is not an exponential distribution with a constant rate. Under suitable continuity, the finite-time first-arrival density is

$$
f_{T_1}(t)=\lambda(t)e^{-\Lambda(t)}.
$$

By count–arrival duality and differentiation,

$$
f_{T_n}(t)
=\lambda(t)e^{-\Lambda(t)}
\frac{\Lambda(t)^{n-1}}{(n-1)!},
\qquad t>0,
$$

at points where the derivative is defined. These are the arrival-density formulas in [Siegrist’s NHPP notes](https://www.randomservices.org/random/poisson/Nonhomogeneous.html).

If $\Lambda(\infty)<\infty$, an arrival may never occur. The displayed finite-time density then integrates to less than 1, with the remaining probability placed at $T_n=\infty$; for example,

$$
P(T_1=\infty)=e^{-\Lambda(\infty)}.
$$

### Conditioning in a nonhomogeneous window

Suppose $0<\Lambda(t)<\infty$. For $N(t)=n$, the joint ordered-arrival density is

$$
\boxed{
f_{T_1,\ldots,T_n\mid N(t)}
(u_1,\ldots,u_n\mid n)
=n!\prod_{j=1}^n\frac{\lambda(u_j)}{\Lambda(t)},
\quad 0<u_1<\cdots<u_n<t.
}
$$

One derivation uses the mixed density–mass $e^{-\Lambda(t)}\prod_j\lambda(u_j)$ and divides by $P(N(t)=n)$. This is the time-varying analogue of the homogeneous joint-to-conditional calculation.

Thus the unordered points have iid density and CDF

$$
g_t(u)=\frac{\lambda(u)}{\Lambda(t)}\mathbf1_{\{0<u<t\}},
\qquad
G_t(u)=\frac{\Lambda(u)}{\Lambda(t)}
\quad(0\le u\le t).
$$

They are uniform in ordinary time only when the intensity is constant almost everywhere on the window. The conditional count allocation becomes

$$
N(s)\mid\{N(t)=n\}
\sim\mathrm{Binomial}
\left(n,\frac{\Lambda(s)}{\Lambda(t)}\right).
$$

In particular,

$$
\frac{\Lambda(T_k)}{\Lambda(t)}
\mid\{N(t)=n\}
\sim\mathrm{Beta}(k,n+1-k).
$$

The beta law applies to transformed time. It does not imply $E[T_k\mid N(t)=n]$ is obtained by applying $\Lambda^{-1}$ to the transformed-time mean.

### Time change and time-dependent thinning

If $M$ is an HPP of rate 1, then

$$
N(t)=M(\Lambda(t))
$$

constructs an NHPP with integrated intensity $\Lambda$ ([Gallager’s nonlinear-time construction](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.04:_Non-homogeneous_Poisson_Processes)). When $\Lambda$ is continuous, strictly increasing, and unbounded,

$$
\Lambda(T_n)\sim\mathrm{Gamma}(n,1),
$$

and transformed gaps are iid $\mathrm{Exp}(1)$. With flat portions one uses an appropriate generalized inverse; unboundedness avoids truncating the transformed process at a finite intensity horizon.

For independent retention with deterministic probability $p(u)$ at time $u$, retained and rejected streams are independent NHPPs with rates

$$
\lambda_R(u)=\lambda(u)p(u),
\qquad
\lambda_D(u)=\lambda(u)(1-p(u)),
$$

as illustrated by [Gallager’s nonhomogeneous splitting construction](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Discrete_Stochastic_Processes_(Gallager)/02:_Poisson_Processes/2.04:_Non-homogeneous_Poisson_Processes). This requires conditionally independent marking given event times, not a rule based on event rank or past retention outcomes.

## Compound Poisson sums

### A random number of contributions

This section connects process models to total expectation and total variance. Let $N$ be an HPP of rate $\lambda$, and let $Y_1,Y_2,\ldots$ be iid marks, independent of the entire count process, with

$$
E[Y_1]=m,\qquad \mathrm{Var}(Y_1)=v<\infty.
$$

Define

$$
S(t)=\sum_{i=1}^{N(t)}Y_i,
$$

with an empty sum equal to zero. This is a compound Poisson process, the random-sum model presented in [Siegrist’s compound Poisson notes](https://www.randomservices.org/random/poisson/Compound.html).

Given $N(t)=n$,

$$
E[S(t)\mid N(t)=n]=nm,
\qquad
\mathrm{Var}(S(t)\mid N(t)=n)=nv.
$$

Therefore

$$
\boxed{E[S(t)]=\lambda tm,}
$$

$$
\begin{aligned}
\mathrm{Var}(S(t))
&=E[N(t)v]+\mathrm{Var}(N(t)m)\\
&=\lambda tv+\lambda tm^2\\
&=\boxed{\lambda tE[Y_1^2]}.
\end{aligned}
$$

The extra $m^2$ term accounts for randomness in the number of contributions. Using only $\lambda tv$ would ignore this source of variation.

If $M_Y(\theta)=E[e^{\theta Y_1}]$ is finite at the chosen argument, conditioning also gives

$$
M_{S(t)}(\theta)
=E[M_Y(\theta)^{N(t)}]
=\exp\{\lambda t(M_Y(\theta)-1)\}.
$$

Unless the marks have a special form, $S(t)$ is not Poisson-distributed. If $Y_i\equiv1$, the original count process is recovered; if $Y_i\equiv2$, it has jumps of size 2 instead.

## Statistical extensions

### Fixed-window likelihood and the constant-rate MLE

Suppose an HPP is observed over a fixed, nonrandom interval $(0,t]$, with $t>0$ and $N(t)=n$. The count likelihood is

$$
L(\lambda;n)
=e^{-\lambda t}\frac{(\lambda t)^n}{n!}.
$$

Ignoring factors constant in $\lambda$,

$$
\ell(\lambda)=n\log\lambda-\lambda t+C.
$$

For $n>0$,

$$
\ell'(\lambda)=\frac n\lambda-t=0
\quad\Longrightarrow\quad
\boxed{\widehat\lambda=\frac nt}.
$$

The second derivative is $-n/\lambda^2<0$. When $n=0$, the maximum is at $\lambda=0$ if zero is allowed, while the parameter space $\lambda>0$ has only a boundary supremum as $\lambda\downarrow0$.

Exact event times on that fixed window give the mixed density–mass

$$
e^{-\lambda t}\lambda^n
\quad\text{on }0<u_1<\cdots<u_n<t.
$$

Its conditional timing factor $n!/t^n$ is free of $\lambda$, so the total count is sufficient for the homogeneous rate in this fixed-window experiment. This conclusion follows directly from the conditional-uniform theorem, rather than requiring a new timing model.

The estimator has

$$
E[\widehat\lambda]=\lambda,
\qquad
\mathrm{Var}(\widehat\lambda)=\frac{\lambda}{t}.
$$

These statements use fixed observation duration. An experiment that stops at the $n$th event has a different sampling distribution and should not inherit the same unbiasedness claim automatically.

### Deterministic time variation versus a random common rate

An NHPP uses a deterministic rate function and retains independent increments. A different model draws a common random rate $L$ once and then, conditional on $L$, runs an HPP with that rate.

For finite $E[L]$ and $\mathrm{Var}(L)$,

$$
E[N(t)]=tE[L],
$$

$$
\begin{aligned}
\mathrm{Var}(N(t))
&=E[Lt]+\mathrm{Var}(Lt)\\
&=tE[L]+t^2\mathrm{Var}(L).
\end{aligned}
$$

If $A,B$ count disjoint intervals of lengths $a,b$, conditional independence gives

$$
\mathrm{Cov}(A,B)
=\mathrm{Cov}(La,Lb)
=ab\,\mathrm{Var}(L).
$$

Thus a shared nondegenerate random rate creates overdispersion and positive dependence after averaging over the latent rate. “Independent given the rate” is not the same as unconditionally independent.

This is a useful model diagnostic, not a unique explanation of overdispersion. The calculations demonstrate why a deterministic time-varying rate and a random common rate should not be conflated.

## Seven fully worked problems

### Fixed total, early-window count, and conditional dependence

An HPP produces exactly eight events in one hour. Find the probability that exactly three occurred in the first 15 minutes, and find the conditional covariance of the first- and second-quarter-hour counts.

Let $K=N(1/4)$. Given $N(1)=8$,

$$
K\sim\mathrm{Binomial}(8,1/4).
$$

Therefore

$$
P(K=3\mid N(1)=8)
=\binom83(1/4)^3(3/4)^5
=\frac{1701}{8192}
\approx0.207642.
$$

For quarter-hour counts $C_1,\ldots,C_4$,

$$
(C_1,\ldots,C_4)\mid\{N(1)=8\}
\sim\mathrm{Multinomial}(8;1/4,1/4,1/4,1/4).
$$

Thus

$$
E[C_1\mid N(1)=8]=2,
$$

$$
\mathrm{Var}(C_1\mid N(1)=8)=\frac32,
\qquad
\mathrm{Cov}(C_1,C_2\mid N(1)=8)=-\frac12.
$$

The original rate is irrelevant after fixing this total. Multiplying unconditional Poisson probabilities would not answer the conditional question.

### First, third, and last arrivals in a fixed window

Given exactly four events in two hours, find the first-arrival survival probability at 30 minutes, the conditional mean of the third arrival, and the probability that all events occur by 90 minutes.

For the first event,

$$
P(T_1>0.5\mid N(2)=4)
=\left(1-\frac{0.5}{2}\right)^4
=\left(\frac34\right)^4
=\frac{81}{256}
=0.31640625.
$$

For the third event,

$$
\frac{T_3}{2}\mid\{N(2)=4\}\sim\mathrm{Beta}(3,2),
$$

so

$$
E[T_3\mid N(2)=4]=\frac{3\cdot2}{5}=\frac65\text{ hours},
$$

$$
\mathrm{Var}(T_3\mid N(2)=4)
=\frac{3\cdot2\cdot2^2}{5^2\cdot6}
=\frac4{25}\text{ hours}^2.
$$

All events occur by 90 minutes exactly when the maximum is at most $1.5$:

$$
P(T_4\le1.5\mid N(2)=4)
=\left(\frac{1.5}{2}\right)^4
=\frac{81}{256}.
$$

The first and final probabilities happen to be equal by reflection symmetry. Neither calculation treats the first or last chronological arrival as a uniform random variable.

### Thinning, reverse conditioning, and moment checks

Events arrive at rate 10 per hour, and each is independently recorded with probability $p=0.3$. Given that four events were recorded in one hour, find the conditional law and moments of the total event count.

Let $R=R(1)$ and $D=D(1)$. Independent thinning gives

$$
R\sim\mathrm{Poisson}(3),
\qquad
D\sim\mathrm{Poisson}(7),
\qquad R\perp D.
$$

Thus

$$
N(1)\mid\{R=4\}\ \overset d=\ 4+\mathrm{Poisson}(7).
$$

For integers $n\ge4$,

$$
P(N(1)=n\mid R=4)
=e^{-7}\frac{7^{n-4}}{(n-4)!}.
$$

The conditional mean is 11 and the conditional variance is 7. The total is not $\mathrm{Poisson}(11)$; its support excludes 0 through 3.

Check the original mean and variance:

$$
E[E[N(1)\mid R]]=E[R+7]=3+7=10,
$$

$$
E[\mathrm{Var}(N(1)\mid R)]
+\mathrm{Var}(E[N(1)\mid R])
=7+\mathrm{Var}(R)=10.
$$

For comparison, if instead the original total were fixed at 10, the recorded count would be $\mathrm{Binomial}(10,0.3)$. Reversing the conditioning changes the distribution.

### Superposition and a race between event numbers

Two independent streams A and B have rates 2 and 3 per hour. Find the next combined waiting-time distribution, the probability A wins the next event, and $P(T_2^A<T_3^B)$.

The combined stream has rate 5. Thus

$$
W\sim\mathrm{Exp}(5),
\qquad
E[W]=\frac15\text{ hour},
\qquad
P(\text{next label is A})=\frac25.
$$

For the numbered-arrival race, consider the first $2+3-1=4$ combined events. A’s second event precedes B’s third if at least two of those four events carry label A.

Therefore

$$
\begin{aligned}
P(T_2^A<T_3^B)
&=P(\mathrm{Binomial}(4,2/5)\ge2)\\
&=1-\left(\frac35\right)^4
-4\left(\frac25\right)\left(\frac35\right)^3\\
&=\frac{328}{625}=0.5248.
\end{aligned}
$$

The waiting times between combined events are irrelevant to this label-order event. Their rate determines the clock, while label probabilities determine which stream reaches its target first.

### A nonhomogeneous rate and conditional timing

Suppose an NHPP has rate

$$
\lambda(u)=2u,\qquad u\ge0,
$$

so $\Lambda(t)=t^2$. Find the law of the count in $(1,2]$, the first-arrival survival function, and the conditional probability that an event lies in the first hour given exactly one event by time 2.

The interval count has parameter

$$
\int_1^2 2u\,du=4-1=3.
$$

Hence

$$
N(2)-N(1)\sim\mathrm{Poisson}(3),
$$

and

$$
P(T_1>t)=e^{-t^2}.
$$

Given $N(2)=1$, the single event has density

$$
g_2(u)=\frac{2u}{4}=\frac u2,
\qquad 0<u<2.
$$

Therefore

$$
P(T_1\le1\mid N(2)=1)
=\int_0^1\frac u2\,du
=\frac14.
$$

It is not $1/2$: the rate increases with time. The conditional mean event time is

$$
E[T_1\mid N(2)=1]
=\int_0^2u\frac u2\,du
=\frac43.
$$

For a residual wait starting at time 1,

$$
P(W_1>u)
=e^{-[(1+u)^2-1]}
=e^{-2u-u^2}.
$$

The transformed time $T_1^2=\Lambda(T_1)$ is $\mathrm{Exp}(1)$, but $T_1$ itself is not exponential. The units implicit in $\lambda(u)=2u$ require the coefficient to carry inverse-squared-time units.

### Compound arrivals and total variance

Claims arrive at rate 2 per year. Each independent claim amount is 100 with probability 0.75 and 300 with probability 0.25; find the mean and variance of the total amount in three years.

For one claim amount $Y$,

$$
m=E[Y]=0.75(100)+0.25(300)=150,
$$

$$
E[Y^2]=0.75(100^2)+0.25(300^2)=30000,
$$

$$
v=30000-150^2=7500.
$$

The number of claims is $N(3)\sim\mathrm{Poisson}(6)$. Conditional moments of the total amount $S(3)$ are

$$
E[S(3)\mid N(3)]=150N(3),
\qquad
\mathrm{Var}(S(3)\mid N(3))=7500N(3).
$$

Hence

$$
E[S(3)]=6(150)=900,
$$

$$
\begin{aligned}
\mathrm{Var}(S(3))
&=6(7500)+6(150)^2\\
&=45000+135000\\
&=180000.
\end{aligned}
$$

The standard deviation is $\sqrt{180000}\approx424.264$ in amount units. The total is supported on nonnegative multiples of 100, not on all real values or on ordinary event counts.

As a second derivation, thinning produces independent counts $A\sim\mathrm{Poisson}(4.5)$ and $B\sim\mathrm{Poisson}(1.5)$ for the two claim sizes. Then

$$
S(3)=100A+300B,
$$

whose variance is $100^2(4.5)+300^2(1.5)=180000$. This agrees with the conditional-variance calculation.

### Why keeping every second event is not thinning

Start with an HPP of rate $\lambda$ and retain arrivals numbered 2, 4, 6, and so on. Find the first retained waiting-time distribution and explain why the retained process is not an HPP of rate $\lambda/2$.

The first retained event occurs at

$$
T_2=X_1+X_2\sim\mathrm{Gamma}(2,\lambda).
$$

Its survival function is

$$
P(T_2>t)=e^{-\lambda t}(1+\lambda t),
$$

rather than $e^{-(\lambda/2)t}$. Its mean and variance are

$$
E[T_2]=\frac2\lambda,
\qquad
\mathrm{Var}(T_2)=\frac2{\lambda^2}.
$$

An exponential variable with the same mean would have variance $4/\lambda^2$. Equal mean gap lengths therefore do not establish equality of processes.

The retained count is

$$
K(t)=\left\lfloor\frac{N(t)}2\right\rfloor,
$$

and

$$
P(K(t)=0)=P(N(t)\le1)=e^{-\lambda t}(1+\lambda t).
$$

This alone contradicts the proposed $\mathrm{Poisson}(\lambda t/2)$ law. Independent Bernoulli decisions, not a matching long-run retention fraction, are what make the thinning theorem work.

## Practice problems

Try to solve each problem from its event description rather than matching keywords to a formula. When a total is fixed, say explicitly which variables have become dependent.

1. **Conditional allocation:** Given $N(3)=6$ for an HPP, find the law and variance of $N(1)$.
2. **A multinomial event:** Given $N(1)=5$, find the probability of counts $(2,1,2)$ in three equal subintervals.
3. **Earliest arrival:** Given $N(1)=3$, find the density and mean of $T_1$.
4. **Conditional CDF duality:** Given $N(2)=5$, compute $P(T_2\le0.5)$ under this conditioning.
5. **Exact epoch versus count:** Find $E[T_2\mid T_4=t]$ and $E[T_2\mid N(t)=4]$, using natural density versions where necessary.
6. **Spacings:** Given $N(2)=3$, find the mean and variance of each of the four window spacings, and the covariance of two distinct spacings.
7. **Unseen events:** An HPP has rate 8 per hour and each event is independently retained with probability $1/4$. Given $R(2)=5$, find the law, mean, and variance of $N(2)$.
8. **Competing streams:** Independent streams have rates 1 and 4. Find the probability A’s third event precedes B’s first.
9. **NHPP allocation:** For $\lambda(u)=3u^2$, find $P(N(1)=1\mid N(2)=4)$ and the unconditional law of $N(2)-N(1)$.
10. **Compound sum:** At rate 3 per hour, marks are iid $\mathrm{Exp}(2)$, independent of event times. Find the mean and variance of their total over two hours.
11. **Fixed-window inference:** In four hours, 18 HPP events are observed. Find the rate MLE and its model-based plug-in standard error.
12. **Latent rate:** A common random rate equals 2 or 6 with probabilities $1/2$ each. Find $E[N(1)]$, $\mathrm{Var}(N(1))$, and the covariance of counts in two disjoint one-hour intervals.
13. **Finite total intensity:** For $\lambda(u)=e^{-u}$, find the probability there is never any event, and explain why the first-arrival finite-time density does not integrate to 1.
14. **Proof task:** Starting from the mixed arrival-density/count formula, derive the NHPP conditional ordered-arrival density and verify its integral is 1.

## Solutions and checks

### Allocation and timing problems

- **Problem 1:** $N(1)\mid\{N(3)=6\}\sim\mathrm{Binomial}(6,1/3)$. Its mean is 2 and its variance is $6(1/3)(2/3)=4/3$.
- **Problem 2:** The probability is $\frac{5!}{2!1!2!}(1/3)^5=\frac{10}{81}\approx0.123457$.
- **Problem 3:** The survival function is $(1-u)^3$ for $0<u<1$. Thus the density is $3(1-u)^2$ on $(0,1)$, and the mean is $1/4$.
- **Problem 4:** The event is at least two of the five points in the first quarter of the two-hour window. Its probability is $1-(3/4)^5-5(1/4)(3/4)^4=\frac{47}{128}=0.3671875$.
- **Problem 5:** The exact-epoch mean is $2t/4=t/2$. The fixed-count mean is $2t/5$; these condition on different information.
- **Problem 6:** Here $n+1=4$ and $t=2$, giving mean $1/2$, variance $3(4)/(16\cdot5)=3/20$, and distinct-spacing covariance $-4/(16\cdot5)=-1/20$.

### Stream and extension problems

- **Problem 7:** The discarded two-hour mean is $8(3/4)(2)=12$. Hence $N(2)\mid\{R(2)=5\}\overset d=5+\mathrm{Poisson}(12)$, with mean 17 and variance 12.
- **Problem 8:** A must produce the first three combined labels, each with probability $1/5$. The probability is $(1/5)^3=1/125=0.008$.
- **Problem 9:** $\Lambda(t)=t^3$, so the allocation probability is $1/8$. The requested probability is $4(1/8)(7/8)^3=343/1024\approx0.334961$, and $N(2)-N(1)\sim\mathrm{Poisson}(7)$.
- **Problem 10:** The mean mark is $1/2$ and its second moment is $2/2^2=1/2$. With mean event count 6, the total mean is 3 and variance is 3.
- **Problem 11:** The MLE is $18/4=4.5$ per hour. The plug-in standard error is $\sqrt{\widehat\lambda/4}=\sqrt{18}/4\approx1.06066$ per hour; it is an estimated sampling standard deviation, not an exact confidence interval.
- **Problem 12:** $E[L]=4$ and $\mathrm{Var}(L)=4$. Thus the one-hour mean is 4, variance is $4+4=8$, and disjoint one-hour counts have covariance 4.
- **Problem 13:** $\Lambda(\infty)=1$, so $P(T_1=\infty)=e^{-1}$. The finite-time density integrates to $1-e^{-1}$, with the missing mass at infinity.

### NHPP normalization proof

For Problem 14, the mixed density–mass on the ordered region is

$$
e^{-\Lambda(t)}\prod_{j=1}^n\lambda(u_j).
$$

Dividing by $e^{-\Lambda(t)}\Lambda(t)^n/n!$ yields

$$
n!\prod_{j=1}^n\frac{\lambda(u_j)}{\Lambda(t)}.
$$

The product is symmetric in the $n$ coordinates. Its integral over the full cube $(0,t)^n$ is

$$
\prod_{j=1}^n\int_0^t\frac{\lambda(u_j)}{\Lambda(t)}\,du_j=1.
$$

The cube splits into $n!$ ordered regions except for the probability-zero tie sets. Therefore the integral of the product over one ordered region is $1/n!$, and the extra $n!$ factor normalizes the ordered density.

This proof requires $0<\Lambda(t)<\infty$. If $\Lambda(t)=0$, only the zero-count event has positive probability and no positive-count conditional density is defined by the ratio.

## R verification and simulation

### Deterministic distribution checks

The following base-R blocks are designed to run in order in a clean session. The checks compare formulas, integrals, and alternative representations; simulations are reported separately with Monte Carlo uncertainty.

```r
checks_run <- 0L
check <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(
    as.numeric(actual), as.numeric(expected), tolerance = tolerance
  )))
  checks_run <<- checks_run + 1L
}

# Fixed-total binomial law obtained as joint / marginal.
rate <- 4
s <- 0.25
t <- 1
n <- 8
k <- 0:n
ratio <- dpois(k, rate * s) * dpois(n - k, rate * (t - s)) /
  dpois(n, rate * t)
check(ratio, dbinom(k, size = n, prob = s/t))
check(sum(ratio), 1)
check(sum(k * ratio), n * s/t)
check(sum((k - n*s/t)^2 * ratio), n*s/t*(1-s/t))
check(dbinom(3, 8, 1/4), 1701/8192)
check(-8 * (1/4) * (1/4), -1/2)
check(dmultinom(c(2, 1, 2), prob = rep(1/3, 3)), 10/81)

# Conditional arrival times: beta density and binomial CDF duality.
check((1 - 0.5/2)^4, 81/256)
check(pbeta(1.5/2, shape1 = 4, shape2 = 1), 81/256)
check(2 * 3/5, 6/5)
check(2^2 * 3 * 2 / (5^2 * 6), 4/25)
u <- 0.5
check(pbeta(u/2, shape1 = 2, shape2 = 4),
      pbinom(1, size = 5, prob = u/2, lower.tail = FALSE))
check(pbeta(u/2, shape1 = 2, shape2 = 4), 47/128)
check(integrate(function(x) dbeta(x/2, 3, 2) / 2,
                0, 2)$value, 1)
check(integrate(function(x) x * dbeta(x/2, 3, 2) / 2,
                0, 2)$value, 6/5)

# Triangular ordered density for exactly two arrivals on (0, 2].
ordered_mass <- integrate(function(v) {
  vapply(v, function(vv) {
    integrate(function(u) rep(2/2^2, length(u)), 0, vv)$value
  }, numeric(1))
}, 0, 2)$value
check(ordered_mass, 1)

# Thinning: joint factorization and reverse conditional PMF.
rate <- 10
p <- 0.3
q <- 1 - p
r <- 4
d <- 0:100
joint <- dpois(r + d, rate) * dbinom(r, size = r + d, prob = p)
check(joint, dpois(r, rate*p) * dpois(d, rate*q))
conditional <- joint / dpois(r, rate*p)
check(conditional, dpois(d, rate*q))
check(sum(conditional), 1)
check(sum((r + d) * conditional), 11)
check(sum((r + d - 11)^2 * conditional), 7)
check(7 + 3, 10)

# Races: label-order probability versus integration of gamma laws.
race <- pbinom(1, size = 4, prob = 2/5, lower.tail = FALSE)
check(race, 328/625)
race_integral <- integrate(function(x) {
  dgamma(x, shape = 2, rate = 2) *
    pgamma(x, shape = 3, rate = 3, lower.tail = FALSE)
}, 0, Inf)$value
check(race_integral, race)
check(integrate(function(x) 2 * exp(-5*x), 0, Inf)$value, 2/5)

# NHPP with lambda(u) = 2u, Lambda(t) = t^2.
Lambda <- function(t) t^2
check(Lambda(2) - Lambda(1), 3)
check(integrate(function(u) u/2, 0, 2)$value, 1)
check(integrate(function(u) u/2, 0, 1)$value, 1/4)
check(integrate(function(u) u * u/2, 0, 2)$value, 4/3)
check(exp(-(Lambda(1.5) - Lambda(1))), exp(-1.25))

# Spacing moments and compound-sum calculations.
check(3 * 2^2 / (4^2 * 5), 3/20)
check(-2^2 / (4^2 * 5), -1/20)
mark_mean <- 0.75 * 100 + 0.25 * 300
mark_second <- 0.75 * 100^2 + 0.25 * 300^2
check(mark_mean, 150)
check(mark_second - mark_mean^2, 7500)
check(6 * mark_mean, 900)
check(6 * mark_second, 180000)
check(100^2 * 4.5 + 300^2 * 1.5, 180000)
check(pgamma(1, shape = 2, rate = 4, lower.tail = FALSE),
      ppois(1, lambda = 4))

cat("Passed", checks_run, "deterministic checks.\n")
```

The discarded-count sums truncate an infinite tail, with omitted mass below tolerance for mean 7 and the chosen range. The exact equality claims in the mathematical derivations do not rely on that numerical truncation.

### Simulating an HPP by count plus uniform locations

The conditional-uniform theorem gives an exact finite-window algorithm: draw the total count, then draw and sort that many uniform locations. It avoids both time discretization and guessing how many exponential gaps to generate.

```r
simulate_hpp_uniform <- function(rate, horizon) {
  stopifnot(length(rate) == 1L, is.finite(rate), rate >= 0,
            length(horizon) == 1L, is.finite(horizon), horizon >= 0)
  n <- rpois(1, lambda = rate * horizon)
  if (n == 0L) return(numeric(0))
  sort(runif(n, min = 0, max = horizon))
}

set.seed(2544)
times <- simulate_hpp_uniform(rate = 4, horizon = 2)
labels <- rbinom(length(times), size = 1, prob = 0.3)
retained_times <- times[labels == 1]
discarded_times <- times[labels == 0]
stopifnot(length(retained_times) + length(discarded_times) == length(times))
print(list(all = times, retained = retained_times, discarded = discarded_times))
```

This construction simulates the joint event pattern on the chosen finite window, not only the total count. If several adjacent windows are simulated independently with the correct rates, their union also has the appropriate HPP law.

### Conditional order statistics and spacings

The next experiment fixes the total at four points in a two-hour window. Sorting independent uniforms produces dependent chronological arrival times, and differencing gives five dependent window spacings.

```r
set.seed(3544)
B <- 100000L
n <- 4L
horizon <- 2

U <- matrix(runif(B * n, 0, horizon), nrow = B, ncol = n)
ordered <- t(apply(U, 1, sort))
spacings <- cbind(ordered, horizon) - cbind(0, ordered)
stopifnot(max(abs(rowSums(spacings) - horizon)) < 1e-12)

report_mean <- function(x, target) {
  c(estimate = mean(x), target = target,
    Monte_Carlo_SE = sd(x) / sqrt(length(x)))
}

print(report_mean(ordered[, 3], 6/5))
print(report_mean(as.numeric(ordered[, 1] > 0.5), 81/256))
print(c(empirical_variance_T3 = var(ordered[, 3]), target = 4/25))
print(c(empirical_spacing_mean = mean(spacings[, 1]), target = 2/5))
print(c(empirical_spacing_covariance = cov(spacings[, 1], spacings[, 2]),
        target = -2^2 / (5^2 * 6)))
```

The covariance target is for spacings, not ordered arrival epochs. Fixing the row sum of the spacings does not make those spacings independent uniforms or independent exponentials.

### Thinning and compound sums by two representations

Here we simulate a total count and then a conditional binomial split. The resulting retained and discarded counts should be unconditionally independent even though they are constrained to sum to the generated total on each replication.

```r
set.seed(4544)
B <- 200000L
total <- rpois(B, lambda = 10)
recorded <- rbinom(B, size = total, prob = 0.3)
discarded <- total - recorded

print(report_mean(recorded, 3))
print(report_mean(discarded, 7))
print(c(empirical_covariance = cov(recorded, discarded), target = 0))

selected <- total[recorded == 4]
stopifnot(length(selected) > 1L)
print(report_mean(selected, 11))
print(c(selected_replicates = length(selected),
        empirical_conditional_variance = var(selected), target = 7))

# Compound claim amount using independent thinned count streams.
small <- rpois(B, lambda = 4.5)
large <- rpois(B, lambda = 1.5)
amount <- 100 * small + 300 * large
print(report_mean(amount, 900))
print(c(empirical_amount_variance = var(amount), target = 180000))
```

The empirical covariance and variance estimates fluctuate from run to run. The reported Monte Carlo standard errors apply to means and indicator-probability estimates, not automatically to those covariance and variance estimates.

### NHPP simulation by inverse integrated intensity

For $\lambda(u)=2u$, the integrated intensity is $\Lambda(u)=u^2$. Conditional on $N(t)=n$, a location has CDF $(u/t)^2$, so its inverse-CDF draw is $t\sqrt V$ with $V\sim\mathrm{Uniform}(0,1)$.

```r
simulate_nhpp_2t <- function(horizon) {
  stopifnot(length(horizon) == 1L, is.finite(horizon), horizon >= 0)
  n <- rpois(1, lambda = horizon^2)
  if (n == 0L) return(numeric(0))
  sort(horizon * sqrt(runif(n)))
}

set.seed(5544)
print(simulate_nhpp_2t(horizon = 2))

# One point conditional on N(2) = 1: not uniform in ordinary time.
B <- 200000L
conditional_location <- 2 * sqrt(runif(B))
print(report_mean(conditional_location, 4/3))
print(report_mean(as.numeric(conditional_location <= 1), 1/4))

# Unconditional first-arrival time: T1^2 ~ Exp(1).
T1 <- sqrt(rexp(B, rate = 1))
print(report_mean(as.numeric(T1 > 1), exp(-1)))
```

Uniform ordinary-time locations would simulate the wrong process here. The transformation supplies a direct computational check of the integrated-intensity reasoning.

### Verification record

All five R blocks were extracted from this Markdown file and run together in a clean R 4.5.2 session; all 38 deterministic numerical checks passed, and all simulation examples completed without errors. Separately, all 141 display equations and 285 inline math expressions passed a strict KaTeX syntax check, and the Markdown parser recognized the same 426 math expressions.

The document uses `$...$` for inline mathematics and separate-line `$$...$$` display blocks, with `\mathrm` for named distributions and operators. These are syntax and implementation checks, not a claim that this file was published to GitHub or tested inside GitHub's own renderer.

## Common errors and lecture preparation

### Errors to catch before submitting

- **Poisson after fixing the total:** A subwindow count conditional on a fixed total is binomial, not its unconditional Poisson law.
- **Independent after fixing the total:** Partition counts and window spacings acquire fixed-sum dependence.
- **“Arrival times are iid uniform”:** The unordered points are iid uniform in an HPP window; the chronological times are their dependent order statistics.
- **Missing the no-further-arrival factor:** The first-$n$-arrival density alone does not encode $N(t)=n$.
- **Wrong order-statistic denominator:** Given $N(t)=n$, $E[T_k]=kt/(n+1)$, not $kt/n$.
- **Confusing two conditionings:** Fixing $T_n=t$ is not the same as fixing $N(t)=n$.
- **Missing the shift:** $N(t)\mid R(t)=r$ starts at $r$; the unobserved rejected count starts at zero.
- **Assuming independent streams:** Superposition requires independence, not merely two Poisson marginal laws.
- **Calling rank-based selection thinning:** Keeping every second event gives gamma gaps, not exponential gaps.
- **Using an endpoint rate:** NHPP interval counts require $\int_s^t\lambda(u)\,du$.
- **Dropping count variation:** A compound-sum variance uses $\lambda tE[Y^2]$, not only $\lambda t\mathrm{Var}(Y)$.
- **Averaging through a nonlinear inverse:** In general, $E[\Lambda^{-1}(Z)]\ne\Lambda^{-1}(E[Z])$.
- **Confusing rate uncertainty with time variation:** A common latent rate can create dependence between disjoint intervals.

### Questions to bring to lecture

- **Core scope:** Are nonhomogeneous and compound Poisson models included in Part II or deferred?
- **Derivations:** Should you reproduce the conditional joint-density proof, or only use the order-statistic theorem?
- **Spacings:** Is the Dirichlet terminology required, or is the constant-density simplex calculation sufficient?
- **Process independence:** Is a finite-partition proof expected for thinning and superposition?
- **Observation scheme:** Are rate-inference questions based on a fixed time window or stopping after a fixed event number?
- **Computation:** Is R expected only for checking, or will simulation algorithms be examinable?

## Compact reference sheet

### Homogeneous conditioning

Assume $0<s<t$ and $N$ is an HPP of rate $\lambda$. All positive-count conditional statements below use valid indices.

$$
N(s)\mid\{N(t)=n\}
\sim\mathrm{Binomial}(n,s/t).
$$

$$
(C_1,\ldots,C_m)\mid\{N(t)=n\}
\sim\mathrm{Multinomial}(n;\Delta_1/t,\ldots,\Delta_m/t).
$$

$$
f_{T_1,\ldots,T_n\mid N(t)}
(\boldsymbol u\mid n)
=\frac{n!}{t^n}\mathbf1_{\{0<u_1<\cdots<u_n<t\}}.
$$

$$
\frac{T_k}{t}\mid\{N(t)=n\}
\sim\mathrm{Beta}(k,n+1-k),
\qquad
E[T_k\mid N(t)=n]=\frac{kt}{n+1}.
$$

$$
\mathrm{Var}(T_k\mid N(t)=n)
=\frac{k(n+1-k)t^2}{(n+1)^2(n+2)}.
$$

$$
P(T_1>u\mid N(t)=n)=(1-u/t)^n,
\qquad
P(T_n\le u\mid N(t)=n)=(u/t)^n.
$$

### Splitting, combining, and random sums

For independent marking and independent input streams, respectively,

$$
N_R\sim\mathrm{HPP}(\lambda p),
\qquad
N_D\sim\mathrm{HPP}(\lambda(1-p)),
\qquad N_R\perp N_D.
$$

$$
\sum_jN_j\sim\mathrm{HPP}\left(\sum_j\lambda_j\right),
\qquad
P(\text{next label}=j)=\frac{\lambda_j}{\sum_i\lambda_i}.
$$

$$
N(t)\mid\{N_R(t)=r\}
\ \overset d=\ r+\mathrm{Poisson}(\lambda(1-p)t).
$$

$$
E[S(t)]=\lambda tE[Y],
\qquad
\mathrm{Var}(S(t))=\lambda tE[Y^2].
$$

### Time-varying rate

For deterministic locally integrable $\lambda(u)\ge0$, let $\Lambda(t)=\int_0^t\lambda(u)\,du$. For conditioning on a positive total, assume $0<\Lambda(t)<\infty$.

$$
N(t)-N(s)\sim\mathrm{Poisson}(\Lambda(t)-\Lambda(s)).
$$

$$
P(W_s>u)
=e^{-[\Lambda(s+u)-\Lambda(s)]}.
$$

$$
N(s)\mid\{N(t)=n\}
\sim\mathrm{Binomial}\left(n,\frac{\Lambda(s)}{\Lambda(t)}\right).
$$

$$
g_t(u)=\frac{\lambda(u)}{\Lambda(t)}\mathbf{1}_{\{0<u<t\}}
$$

is the density of an unordered point in the fixed-count window. Chronological points are its order statistics, not independent draws in their chronological order.

The final diagnostic question is: **What is fixed, what remains random, and which independence statements still hold after conditioning?** Answering it usually determines both the distribution and the correct support.
