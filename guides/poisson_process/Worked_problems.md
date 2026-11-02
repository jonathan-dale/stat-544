
# Problem 1. Manufacturing Defects 
A factory manufactures electronic sensors used in medical devices. Each sensor independently has a probability of 0.002 of being defective. A quality-control engineer examines 1,500 sensors produced during one shift.  
(a) Let X denote the number of defective sensors. Identify the exact probability distribution of X.  
(b) Use a Poisson distribution to approximate P(X = 4).  
(c) Approximate P(X ≥ 5).  
(d) Calculate the exact binomial probabilities for (b) and (c), and compare them with the Poisson approximations.  


The exact model is **binomial**, and the Poisson approximation uses $\lambda=np=3$. Here are the answers first; then we’ll unpack why these distributions apply.

| Part | Requested result | Answer |
|---|---|---:|
| (a) | Exact distribution | $X\sim\mathrm{Binomial}(1500,0.002)$ |
| (b) | Poisson approximation to $P(X=4)$ | $0.168031$ |
| (c) | Poisson approximation to $P(X\geq5)$ | $0.184737$ |
| (d) | Exact binomial $P(X=4)$ | $0.168200$ |
| (d) | Exact binomial $P(X\geq5)$ | $0.184568$ |

The central idea is that a Poisson distribution approximates a binomial count when there are many independent opportunities for an event, but the event is rare on each opportunity. The approximation matches the expected count, $\lambda=np$. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

## (a) Why the exact distribution is binomial

For each sensor, define an indicator random variable:

$$
I_i=
\begin{cases}
1,&\text{sensor }i\text{ is defective},\\
0,&\text{sensor }i\text{ is not defective}.
\end{cases}
$$

Each $I_i$ is a Bernoulli random variable with parameter $p=0.002$. The total number of defective sensors is

$$
X=I_1+I_2+\cdots+I_{1500}.
$$

A binomial random variable counts successes in a fixed number of independent Bernoulli trials with the same success probability. Here, “success” is simply the outcome we are counting—being defective—not a desirable outcome. Independence and a constant probability are essential to this model. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

Thus,

$$
\boxed{X\sim\mathrm{Binomial}(n=1500,p=0.002).}
$$

Its probability mass function is

$$
P(X=k)=\binom{1500}{k}(0.002)^k(0.998)^{1500-k},
\qquad k=0,1,\ldots,1500.
$$

To understand that formula, break it into two pieces:

- $(0.002)^k(0.998)^{1500-k}$ is the probability of one particular arrangement of $k$ defective and $1500-k$ nondefective sensors, by independence.
- $\binom{1500}{k}$ counts how many different arrangements have exactly $k$ defective sensors.

The binomial mean and variance are $np$ and $np(1-p)$, respectively.  In this problem, [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson24)

$$
E[X]=1500(0.002)=3,
$$

$$
\mathrm{Var}(X)=1500(0.002)(0.998)=2.994.
$$

So we expect three defective sensors per shift on average—not necessarily exactly three in any particular shift.

## Why a Poisson approximation makes sense

The relevant limiting result is

$$
n\to\infty,\qquad p\to0,\qquad np\to\lambda
\quad\Longrightarrow\quad
\mathrm{Binomial}(n,p)\text{ approaches }\mathrm{Poisson}(\lambda).
$$

In words: as the number of independent opportunities grows and the event probability shrinks, with the expected count staying finite, the binomial probabilities approach Poisson probabilities. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

Here, $n=1500$ is large and $p=0.002$ is small, so we use a separate approximating variable

$$
Y\sim\mathrm{Poisson}(3).
$$

Writing $Y$ separately helps keep the distinction clear: $X$ is exactly binomial; $Y$ provides approximate probabilities for $X$.

The Poisson probability mass function is

$$
P(Y=k)=\frac{e^{-\lambda}\lambda^k}{k!},
\qquad k=0,1,2,\ldots.
$$

Unlike the binomial formula, it depends only on the expected count $\lambda$, rather than separately on $n$ and $p$. [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

### Where the Poisson formula comes from

For a fixed count $k$, substitute $p=\lambda/n$ into the binomial formula:

$$
\begin{aligned}
P(X=k)
&=\binom{n}{k}
\left(\frac{\lambda}{n}\right)^k
\left(1-\frac{\lambda}{n}\right)^{n-k}\\
&=
\frac{\lambda^k}{k!}
\left[\prod_{j=0}^{k-1}\left(1-\frac{j}{n}\right)\right]
\left(1-\frac{\lambda}{n}\right)^n
\left(1-\frac{\lambda}{n}\right)^{-k}.
\end{aligned}
$$

As $n\to\infty$, the three factors after $\lambda^k/k!$ approach $1$, $e^{-\lambda}$, and $1$, respectively. That leaves

$$
P(X=k)\longrightarrow \frac{e^{-\lambda}\lambda^k}{k!}.
$$

This is the mathematical basis of the approximation—not just a rule about “large $n$, small $p$.” [online.stat.psu](https://online.stat.psu.edu/stat504/Lesson01)

## (b) Approximate $P(X=4)$

Using $Y\sim\mathrm{Poisson}(3)$,

$$
\begin{aligned}
P(X=4)
&\approx P(Y=4)\\
&=\frac{e^{-3}3^4}{4!}\\
&=\frac{81e^{-3}}{24}\\
&\approx \boxed{0.168031}.
\end{aligned}
$$

Thus, the approximate probability of exactly four defective sensors is $16.8031\%$.

## (c) Approximate $P(X\geq5)$

“At least five” includes $5,6,7,\ldots$. Rather than adding all those probabilities, use the complement:

$$
P(Y\geq5)=1-P(Y\leq4).
$$

The endpoint matters: the complement of “at least five” is “at most four,” not “at most five.”

Therefore,

$$
\begin{aligned}
P(X\geq5)
&\approx 1-\sum_{k=0}^{4}\frac{e^{-3}3^k}{k!}\\
&=1-e^{-3}
\left(
1+3+\frac{9}{2}+\frac{27}{6}+\frac{81}{24}
\right)\\
&=1-16.375e^{-3}\\
&\approx \boxed{0.184737}.
\end{aligned}
$$

So the approximate probability of at least five defective sensors is $18.4737\%$.

## (d) Exact binomial probabilities and comparison

“Exact binomial” means using the original binomial model rather than approximating it with another distribution. The displayed decimals are still rounded.

For exactly four defects,

$$
\begin{aligned}
P(X=4)
&=\binom{1500}{4}(0.002)^4(0.998)^{1496}\\
&\approx \boxed{0.168200}.
\end{aligned}
$$

For at least five defects, use the same complement strategy:

$$
\begin{aligned}
P(X\geq5)
&=1-P(X\leq4)\\
&=1-\sum_{k=0}^{4}
\binom{1500}{k}(0.002)^k(0.998)^{1500-k}\\
&\approx \boxed{0.184568}.
\end{aligned}
$$

The numerical comparison is:

| Event | Exact binomial | Poisson approximation | Absolute error | Relative error |
|---|---:|---:|---:|---:|
| $X=4$ | 0.168200 | 0.168031 | 0.000168 | 0.1000% |
| $X\geq5$ | 0.184568 | 0.184737 | 0.000168 | 0.0912% |

Here,

$$
\text{relative error}
=
\frac{|\text{approximation}-\text{exact}|}{\text{exact}}\times100\%.
$$

The Poisson approximation slightly underestimates $P(X=4)$ and slightly overestimates $P(X\geq5)$. Both errors are only about $0.0168$ percentage points. This close agreement is consistent with the large-$n$, small-$p$ setting in which the approximation works well. [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson12)

### Checking the results in R

These commands reproduce all four probabilities:

```r
n <- 1500
p <- 0.002
lambda <- n * p

# Poisson approximations
dpois(4, lambda = lambda)
ppois(4, lambda = lambda, lower.tail = FALSE)

# Exact binomial probabilities
dbinom(4, size = n, prob = p)
pbinom(4, size = n, prob = p, lower.tail = FALSE)
```

For both tail calculations, `lower.tail = FALSE` with a cutoff of $4$ gives $P(X>4)=P(X\geq5)$.

The reusable pattern is: identify the exact count model first, set $\lambda=np$ for the Poisson approximation, and translate the event carefully before calculating—especially the cutoff in a tail probability.


# Problem 2. Rare Adverse Events in a Clinical Study
A clinical study enrolls 600 patients to evaluate the safety of a treatment. Patients are divided into three independent groups according to baseline health conditions. Each patient can experience at most one adverse event in the observation period; indicators are independent.
| Patient group | Number of patients | Adverse-event probability |
|---|---|---:|
| Low risk | 300 | 0.002 |
| Moderate risk | 200 | 0.005 |
| High risk | 100 | 0.010 |


(a) Let X denote the number of patients experiencing an event. Explain why X is not binomial.  
(b) Approximate the distribution of X using a Poisson distribution and calculate its parameter.  
(c) Approximate P(X = 0).  
(d) Approximate P(X ≥ 3).  
(e) Compute Σ pᵢ², the bound presented in the lecture, and interpret what it tells us.  

The main change from Problem 1 is that patients now have different event probabilities. The total is therefore not binomial, but a **Poisson** approximation still works: its parameter is the sum of all the individual event probabilities. This is the setting covered by Le Cam’s Poisson approximation theorem. [en.wikipedia](https://en.wikipedia.org/wiki/Le_Cam's_theorem)

Here are the results:

| Part | Answer |
|---|---|
| (a) | $X$ is not binomial because the event probability is not the same for every patient. Its exact distribution is called Poisson-binomial. |
| (b) | $X\approx\mathrm{Poisson}(2.6)$. |
| (c) | $P(X=0)\approx 0.074274$. |
| (d) | $P(X\geq3)\approx 0.481570$. |
| (e) | $\sum_{i=1}^{600}p_i^2=0.0162$, giving an absolute probability-error bound of $0.0162$, or $1.62$ percentage points, for any event involving the total count. |

Let’s build the reasoning behind each part.

## (a) Why $X$ is not binomial

Define an indicator for patient $i$:

$$
I_i=
\begin{cases}
1,&\text{patient }i\text{ experiences an adverse event},\\
0,&\text{otherwise}.
\end{cases}
$$

Then

$$
X=\sum_{i=1}^{600} I_i.
$$

Each $I_i$ is Bernoulli with event probability $p_i$. However, those probabilities differ:

$$
p_i=
\begin{cases}
0.002,&\text{for a low-risk patient},\\
0.005,&\text{for a moderate-risk patient},\\
0.010,&\text{for a high-risk patient}.
\end{cases}
$$

An ordinary binomial model requires independent trials with the same event probability. A sum of independent Bernoulli variables with potentially different probabilities instead has a Poisson-binomial distribution. Despite its name, that is an exact distribution, not a Poisson approximation. [en.wikipedia](https://en.wikipedia.org/wiki/Le_Cam's_theorem)

Another useful way to represent the exact model is to define the three group counts:

$$
X_L\sim\mathrm{Binomial}(300,0.002),
$$

$$
X_M\sim\mathrm{Binomial}(200,0.005),
$$

$$
X_H\sim\mathrm{Binomial}(100,0.010).
$$

Then

$$
X=X_L+X_M+X_H,
$$

with the three counts independent. Each group count is binomial, but their total is not an ordinary binomial because the probabilities differ across groups.

## (b) Find the Poisson parameter

The generalized rare-event approximation is

$$
X=\sum_i I_i
\quad\approx\quad
Y\sim\mathrm{Poisson}(\lambda),
\qquad
\lambda=\sum_i p_i,
$$

when the indicators are independent and their individual probabilities are small. Unlike a binomial model, this approximation does not require equal probabilities. [en.wikipedia](https://en.wikipedia.org/wiki/Le_Cam's_theorem)

The reason for summing the probabilities is that we match the expected count:

$$
E[X]
=
E\left[\sum_i I_i\right]
=
\sum_i E[I_i]
=
\sum_i p_i.
$$

Calculate each group’s contribution:

| Group | Expected number experiencing an event |
|---|---:|
| Low risk | $300(0.002)=0.6$ |
| Moderate risk | $200(0.005)=1.0$ |
| High risk | $100(0.010)=1.0$ |
| Total | $2.6$ |

Therefore,

$$
\boxed{\lambda=0.6+1.0+1.0=2.6}
$$

and we use

$$
\boxed{Y\sim\mathrm{Poisson}(2.6)}
$$

to approximate $X$.

Connecting this to your earlier question: $2.6$ is the expected count, not the count whose probability we are asking for. It stays $2.6$ in both parts (c) and (d).

## (c) Approximate $P(X=0)$

The Poisson probability formula is

$$
P(Y=k)=\frac{e^{-\lambda}\lambda^k}{k!}.
$$

Here, $\lambda=2.6$ and the requested count is $k=0$. Thus,

$$
\begin{aligned}
P(X=0)
&\approx P(Y=0)\\
&=\frac{e^{-2.6}(2.6)^0}{0!}\\
&=e^{-2.6}\\
&\approx \boxed{0.074274}.
\end{aligned}
$$

We used $(2.6)^0=1$ and $0!=1$. The approximate probability that no patient experiences an adverse event is $7.43\%$.

## (d) Approximate $P(X\geq3)$

Use the complement. The opposite of “at least three” is “zero, one, or two”:

$$
P(X\geq3)
\approx
1-P(Y\leq2).
$$

Therefore,

$$
\begin{aligned}
P(X\geq3)
&\approx 1-\bigl[P(Y=0)+P(Y=1)+P(Y=2)\bigr]\\
&=1-e^{-2.6}
\left(
1+2.6+\frac{(2.6)^2}{2!}
\right)\\
&=1-e^{-2.6}(1+2.6+3.38)\\
&=1-6.98e^{-2.6}\\
&\approx \boxed{0.481570}.
\end{aligned}
$$

So the approximate probability of at least three patients experiencing an event is $48.16\%$.

## (e) Calculate and interpret $\sum_i p_i^2$

The index $i$ runs over all 600 patients, not just the three groups. Consequently, each group’s squared probability must be multiplied by its number of patients:

$$
\begin{aligned}
\sum_{i=1}^{600}p_i^2
&=300(0.002)^2
+200(0.005)^2
+100(0.010)^2\\
&=300(0.000004)
+200(0.000025)
+100(0.0001)\\
&=0.0012+0.0050+0.0100\\
&=\boxed{0.0162}.
\end{aligned}
$$

Notice that we square each patient’s probability before summing. This is not $\left(\sum_i p_i\right)^2$.

### What the bound guarantees

Le Cam’s inequality implies that, for any set of possible counts $A$,

$$
\left|P(X\in A)-P(Y\in A)\right|
\leq \sum_i p_i^2
=0.0162.
$$

This is a bound on total variation distance: the largest absolute difference between the two models’ probabilities for the same event. [en.wikipedia](https://en.wikipedia.org/wiki/Le_Cam's_theorem)

For example, it guarantees both

$$
\left|P(X=0)-0.0742736\right|\leq0.0162
$$

and

$$
\left|P(X\geq3)-0.4815704\right|\leq0.0162.
$$

Thus, the approximation error is at most $1.62$ percentage points for either probability. This is an absolute error bound, not a $1.62\%$ relative error, and it is a worst-case guarantee rather than a claim that the actual error equals $0.0162$.

### Why you might see a factor of two

Some statements of Le Cam’s inequality use the sum of absolute errors across all individual counts:

$$
\sum_{k=0}^{\infty}
\left|P(X=k)-P(Y=k)\right|
\leq 2\sum_i p_i^2
=0.0324.
$$

Total variation distance is half of that sum, so the corresponding bound for any single event remains $0.0162$. The two formulations are consistent. [en.wikipedia](https://en.wikipedia.org/wiki/Le_Cam's_theorem)

The useful lesson is that independence and small individual probabilities allow a Poisson approximation even when the probabilities differ. The parameter $\sum_i p_i=2.6$ sets the expected count, while $\sum_i p_i^2=0.0162$ controls the approximation error.