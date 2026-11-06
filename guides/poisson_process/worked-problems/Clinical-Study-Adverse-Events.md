# Rare Adverse Events in a Clinical Study
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

## Calculate and interpret $\sum_i p_i^2$: part 

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