# Conditional Distributions: A theory-first lesson

#### Problems:
1. It is well-known that hypertension (high blood pressure) and type II diabetes are major risk factors for cardiovascular diseases. Suppose in a group of patients with both hypertension and diabetes, a bivariate normal distribution is used to model the diastolic blood pressure (B) and fasting plasma glucose (G) with the following statistic:   
$𝐸(𝐵)=80\ mmHg$,  
$𝑆𝐷(𝐵)=7.5 mmHg$,  
$𝐸(𝐺)=117 mg/dL$,  
$𝑆𝐷(𝐺)=24 mg/dL$.  
Also assume that the correlation between the diastolic blood pressure and fasting plasma glucose is 0.3.

    A. Estimate the proportion of patients whose diastolic blood pressure at least 80 AND the fasting glucose is at most 126.

    B. Suppose your colleague ignored the correlation and treated the diabolic blood pressure and the fasting glucose as two independent variables. Then they estimated the proportion of patients whose diastolic blood pressure at least 80 AND the fasting glucose is at most 126. What would be their estimate ? How would you explain the discrepancy between your and their answers?

2. Suppose the lifetime of a certain type of light bulb can be modeled by a gamma distribution with shape parameter 𝛼=2 and rate parameter 𝜆=1/500.

    A. Find the (theoretical) mean and standard deviation of the distribution.  
    B. Generate 1000 lifetimes and report the summary statistics.   
    C. Select lifetimes that are at least 150 hours and reduce that amount from each remaining observation. Compare the summary statistics with that in (b).  
    D. Overlay the histograms using data in (b) and (c). Do they match well? Explain.  

See [problems.pdf](./problems.pdf) for fully worked solutions.

## Lesson
The idea connecting the two problems is **conditioning**: learning something about an observation changes which outcomes you should consider and how likely those outcomes are.

- In the patient problem, knowing blood pressure is at least 80 changes the distribution of glucose because the variables are correlated.
- In the bulb problem, knowing a bulb has survived 150 hours changes the distribution of its lifetime. Subtracting 150 then measures its remaining life.

This lesson builds those connections step by step. If you have limited time, prioritize Sections 2, 4, and 5; those contain the main reasoning behind the homework.

## 1. The foundation: What a probability distribution actually tells you

### A random variable is a numerical outcome

For these problems:

- $B$: a randomly selected patient’s diastolic blood pressure.
- $G$: that patient’s fasting glucose.
- $X$: a randomly selected bulb’s total lifetime.

A distribution describes how probability is spread across the possible values of a random variable.

The distribution is the model. Your simulated observations are a sample from that model.

### Density is not probability

For a continuous random variable, a probability density function $f(x)$ describes where outcomes are concentrated. Probability is an area under the density:

$$
P(a\leq X\leq b)=\int_a^b f(x)\ dx
$$

For example:

$$
P(X\geq150)=\int_{150}^{\infty}f(x)\ dx
$$

A density’s height at one point is not the probability of that exact outcome. In a continuous model, individual points have probability zero. Consequently, using $>$ rather than $\geq$ does not change the theoretical probabilities in this homework.

### Three descriptions of probability you will encounter

| Description | Notation | Meaning |
|---|---|---|
| Density | $f(x)$ | How concentrated outcomes are near $x$ |
| Cumulative distribution function | $F(x)=P(X\leq x)$ | Probability accumulated up to $x$ |
| Survival function | $S(x)=P(X>x)=1-F(x)$ | Probability of exceeding $x$ |

The survival function is particularly useful for lifetimes: $S(150)$ is the proportion of bulbs expected to last beyond 150 hours.

## 2. Patient problem: Understanding joint probability

### 2.1 One patient has two measurements

Blood pressure and glucose are not two unrelated lists here. They are two measurements on the same person.

A joint distribution describes the pairs:

$$
(B,G)
$$

The problem assumes a bivariate normal distribution: a **two-dimensional normal model**. Its parameters are two means, two variances, and a correlation. A useful property is that each individual variable has a normal marginal distribution. [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson21.html)

Here:

$$
B\sim N(80,7.5^2),
\qquad
G\sim N(117,24^2),
\qquad
\rho=0.3
$$

I am using $N(\mu,\sigma^2)$, so the second argument is variance, not standard deviation.

### 2.2 Translate the words into an event before calculating

The question asks for:

> Blood pressure at least 80 AND glucose at most 126.

Translate each phrase:

$$
B\geq80,
\qquad
G\leq126
$$

“AND” means both must hold:

$$
P(B\geq80,\ G\leq126)
$$

Imagine a scatterplot with blood pressure on the horizontal axis and glucose on the vertical axis:

1. Draw a vertical line at $B=80$.
2. Keep the region to its right.
3. Draw a horizontal line at $G=126$.
4. Keep the region below it.

The desired probability is the amount of the joint distribution in that overlapping region.

This is why calculating only a blood-pressure probability or only a glucose probability is insufficient.

### 2.3 Why you cannot automatically multiply

The general multiplication rule is:

$$
P(A\cap C)=P(A)P(C\mid A)
$$

Here, that becomes:

$$
P(B\geq80,\ G\leq126)=P(B\geq80)P(G\leq126\mid B\geq80)
$$

Read it as:

> First find the proportion meeting the blood-pressure condition. Then, among those people, find the proportion meeting the glucose condition.

Independence would allow you to replace the conditional probability with the marginal probability:

$$
P(G\leq 126\mid B\geq 80)=P(G\leq 126)
$$

But that is an assumption, not a universal rule. In a bivariate normal model, nonzero correlation rules out independence; zero correlation and independence are equivalent for this particular family. That equivalence does not hold for arbitrary distributions. [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson21.html)

### 2.4 What positive correlation changes

Positive correlation means higher blood pressure tends to occur with higher glucose—not for every person, but as a statistical tendency.

For a bivariate normal distribution:

$$
E(G\mid B=b)=117+0.3\frac{24}{7.5}(b-80)
$$

Thus:

$$
E(G\mid B=b)=117+0.96(b-80)
$$

For example:

- At $(B=80)$, the conditional mean glucose is $117$.
- At $(B=90)$, the conditional mean glucose is $126.6$.

The conditional variance is:

$$
\mathrm{Var}(G\mid B=b)=24^2(1-0.3^2).
$$

So the normal model explicitly tells us how the glucose distribution changes when blood pressure is known. [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson21.html)

Notice an important distinction:

$$
G\mid B=80
\quad\text{and}\quad
G\mid B\geq80
$$

are different distributions. The first fixes blood pressure at one value; the second combines all blood pressures above a threshold.

### 2.5 The independent calculation

Because 80 is the mean of a normal distribution:

$$
P(B\geq80)=0.5
$$

For glucose, standardize:

$$
Z=\frac{G-117}{24}
$$

At $(G=126)$:

$$
z=\frac{126-117}{24}=0.375
$$

Therefore:

$$
P(G\leq126)=\Phi(0.375)\approx0.64617
$$

where $(\Phi)$ is the standard normal CDF.

If independence were assumed:

$$
P(B\geq80,\;G\leq126)=0.5(0.64617)\approx0.32308
$$

That is approximately $32.31\%$.

Standardization changes the units, not the dependence: standardizing $B$ and $G$ does not make them independent. Their correlation remains $0.3$. [probabilitycourse](https://www.probabilitycourse.com/chapter5/5_3_2_bivariate_normal_dist.php)

### 2.6 The correlated answer and its interpretation

Our earlier joint calculation gave:

$$
P(B\geq80,\;G\leq126)\approx0.27799.
$$

That is approximately \(27.80\%\).

Consequently:

$$
P(G\leq126\mid B\geq80)=\frac{0.27799}{0.5}\approx0.55597
$$

Compare:

| Group | Probability glucose is at most 126 |
|---|---:|
| All patients | $64.62\%$ |
| Patients with blood pressure at least 80 | $55.60\%$ |

The independence calculation overestimates the answer because it uses the overall glucose probability instead of the probability within the selected blood-pressure group.

A concise explanation for this problem is:

> Because blood pressure and glucose are positively correlated, selecting patients with blood pressure at least 80 favors higher glucose values. Therefore, glucose at most 126 is less likely within that selected group than in the full population. Ignoring correlation overestimates the joint probability.

## 3. Bulb problem: What the gamma distribution means

### 3.1 Why gamma is a plausible lifetime model

A gamma distribution models positive quantities and can accommodate right-skewed outcomes: many moderate values and a smaller number of large values.

For integer shape $\alpha$, it also has a useful interpretation as the sum of $\alpha$ independent exponential waiting times with the same rate. [distribution-explorer.github](https://distribution-explorer.github.io/continuous/gamma.html)

For your parameters:

$$
X\sim\mathrm{Gamma}\left(\alpha=2,\lambda=\frac1{500}\right).
$$

One mathematical representation is:

$$
X=T_1+T_2,
$$

where $T_1$ and $T_2$ are independent exponential random variables, each with mean 500 hours.

This is an interpretation of the distribution—not evidence that the actual bulb physically fails through two such stages.

### 3.2 Rate and scale are reciprocals

Your rate is:

$$
\lambda=\frac1{500}=0.002\text{ per hour}
$$

The equivalent scale is:

$$
\theta=\frac1\lambda=500\text{ hours}
$$

The formulas are:

| Quantity | Shape–rate form | Shape–scale form |
|---|---|---|
| Mean | $\alpha/\lambda$ | $\alpha\theta$ |
| Variance | $\alpha/\lambda^2$ | $\alpha\theta^2$ |
| Standard deviation | $\sqrt{\alpha}/\lambda$ | $\sqrt{\alpha}\theta$ |

Both parameterizations describe the same distribution. [distribution-explorer.github](https://distribution-explorer.github.io/continuous/gamma.html)

Substituting your values:

$$
E(X)=1000\text{ hours},
$$

$$
\mathrm{Var}(X)=500{,}000\text{ hours}^2,
$$

$$
SD(X)=500\sqrt2\approx707.11\text{ hours}.
$$

Variance has squared units; standard deviation has the original units. That is one reason SD is easier to interpret.

## 4. Simulation: Why your sample does not equal the theoretical model

In part B, you generated:

$$
X_1,\ldots,X_{1000}
$$

as independent observations from the gamma distribution.

The theoretical mean is:

$$
\mu=E(X)=1000.
$$

The sample mean is:

$$
\bar X=\frac1{1000}\sum_{i=1}^{1000}X_i.
$$

Your sample mean was $965.146$. That does not contradict the model. A random sample fluctuates.

### 4.1 Spread of observations versus uncertainty in their average

These are different:

- $SD(X)$: how much individual bulb lifetimes vary.
- $SD(\bar X)$: how much sample means vary across repeated samples.

For independent observations:

$$
SD(\bar X)=\frac{SD(X)}{\sqrt n}.
$$

Here:

$$
SD(\bar X)=\frac{707.11}{\sqrt{1000}}\approx22.36\text{ hours}.
$$

Your mean is about $34.85$ hours below the theoretical mean, or roughly $1.56$ standard errors. That is not an unusual discrepancy.

Increasing the sample size makes the sample mean more stable. It does not make individual bulb lifetimes less variable.

### 4.2 How to read your summary

Your original sample had:

| Statistic | Value |
|---|---:|
| Minimum | 6.295 |
| First quartile | 452.561 |
| Median | 820.899 |
| Mean | 965.146 |
| Third quartile | 1321.466 |
| Maximum | 4873.299 |

Interpretations:

- About 25% of observed lifetimes are below the first quartile.
- About 50% are below the median.
- About 75% are below the third quartile.
- The mean uses the magnitudes of all observations, so large lifetimes pull it upward.

The mean exceeding the median is consistent with a right-skewed sample.

Your observed minimum and maximum are sample extremes, not the endpoints of the gamma distribution. The model allows positive lifetimes beyond your observed range. [distribution-explorer.github](https://distribution-explorer.github.io/continuous/gamma.html)

## 5. Part C: Selecting survivors and measuring remaining life

This is the central conceptual step in the bulb problem.

### 5.1 Two operations happen—not one

The instruction is:

> Select lifetimes at least 150 hours and reduce that amount from each remaining observation.

That means:

1. Keep observations satisfying $X\geq150$.
2. Subtract 150 from each retained observation.

The resulting random variable is described by:

$$
R=(X-150\mid X\geq150).
$$

In words:

> Remaining lifetime, given that the bulb survived to age 150.

For example:

| Original lifetime | Included? | Remaining lifetime |
|---|---|---:|
| 100 | No | — |
| 200 | Yes | 50 |
| 1000 | Yes | 850 |

### 5.2 Why the new mean is not $1000-150$

If you subtracted 150 from every observation without selecting, then:

$$
E(X-150)=E(X)-150=850.
$$

But that is not the task.

The task first removes short-lived bulbs. Consequently:

$$
E(R)=E(X\mid X\geq150)-150.
$$

The conditional mean $E(X\mid X\geq150)$ is larger than the original mean because the smallest lifetimes have been excluded.

That selection partly offsets the subsequent subtraction.

For this particular gamma model, the theoretical remaining-life mean is:

$$
E(R)\approx884.62\text{ hours}.
$$

Thus:

$$
E(X\mid X\geq150)\approx1034.62\text{ hours},
$$

and:

$$
1034.62-150=884.62.
$$

These are population benchmarks. Your simulated remaining-life mean will differ.

### 5.3 A short derivation of the new distribution

Let $t$ be the age already survived, and let $r\geq0$ be an additional amount of time.

Then:

$$
P(R>r)=P(X>t+r\mid X>t).
$$

Using conditional probability:

$$
P(R>r)=\frac{P(X>t+r)}{P(X>t)}=\frac{S_X(t+r)}{S_X(t)}.
$$

This equation is the general recipe for a remaining-lifetime distribution.

For a shape-2 gamma distribution, integrating its density gives:

$$
S_X(x)=e^{-\lambda x}(1+\lambda x).
$$

Therefore:

$$
P(R>r)=e^{-\lambda r}\frac{1+\lambda t+\lambda r}{1+\lambda t}.
$$

At $t=150$, $\lambda t=0.3$:

$$
P(R>r)=e^{-r/500}\left(1+\frac{r}{650}\right).
$$

Compare this with the original survival function:

$$
P(X>r)=e^{-r/500}\left(1+\frac{r}{500}\right).
$$

They are not the same. This directly establishes that the remaining-life distribution differs from the original gamma distribution. The derivation uses the shape-2 gamma density supplied by the model. [distribution-explorer.github](https://distribution-explorer.github.io/continuous/gamma.html)

### 5.4 How much should selection change the sample size?

The survival probability is:

$$
P(X\geq150)
=
e^{-0.3}(1.3)
\approx0.96306.
$$

So approximately \(96.31\%\) of bulbs survive to 150 hours.

Out of 1000 simulated lifetimes, the expected number retained is:

$$
1000(0.96306)\approx963.
$$

The actual count is random; it need not equal 963.

### 5.5 What happens to the standard deviation?

Subtracting a constant does not change variance:

$$
\mathrm{Var}(Y-c)=\mathrm{Var}(Y).
$$

But selection can change variance.

Therefore:

$$
SD(R)=SD(X\mid X\geq150),
$$

not necessarily $SD(X)$.

For this model, the theoretical comparison is:

| Quantity | Original lifetime $X$ | Remaining life $R$ |
|---|---:|---:|
| Mean | 1000 hours | 884.62 hours |
| Standard deviation | 707.11 hours | 697.63 hours |
| Distribution represented | All bulbs | Bulbs that survived 150 hours, with elapsed time removed |

The remaining-life values come from the conditional distribution derived above. Do not substitute these theoretical values for the sample statistics requested in part C.

## 6. Part D: Should the histograms match?

### 6.1 They should be similar, but not identical

Both distributions are right-skewed, and only about $3.69\%$ of observations are excluded. Therefore, substantial overlap is reasonable.

But they represent different quantities:

- B: total lifetimes of new bulbs.
- C: remaining lifetimes of bulbs that have survived 150 hours.

The conditional-distribution calculation shows that their population distributions differ. Even with an extremely large simulation, the difference would not disappear.

### 6.2 The missing concept: Memorylessness

A distribution is memoryless when:

$$
P(X>t+r\mid X>t)=P(X>r).
$$

In words:

> Having already survived $t$ hours does not change the distribution of the additional lifetime.

The exponential distribution has this property. A gamma distribution with shape 1 is exponential; your gamma distribution has shape 2, so this special case does not apply. [distribution-explorer.github](https://distribution-explorer.github.io/continuous/gamma.html)

You can see the distinction algebraically. For an exponential distribution:

$$
S(x)=e^{-\lambda x}.
$$

Then:

$$
\begin{aligned}
\frac{S(t+r)}{S(t)}
=
\frac{e^{-\lambda(t+r)}}{e^{-\lambda t}}
=
e^{-\lambda r}
=
S(r).
\end{aligned}
$$

For your shape-2 gamma model, the extra factor involving \(t\) does not cancel. That is why the original and remaining-life distributions differ.

### 6.3 Separate distribution differences from plotting differences

A histogram is an approximation based on data and bin choices.

For a fair comparison:

- Use the same bin boundaries.
- Use the same horizontal scale.
- Prefer density rather than counts because the sample sizes differ.

On a density histogram:

$$
\text{bar height}
=
\frac{\text{count in bin}}
{\text{sample size}\times\text{bin width}}.
$$

Thus:

$$
\text{bar area}
=
\frac{\text{count in bin}}{\text{sample size}},
$$

so area represents relative frequency. Each histogram has total area 1. [rdocumentation](https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/hist)

A suitable explanation for D is:

> The histograms overlap substantially, but they do not match exactly. Part C represents remaining lifetimes conditional on survival to 150 hours, rather than original lifetimes. The gamma distribution with shape 2 is not memoryless, so these distributions differ. Random sampling and histogram binning also contribute to visible differences.

## 7. A reusable problem-solving checklist

Before opening R, answer these questions.

### Step 1: What does the random variable measure?

For example:

$$
X=\text{total bulb lifetime}.
$$

Do not confuse total lifetime with remaining lifetime.

### Step 2: What distribution and parameterization are given?

For example:

$$
X\sim\mathrm{Gamma}(\text{shape}=2,\text{rate}=1/500).
$$

Write “rate” or “scale” explicitly.

### Step 3: Translate the question into symbols

Examples:

$$
P(B\geq80,\;G\leq126),
$$

$$
E(X),
$$

$$
(X-150)\mid X\geq150.
$$

### Step 4: Is there dependence or selection?

Look for words such as:

- AND
- given
- among
- survived
- select
- at least

“AND” suggests a joint event. “Given,” “among,” and survivor selection suggest conditioning.

### Step 5: Choose the mathematical rule

| Situation | Rule |
|---|---|
| Both events occur | $P(A\cap C)=P(A)P(C\mid A)$ |
| Independent events | $P(A\cap C)=P(A)P(C)$ |
| Standardize a normal variable | $Z=(X-\mu)/\sigma$ |
| Subtract a constant | $E(X-c)=E(X)-c$ |
| Spread after subtracting a constant | $\mathrm{Var}(X-c)=\mathrm{Var}(X)$ |
| Remaining life after survival | $P(R>r)=S(t+r)/S(t)$ |

### Step 6: Check whether the answer makes sense

For this homework:

- Positive correlation should reduce the “higher blood pressure, lower glucose” joint probability relative to independence.
- Gamma mean and SD must be positive and have units of hours.
- Selecting survivors reduces sample size.
- Subtracting elapsed time cannot produce negative remaining lifetimes after the selection.
- Remaining-life statistics need not equal original statistics minus 150.
- A simulation should approximate theoretical quantities, not reproduce them exactly.

## 8. Where to read without getting lost in the books

Because your syllabus does not specify editions, section titles and index terms are safer than page numbers.

Use Ross or Rice to look up only the topic you need:

| Homework concept | Search in the contents or index |
|---|---|
| Probability as area | Continuous random variables; probability density |
| Blood pressure and glucose together | Joint distributions; bivariate normal |
| Multiplying probabilities | Independence; conditional probability |
| How glucose changes with blood pressure | Conditional distributions; conditional expectation |
| Converting correlation to covariance | Covariance; correlation coefficient |
| Gamma mean and SD | Gamma distribution; expectation; variance |
| Why simulated statistics fluctuate | Sampling distributions; law of large numbers |
| Remaining bulb lifetime | Exponential distribution; memoryless property; survival function |

For short, freely available explanations, Penn State’s lesson on bivariate normal distributions covers the patient problem’s conditional-normal formulas and independence result. The Probability Distribution Explorer gives a compact gamma explanation, including the waiting-time interpretation and mean/variance formulas. [online.stat.psu](https://online.stat.psu.edu/stat414/Lesson21.html)

## 9. Quick self-check

Try answering these without calculation.

1. Why does “AND” not automatically mean multiply the two marginal probabilities?
2. If you subtract 150 from every observation, what happens to its SD?
3. Why can the SD change in part C anyway?
4. Would the remaining-life histogram match the original population distribution if the lifetime model were exponential?
5. Why is a sample mean of 965 compatible with a theoretical mean of 1000?

### Answers

1. Multiplying marginal probabilities requires independence. Otherwise, one factor must be conditional.
2. It stays the same.
3. Part C also removes observations, changing the group being summarized.
4. Yes, at the population-distribution level, because exponential lifetimes are memoryless. Finite-sample histograms would still differ randomly. [dylanspicker.quarto](https://dylanspicker.quarto.pub/the-memoryless-property-of-the-exponential-distribution/)
5. A sample mean varies from sample to sample; it is an estimate, not the population mean itself.