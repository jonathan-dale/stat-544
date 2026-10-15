### A three-layer method
For every distribution, theorem, or R function, connect these three layers.

| Layer | Question to answer | Example: **Binomial** |
| --- | --- | --- |
| Story/process | What real mechanism could generate this? | n independent requests, each succeeds/fails                                                                                           |
| Math/model | What are the assumptions and key properties? | $$X∼Binomial(n,p),  \\E[X]=np, \\Var(X)=np(1−p)$$ |
| R/experiment | How do I calculate and simulate it? | dbinom, pbinom, qbinom, rbinom |



# Distribution card

> Use this checklist as a guide to help determine what distribution the data is most likely to follow.

**Support**: What values are possible?

**Parameters**: What do they mean operationally?

**Generative story**: What process produces it?

**PMF/PDF**: What is the probability model?

**CDF**: What probability question does it answer?

**Mean and variance**: What are the center and spread?

**Assumptions**: What would make this inappropriate?

**R mapping**: `d*`, `p*`, `q*`, and `r*`

**One simulation**: Verify theory numerically and visually.





### A three-layer method
For every distribution, theorem, or R function, connect these three layers.

| Layer | Question to answer | Example: **Binomial** |
| --- | --- | --- |
| Story/process | What real mechanism could generate this? | n independent requests, each succeeds/fails                                                                                           |
| Math/model | What are the assumptions and key properties? | $$X∼Binomial(n,p), \\ E[X]=np, \\ Var(X)=np(1−p)$$ |
| R/experiment | How do I calculate and simulate it? | dbinom, pbinom, qbinom, rbinom |
