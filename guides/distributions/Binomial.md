# Binomial Distribution

```{r}
# Let:
# n = num trials
# p = prob of success
# q = (1-p)
# s = success
# f = failure
# 
# P(X = s) = n(choose)s * p^s * q^f
# 
# Binomial probabilities are calculated in R using dbinom() and pbinom()
# dbinom(x=successes, n=number, p=prob)
```

```{r}
# prob X equals 5
dbinom(5, 10, 0.35) # 5 successes out of 10 trials, with p = 0.35

# prob X is less than or equal to 5
pbinom(5, 10, 0.35)

# prob X is greater than 5
pbinom(5, 10, 0.35, lower.tail = FALSE)
```


```
###########
# Example #
###########
# A clinical trial with 8 patients, prob of a good outcome is 0.15. Simulate the trial with 500 replicates.
# Run the code a few times and see that often less than 3 patients respond well to the drug.
# In very few cases, we see 5-6 good responses, this is very unlikely though it can happen.
```

```{r}
n = 8
prob = .15
trials = 500

set.seed(32)


responses.replicates = vector("numeric", n)
number.responses = vector("numeric", trials)

for (k in 1:trials){
  number.responses.replicates = sample(c(0,1), size = n,
                              prob = c(1 - prob, prob),
                              replace = TRUE)
  number.responses[k] = sum(number.responses.replicates)
}
table(number.responses)
barplot(table(number.responses))


```



