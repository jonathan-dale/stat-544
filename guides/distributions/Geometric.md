## The geometric random variable counts failures before the first success.

Example:

In a biology lab, students are checking the eye color for a large number of fruit flies. Suppose that the probability of white eyes is 1/4 and the probability of red eyes is 3/4.

a. What is the probability that the first fruit fly to have white eyes is the fifth fly checked?

b. What is the probability that students need to check at most six fruit flies to observe the first white-eyed specimen?

Solution: 

Part A. 
Here failures are the red eye flies and success is the white eye fly. If Y is the number of red flies before the first white fly then Y ~ Goem(p = 1/4).

If the first white eyed fly is the fifth checked:
 - the first 4 flies are red (failures)
 - the 5th fly is white (success)

$$
P(X = 5) = (3 choose 4)^4 (1 choose 4) = 81/1024 = 0.7910156.
$$


Part B. 
The question states that "at most six flies before the first white" means that the first white-eyed fly appers on check 1 - 6.

P(X <= 6) = 1 - P(X > 6) = 1 - (3/4)^6 = 3367/4096 ~ .8220

So about 82% chance that the researcher needs to check at least 6 flies before success.


```{r}
p_white = 1/4
p_red = 3/4

# a. First white-eyed fly is the fifth checked:
# 4 red flies occur before the first white fly
part_a = dgeom(4, prob = p_white)

## Or calculate it directly
p_red^4 * p_white

# b. First white-eyed fly appears within the first 6 checks:
# At most 5 red flies occur before the first white fly
part_b = pgeom(5, prob = p_white)

part_a
part_a * 100

part_b
part_b * 100
```


