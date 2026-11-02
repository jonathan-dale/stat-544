# Convergences 


#### Modes of convergences
- Convergence in Distribution
- Convergence in Probability
- Almost Sure Convergence

Probability theory defines several ways a sequence of random variables can converge, ordered from strongest to weakest:
1. Almost Sure Convergence $X_n \xrightarrow{a.s.} X$:  
    The sequence of random variables equals the limit value with probability 1 as $n \to \infty$. (Used in the Strong Law of Large Numbers).
2. Convergence in Probability $X_n \xrightarrow{p} X$:  
    The chance that $X_{n}$ deviates from $X$ by more than a tiny amount $\epsilon$ shrinks to zero. (Used in the Weak Law of Large Numbers).
3. Convergence in Law / Distribution $X_n \xrightarrow{d} X$:  
    The weakest major form of convergence; only the probability distributions need to match up in the limit.

Key Hierarchy Rules:

• Almost sure convergence $\implies $ Convergence in probability $\implies $ Convergence in distribution.
• The reverse implications are not true in general.
• Special exception: If a sequence converges in distribution to a constant $c$, then it also converges in probability to $c$ $(X_n \xrightarrow{d} c \implies X_n \xrightarrow{p} c$).





### Convergence in Law/Distribution:

• Definition: A sequence of random variables $X_1, X_2, \dots, X_n$ converges in law (or in distribution) to a random variable $X$ if the sequence of their cumulative distribution functions (CDFs) $F_{X_n}(x)$ converges pointwise to the CDF $F_X(x)$ of the limit distribution at every single point $x$ where $F_{X}$ is continuous.
• Notation: Written as $X_n \xrightarrow{d} X$ or $X_n \xrightarrow{law} X$.
• Core Idea: You do not care if the individual random variables track each other closely; you only care that their overall shape or probability distribution settles down and looks like the target distribution as $n \to \infty$.






Lesson 1: Convergence in Probability $X_n \xrightarrow{p} X$

With a regular sequence like $\frac{1}{n}$, the numbers are fixed. With random variables, the sequence items $X_1, X_2, X_3, \dots$ yield different values every time you run the experiment.
Because we cannot guarantee a random variable will never throw a wild, unpredictable value, we change the game. Instead of demanding that the values themselves get perfectly trapped inside the $\epsilon$ margin, we demand that the probability of the variable escaping the margin shrinks to zero.

The Formal Definition

A sequence of random variables $X_{n}$ converges in probability to $X$ if for any tiny error margin $\epsilon > 0$ :

$$\lim _{n\rightarrow \infty }P(|{}X_{n}-X|{}\ge \epsilon )=0$$

Breaking It Down (Using Your Anchor)

Think of this as a regular sequence of numbers—but the sequence we are tracking is a sequence of probabilities.
1. The Error Box: You pick a tiny margin $\epsilon$ (say, $0.01$). You look at the probability that $X_{n}$ lands outside this safety box away from $X$.
2. The Shrinking Sequence: Let's call the probability of escaping the box $p_n = P(\vert{}X_n - X\vert{} \geq \epsilon)$.
3. The Trap: As \(n\) grows, this probability sequence $p_{n}$ acts exactly like $\frac{1}{n}$. It shrinks, and shrinks, and settles down to exactly $0$.
The Intuition: As $n$ gets larger, it becomes mathematically near-impossible for $X_{n}$ to freak out and land far away from $X$.