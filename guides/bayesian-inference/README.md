# Bayesian Inference

[Bayesian Inference](Bayesian_Inference.md) is a single-part STAT 544 study guide that builds from probability and statistics prerequisites to Bayes updates, beta-binomial conjugacy, posterior summaries, credible intervals, and loss-based estimation. It follows the supplied Fall 2026 lecture slides and includes explicit derivations, seven worked examples, seventeen practice problems with solutions, and a compact formula reference.

Clearly labeled extensions connect the lecture to posterior prediction, Poisson-gamma and normal-normal inference, prior sensitivity, improper priors, parameter transformations, and point-null comparisons. The explanations assume graduate-level goals without assuming prior undergraduate statistics coursework.

## R verification

The base-R code in the guide is also provided in [Bayesian_Inference_checks.R](Bayesian_Inference_checks.R). Run the complete checks and simulations from this directory with:

```sh
Rscript --vanilla Bayesian_Inference_checks.R
```

No external R packages are required. Deterministic checks compare formulas, numerical integrals, and alternate representations; simulations report Monte Carlo standard errors rather than asserting exact random equality.

## Verification record

Verified on November 8, 2026: all 113 deterministic assertions passed in a clean R 4.5.2 session, and the separate simulations completed with seed 544 and 200,000 draws per model. All mathematical expressions passed KaTeX syntax checks, local navigation and companion-file links resolved, and code fences were balanced.

The companion script matches the five R blocks in the guide exactly, apart from explanatory script headers. Syntax and numerical checks are not a proof of model appropriateness, and GitHub's own renderer was not tested.
