# Bayesian Hierarchical Hypothesis Testing (BHHT)


The function `Bayes_HHT()` takes posterior samples from a Bayesian variable selection model (`B`, effects in columns, samples in rows), a hierarchical clustering of the predictors (i.e., the columns of `X` in `y=Xb+e`), and a Bayesian FDR threshold (`alpha`) and returns credible sets, constiting on predictors that are jointly associated with the outcome. Further details about the methodology can be found in (Samaddar, Maiti, and de los Campos, 2023)[]. The following example demonstrate the use of the function.


# Reference


