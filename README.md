# Bayesian Hierarchical Hypothesis Testing (BHHT)


The function `Bayes_HHT()` takes posterior samples from a Bayesian variable selection model (`B`, effects in columns, samples in rows), a hierarchical clustering of the predictors (i.e., the columns of `X` in `y=Xb+e`), and a Bayesian FDR threshold (`alpha`) and returns credible sets, consisting of predictors that are jointly associated with the outcome. Further details about the methodology can be found in [Samaddar, Maiti, and de los Campos, 2023](). The following example demonstrates the use of the function.


# Examples

The following examples (used in [Samaddar, Maiti, and de los Campos, 2023]()) illustrate the application of Bayesian Hierarchical Hypothesis Testing.

  - [Simple demonstration](https://github.com/AnirbanSamaddar/Bayes_HHT/blob/main/Mice_example/mice_example.md). This example uses a publicly available mice data set and presents a simple example that illustrates the methodology. 
  - [Simulation 1](https://github.com/AnirbanSamaddar/Bayes_HHT/blob/main/Simulation_1/README.md). This simulation uses simulated genotypes and the script is stand-alone (i.e., it does not require having data subject to access restriction).
  - [Simulation 2](https://github.com/AnirbanSamaddar/Bayes_HHT/blob/main/Simulation_2/README.md). This simulation uses real human genotypes from the UK-Biobank.

