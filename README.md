# Bayesian Hierarchical Hypothesis Testing (BHHT)


The function `Bayes_HHT()` takes posterior samples from a Bayesian variable selection model (`B`, effects in columns, samples in rows), a hierarchical clustering of the predictors (i.e., the columns of `X` in `y=Xb+e`), and a Bayesian FDR threshold (`alpha`) and returns credible sets, constiting on predictors that are jointly associated with the outcome. Further details about the methodology can be found in (Samaddar, Maiti, and de los Campos, 2023)[]. The following example demonstrate the use of the function.


# Reference

# Work order

~~- Add the demonstrating example~~
~~- Change the test function~~
~~- Make necessary changes to the example~~
~~- Add Simulation 1 and Simulation 2 (using mice)~~
- Create readme or markdown files and demonstrate one setting for both simulation 1 and 2
- Fill the texts on the landing page
- Modify the output of Bayes_HHT
   
# Bayes_HHT
Implementation of Bayesian Hierarchical Hypothesis Testing 

In the simualation study, we have run each setting for 1000 MC reps. The following code is for demonstration where we run 10 MC reps for the setting -- ```n=10K S=20 r=0.9 h2=1.25%```.
```
$cd src
$core=1
$Rscript Simulation_study.R ${core}
```
The variable ```core``` defines the setting in ```Simulation_study.R``` file. Changing ```core=2``` will run another predefined setting where ```n=50K```. Please change the parameters in the lines 49-54 in ```Simulation_study.R``` file to run different settings. The output from the simulation study setting will be stored in the folder ```~\output```. 

To generate plots similar to Figure 3, please run the following code.
```
$cd src
$resolution='c(1,5)'
$core=1
$Rscript Make_plot.R ${resolution} ${core}
```
Note that to generate the power-FDR plot like Figure 3 one needs to fix the resolution. In the above code, change the array variable ```resolution``` to generate Figure 3 for different resolutions. The output will be stored in the folder ```$HOME'/figures/Setting_'${core}```.
