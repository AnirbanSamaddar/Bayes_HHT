# Bayes_HHT
Implementation of Bayesian Hierarchical Hypothesis Testing 

# Replication
To replicate plots from Motivation, run the following code snippet.
```
$cd src
$Rscript Motivation.R
```
The output will be stored in the folder ```~\figures\Motivation```.
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
