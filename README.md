## Setup
Before running any code, go to *Analysis\bads_master* and run **install.m** to install BADS. 

## Figures
- **Manuscript_AllPlots.m** creates all figures for the manuscript. It contains some functions itself, and additionally calls **Manuscript_UJoint_RespDistrVisualization.m**, **Manuscript_UJoint_RespDistrVisualization_semiparam.m**, **Manuscript_AllFits_RespDistrVisualization_resc.m**, **Manuscript_AllFits_RespDistrVisualization_semiparamInsp_resc.m** to create model fit response distribution plots. These secondary functions in turn call **Manuscript_UnimodalFits_Visualization.m**, **Manuscript_BimodalCFits_Visualization_resc.m**, **Manuscript_BimodalAVFits_Visualization_resc.m**, which can be universally used for parametric, nonparamIndv, and nonparamInsp model fits. 

## Model fits
- **fit_UJointModel_parametric.m**, **fit_UJointModel_semiparam.m**, **fit_AllDataModel_semiparamInsp_resc.m**, **fit_AllDataModel_parametric_resc.m** are the improved model fits code. They fit parametric models on UV+UA data, the semiparametric model on UV+UA data, the semiparametric-inspired models on all the data, and parametric models on all the data respectively. 
- The NLL code for parametric models are is **NLLfun_UAV_parametric.m**, **NLLfun_BC_parametric.m**, and **NLLfun_BAV_parametric.m**. They each take in a set of parameters and data, and return either their summed NLL (for model fitting), each trial's likelihood (never used; can be an option for visualizing UV data and model fits), or generatives posterior predictive samples using the dataset's stimulus and input parameters (for all model fit visualizations). 
- The NLL code for the semiparametric model is **NLLfun_UAV_semiparam.m**. Its functions are similar to above.
- The NLL code for semiparametric inspired models are **NLLfun_UAV_semiparamInsp.m**, **NLLfun_BC_semiparamInsp.m**, and **NLLfun_BAV_semiparamInsp.m**. Its functions are similar to above.

- The data files they save follow a new naming convention that is more organized:
```

semiparam on UV+UA data: 
time reserve 10:00:00
NLLfun_UAV_semiparam.m
"fittedparams_UJoint_semiparamIndv_Arescalefree_" + iter;

semiparamInsp parametric on all data: 
time reserve 10:00:00
???, NLLfun_BC_parametric.m, NLLfun_BAV_parametric.m
'fittedparams_All_UBresc_semiparamInspired_'+causal_inf_strategy+"_lapse"+lapse_type+"_Arescalefree_" + iter;

parametric on all data: 
time reserve 23:00:00
NLLfun_UAV_semiparamInsp.m, NLLfun_BC_semiparamInsp.m, and NLLfun_BAV_semiparamInsp.m
'fittedparams_All_UBresc_'+hetero_type+"-"+prior_type+"-"+causal_inf_strategy+"_lapse"+lapse_type+"_Arescalefree_" + iter;
```
