## Setup
Before running any code, go to *Analysis\bads_master* and run **install.m** to install BADS. 

## Data
- The actual datafiles used for ModelFits and Plots is only the following four files: *data_stratified_UV.mat*, *data_stratified_UA.mat*, *BC_data.mat*, *BV_data.mat*. They store the 15 human subject's data on UV. UA, BC, BA+BV trials respectively. Each *.mat* file is a cell(1,15) corresponding to the 15 subjects. 
  - For *data_stratified_UV.mat*, *data_stratified_UA.mat*, each of the 15 subjects is additionally stratfied to cell(1,7) on a second level, stratified into 7 stimulus bins by true stimulus location (useful for visualizations purposes). 
  - For *BC_data.mat*, the second level is just a matrix containing trial information. Different rows are trials. Col 1 is visual reliability level, Col 2 is the auditory stimulus location s_A, Col 3 is the visual stimulus location s_V, Col 4 is the subject's response r_C (1 is "same", 2 is "different").
  - For *BV_data.mat*, the second level is also just a matrix containing trial information. Different rows are trials. Col 1 is visual reliability level, Col 2 is the auditory stimulus location s_A, Col 3 is the visual stimulus location s_V, Col 4 is the subject's response r_S, and Col 5 denotes whether the response is a location estimate for the visual (1) or auditory (2) stimulus.

## ModelFits
- **fit_UJointModel_parametric.m**, **fit_UJointModel_semiparam.m**, **fit_AllDataModel_semiparamInsp_resc.m**, **fit_AllDataModel_parametric_resc.m** are the improved model fits code. They fit parametric models on UV+UA data, the semiparametric model on UV+UA data, the semiparametric-inspired models on all the data, and parametric models on all the data respectively. 
- The NLL code for parametric models are is **NLLfun_UAV_parametric.m**, **NLLfun_BC_parametric.m**, and **NLLfun_BAV_parametric.m**. They each take in a set of parameters and data, and return either their summed NLL (for model fitting), each trial's likelihood (never used; can be an option for visualizing UV data and model fits), or generatives posterior predictive samples using the dataset's stimulus and input parameters (for all model fit visualizations). 
- The NLL code for the semiparametric model is **NLLfun_UAV_semiparam.m**. Its functions are similar to above.
- The NLL code for semiparametric inspired models are **NLLfun_UAV_semiparamInsp.m**, **NLLfun_BC_semiparamInsp.m**, and **NLLfun_BAV_semiparamInsp.m**. Its functions are similar to above.

The data files they save follow the naming conventions below:
```
semiparam on UV+UA data: 
time reserve 10:00:00
NLLfun_UAV_semiparam.m
"fittedparams_UJoint_semiparamIndv_Arescalefree_" + iter;

semiparamInsp parametric on all data: 
time reserve 10:00:00
NLLfun_UAV_parametric.m, NLLfun_BC_parametric.m, NLLfun_BAV_parametric.m
'fittedparams_All_UBresc_semiparamInspired_'+causal_inf_strategy+"_lapse"+lapse_type+"_Arescalefree_" + iter;

parametric on all data: 
time reserve 23:00:00
NLLfun_UAV_semiparamInsp.m, NLLfun_BC_semiparamInsp.m, and NLLfun_BAV_semiparamInsp.m
'fittedparams_All_UBresc_'+hetero_type+"-"+prior_type+"-"+causal_inf_strategy+"_lapse"+lapse_type+"_Arescalefree_" + iter;
```

## Plots
- **Manuscript_AllPlots.m** creates all figures for the manuscript. It contains some functions itself, and additionally calls **Manuscript_UJoint_RespDistrVisualization.m**, **Manuscript_UJoint_RespDistrVisualization_semiparam.m**, **Manuscript_AllFits_RespDistrVisualization_resc.m**, **Manuscript_AllFits_RespDistrVisualization_semiparamInsp_resc.m** to create model fit response distribution plots. These secondary functions in turn call **Manuscript_UnimodalFits_Visualization.m**, **Manuscript_BimodalCFits_Visualization_resc.m**, **Manuscript_BimodalAVFits_Visualization_resc.m**, which can be universally used for parametric, nonparamIndv, and nonparamInsp model fits. 
