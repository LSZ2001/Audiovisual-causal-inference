## Setup
Before running any code, go to *Analysis\bads_master* and run **install.m** to install BADS. 

For description of each folder's contents, see below.

## Data
- The actual datafiles used for ModelFits and Plots is only the following four files: *data_stratified_UA.mat*, *data_stratified_UV.mat*, *BC_data.mat*, *BAV_data.mat*. They store the 15 human subject's data on UA, UV, BC, BA+BV trials respectively. Each *.mat* file is a cell(1,15) corresponding to the 15 subjects. 
  - For *data_stratified_UV.mat*, *data_stratified_UA.mat*, each of the 15 subjects is additionally stratfied to cell(1,7), corresponding to the 7 stimulus bins by true stimulus location (useful for visualizations purposes). For either file, different rows are trials. Col 1 is the stimulus location s, Col 2 is the subject's response r_S. For UV data, there is an additional Col 3 denoting the visual reliability level (1 corresponding to highest reliability, 3 lowest). To unstratify by the 7 stimulus bins (removing the second-level cell(1,7)'s), use **data_stratified_to_data.m** (which is often used in the model fitting code). 
  - For *BC_data.mat*, the second level is just a matrix containing trial information. Different rows are trials. Col 1 is visual reliability level, Col 2 is the auditory stimulus location s_A, Col 3 is the visual stimulus location s_V, Col 4 is the subject's response r_C (1 is "same", 2 is "different").
  - For *BAV_data.mat*, the second level is also just a matrix containing trial information. Different rows are trials. Col 1 is visual reliability level, Col 2 is the auditory stimulus location s_A, Col 3 is the visual stimulus location s_V, Col 4 is the subject's response r_S, and Col 5 denotes whether the response is a location estimate for the visual (1) or auditory (2) stimulus.

### How the four data files are generated
 - *alldata.mat* is the raw data files, containing data from all 5 tasks. 
 - **UA_data_visualization.m**, **UV_data_visualization.m**, **BC_data_visualization.m**, **BAV_data_visualization.m** operate on *alldata.mat* to generate *data_stratified_UA.mat*, *data_stratified_UV.mat*, *BC_data.mat*, *BAV_data.mat* respectively.


## ModelFits
This folder only saves the fitted params for each model. For how these models are fitted (Note: the *.m* code are all in the *Analysis* folder instead), see descriptions below.
- In the *Analysis* folder, **fit_UJointModel_parametric.m**, **fit_UJointModel_semiparam.m**, **fit_AllDataModel_semiparamInsp_resc.m**, **fit_AllDataModel_parametric_resc.m** are the improved model fits code. They fit parametric models on UV+UA data, the semiparametric model on UV+UA data, the semiparametric-inspired models on all the data, and parametric models on all the data respectively. The fitted parameters are saved in the *ModelFits* folder.
  - **fit_AllDataModel_semiparamInsp_resc.m** requires the saved fitted parameters from **fit_UJointModel_semiparam.m**.
- In the *Analysis* folder, the NLL code for parametric models are is **NLLfun_UAV_parametric.m**, **NLLfun_BC_parametric.m**, and **NLLfun_BAV_parametric.m**. They each take in a set of parameters and data, and return either their summed NLL (for model fitting), each trial's likelihood (never used; can be an option for visualizing UV data and model fits), or generatives posterior predictive samples using the dataset's stimulus and input parameters (for all model fit visualizations). 
- In the *Analysis* folder, the NLL code for the semiparametric model is **NLLfun_UAV_semiparam.m**. Its functions are similar to above. The model parameters are explained below:
```
% For each subject, conduct a separate fit with the following 40 parameters:
% theta = [sigma_fun_vis_pivot_params, sigma_fun_aud_pivot_params, prior_pivot_params, scale_rel_med, scale_rel_low, lapse, sigma_motor, rescale_aud]

s_pivot = [0,0.1,0.3,1,2,4,6,8,10,15,20,45]; 
num_pivots = length(s_pivot); % = 12
sigma_fun_vis_rel_high_pivots = cumsum(theta(1:num_pivots)); % cumsum of sigma_fun_vis_pivot_params -> monotonically increasing.
sigma_fun_aud_pivots = cumsum(theta((num_pivots+1):(2*num_pivots))); % cumsum of sigma_fun_aud_pivot_params
sigma_fun_vis_rel_high_pivots = [fliplr(sigma_fun_vis_rel_high_pivots(2:end)), sigma_fun_vis_rel_high_pivots]; % assume symmetry across s=0
sigma_fun_aud_pivots = [fliplr(sigma_fun_aud_pivots(2:end)), sigma_fun_aud_pivots];
prior_pivots = cumsum([1,theta((2*num_pivots+1):(3*num_pivots-1))]); % monotonically decreasing
prior_pivots = [fliplr(prior_pivots(2:end)), prior_pivots]; % the pivot at s=0 is always 1, to reduce 1 df. This works due to later normalization in s_hat_PM.

sigma_fun_vis = @(s)  min([repmat(45,length(s),1) , interp1(s_pivot_full, exp(sigma_fun_vis_rel_high_pivots), s, 'pchip')], [], 2);       
sigma_fun_aud = @(s)  min([repmat(45,length(s),1) , interp1(s_pivot_full, exp(sigma_fun_aud_pivots), s, 'pchip')], [], 2);    
p_s = exp(interp1(s_pivot_full, (prior_pivots), s_grid_integrate, 'pchip'));
log_prior = log(p_s);
```
- In the *Analysis* folder, the NLL code for semiparametric inspired models are **NLLfun_UAV_semiparamInsp.m**, **NLLfun_BC_semiparamInsp.m**, and **NLLfun_BAV_semiparamInsp.m**. Its functions are similar to above.

The data files they save follow the naming conventions below. These saved *.mat* files are the only contents in the *ModelFits* folder.
```
parametric on UA+UV data: 
Not parallelized -- each fit takes around 6-10 hours.
bads.m optimizer
NLLfun_UAV_parametric.m;
'fittedparams_UJoint_'+hetero_type+"-"+prior_type + "_rescale"+rescale + "_lapse"+lapse_type + iter;

semiparam on UV+UA data: 
time reserved per iter 10:00:00
cmaes.m optimizer
NLLfun_UAV_semiparam.m
"fittedparams_UJoint_Semiparam_rescalefree_lapseUniform" + iter;

semiparamInsp parametric on all data: 
time reserved per iter 10:00:00
bads.m optimizer
NLLfun_UAV_semiparamInsp.m, NLLfun_BC_semiparamInsp.m, and NLLfun_BAV_semiparamInsp.m
'fittedparams_All_UBresc_SemiparamInspired_'+causal_inf_strategy+"_rescalefree_lapseUniform" + iter;

parametric on all data: 
time reserved per iter 23:00:00
bads.m optimizer
NLLfun_UAV_parametric.m, NLLfun_BC_parametric.m, NLLfun_BAV_parametric.m
'fittedparams_All_UBresc_'+hetero_type+"-"+prior_type+"-"+causal_inf_strategy+"_rescalefree_lapseUniform" + iter;
```

## Plots
- **Manuscript_AllPlots.m** creates all figures for the manuscript. It contains some functions itself, and additionally calls **Manuscript_UJoint_RespDistrVisualization.m**, **Manuscript_UJoint_RespDistrVisualization_semiparam.m**, **Manuscript_AllFits_RespDistrVisualization_resc.m**, **Manuscript_AllFits_RespDistrVisualization_semiparamInsp_resc.m** to create model fit response distribution plots. These secondary functions in turn call **Manuscript_UnimodalFits_Visualization.m**, **Manuscript_BimodalCFits_Visualization_resc.m**, **Manuscript_BimodalAVFits_Visualization_resc.m**, which can be universally used for parametric, nonparamIndv, and nonparamInsp model fits. 
