# AudioVisual-causal-inference

AudioVisual Perception Modeling with Dr. Luigi Acerbi and Dr. Wei Ji Ma. 
- **Bold** words are .m filenames.
- *Italicized* words are .mat filenames or directory folder names.

Before running any code, install BADS at the following link for model fitting:
  - Acerbi, L. & Ma, W. J. (2017). Practical Bayesian Optimization for Model Fitting with Bayesian Adaptive Direct Search. In Advances in Neural Information Processing Systems 30, pages 1834-1844. https://github.com/acerbilab/bads

All code assumes that the current directory of Matlab is the main folder.

## Main folder
- **fast_fit_visualzize.m** is an fast example code for fitting and visualizing models.  It serves as a pipeline example for what the other two pieces of code above do. It fits the Exp-GaussianLaplace model with free audiovisual rescale and uniform lapse on UV+UA data (only 1 init per subject), and visualizes the fits. The resulting model fits are saved in a separate subfolder *fast_modefit_example*, to avoid confusion and overwriting of the proper model fits saved in the *modelfits* subfolder. 

- **fit_models.m** contains example model fits for one parametric model on UA+UV data, the semiparametric model on UA+UV data, one semiparamInspired model on all data, and one parametric model on all data. The number of iters per subject in the code is used to generate the saved *.mat* files (model fits) in the *modelfits* folder.
- **manuscript_allplots.m** creates all figures for the manuscript, based on the datafiles and the saved model fits.

## *data* subfolder
- *alldata.mat* is the raw datafile.
- **parse_data.m** creates the datafiles needed *data_stratified_UA.mat*, *data_stratified_UV.mat*, *BC_data.mat*, *BAV_data.mat*, from the raw datafile *alldata.mat*.
- The actual datafiles used for ModelFits and Plots is only the following four files: *data_stratified_UA.mat*, *data_stratified_UV.mat*, *BC_data.mat*, *BAV_data.mat*. They store the 15 human subject's data on UA, UV, BC, BA+BV trials respectively. Each *.mat* file is a cell(1,15) corresponding to the 15 subjects. 
  - For *data_stratified_UV.mat*, *data_stratified_UA.mat*, each of the 15 subjects is additionally stratfied to cell(1,7), corresponding to the 7 stimulus bins by true stimulus location (useful for visualizations purposes). For either file, different rows are trials. Col 1 is the stimulus location s, Col 2 is the subject's response r_S. For UV data, there is an additional Col 3 denoting the visual reliability level (1 corresponding to highest reliability, 3 lowest). To unstratify by the 7 stimulus bins (removing the second-level cell(1,7)'s), use **data_stratified_to_data.m** (which is often used in the model fitting code). 
  - For *BC_data.mat*, the second level is just a matrix containing trial information. Different rows are trials. Col 1 is visual reliability level, Col 2 is the auditory stimulus location s_A, Col 3 is the visual stimulus location s_V, Col 4 is the subject's response r_C (1 is "same", 2 is "different").
  - For *BAV_data.mat*, the second level is also just a matrix containing trial information. Different rows are trials. Col 1 is visual reliability level, Col 2 is the auditory stimulus location s_A, Col 3 is the visual stimulus location s_V, Col 4 is the subject's response r_S, and Col 5 denotes whether the response is a location estimate for the visual (1) or auditory (2) stimulus.

## *modelfits* subfolder
This folder only contains the saved fitted params for each model. For specifics on the filenames, number of iterations per subject, and reserved time on the HPC for these model fits, see below:
```
parametric on UA+UV data: 
Not parallelized for HPC -- each fit (all 15 subjects, 10 inits each) takes around 6-10 hours.
bads.m optimizer
NLLfun_UAV_parametric.m;
'fittedparams_UJoint_'+hetero_type+"-"+prior_type + "_rescale"+rescale + "_lapse"+lapse_type + iter;
Note: it allows the fitting of models with truncated Gaussian lapse (additional parameter sigma_{lapse} involved).

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

## *analysis* subfolder
### Model fitting code 
- **fit_ujointmodel_parametric.m**, **fit_ujointmodel_semiparam.m**, **fit_alldatamodel_semiparaminsp_resc.m**, **fit_alldatamodel_parametric_resc.m** are the model fits code being called by **fit_models.m**. They fit parametric models on UV+UA data, the semiparametric model on UV+UA data, the semiparamInspired models on all the data, and parametric models on all the data respectively. The fitted parameters are saved in the *modelfits* folder.
  - **fit_alldatamodel_semiparaminsp_resc.m** requires the saved fitted parameters from **fit_ujointmodel_semiparam.m**, called *fittedparams_UJoint_Semiparam_rescalefree_lapseUniform.mat*.
- The NLL code for parametric models are is **nllfun_uav_parametric.m**, **nllfun_bc_parametric.m**, and **nllfun_bav_parametric.m**. They each take in a set of parameters and data, and return either their summed NLL (for model fitting), each trial's likelihood (never used; can be an option for visualizing UV data and model fits), or generatives posterior predictive samples using the dataset's stimulus and input parameters (for all model fit visualizations). 
- The NLL code for the semiparametric model is **nllfun_uav_semiparam.m**. Its functions are similar to above. The model parameters are explained below:
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
- The NLL code for semiparamInspired models are **nllfun_uav_semiparaminsp.m**, **nllfun_bc_semiparaminsp.m**, and **nllfun_bav_semiparaminsp.m**. Its functions are similar to above.

### Model visualization code
**manuscript_allplots** in the main folder calls functions in this subfolder. It calls **manuscript_ujoint_respdistrvisualization.m**, **manuscript_ujoint_respdistrvisualization_semiparam.m**, **manuscript_allfits_respdistrvisualization_resc.m**, **manuscript_allfits_respdistrvisualization_semiparaminsp_resc.m** to create model fit response distribution plots. These 4 secondary functions in turn call **manuscript_unimodalfits_visualization.m**, **manuscript_bimodalcFits_visualization_resc.m**, **manuscript_bimodalavFits_visualization_resc.m**, which can be universally used for parametric, semiparam, and semiparamInsp model fits. 

## *plots* subfolder
This folder only saves the manuscript figures created by **manuscript_allplots**. 

## *utils* subfolder
#### Helper functions
- **heterotype_to_sigmafun.m** is for parametric models only. It generates symbolic functions for the sensory noise function sigma(s) from the function family names "constant" or "exp".
- **complete_thetaua_for_ujointFits.m** is for parametric model's NLL on UA data only. It completes the UA parameters by joining UA unique parameters and UA+UV shared parameters. 
- **merge_ujoint_badsbounds.m** and **sigmafun_badsbounds_comprehensive.m** output the parameter bounds for BADS model fitting, according to the parametric model specifications.
- **qtrapz.m** is code for trapezoidal numerical integration, used for computing posterior means. 
  - Acerbi, L. (2019). An Exploration of Acquisition and Mean Functions in Variational Bayesian Monte Carlo. In Proc. Machine Learning Research 96: 1-10. 1st Symposium on Advances in Approximate Bayesian Inference, Montréal, Canada. https://github.com/acerbilab/vbmc

#### Optimizer-related
- **cmaes.m** is a numerical optimizer used for Semiparametric model fitting.
  - N. Hansen, S. D. Müller and P. Koumoutsakos, "Reducing the Time Complexity of the Derandomized Evolution Strategy with Covariance Matrix Adaptation (CMA-ES)," in Evolutionary Computation, vol. 11, no. 1, pp. 1-18, March 2003, doi: 10.1162/106365603321828970. https://cma-es.github.io/
- **cmaes_modded** is a set of CMAES options developed in: 
  - Acerbi, L. (2019). An Exploration of Acquisition and Mean Functions in Variational Bayesian Monte Carlo. In Proc. Machine Learning Research 96: 1-10. 1st Symposium on Advances in Approximate Bayesian Inference, Montréal, Canada. https://github.com/acerbilab/vbmc

#### Visualization aesthetics
- **colorbrewer.m** contribute to visualization aesthetics. 
  - Stephen23 (2023). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.3), GitHub. Retrieved May 31, 2023.
 
#### Truncated-Gaussian sampling
- **trandn.m** generates posterior predictive samples out of the truncated Gaussian lapse distribution. Used for visualization for truncated-Gaussian lapse models. 
  - Zdravko Botev (2023). Truncated Normal Generator (https://www.mathworks.com/matlabcentral/fileexchange/53180-truncated-normal-generator), MATLAB Central File Exchange. Retrieved May 31, 2023.
