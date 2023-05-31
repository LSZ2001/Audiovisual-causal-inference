## Setup
Before running any code, go to Analysis\bads_master and run **install.m** to install BADS. 

## Figures
- **Manuscript_AllPlots.m** creates all figures for the manuscript. It contains some functions itself, and additionally calls **Manuscript_UJoint_RespDistrVisualization.m**, **Manuscript_UJoint_RespDistrVisualization_semiparam.m**, **Manuscript_AllFits_RespDistrVisualization_resc.m**, **Manuscript_AllFits_RespDistrVisualization_semiparamInsp_resc.m** to create model fit response distribution plots. These secondary functions in turn call **Manuscript_UnimodalFits_Visualization.m**, **Manuscript_BimodalCFits_Visualization_resc.m**, **Manuscript_BimodalAVFits_Visualization_resc.m**, which can be universally used for parametric, nonparamIndv, and nonparamInsp model fits. 
## Model fits
- **fit_UJointModel_parametric.m**, **fit_UJointModel_semiparam.m**, **fit_AllDataModel_semiparamInsp_resc.m**, **fit_AllDataModel_parametric_resc.m** are the improved model fits code. They are largely identical to **UJointFits_PMPMmufree_Heterosk_GPC.m**, **fitAllData_nonparam_indv_dir_CMAES_GPC_flat.m**, **fit_AllDataModel_jobArray_newresc_nonparaminspired.m**, **fit_AllDataModel_jobArray_newresc.m** respectively, but renamed for future data-sharing. 
- The NLL code for parametric (previously called "nonparametric inspired" fits) models are is **???**, **NLLfun_BC_parametric.m**, and **NLLfun_BAV_parametric.m**, which are identical to **midpoint_postmean_NLLfun_comprehensive.m**, **midpoint_postmean_NLLfun_comprehensive_BAV_newresc.m**, and **midpoint_postmean_NLLfun_comprehensive_BAV_newresc.m** respectively but renamed properly
- The NLL code for the semiparametric (previously called "nonparametric") model is **NLLfun_UAV_semiparam.m**, which is identical to **midpoint_NLLfun_comprehensive_nonparam_mixture_indv_CMAES.m** but renamed properly.
- The NLL code for semiparametric inspired (previously called "nonparametric inspired") models are **NLLfun_UAV_semiparamInsp.m**, **NLLfun_BC_semiparamInsp.m**, and **NLLfun_BAV_semiparamInsp.m**, which are identical to **NLLfun_UAV_nonparamInspired.m**, **NLLfun_BC_nonparamInspired.m**, and **NLLfun_BAV_nonparamInspired.m** respectively but renamed properly.

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
