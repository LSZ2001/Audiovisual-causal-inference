clear all; close all;
% cd('C:\Users\liu_s\Audiovisual-causal-inference\analysis')
model_path = "..\modelfits\";
data_path = "..\data\";
num_subjects = 15;

%% Parametric models on UA+UV data
prior_type = "SingleGaussian"; % Can also be "SingleGaussian", "TwoGaussiansBothFixedZero"
hetero_type = "constant"; % Can also be "constant";
lapse_type = "Uniform"; % Can also be "Gaussian";
rescale_aud = "1"; % Can also be "4/3" or "free"l
num_inits_persubj = 10;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_ujointmodel_parametric(iter,prior_type,hetero_type, lapse_type, rescale_aud)
end

%% Semiparam models on UA+UV data
num_inits_persubj = 81;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_ujointmodel_semiparam(iter);
end

%% SemiparamInsp models on all data
causal_inf_strategy = "ProbMatching"; % Can also either be "ModelSelection" or "ModelAveraging"
num_inits_persubj = 100;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_alldatamodel_semiparaminsp_resc(iter,causal_inf_strategy)
end

%% Parametric models on all data
prior_type = "GaussianLaplaceBothFixedZero"; % Can also be "SingleGaussian", "TwoGaussiansBothFixedZero"
hetero_type = "exp"; % Can also be "constant";
causal_inf_strategy = "ProbMatching"; % Can also either be "ModelSelection" or "ModelAveraging"
num_inits_persubj = 10;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_alldatamodel_parametric_resc(iter,prior_type,hetero_type,causal_inf_strategy)
end
