function [out_struct] = fit_alldatamodel_semiparaminsp_resc(iter,causal_inf_strategy, num_inits, data_path, model_path, model_path_temp)
if(nargin==0)
    iter=1;
    causal_inf_strategy = "ModelAveraging"; 
    num_inits = 100; data_path = "data\"; model_path = "modelfits\";  model_path_temp = "modelfits\temp";
elseif(nargin==1)
    causal_inf_strategy = "ModelAveraging"; 
    num_inits = 100; data_path = "data\"; model_path = "modelfits\"; model_path_temp = "modelfits\temp";
elseif(nargin==2)
    num_inits = 100; data_path = "data\"; model_path = "modelfits\"; model_path_temp = "modelfits\temp";
elseif(nargin==3)
    data_path = "data\"; model_path = "modelfits\"; model_path_temp = "modelfits\temp";
elseif(nargin==4)
    model_path = "modelfits\"; model_path_temp = "modelfits\temp";
elseif(nargin==5)
    model_path_temp = "modelfits\temp";
end
lapse_type="Uniform"; rescale_aud = "free";

causal_inf_strategy
lapse_type
rescale_aud

rng('default')
rng(iter)

PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
consider_lapse=true; % fit a lapse parameter.

%% Load the fake data for param recov. 
load(data_path+"bav_data.mat")
load(data_path+"bc_data.mat")
load(data_path+"data_stratified_uv.mat");
load(data_path+"data_stratified_ua.mat");
num_subjects = length(data_stratified_UA);
data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
data_UA = data_stratified_to_data(data_stratified_UA, false, false);
UAV_data = cell(1,num_subjects);
Gaussian_lapse_SDs = zeros(1,num_subjects); % If assume Gaussian Lapse, then its SD is just the SD across all UAV trials for that subject. 
for i=1:num_subjects
    data_UA{i}(:,3) = 4;
    UAV_data{i} = [data_UV{i}; data_UA{i}];
    Gaussian_lapse_SDs(i) = std([UAV_data{i}(:,2); BAV_data{i}(:,4)]);
end

%% Load data corresponding to subj specified by iter.
% Get which subject and init to fit.
subjidxs = repelem(1:num_subjects,num_inits);
rand_inits = repmat(1:num_inits, 1,num_subjects)';
subjidx = subjidxs(iter);
init=rand_inits(iter);

%% Load the nonparametrically fitted params, as base shapes for sigma(s) vis, sigma(s) aud, p(s). 
load(model_path + "fittedparams_UJoint_Semiparam_rescalefree_lapseUniform");
theta_fitted_nonparam = theta_fitted;
theta_nonparam = theta_fitted_nonparam(subjidx,:);
F_vals_nonparam = F_vals;
clear T_Ends Theta_fitted theta_inits F_vals theta_fitted

% Extract the 3 function shapes.
s_pivot = [0,0.1,0.3,1,2,4,6,8,10,15,20,45];
s_pivot_full = [-fliplr(s_pivot(2:end)), s_pivot];
num_pivots = length(s_pivot);
sigma_fun_vis_rel_high_pivots = cumsum(theta_nonparam(1:num_pivots));
sigma_fun_aud_pivots = cumsum(theta_nonparam((num_pivots+1):(2*num_pivots)));
sigma_fun_vis_rel_high_pivots = [fliplr(sigma_fun_vis_rel_high_pivots(2:end)), sigma_fun_vis_rel_high_pivots];
sigma_fun_aud_pivots = [fliplr(sigma_fun_aud_pivots(2:end)), sigma_fun_aud_pivots];
prior_pivots = cumsum([1,theta_nonparam((2*num_pivots+1):(3*num_pivots-1))]);
prior_pivots = [fliplr(prior_pivots(2:end)), prior_pivots];

% Get sigma(s) as function handles
sigma_fun_vis_rel_high = @(s)  min([repmat(45,length(s),1) , interp1(s_pivot_full, exp(sigma_fun_vis_rel_high_pivots), s, 'pchip')], [], 2);       
sigma_fun_aud = @(s) min([repmat(45,length(s),1) , interp1(s_pivot_full, exp(sigma_fun_aud_pivots), s, 'pchip')], [], 2);       
sigma_fun_all = {sigma_fun_vis_rel_high;sigma_fun_vis_rel_high;sigma_fun_vis_rel_high;sigma_fun_aud};

% Get unnormalized prior as vector over s_grid_integrate.
a=PMIntegrationParams(1); b=PMIntegrationParams(2); num_bins = PMIntegrationParams(3);
binedges = linspace(a,b,num_bins+1);
binwidth = binedges(2)-binedges(1);
s_grid_integrate = (binedges(1:end-1) + binwidth/2)';
p_s = exp(interp1(s_pivot_full, (prior_pivots), s_grid_integrate, 'pchip'));
%log_prior = log(p_s);

%% Specify modeling assumptions
ModelComponents.CausalInfStrategy = causal_inf_strategy;
ModelComponents.SGridIntegrate=s_grid_integrate; % Use this directly! Do not recreate it from [a,b,num_bins] in NLL function!
ModelComponents.PriorUnnormalized=p_s; % Because we need to apply power, so need to normalize prior.^UB_shared_prior_power IN NLL code before using prior for anything!
ModelComponents.PriorAssumption="independent";
ModelComponents.SigmaFuns=sigma_fun_all;
ModelComponents.LapseRange = [-45,45];
ModelComponents.RescaleVis="1";
ModelComponents.RescaleAud=rescale_aud;

%% LB, UB of free parameters
% [scale_rel_med, scale_rel_low, lapse, sig_motor, rescale_aud, p_same,
% Bimod_rescale_vis, Bimod_rescale_aud, UB_shared_prior_power];

LB = [0.7, 0.7, 0, 0.25, 1/3, 0, 0.5, 0.5, 0.1];
PLB = [1, 1, 0.01, 0.5, 2/3, 0, 0.8, 0.8, 0.5];
PUB = [5, 5, 0.2, 1.5, 1.5, 1, 3, 3, 2];
UB = [20, 20, 1, 2, 3, 1, 5, 5, 10];

% IBS/BADS settings
options = bads('defaults');
options.SpecifyTargetNoise = false;  
options.UncertaintyHandling = false;  


%% Test code for new NLL code.
% nonparam_NLL_BADS = zeros(1,15);
% for subj = 1:15
%     subj
%     theta_nonparam = theta_fitted_nonparam(subj,:);
%     
%     % Extract the 3 function shapes.
%     sigma_fun_vis_rel_high_pivots = cumsum(theta_nonparam(1:num_pivots));
%     sigma_fun_aud_pivots = cumsum(theta_nonparam((num_pivots+1):(2*num_pivots)));
%     sigma_fun_vis_rel_high_pivots = [fliplr(sigma_fun_vis_rel_high_pivots(2:end)), sigma_fun_vis_rel_high_pivots];
%     sigma_fun_aud_pivots = [fliplr(sigma_fun_aud_pivots(2:end)), sigma_fun_aud_pivots];
%     prior_pivots = cumsum([1,theta_nonparam((2*num_pivots+1):(3*num_pivots-1))]);
%     prior_pivots = [fliplr(prior_pivots(2:end)), prior_pivots];
% 
%     % Get sigma(s) as function handles
%     sigma_fun_vis_rel_high = @(s)  min([repmat(45,length(s),1) , interp1(s_pivot_full, exp(sigma_fun_vis_rel_high_pivots), s, 'pchip')], [], 2);       
%     sigma_fun_aud = @(s) min([repmat(45,length(s),1) , interp1(s_pivot_full, exp(sigma_fun_aud_pivots), s, 'pchip')], [], 2);       
%     sigma_fun_all = {sigma_fun_vis_rel_high;sigma_fun_vis_rel_high;sigma_fun_vis_rel_high;sigma_fun_aud};
%     
%     % Get normalized prior as vector over s_grid_integrate.
%     a=PMIntegrationParams(1); b=PMIntegrationParams(2); num_bins = PMIntegrationParams(3);
%     binedges = linspace(a,b,num_bins+1);
%     binwidth = binedges(2)-binedges(1);
%     s_grid_integrate = (binedges(1:end-1) + binwidth/2)';
%     p_s = exp(interp1(s_pivot_full, (prior_pivots), s_grid_integrate, 'pchip'));
%     
%     % Specify modeling assumptions
%     ModelComponents1.CausalInfStrategy = causal_inf_strategy;
%     ModelComponents1.SGridIntegrate=s_grid_integrate; % Use this directly! Do not recreate it from [a,b,num_bins] in NLL function!
%     ModelComponents1.PriorUnnormalized=p_s; % Because we need to apply power, so need to normalize prior.^UB_shared_prior_power IN NLL code before using prior for anything!
%     ModelComponents1.PriorAssumption="independent";
%     ModelComponents1.SigmaFuns=sigma_fun_all;
%     ModelComponents1.LapseRange = [-45,45];
%     ModelComponents1.RescaleVis="1";
%     ModelComponents1.RescaleAud="free";
%     
%     theta0 = [theta_fitted_nonparam(subj,(end-4):end),0.5,1,1,1];
%     UAV_data_subj = UAV_data{subj};
%     BC_data_subj = BC_data{subj};
%     BAV_data_subj = BAV_data{subj};
%     nonparam_NLL_BADS(subj) = NLLfun_UAV_nonparamInspired(theta0([1:5,end]), UAV_data_subj, ModelComponents1, false, false, consider_lapse);    
%     %nonparam_NLL_BADS(subj) = NLLfun_BC_UBresc_nonparamInspired(theta0, BC_data_subj, ModelComponents1, false, false, consider_lapse);    
%     %nonparam_NLL_BADS(subj) = NLLfun_BAV_UBresc_nonparamInspired(theta0, BAV_data_subj, ModelComponents1, false, false, consider_lapse);    
% 
% end
% %save("nonparam_81inits_BADS_NLL",'nonparam_NLL_BADS') 
% %save("nonparam_BC_BADS_NLL",'nonparam_NLL_BADS')
% %save("nonparam_BAV_BADS_NLL",'nonparam_NLL_BADS') 
% 
% [F_min_val_cmaes, min_idx_cmaes] = min(F_vals_nonparam,[],1);
% 
% figure
% hold on
% plot(1:num_subjects, F_min_val_cmaes-nonparam_NLL_BADS, ".-")
% plot(1:num_subjects, zeros(1,15), "k--")
% xlim([0.5,15.5])
% xlabel("Subject")
% ylabel("NLL (bads eval of cmaes best - cmaes best)")
% title("UV+UA")

%% Fit params
[subjidx, init]

% Random initialization
theta0 = rand(size(LB)).*(PUB-PLB) + PLB;

theta0

% BAV data
BAV_data_subj = BAV_data{subjidx};
% S_V_BAV = BAV_data_subj(:,[1,3]); % include reliability level info for each trial.
% S_A_BAV = BAV_data_subj(:,2);
% R_BAV = BAV_data_subj(:,[4,5]); % R is just C in this task type.

% BC data
BC_data_subj = BC_data{subjidx};
% S_V_BC = BC_data_subj(:,[1,3]); % include reliability level info for each trial.
% S_A_BC = BC_data_subj(:,2);
% R_BC = BC_data_subj(:,4); % R is just C in this task type.

% UA, UV data
UAV_data_subj = UAV_data{subjidx};

nllfun = @(theta) nllfun_bav_ubresc_semiparaminsp(theta, BAV_data_subj, ModelComponents, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx))...
+ nllfun_bc_ubresc_semiparaminsp(theta, BC_data_subj, ModelComponents, false, false, consider_lapse)...
+ nllfun_uav_semiparaminsp(theta([1:5,end]), UAV_data_subj, ModelComponents, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx));

tStart = tic;
[Theta_fitted_cell, F_vals_cell] = bads(nllfun,theta0,LB,UB,PLB,PUB,[],options);
T_Ends_cell = toc(tStart); 

out_struct.iter = iter;
out_struct.subjidx = subjidx;
out_struct.init = init;
out_struct.theta_init = theta0;
out_struct.theta_fitted = Theta_fitted_cell;
out_struct.F_vals = F_vals_cell;
out_struct.T_Ends = T_Ends_cell;

%% Save fitted parameters.
filename_basis = 'fittedparams_All_UBresc_SemiparamInspired_'+ModelComponents.CausalInfStrategy;
rescale = ModelComponents.RescaleAud;
if(rescale=="free")
    filename = filename_basis + "_rescalefree";
elseif(rescale=="1")
    filename = filename_basis + "_rescale1";
elseif(rescale=="4/3")
    filename = filename_basis + "_rescale4over3";
end
filename = filename + "_lapse"+lapse_type;
save(model_path_temp+filename+"_"+num2str(iter), 'out_struct', 'ModelComponents');
end
