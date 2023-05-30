model_path = "..\ModelFits\";
data_path = "..\Data\";

prior_type = "GaussianLaplaceBothFixedZero";
hetero_type = "exp";
aud_rescale = "free";
lapse_type = "Gaussian";

%%
s_a_range = -15:5:15;
s_v_range = -20:1:20;

load(data_path+"data_stratified_UV.mat");
data_stratified_UV = data_stratified;
load(data_path+"data_stratified_UA.mat");
data_stratified_UA = data_stratified;
% Unstratify data. 
data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
data_UA = data_stratified_to_data(data_stratified_UA, false, false);
num_subjects = length(data_UA);
UAV_data = cell(1,num_subjects);
Gaussian_lapse_SDs = zeros(1,num_subjects); % If assume Gaussian Lapse, then its SD is just the SD across all UAV trials for that subject. 
for i=1:num_subjects
    data_UA{i}(:,3) = 4;
    UAV_data{i} = [data_UV{i}; data_UA{i}];
    Gaussian_lapse_SDs(i) = std(UAV_data{i}(:,2));
end

%%
model_type="PM"; % two gaussians components of the prior both centered at 0. 
PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
consider_lapse=true; % fit a lapse parameter.
randomize_trials=false;

% UV data model components
ModelComponents_UV.PriorType=prior_type;
ModelComponents_UV.SensoryNoise=hetero_type;
ModelComponents_UV.LapseRange = [-45,45];
ModelComponents_UV.IsFixedPriorMean=(model_type=="PM");
ModelComponents_UV.Rescale="1";
ModelComponents_UV.MotorNoise="Gaussian";
ModelComponents_UV.IsMLModel=(model_type=="ML");
ModelComponents_UV.NumReliabilityLevels = 3; % Get number of reliability levels
ModelComponents_UV.LapseType = lapse_type;
ModelComponents_UV.PMIntegrationParams = PMIntegrationParams;
ModelComponents_UV.GaussianLapseSDs = Gaussian_lapse_SDs;
sigma_fun = heterotype_to_sigmafun(ModelComponents_UV.SensoryNoise);
[LB_UV, UB_UV, PLB_UV, PUB_UV] = sigmafun_BADSbounds_comprehensive(ModelComponents_UV);
num_UV_params = length(LB_UV);


% UA data model components -> same prior/noise model/lapse/motor noise as UV,
% but different noise params sigma0, k1, k2 and include a rescale.
ModelComponents_UA = ModelComponents_UV;
ModelComponents_UA.Rescale=aud_rescale;
ModelComponents_UA.NumReliabilityLevels = 1; % Get number of reliability levels
[LB_UA, UB_UA, PLB_UA, PUB_UA] = sigmafun_BADSbounds_comprehensive(ModelComponents_UA);

% Merge the BADS bounds of UV and UA to one vector, UV in front.
[LB, UB, PLB, PUB, UA_param_keep_idx] = merge_UJoint_BADsbounds(LB_UV,UB_UV,PLB_UV,PUB_UV,LB_UA,UB_UA,PLB_UA,PUB_UA,ModelComponents_UA);


num_UA_uniq_params = length(LB) - length(LB_UA);

%% If Gaussian lapse, add a parameter for the Gaussian lapse distribution's width.
if(lapse_type=="Gaussian")
    LB = [LB, 1];
    PLB = [PLB, 5];
    PUB = [PUB, 20];
    UB = [UB, 90];
end

%% IBS/BADS settings
options = bads('defaults');
options.SpecifyTargetNoise = false;  
options.UncertaintyHandling = false;  

%% Fit params for 15 squarequad PM observers
num_inits=10; num_params = length(LB); % (sigma0, k, sigma_s, lapse, rescale).
Theta_fitted = zeros(num_subjects, num_inits, num_params); % Store run times for each random initialization
F_vals = zeros(num_subjects, num_inits); % Store Log Likelihood for each fitted theta. 
T_Ends = zeros(num_subjects, num_inits); % Store Log Likelihood for each fitted theta.

%% Fit params
theta_inits = zeros(num_subjects, num_inits, num_params);
for subjidx=1:num_subjects
    subjidx
    %theta_true = [sigma0_fakesubjs(subjidx), k_fakesubjs(subjidx), sigma_s_fakesubjs(subjidx)];
    data_subj_UV = data_UV{subjidx};
    data_subj_UA = data_UA{subjidx};
    S_UV = data_subj_UV(:,[1,3]); % include reliability level info for each trial.
    S_UA = data_subj_UA(:,1);
    R_UV = data_subj_UV(:,2); % add in rescale
    R_UA = data_subj_UA(:,2); % add in rescale
    
    switch lapse_type
        case "Uniform"
            nllfun = @(theta) NLLfun_UAV_parametric(ModelComponents_UV, theta(1:num_UV_params),R_UV,S_UV, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx)) ...
                    + NLLfun_UAV_parametric(ModelComponents_UA, complete_thetaUA_for_UJointFits(theta, UA_param_keep_idx, ModelComponents_UV.Rescale=="free"),R_UA,S_UA, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx));
        case "Gaussian"
            nllfun = @(theta) NLLfun_UAV_parametric(ModelComponents_UV, theta([1:num_UV_params, length(theta)]),R_UV,S_UV, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx)) ...
                    + NLLfun_UAV_parametric(ModelComponents_UA, [complete_thetaUA_for_UJointFits(theta(1:(end-1)), UA_param_keep_idx, ModelComponents_UV.Rescale=="free"), theta(end)],R_UA,S_UA, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx));
    end
    
    parfor i=1:num_inits
        [subjidx, i]
        theta0 = rand(size(LB)).*(PUB-PLB) + PLB;
        theta_inits(subjidx,i,:) = theta0;
        
        tStart2 = tic;
        [Theta_fitted(subjidx,i,:), F_vals(subjidx,i)] = bads(nllfun,theta0,LB,UB,PLB,PUB,[],options);
        T_Ends(subjidx,i) = toc(tStart2);
    end
end

m = mean(T_Ends(:))
sem = std(T_Ends(:)) ./ sqrt(num_subjects*num_inits);
[min_val, min_idx] = min(F_vals,[],2);
theta_fitted = zeros(num_subjects, num_params);
NLL_fitted = zeros(num_subjects,1);
for i=1:num_subjects
    theta_fitted(i,:) = Theta_fitted(i,min_idx(i),:);
    NLL_fitted(i) = F_vals(i, min_idx(i));
end


% Find best fitted params for each subject, across random initializations.
close all;
%cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Analysis')
filename_basis = 'fittedparams_UJoint_'+hetero_type+"-"+prior_type;
aud_rescale = ModelComponents_UA.Rescale;
if(aud_rescale=="free")
    filename = filename_basis + "_rescalefree";
elseif(aud_rescale=="1")
    filename = filename_basis + "_rescale1";
elseif(aud_rescale=="4/3")
    filename = filename_basis + "_rescale4over3";
end
filename = filename+"_lapse"+lapse_type;
theta_fitted = zeros(num_subjects, num_params);
[min_val, min_idx] = min(F_vals, [], 2);

for i=1:num_subjects
    theta_ibs_subj = squeeze(Theta_fitted(i,min_idx(i),:)); 
    theta_ibs_subj(theta_ibs_subj==0) = NaN;
    theta_fitted(i,:) = theta_ibs_subj;
end
save(model_path + filename, 'T_Ends', 'Theta_fitted', 'F_vals', 'theta_fitted','theta_inits');
