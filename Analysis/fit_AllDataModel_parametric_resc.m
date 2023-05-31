function [out_struct] = fit_alldatamodel_parametric_resc(iter,prior_type,hetero_type,causal_inf_strategy, num_inits, data_path, model_path_temp)
    if(nargin==0)
        iter=1; prior_type = "SingleGaussian"; hetero_type="constant"; causal_inf_strategy="ModelAveraging";
        num_inits = 10; data_path = "..\data\"; model_path_temp = "..\modelfits\temp\";
    elseif(nargin==1)
         prior_type = "SingleGaussian"; hetero_type="constant"; causal_inf_strategy="ModelAveraging";
         num_inits = 10; data_path = "..\data\"; model_path_temp = "..\modelfits\temp\";
    elseif(nargin==2)
         hetero_type="constant"; causal_inf_strategy="ModelAveraging";
         num_inits = 10; data_path = "..\data\"; model_path_temp = "..\modelfits\temp\";
    elseif (nargin==3)
        causal_inf_strategy="ModelAveraging";
        num_inits = 10; data_path = "..\data\"; model_path_temp = "..\modelfits\temp\";
    elseif(nargin==4)
        num_inits = 10; data_path = "..\data\"; model_path_temp = "..\modelfits\temp\";
    elseif(nargin==5)
        data_path = "..\data\"; model_path_temp = "..\modelfits\temp\";
    elseif(nargin==6)
        model_path_temp = "..\modelfits\temp\";
    end
    lapse_type="Uniform"; rescale_aud = "free";
    
    prior_type
    hetero_type
    causal_inf_strategy
    lapse_type
    rescale_aud
    prior_assumption = "independent" % independent or empirical.

    s_a_range = -15:5:15;
    s_v_range = -20:1:20;

    model_type="PM"; % two gaussians components of the prior both centered at 0. 
    PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
    consider_lapse=true; % fit a lapse parameter. 

    % Load the fake data for param recov. 
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

    % UV data model components
    ModelComponents_V.CausalInfStrategy = causal_inf_strategy;
    ModelComponents_V.PriorType=prior_type;
    ModelComponents_V.PriorAssumption=prior_assumption;
    ModelComponents_V.SensoryNoise=hetero_type;
    ModelComponents_V.LapseRange = [-45,45];
    ModelComponents_V.IsFixedPriorMean=(model_type=="PM");
    ModelComponents_V.Rescale="1";
    ModelComponents_V.MotorNoise="Gaussian";
    ModelComponents_V.IsMLModel=(model_type=="ML");
    ModelComponents_V.NumReliabilityLevels = 3; % Get number of reliability levels
    ModelComponents_V.LapseType = lapse_type;
    ModelComponents_V.GaussianLapseSDs = Gaussian_lapse_SDs;
    ModelComponents_V.PMIntegrationParams = PMIntegrationParams;

    sigma_fun = heterotype_to_sigmafun(ModelComponents_V.SensoryNoise);
    [LB_V, UB_V, PLB_V, PUB_V] = sigmafun_badsbounds_comprehensive(ModelComponents_V);
    num_V_params = length(LB_V);


    % UA data model components -> same prior/noise model/lapse/motor noise as UV,
    % but different noise params sigma0, k1, k2 and include a rescale.
    ModelComponents_A = ModelComponents_V;
    ModelComponents_A.Rescale = rescale_aud;
    ModelComponents_A.NumReliabilityLevels = 1; % Get number of reliability levels
    [LB_A, UB_A, PLB_A, PUB_A] = sigmafun_badsbounds_comprehensive(ModelComponents_A);

    % Merge the BADS bounds of UV and UA to one vector, UV in front.
    [LB, UB, PLB, PUB, A_param_keep_idx] = merge_ujoint_badsbounds(LB_V,UB_V,PLB_V,PUB_V,LB_A,UB_A,PLB_A,PUB_A,ModelComponents_A);

    % Add in final parameter constraints -- p_same, which is just Pr[C=1], and
    % then the unimodal-bimodal rescale. 
    LB = [LB, 0, 0.7, 0.7];
    PLB = [PLB, 0, 1, 1];
    PUB = [PUB, 1, 2, 2];
    UB = [UB, 1, 3, 3];
    num_A_uniq_params = length(LB) - length(LB_A) - 2;


    % IBS/BADS settings
    options = bads('defaults');
    options.SpecifyTargetNoise = false;  
    options.UncertaintyHandling = false;  

    %% Fit params for 15 squarequad PM observers
    num_params = length(LB); % (sigma0, k, sigma_s, lapse, rescale).

    %% Fit params
    dx_max = 0.5;
    subjidxs = repelem(1:num_subjects,num_inits);
    rand_inits = repmat(1:num_inits, 1,num_subjects)';

    %% iter is provided in .sh file (slurm job array)
        subjidx = subjidxs(iter);
        i=rand_inits(iter);   
        [subjidx, i]

        % BAV data
        BAV_data_subj = BAV_data{subjidx};
        S_V_BAV = BAV_data_subj(:,[1,3]); % include reliability level info for each trial.
        S_A_BAV = BAV_data_subj(:,2);
        R_BAV = BAV_data_subj(:,[4,5]); % R is just C in this task type.

        % BC data
        BC_data_subj = BC_data{subjidx};
        S_V_BC = BC_data_subj(:,[1,3]); % include reliability level info for each trial.
        S_A_BC = BC_data_subj(:,2);
        R_BC = BC_data_subj(:,4); % R is just C in this task type.

        % UA, UV data
        data_subj_UV = data_UV{subjidx};
        data_subj_UA = data_UA{subjidx};
        S_UV = data_subj_UV(:,[1,3]); % include reliability level info for each trial.
        S_UA = data_subj_UA(:,1);
        R_UV = data_subj_UV(:,2); % add in rescale
        R_UA = data_subj_UA(:,2); % add in rescale
        nllfun = @(theta) nllfun_bav_parametric(ModelComponents_V, ModelComponents_A, theta,R_BAV,S_V_BAV,S_A_BAV, PMIntegrationParams, false, false, consider_lapse, dx_max, ModelComponents_V.CausalInfStrategy, lapse_type, Gaussian_lapse_SDs(subjidx))...
        + nllfun_bc_parametric(ModelComponents_V, ModelComponents_A, theta,R_BC,S_V_BC,S_A_BC, false, false, consider_lapse, ModelComponents_V.CausalInfStrategy)...
        + nllfun_uav_parametric(ModelComponents_V, theta(1:num_V_params),R_UV,S_UV, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx))...
        + nllfun_uav_parametric(ModelComponents_A, complete_thetaua_for_ujointfits(theta(1:(end-3)), A_param_keep_idx, ModelComponents_V.Rescale=="free"),R_UA,S_UA, false, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx));

        theta0 = rand(size(LB)).*(PUB-PLB) + PLB;
        tStart2 = tic;
        [Theta_fitted_cell, F_vals_cell] = bads(nllfun,theta0,LB,UB,PLB,PUB,[],options);
        T_Ends_cell = toc(tStart2); 
 
    model_spec.prior_type = prior_type;
    model_spec.prior_assumption = prior_assumption; % independent or empirical.
    model_spec.hetero_type = hetero_type;
    model_spec.causal_inf_strategy = causal_inf_strategy; % ModelAveraging, ModelSelection, or ProbMatching.
    model_spec.Arescale = ModelComponents_A.Rescale;
    model_spec.LapseType = lapse_type;

    out_struct.iter = iter;
    out_struct.subjidx = subjidx;
    out_struct.model_specs = model_spec;
    out_struct.randinit = i;
    out_struct.theta_init = theta0;
    out_struct.theta_fitted = Theta_fitted_cell;
    out_struct.F_vals = F_vals_cell;
    out_struct.T_Ends = T_Ends_cell;

    %% Find best fitted params for each subject, across random initializations.
    close all;
    filename_basis = 'fittedparams_All_UBresc_'+hetero_type+"-"+prior_type+"-"+causal_inf_strategy;
    rescale = ModelComponents_A.Rescale;
    if(rescale=="free")
        filename = filename_basis + "_rescalefree";
    elseif(rescale=="1")
        filename = filename_basis + "rescale1";
    elseif(rescale=="4/3")
        filename = filename_basis + "_rescale4over3";
    end
    filename = filename + "_lapse"+lapse_type;
    save(model_path_temp+filename+"_"+num2str(iter), 'out_struct','model_spec');
end
