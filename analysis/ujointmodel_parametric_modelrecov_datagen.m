function [fakedata_UV, fakedata_UA] = ujointmodel_parametric_modelrecov_datagen(iter,prior_type,hetero_type,lapse_type, rescale_aud, num_inits, data_path, model_path_temp)
    if(nargin==0)
        iter=1; prior_type = "SingleGaussian"; hetero_type="constant"; lapse_type="Uniform"; rescale_aud = "free";
        num_inits = 10; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==1)
        prior_type = "SingleGaussian"; hetero_type="constant"; lapse_type="Uniform"; rescale_aud = "free";
        num_inits = 10; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==2)
        hetero_type="constant"; lapse_type="Uniform"; rescale_aud = "free";
        num_inits = 10; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif (nargin==3)
        lapse_type="Uniform"; rescale_aud = "free";
        num_inits = 10; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==4)
        rescale_aud = "free";
        num_inits = 10; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==5)
        num_inits = 10; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==6)
        data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==7)
        model_path_temp = "modelfits\temp\";
    end
    model_path = "modelfits\";
    modelrecovdata_path = "data\model_recov\";

    prior_type
    hetero_type
    lapse_type
    rescale_aud

    s_a_range = -15:5:15;
    s_v_range = -20:1:20;

    model_type="PM"; % two gaussians components of the prior both centered at 0. 
    PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
    consider_lapse=true; % fit a lapse parameter. 

    load(data_path+"data_stratified_uv.mat");
    load(data_path+"data_stratified_ua.mat");
    data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
    data_UA = data_stratified_to_data(data_stratified_UA, false, false);

    % Create fake data structures to hold data generated from fitted models
    fakedata_UV = data_UV;
    fakedata_UA = data_UA;

    num_subjects = length(data_stratified_UA);
    UAV_data = cell(1,num_subjects);
    Gaussian_lapse_SDs = zeros(1,num_subjects); % If assume Gaussian Lapse, then its SD is just the SD across all UAV trials for that subject. 
    for i=1:num_subjects
        data_UA{i}(:,3) = 4;
        UAV_data{i} = [data_UV{i}; data_UA{i}];
        Gaussian_lapse_SDs(i) = std(UAV_data{i}(:,2));
    end

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
    ModelComponents_UV.GaussianLapseSDs = Gaussian_lapse_SDs;
    ModelComponents_UV.PMIntegrationParams = PMIntegrationParams;

    sigma_fun = heterotype_to_sigmafun(ModelComponents_UV.SensoryNoise);
    [LB_UV, UB_UV, PLB_UV, PUB_UV] = sigmafun_badsbounds_comprehensive(ModelComponents_UV);
    num_UV_params = length(LB_UV);


    % UA data model components -> same prior/noise model/lapse/motor noise as UV,
    % but different noise params sigma0, k1, k2 and include a rescale.
    ModelComponents_UA = ModelComponents_UV;
    ModelComponents_UA.Rescale = rescale_aud;
    ModelComponents_UA.NumReliabilityLevels = 1; % Get number of reliability levels
    [LB_UA, UB_UA, PLB_UA, PUB_UA] = sigmafun_badsbounds_comprehensive(ModelComponents_UA);

    % Merge the BADS bounds of UV and UA to one vector, UV in front.
    [LB, UB, PLB, PUB, UA_param_keep_idx] = merge_ujoint_badsbounds(LB_UV,UB_UV,PLB_UV,PUB_UV,LB_UA,UB_UA,PLB_UA,PUB_UA,ModelComponents_UA);

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
    num_params = length(LB); % (sigma0, k, sigma_s, lapse, rescale).

    %% Fit params
    dx_max = 0.5;
    subjidxs = repelem(1:num_subjects,num_inits);
    rand_inits = repmat(1:num_inits, 1,num_subjects)';

    %% iter is provided in .sh file (slurm job array)
      % Load fitted human parameters
    filename = 'fittedparams_UJoint_'+hetero_type+"-"+prior_type+"_rescale"+rescale_aud+"_lapse"+lapse_type;
    load(model_path+filename);

    for subjidx =1:num_subjects
        [subjidx]
    
        % UA, UV data
        data_subj_UV = data_UV{subjidx};
        data_subj_UA = data_UA{subjidx};
        S_UV = data_subj_UV(:,[1,3]); % include reliability level info for each trial.
        S_UA = data_subj_UA(:,1);
        R_UV = data_subj_UV(:,2); % add in rescale
        R_UA = data_subj_UA(:,2); % add in rescale
        
        theta = theta_fitted(subjidx,:);
        switch lapse_type
            case "Uniform"
                % UV
                out_UV = nllfun_uav_parametric(ModelComponents_UV, theta(1:num_UV_params),R_UV,S_UV, true, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx));
                % UA
                out_UA = nllfun_uav_parametric(ModelComponents_UA, complete_thetaua_for_ujointfits(theta, UA_param_keep_idx, ModelComponents_UV.Rescale=="free"),R_UA,S_UA, true, false, consider_lapse, lapse_type, Gaussian_lapse_SDs(subjidx));
    
                out_UV = out_UV(:,1);
                out_UA = out_UA(:,1);
                fakedata_UV{subjidx}(:,2) = out_UV;
                fakedata_UA{subjidx}(:,2) = out_UA;
        end
    end

    
    save(modelrecovdata_path+filename+"_modelrecovdata", 'fakedata_UV', 'fakedata_UA', 'theta_fitted', 'Theta_fitted', 'F_vals','T_Ends');
end
