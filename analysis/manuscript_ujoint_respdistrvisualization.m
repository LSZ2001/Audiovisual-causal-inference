% Same as ....oversubj2, but plotting stuff in subplots.

% This function uses "data_stratified_UA.mat" generated by
% "UnimodalRawDataPlots_Aud_Function.m", and "fitted_params_PM_UA.mat"
% generated by "UnimodalFits_ML_PM_Homosk.m".
function [] = manuscript_ujoint_respdistrvisualization(prior_type, hetero_type, rescale_text, fontsize, figspecs, model_path, plot_lapse, lapse_type, plot_individual, connect_human_errorbars)
    if(nargin==0)
        figsize = get(0, 'ScreenSize');fontsize = 9; figspecs = [0 0 figsize(4)*0.9*0.65 figsize(4)*0.9];
        model_path = "modelfits\";
        plot_lapse = true; lapse_type = "Uniform";  plot_individual=false; connect_human_errorbars=false;
    elseif(nargin==6) % Default is plot the lapse component in fitted response distributions
        plot_lapse = true; lapse_type = "Uniform";  plot_individual=false; connect_human_errorbars=false;
    elseif(nargin==7)
        lapse_type = "Uniform";  plot_individual=false; connect_human_errorbars=false;
    elseif(nargin==8)
        plot_individual=false; connect_human_errorbars=false;
    elseif(nargin==9)
        connect_human_errorbars=false;
    end
    data_path = "data\";
    not_plot_model_ribbons = (prior_type=="NaN");

    
    % For UAV plots, in case screen too small, need to shrink size. 
    set(0,'units','pixels');
    screensize = get(0, 'ScreenSize');
    screen_max_height_ratio = screensize(4)/figspecs(4);
    
    s_a_range = -15:5:15;
    s_v_range = linspace(-20,20,8);
    s_v_range = (s_v_range(2:end) - s_v_range(1:(end-1)))./2 + s_v_range(1:(end-1));
    load(data_path+"bav_data.mat")
    load(data_path+"bc_data.mat")
    load(data_path+"data_stratified_uv.mat");
    load(data_path+"data_stratified_ua.mat");
    % Unstratify data. 
    data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
    data_UA = data_stratified_to_data(data_stratified_UA, false, false);
    num_subjects = length(data_UA);
    UAV_data = cell(1,num_subjects);
    Gaussian_lapse_SDs = zeros(1,num_subjects); % If assume Gaussian Lapse, then its SD is just the SD across all UAV trials for that subject. 
    %Gaussian_lapse_SDs_bimod = zeros(1,num_subjects);
    for i=1:num_subjects
        data_UA{i}(:,3) = 4;
        UAV_data{i} = [data_UV{i}; data_UA{i}];
        Gaussian_lapse_SDs(i) = std(UAV_data{i}(:,2));
        %Gaussian_lapse_SDs(i) = std([UAV_data{i}(:,2); BAV_data{i}(:,4)]);
    end

    %%
    model_type="PM"; % two gaussians components of the prior both centered at 0. 
    PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
    randomize_trials=false;
    mu = 0; % prior distribution mean
    [~, num_subjects] = size(data_stratified_UA);

    % UV data model components
    ModelComponents_UV.PriorType=prior_type;
    ModelComponents_UV.SensoryNoise=hetero_type;
    ModelComponents_UV.LapseRange = [-45,45];
    ModelComponents_UV.IsFixedPriorMean=(model_type=="PM");
    ModelComponents_UV.Rescale="1";
    ModelComponents_UV.MotorNoise="Gaussian";
    ModelComponents_UV.IsMLModel=(model_type=="ML");
    ModelComponents_UV.NumReliabilityLevels = 3; % Get number of reliability levels
    ModelComponents_UV.PMIntegrationParams = PMIntegrationParams; 

    ModelComponents_UA = ModelComponents_UV;
    ModelComponents_UA.Rescale=rescale_text;
    ModelComponents_UA.NumReliabilityLevels = 1; % Get number of reliability levels


    %% Load fitted params
    if(not_plot_model_ribbons) % Only plot data visualizations
        fitted_params_PM_UV = NaN(num_subjects,1);
        fitted_params_PM_UA = NaN(num_subjects,1);
        
    else % Plot data and one model's predictions

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
        load(model_path + filename)

    %     theta_fitted_old = theta_fitted;
    %     [min_val_old, min_idx_old] = min(F_vals,[],2);

    % if(lapse_type == "Gaussian")
    %     load("fittedparams_UJoint_exp-GaussianLaplaceBothFixedZero_lapseGaussian_rescalefree.mat");
    % end
    %   load("fittedparams_UJoint_exp-GaussianLaplaceBothFixedZero_lapseGaussian_rescalefree.mat");
    % %     %load("fittedparams_UJoint_GaussianLaplaceBothFixedZeroPMexpdnoiseGaussian_rescalefree_lapseGaussian.mat");
    % %     load("fittedparams_UJoint_GaussianLaplaceBothFixedZeroPMexpdnoiseGaussianlapseGaussian_rescalefree.mat")
    %     figure
    %     plot(theta_fitted_old(:,7), theta_fitted(:,7),"k*")
    %     hold on
    %     plot([0,0.045],[0,0.045],"b--")
    %     xlabel("lapse rate of fitted model assuming Uniform lapse")
    %     ylabel("lapse rate of fitted model assuming Gaussian lapse")

        [~,~,num_params]=size(Theta_fitted);
        [min_val, min_idx] = min(F_vals,[],2);
        theta_fitted = zeros(num_subjects, num_params);
        NLL_fitted = zeros(num_subjects,1);
        for i=1:num_subjects
            theta_fitted(i,:) = Theta_fitted(i,min_idx(i),:);
            NLL_fitted(i) = F_vals(i, min_idx(i));
        end
        fitted_params_PM = theta_fitted;
        nll_PM = F_vals;

    %     %% Test
    %     NLLs = [sum(min_val_old); sum(min_val)];
    %     AICs = [2*(sum(min_val_old) + 14); 2*(sum(min_val) + 15)];
    %     load("NumTrials_allsubjects")
    %     n_data_UAV = n_data_bytasktype(4,:);
    %     BICs_subj = [min_val_old'; min_val'];
    %     BICs_subj = 2.*BICs_subj + [14;15].*n_data_UAV;
    %     BICs = sum(BICs_subj,2);
    %     modelselection = [NLLs, AICs, BICs];


        %% Plots
        [LB_UA, UB_UA, PLB_UA, PUB_UA] = sigmafun_badsbounds_comprehensive(ModelComponents_UA);
        [LB_UV, UB_UV, PLB_UV, PUB_UV] = sigmafun_badsbounds_comprehensive(ModelComponents_UV);

        [~,~,~,~,UA_param_keep_idx] = merge_ujoint_badsbounds(LB_UV,UB_UV,PLB_UV,PUB_UV,LB_UA,UB_UA,PLB_UA,PUB_UA, ModelComponents_UA);

        %% Get fitted params for UV and UA separately.
        switch lapse_type
            case "Uniform"
                fitted_params_PM_UV = fitted_params_PM(:, 1:(end-length(UA_param_keep_idx)));
                fitted_params_PM_UA = [];
                for i=1:num_subjects
                    fitted_params_PM_UA = [fitted_params_PM_UA; complete_thetaua_for_ujointfits(theta_fitted(i,:), UA_param_keep_idx, ModelComponents_UV.Rescale=="free")];
                end
            case "Gaussian"
                fitted_params_PM_UV = fitted_params_PM(:, [1:(end-length(UA_param_keep_idx)-1),end]);
                fitted_params_PM_UA = [];
                for i=1:num_subjects
                    fitted_params_PM_UA = [fitted_params_PM_UA; [complete_thetaua_for_ujointfits(theta_fitted(i,1:(end-1)), UA_param_keep_idx, ModelComponents_UV.Rescale=="free")], theta_fitted(i,end)];
                end
        end
    end

    %%
    is_bimodal=false; is_save=false; 

    %colors_orig = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840 ];
    colors = brewermap(length(s_v_range)+1,"Set1");
    colors = colors([1,2,3,4,5,7,8],:);
    mu=0;

    save_suffix = "_";
    if(is_bimodal)
        save_suffix = save_suffix + "BA";  
    else
        save_suffix = save_suffix + "UA";
    end

    if(plot_individual)
        fontsize_UAV = fontsize;
        figspecs_UAV = figspecs;
    else
        fontsize_UAV = fontsize*screen_max_height_ratio;
        figspecs_UAV = figspecs.*screen_max_height_ratio;
    end
    %%
    UV_use_pred_samples = true;
    UA_use_pred_samples = true; % can be false
    model_family = "parametric";
    manuscript_unimodalfits_visualization(data_stratified_UV, fitted_params_PM_UV, ModelComponents_UV, true, colors, s_v_range, model_family, plot_lapse, UV_use_pred_samples, fontsize_UAV, figspecs_UAV, lapse_type, Gaussian_lapse_SDs, plot_individual, connect_human_errorbars)
    manuscript_unimodalfits_visualization(data_stratified_UA, fitted_params_PM_UA, ModelComponents_UA, false, colors, s_a_range, model_family, plot_lapse, UA_use_pred_samples, fontsize_UAV, figspecs_UAV, lapse_type, Gaussian_lapse_SDs, plot_individual, connect_human_errorbars)

end
