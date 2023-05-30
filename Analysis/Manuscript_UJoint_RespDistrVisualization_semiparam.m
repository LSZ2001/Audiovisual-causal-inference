% This code does model recovery. It is same as 3_withlapse, but it is now
% compatible with a fake dataset with lapse rates.
function [] = Manuscript_UJoint_RespDistrVisualization_semiparam(fontsize, figspecs, model_path, plot_lapse, lapse_type)
    if(nargin==0)
        figsize = get(0, 'ScreenSize');fontsize = 9; figspecs = [0 0 figsize(4)*0.9*0.65 figsize(4)*0.9];
        model_path = "ModelFits\"
        plot_lapse = true; lapse_type = "Uniform";
    elseif(nargin==3) % Default is plot the lapse component in fitted response distributions
        plot_lapse = true; lapse_type = "Uniform";
    elseif(nargin==4)
        lapse_type = "Uniform";
    end
    
    %cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Analysis')
    s_a_range = -15:5:15;
    s_v_range = linspace(-20,20,8);
    s_v_range = (s_v_range(2:end) - s_v_range(1:(end-1)))./2 + s_v_range(1:(end-1));

    load("BAV_data.mat")
    load("BC_data.mat")
    load("data_stratified_UV.mat");
    data_stratified_UV = data_stratified;
    load("data_stratified_UA.mat");
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
        %Gaussian_lapse_SDs(i) = std([UAV_data{i}(:,2); BAV_data{i}(:,4)]);
    end
    
    %%
    model_type="PM"; % two gaussians components of the prior both centered at 0. 
    PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
    consider_lapse=true; % fit a lapse parameter.
    randomize_trials=false;
    mu = 0; % prior distribution mean
    [~, num_subjects] = size(data_stratified_UV);

    % UV data model components
    ModelComponents.SPivot = [0,0.1,0.3,1,2,4,6,8,10,15,20,45];
    ModelComponents.LapseRange = [-45,45];
    ModelComponents.RescaleVis="free"; % fit k_vis
    ModelComponents.RescaleAud="free"; % fit k_aud
    ModelComponents.NumReliabilityLevels = 3; % Get number of visual reliability levels
    ModelComponents.MotorNoise = "Gaussian";

    num_pivots = length(ModelComponents.SPivot)-1;
    
    %% Load fitted params
    load(model_path+"fittedparams_UJoint_Semiparam_rescalefree_lapseUniform.mat");
    fitted_params_PM = theta_fitted;

    %% Plots
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

    UV_use_pred_samples = true;
    UA_use_pred_samples = false;
    model_type="PM";
    
    %% Plot

%     Manuscript_UnimodalFits_Heterosk_Visualization_nonparamIndv(data_stratified_UV, fitted_params_PM, ModelComponents, PMIntegrationParams, true, colors, s_v_range, model_type, is_save, "_UV", plot_lapse, UV_use_pred_samples, fontsize)
%     Manuscript_UnimodalFits_Heterosk_Visualization_nonparamIndv(data_stratified_UA, fitted_params_PM, ModelComponents, PMIntegrationParams, false, colors, s_a_range, model_type, is_save, "_UA", plot_lapse, UA_use_pred_samples, fontsize)

    UV_use_pred_samples = true;
    UA_use_pred_samples = true; % can be false;
    model_family="semiparam";
    Manuscript_UnimodalFits_Visualization(data_stratified_UV, fitted_params_PM, ModelComponents, true, colors, s_v_range, model_family, plot_lapse, UV_use_pred_samples, fontsize, figspecs, lapse_type, Gaussian_lapse_SDs)
    Manuscript_UnimodalFits_Visualization(data_stratified_UA, fitted_params_PM, ModelComponents, false, colors, s_a_range, model_family, plot_lapse, UA_use_pred_samples, fontsize, NaN, lapse_type, Gaussian_lapse_SDs)

end
