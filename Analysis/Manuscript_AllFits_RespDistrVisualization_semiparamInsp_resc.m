function [] = Manuscript_AllFits_RespDistrVisualization_semiparamInsp_resc(causal_inf_strategy, fontsize, figspecs, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type)
    model_family="semiparamInsp";
    if(nargin==0)
        causal_inf_strategy = "ProbMatching";
        fontsize=9;
        figsize = get(0, 'ScreenSize');
        figspecs = [0 0 figsize(4)*4/3 figsize(4)];
        model_path = "ModelFits\";
        plot_lapse = true; lapse_type = "Uniform";
    elseif(nargin==7)
        % plot the lapse component in fitted response distributions
        plot_lapse = true; lapse_type = "Uniform";
    elseif(nargin==8)
        lapse_type = "Uniform";
    end
    
    %%
    s_a_range = -15:5:15;
    s_v_range = linspace(-20,20,8);
    s_v_range = (s_v_range(2:end) - s_v_range(1:(end-1)))./2 + s_v_range(1:(end-1));
    colors = brewermap(length(s_v_range)+1,"Set1");
    colors = colors([1,2,3,4,5,7,8],:);

    %% Load data files.
    load("BAV_data.mat")
    load("BC_data.mat")
    load("data_stratified_UV.mat");
    num_subjects = length(BAV_data);
    load("data_stratified_UA.mat");
    data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
    data_UA = data_stratified_to_data(data_stratified_UA, false, false);
    num_subjects = length(data_UA);
    UAV_data = cell(1,num_subjects);
    Gaussian_lapse_SDs = zeros(1,num_subjects); % If assume Gaussian Lapse, then its SD is just the SD across all UAV trials for that subject. 
    for i=1:num_subjects
        data_UA{i}(:,3) = 4;
        UAV_data{i} = [data_UV{i}; data_UA{i}];
        Gaussian_lapse_SDs(i) = std([UAV_data{i}(:,2); BAV_data{i}(:,4)]);
    end

    %% Load nonparamInspired fitted param0 (and ModelComponents).
    filename = 'fittedparams_All_UBresc_SemiparamInspired_'+causal_inf_strategy+"_rescalefree";
    filename = filename + "_lapse"+lapse_type;    
    load(model_path + filename+'.mat');
    fitted_params_PM = theta_fitted;
    
    PMIntegrationParams = [-45,45,201];
    
    %% Plots
    is_bimodal=false; is_save=false;

    %%
    close all
    model_type="PM";
    %dx_max = 0.5;
    unique_bins_subj = false; % Whether use unique binedges for each subject, or same binedges shared across subjects.

    %% UV, UA

    UV_use_pred_samples = true;
    UA_use_pred_samples = true; % can be false;
    ModelComponents.NumReliabilityLevels=3;
    Manuscript_UnimodalFits_Visualization(data_stratified_UV, fitted_params_PM, ModelComponents, true, colors, s_v_range, model_family, plot_lapse, UV_use_pred_samples, fontsize, figspecs, lapse_type, Gaussian_lapse_SDs)
    ModelComponents.NumReliabilityLevels=1;
    Manuscript_UnimodalFits_Visualization(data_stratified_UA, fitted_params_PM, ModelComponents, false, colors, s_a_range, model_family, plot_lapse, UA_use_pred_samples, fontsize, NaN, lapse_type, Gaussian_lapse_SDs)

    exportgraphics(gcf,figpath+save_name+'-UAV.png','Resolution',png_dpi);
    exportgraphics(gcf,figpath+save_name+'-UAV.pdf',"ContentType","vector");


    %% BC
    Manuscript_BimodalCFits_Visualization_resc(BC_data, fitted_params_PM, ModelComponents, ModelComponents, model_family, plot_lapse, ModelComponents.CausalInfStrategy, fontsize, figspecs)

    exportgraphics(gcf,figpath+save_name+'-BC.png','Resolution',png_dpi);
    exportgraphics(gcf,figpath+save_name+'-BC.pdf',"ContentType","vector");

    %% BA, BV
    Manuscript_BimodalAVFits_Visualization_resc(BAV_data, fitted_params_PM, ModelComponents, ModelComponents, PMIntegrationParams, NaN, model_family, plot_lapse, ModelComponents.CausalInfStrategy, unique_bins_subj, fontsize, figspecs, figpath, save_name, png_dpi, lapse_type, Gaussian_lapse_SDs);

end
