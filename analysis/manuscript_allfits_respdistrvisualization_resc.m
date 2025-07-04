function [] = manuscript_allfits_respdistrvisualization_resc(prior_type, hetero_type, causal_inf_strategy, fontsize, figspecs, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type, plot_individual)
    model_family = "parametric";
    if(nargin==0)
        prior_type = "GaussianLaplaceBothFixedZero";
        hetero_type = "exp";
        causal_inf_strategy = "ModelSelection"; % ModelAveraging, ModelSelection, or ProbMatching.
        fontsize=9;
        figsize = get(0, 'ScreenSize');
        figspecs = [0 0 figsize(4)*4/3 figsize(4)];
        model_path = "modelfits\";
        plot_lapse = true; lapse_type = "Uniform"; plot_individual=false;
    elseif(nargin==9)
        % plot the lapse component in fitted response distributions
        plot_lapse = true; lapse_type = "Uniform"; plot_individual=false;
    elseif(nargin==10)
        lapse_type = "Uniform"; plot_individual=false;
    elseif(nargin==11)
        plot_individual=false;
    end
    data_path = "data\";

    % For UAV plots, in case screen too small, need to shrink size. 
    set(0,'units','pixels');
    screensize = get(0, 'ScreenSize');
    screen_max_height_ratio = screensize(4)/figspecs(4);

    %%
    s_a_range = -15:5:15;
    s_v_range = linspace(-20,20,8);
    s_v_range = (s_v_range(2:end) - s_v_range(1:(end-1)))./2 + s_v_range(1:(end-1));
    load(data_path+"bav_data.mat")
    load(data_path+"bc_data.mat")
    load(data_path+"data_stratified_uv.mat");
    load(data_path+"data_stratified_ua.mat");
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
    

    %%
    model_type="PM"; % two gaussians components of the prior both centered at 0. 
    PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
    [~, num_subjects] = size(BAV_data);
    prior_assumption = "independent"; % independent or empirical.

    % UV data model components
    ModelComponents_UV.CausalInfStrategy = causal_inf_strategy;
    ModelComponents_UV.PriorType=prior_type;
    ModelComponents_UV.PriorAssumption=prior_assumption;
    ModelComponents_UV.SensoryNoise=hetero_type;
    ModelComponents_UV.LapseRange = [-45,45];
    ModelComponents_UV.IsFixedPriorMean=(model_type=="PM");
    ModelComponents_UV.Rescale="1";
    ModelComponents_UV.MotorNoise="Gaussian";
    ModelComponents_UV.IsMLModel=(model_type=="ML");
    ModelComponents_UV.NumReliabilityLevels = 3; % Get number of reliability levels
    ModelComponents_UV.PMIntegrationParams = PMIntegrationParams;
    sigma_fun = heterotype_to_sigmafun(ModelComponents_UV.SensoryNoise);
    [LB_V, UB_V, PLB_V, PUB_V] = sigmafun_badsbounds_comprehensive(ModelComponents_UV);
    LB_V(1) = 1; PLB_V(1) = 1; % since dx_max=0.5, simga0 must be much larger than it.
    num_V_params = length(LB_V);


    % UA data model components -> same prior/noise model/lapse/motor noise as UV,
    % but different noise params sigma0, k1, k2 and include a rescale.
    ModelComponents_UA = ModelComponents_UV;
    ModelComponents_UA.Rescale="free";
    ModelComponents_UA.NumReliabilityLevels = 1; % Get number of reliability levels
    [LB_A, UB_A, PLB_A, PUB_A] = sigmafun_badsbounds_comprehensive(ModelComponents_UA);
    LB_A(1) = 1; PLB_A(1) = 1; % since dx_max=0.5, simga0 must be much larger than it.

    % Merge the BADS bounds of UV and UA to one vector, UV in front.
    [LB, UB, PLB, PUB, A_param_keep_idx] = merge_ujoint_badsbounds(LB_V,UB_V,PLB_V,PUB_V,LB_A,UB_A,PLB_A,PUB_A,ModelComponents_UA);

    % Add in final parameter constraints -- p_same, which is just Pr[C=1].
    LB = [LB, 0];
    LB(8) = 1;
    PLB = [PLB, 0];
    PLB(8) = 1.1;
    PUB = [PUB, 1];
    UB = [UB, 1];

    num_A_uniq_params = length(LB) - length(LB_A) - 1;

    %% Load fitted params
    filename = 'fittedparams_All_UBresc_'+hetero_type+"-"+prior_type+"-"+causal_inf_strategy+"_rescalefree";
    filename = filename + "_lapse"+lapse_type;
    load(model_path + filename);
    
    % load("temp_all_SingleGaussianPMexpdnoriseGaussianProbaMatching_Arescalefree_Subj1to7.mat");
    % 
    % [~,~,num_params]=size(Theta_fitted);
    % [min_val, min_idx] = min(F_vals,[],2);
    % theta_fitted = zeros(num_subjects, num_params);
    % NLL_fitted = zeros(num_subjects,1);
    % for i=1:num_subjects
    %     theta_fitted(i,:) = Theta_fitted(i,min_idx(i),:);
    %     NLL_fitted(i) = F_vals(i, min_idx(i));
    % end
    fitted_params_PM = theta_fitted;
    % sigma0s = fitted_params_PM(:,[1,11]);
    % sigma0s = max(sigma0s, 1);
    % fitted_params_PM(:,[1,11]) = sigma0s;
    theta_fitted = theta_fitted(:,1:(end-3));

    %% Plots
    fitted_params_PM_UV = theta_fitted(:, 1:(end-length(A_param_keep_idx)));
    fitted_params_PM_UA = [];
    for i=1:num_subjects
        fitted_params_PM_UA = [fitted_params_PM_UA; complete_thetaua_for_ujointfits(theta_fitted(i,:), A_param_keep_idx, ModelComponents_UV.Rescale=="free")];
    end

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

    %% Check NLL code
    % fittypes = ["BV", "BA"];
    % response_type=2;
    % fittype=fittypes(response_type);
    % 
    % subjidx=1;
    % BAV_datasubj = BAV_data{subjidx};
    % nll1 = zeros(1,2);
    % fittypes = ["BV", "BA"];
    % for resp_type=1:2
    %     B_datasubj = BAV_datasubj(BAV_datasubj(:,5)==resp_type,:);
    %     R = B_datasubj(:,4);
    %     S_V = B_datasubj(:,[1,3]); S_A = B_datasubj(:,2);
    %     nll1(resp_type) = midpoint_postmean_NLLfun_comprehensive_BAV_old(ModelComponents_V, ModelComponents_A, fitted_params_PM(1,:), R, S_V, S_A, PMIntegrationParams, false, false, true, 0.5, "ModelAveraging", fittypes(resp_type));
    % end
    % 
    % R = BAV_datasubj(:,4);
    % S_V = BAV_datasubj(:,[1,3]); S_A = BAV_datasubj(:,2);
    % R = BAV_datasubj(:,[4,5]);
    % nll2 = midpoint_postmean_NLLfun_comprehensive_BAV(ModelComponents_V, ModelComponents_A, fitted_params_PM(1,:), R, S_V, S_A, PMIntegrationParams, false, false, true, 0.5, "ModelAveraging");
    % 
    % sum(nll1)
    % nll2



    %%
    close all
    model_type="PM";
    dx_max = 0.5;
    unique_bins_subj = false; % Whether use unique binedges for each subject, or same binedges shared across subjects.


    model_family = "parametric";
    
    %% UV, UA
    UV_use_pred_samples = true;
    UA_use_pred_samples = true; % can be false
    if(plot_individual)
        fontsize_UAV = fontsize;
        figspecs_UAV = figspecs;
        manuscript_unimodalfits_visualization(data_stratified_UV, fitted_params_PM_UV, ModelComponents_UV, true, colors, s_v_range, model_family, plot_lapse, UV_use_pred_samples, fontsize_UAV, figspecs_UAV, lapse_type, Gaussian_lapse_SDs, plot_individual)
        manuscript_unimodalfits_visualization(data_stratified_UA, fitted_params_PM_UA, ModelComponents_UA, false, colors, s_a_range, model_family, plot_lapse, UA_use_pred_samples, fontsize_UAV, figspecs_UAV, lapse_type, Gaussian_lapse_SDs, plot_individual)

        figure(1)
        saveas(gca, figpath+save_name+'-UAV_Individualmean.fig')
        exportgraphics(gcf,figpath+save_name+'-UAV_Individualmean.png','Resolution',png_dpi);
        exportgraphics(gcf,figpath+save_name+'-UAV_Individualmean.pdf',"ContentType","vector");
        figure(2)
        saveas(gca, figpath+save_name+'-UAV_IndividualSD.fig')
        exportgraphics(gcf,figpath+save_name+'-UAV_IndividualSD.png','Resolution',png_dpi);
        exportgraphics(gcf,figpath+save_name+'-UAV_IndividualSD.pdf',"ContentType","vector");
    else
        fontsize_UAV = fontsize*screen_max_height_ratio;
        figspecs_UAV = figspecs.*screen_max_height_ratio;
        manuscript_unimodalfits_visualization(data_stratified_UV, fitted_params_PM_UV, ModelComponents_UV, true, colors, s_v_range, model_family, plot_lapse, UV_use_pred_samples, fontsize_UAV, figspecs_UAV, lapse_type, Gaussian_lapse_SDs, plot_individual)
        manuscript_unimodalfits_visualization(data_stratified_UA, fitted_params_PM_UA, ModelComponents_UA, false, colors, s_a_range, model_family, plot_lapse, UA_use_pred_samples, fontsize_UAV, figspecs_UAV, lapse_type, Gaussian_lapse_SDs, plot_individual)

        exportgraphics(gcf,figpath+save_name+'-UAV.png','Resolution',png_dpi);
        exportgraphics(gcf,figpath+save_name+'-UAV.pdf',"ContentType","vector");
    end

    %% BC
    manuscript_bimodalcfits_visualization_resc(BC_data, fitted_params_PM, ModelComponents_UV, ModelComponents_UA, model_family, plot_lapse, ModelComponents_UV.CausalInfStrategy, fontsize, figspecs, plot_individual)
    if(plot_individual)
        saveas(gca, figpath+save_name+'-BC_Individual.fig')
        exportgraphics(gcf,figpath+save_name+'-BC_Individual.png','Resolution',png_dpi);
        exportgraphics(gcf,figpath+save_name+'-BC_Individual.pdf',"ContentType","vector");
    else
        exportgraphics(gcf,figpath+save_name+'-BC.png','Resolution',png_dpi);
        exportgraphics(gcf,figpath+save_name+'-BC.pdf',"ContentType","vector");
    end
    %% BA, BV
    manuscript_bimodalavfits_visualization_resc(BAV_data, fitted_params_PM, ModelComponents_UV, ModelComponents_UA, PMIntegrationParams, dx_max, model_family, plot_lapse, ModelComponents_UV.CausalInfStrategy, unique_bins_subj, fontsize, figspecs, figpath, save_name, png_dpi, lapse_type, Gaussian_lapse_SDs, plot_individual);

end
