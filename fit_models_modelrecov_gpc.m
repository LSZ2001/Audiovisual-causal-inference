
function [] = fit_models_modelrecov_gpc(iteration)
    %cd('audio_visual_perception')
    model_path = "modelfits\";
    model_path_temp = model_path+"temp\";
    data_path = "data\";
    modelrecovdata_path = "data\model_recov\";
    analysis_path = "analysis\";
    num_subjects = 15;
    addpath(analysis_path,data_path,model_path,model_path_temp,"utils\","bads-master");
    
    
    
    prior_types = ["SingleGaussian","GaussianLaplaceBothFixedZero","SingleGaussian","SingleGaussian","GaussianLaplaceBothFixedZero","GaussianLaplaceBothFixedZero"]; % "SingleGaussian", "GaussianLaplaceBothFixedZero", or "TwoGaussiansBothFixedZero"
    hetero_types = ["constant","exp","constant","exp","constant","exp"]; % "constant" or "exp";
    lapse_types = repmat("Uniform",1,length(prior_types)); % "Uniform" or "Gaussian";
    rescale_auds = ["free","free","1","free","free","1"]; % "1", "4/3", or "free";
    % Model parameters
    num_params = [8,14,7,12,10,13];
    
    prior_type_modelrecovdatas = prior_types;
    hetero_type_modelrecovdatas = hetero_types; % "constant" or "exp";
    lapse_type_modelrecovdatas = lapse_types; % "Uniform" or "Gaussian";
    rescale_aud_modelrecovdatas = rescale_auds; % "1", "4/3", or "free";
    
    %% Fit the models to simulated datasets
    dataset_idx = floor(iteration/6)+1;
    model_idx = mod(iteration,6);
    if(model_idx==0)
        model_idx = 6;
        dataset_idx = dataset_idx-1;
    end
    dataset_idx
    model_idx
    prior_type = prior_types(model_idx);
    hetero_type = hetero_types(model_idx);
    lapse_type = lapse_types(model_idx);
    rescale_aud = rescale_auds(model_idx);

    num_inits_persubj = 10; % Number of random initializations per subject.
    iters = 1:(num_inits_persubj * num_subjects);
    prior_type_modelrecovdata = prior_type_modelrecovdatas(dataset_idx); % "SingleGaussian", "GaussianLaplaceBothFixedZero", or "TwoGaussiansBothFixedZero"
    hetero_type_modelrecovdata = hetero_type_modelrecovdatas(dataset_idx); % "constant" or "exp";
    lapse_type_modelrecovdata = lapse_type_modelrecovdatas(dataset_idx); % "Uniform" or "Gaussian";
    rescale_aud_modelrecovdata = rescale_aud_modelrecovdatas(dataset_idx); % "1", "4/3", or "free";

    parfor iter = iters
        fit_ujointmodel_parametric_modelrecov(iter,prior_type,hetero_type, lapse_type, rescale_aud, ...
            prior_type_modelrecovdata, hetero_type_modelrecovdata, lapse_type_modelrecovdata, rescale_aud_modelrecovdata,...
            num_inits_persubj, modelrecovdata_path, model_path_temp)
    end
    merge_parametric_ujointfits_files(prior_type,hetero_type,lapse_type,rescale_aud,...
        prior_type_modelrecovdata, hetero_type_modelrecovdata, lapse_type_modelrecovdata, rescale_aud_modelrecovdata,...
        num_inits_persubj, model_path, model_path_temp)
end

%% Helper functions to merge different subjects/inits for the same model.
function [] = merge_parametric_ujointfits_files(prior_type,hetero_type,lapse_type,rescale_aud, prior_type_modelrecovdata, hetero_type_modelrecovdata, lapse_type_modelrecovdata, rescale_aud_modelrecovdata, num_inits_persubj, model_path, model_path_temp)
    filename = 'fittedparams_UJoint_'+hetero_type+"-"+prior_type+"_rescale"+rescale_aud+"_lapse"+lapse_type;
    datafilename = hetero_type_modelrecovdata+"-"+prior_type_modelrecovdata+"_rescale"+rescale_aud_modelrecovdata+"_lapse"+lapse_type_modelrecovdata;
    filename_final = filename + "__" + datafilename;
    load(model_path_temp+filename_final+"_1.mat");
    num_subjects=15; num_params = length(out_struct.theta_fitted);

    Theta_fitted = zeros(num_subjects, num_inits_persubj, num_params);
    T_Ends = zeros(num_subjects, num_inits_persubj);
    F_vals = zeros(num_subjects, num_inits_persubj);
    
    for i=1:(num_subjects.*num_inits_persubj)
        subjidx = ceil(i/num_inits_persubj);
        randinit = mod(i,num_inits_persubj);
        if(randinit==0)
            randinit=num_inits_persubj;
        end
        [subjidx, randinit]
        load(model_path_temp+filename_final+"_"+num2str(i));
        T_Ends(subjidx, randinit) = out_struct.T_Ends;
        F_vals(subjidx, randinit) = out_struct.F_vals;
        Theta_fitted(subjidx, randinit,:) = out_struct.theta_fitted;
    end

    [min_val, min_idx] = min(F_vals,[],2);
    theta_fitted = zeros(num_subjects, num_params);
    NLL_fitted = zeros(num_subjects,1);
    for i=1:num_subjects
        theta_fitted(i,:) = Theta_fitted(i,min_idx(i),:);
        NLL_fitted(i) = F_vals(i, min_idx(i));
    end
    
    save(model_path+filename_final, "Theta_fitted", "T_Ends","F_vals", "theta_fitted");
end


%% Model recovery AIC/BIC
function [UnimodalData_ModelComparison_FinalTables] = unimodaldata_modelcomparison_visualize(priors, noises, rescales, model_types, num_params, is_plot, fontsize, model_path, data_path, plot_BIC_only)
    num_models = length(priors);
    num_subjects = 15;

    load(data_path+"data_stratified_UV.mat");
    load(data_path+"data_stratified_UA.mat");
    data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
    data_UA = data_stratified_to_data(data_stratified_UA, false, false);
    n_data = zeros(1,15);
    for subjidx=1:num_subjects
        n_data(subjidx) = length(data_UA{subjidx}) +length(data_UV{subjidx}); 
    end
    % AIC, BIC
    NLLs = zeros(num_models, num_subjects);
    AICs = zeros(num_models, num_subjects);
    BICs = zeros(num_models, num_subjects);
    
    % load nonparam indv UJoint fits
    for model=1:(num_models)
        if(model_types(model)=="semiparametric")
            load(model_path + "fittedparams_UJoint_Semiparam_rescalefree_lapseUniform.mat")
            [min_val, min_idx] = min(F_vals,[],1);
            [num_inits,num_params_nonparam]=size(theta_fitted);
            NLLs(model,:) = min_val;
            AICs(model,:) = 2.*min_val + 2.* num_params(model);
            BICs(model,:) = 2.*min_val + num_params(model).*log(n_data);

        else
            load(model_path + "fittedparams_UJoint_"+noises(model)+"-"+priors(model)+"_rescale"+rescales(model)+"_lapseUniform.mat")
            [min_val, min_idx] = min(F_vals,[],2);
            NLLs(model,:) = min_val';
            AICs(model,:) = 2.*min_val' + 2.* num_params(model);
            BICs(model,:) = 2.*min_val' + num_params(model).*log(n_data);
        end
    end
    NLL_sum = sum(NLLs,2);
    AIC_sum = sum(AICs,2);
    BIC_sum = sum(BICs,2);

    NLL_sumvalues_diff = NLL_sum - min(NLL_sum);
    [~,min_NLL_model] = min(NLL_sum);
    AIC_sumvalues_diff = AIC_sum - min(AIC_sum);
    [~,min_AIC_model] = min(AIC_sum);
    BIC_sumvalues_diff = BIC_sum - min(BIC_sum);
    [~,min_BIC_model] = min(BIC_sum);

    % AIC/BIC bootstrapping
    num_bootstrap_samps = 100000;
    NLL_sum_bootstraps = zeros(num_models,num_bootstrap_samps);
    AIC_sum_bootstraps = zeros(num_models,num_bootstrap_samps);
    BIC_sum_bootstraps = zeros(num_models,num_bootstrap_samps);
    rng("default")
    rng(0)
    for samp = 1:num_bootstrap_samps
        sampled_subj = datasample(1:num_subjects,num_subjects); 
        for model=1:(num_models)
            NLL_sum_bootstraps(model, samp) = sum(NLLs(model, sampled_subj));
            AIC_sum_bootstraps(model, samp) = sum(AICs(model, sampled_subj));
            BIC_sum_bootstraps(model, samp) = sum(BICs(model, sampled_subj));
        end
    end
    NLL_bootstraps_diff = NLL_sum_bootstraps - NLL_sum_bootstraps(min_NLL_model,:);
    AIC_bootstraps_diff = AIC_sum_bootstraps - AIC_sum_bootstraps(min_AIC_model,:);
    BIC_bootstraps_diff = BIC_sum_bootstraps - BIC_sum_bootstraps(min_BIC_model,:);

    NLL_bootstraps_errorbars = prctile(NLL_bootstraps_diff,[2.5, 97.5], 2);
    AIC_bootstraps_errorbars = prctile(AIC_bootstraps_diff,[2.5, 97.5], 2);
    BIC_bootstraps_errorbars = prctile(BIC_bootstraps_diff,[2.5, 97,5], 2);

    % Plot
    bootstraps_errorbars_allstats = {NLL_bootstraps_errorbars, AIC_bootstraps_errorbars, BIC_bootstraps_errorbars};
    allstats = {NLL_sumvalues_diff, AIC_sumvalues_diff, BIC_sumvalues_diff};
    stat_names = ["\DeltaNLL", "\DeltaAIC", "\DeltaBIC"];
    
    if(is_plot)
        if(plot_BIC_only)
            statistics = 3;
        else
            statistics = 1:length(stat_names)
        end
        tiledlayout(length(statistics),1, 'TileSpacing', 'tight','Padding', 'none')
        for statistic=statistics
            nexttile
            set(gca,'TickDir','out');
            bootstraps_errorbars = bootstraps_errorbars_allstats{statistic};
            mean_stat = allstats{statistic};
                hold on
                %for model=[1,2,3]
                cats = model_types;
                cats = insertBefore(cats,"_no","\");
                C = categorical(cats);
                C = reordercats(C,cellstr(C)');
                barh(0:(length(C)-1),mean_stat,'FaceColor','k', 'FaceAlpha',0.2)
                errorbar(mean_stat,0:(length(C)-1),mean_stat-squeeze(bootstraps_errorbars(:, 1)),squeeze(bootstraps_errorbars(:, 2))-mean_stat,'horizontal', 'k.')
                yticks(0:(length(C)-1))
                C_capitalized = C;
                for c=1:length(C)
                    cat_char = char(string(C(c)));
                    C_capitalized(c) = convertCharsToStrings([upper(cat_char(1)), cat_char(2:end)]);
                end
                yticklabels(C_capitalized)
                ylim([0-0.5, length(C)-0.5])
                set(gca,'YDir','reverse')
                %end
                if(statistic~=length(stat_names))
                    xticklabels([]);
                    set(gca,'xtick',[])
                end
                set(gca,'FontSize',9)
                xlabel(stat_names(statistic))
                %set(gca,'xticklabel',["diff","max","ent"].')
        end
    end
    % Create a .mat file with the delta NLL, AIC, and BIC tables
    UnimodalData_ModelComparison_FinalTables = cell(1,3); 
    for stat=1:3
        final_table = zeros(length(allstats{1}),3);
        final_table(:,1) = bootstraps_errorbars_allstats{stat}(:,1); % 2.5% percentile
        final_table(:,2) = allstats{stat}; % Sum
        final_table(:,3) = bootstraps_errorbars_allstats{stat}(:,2); % 97.5% percentile
        UnimodalData_ModelComparison_FinalTables{stat} = final_table;
    end
end
