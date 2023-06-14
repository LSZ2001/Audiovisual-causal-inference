clear all; close all;
cd('C:\Users\liu_s\Audiovisual-causal-inference')

figsize = get(0, 'ScreenSize');
figformat = "svg";
figpath = "plots\";
fontsize=9;
png_dpi = 500;
plot_lapse = true;
lapse_type = "Uniform";
model_path = "modelfits\";
data_path = "data\";
analysis_path = "analysis\";

% Unimodal data ModelComparison
priors = ["","GaussianLaplaceBothFixedZero","GaussianLaplaceBothFixedZero","GaussianLaplaceBothFixedZero", "SingleGaussian", "GaussianLaplaceBothFixedZero", "TwoGaussiansBothFixedZero","SingleGaussian","SingleGaussian","SingleGaussian"];
noises = ["","exp", "exp", "exp", "exp", "constant", "exp", "constant","constant","constant"];
rescales = ["","free", "4over3", "1", "free","free","free","free","4over3","1"];
model_types = ["semi-parametric","exp-GaussianLaplace", "exp-GaussianLaplace\_4/3","exp-GaussianLaplace\_1", "exp-SingleGaussian", "const-GaussianLaplace","exp-TwoGaussians","const-SingleGaussian","const-SingleGaussian\_4/3","const-SingleGaussian\_1"];
num_params = [40,14,13,13,12,10,14,8,7,7];

%% Vanila 3 models on unimodal data only
figure('Position', [0 0 figsize(4)*0.9*0.65 figsize(4)*0.9]);
set(gcf, 'Color', 'w')
UnimodalData_ModelComparison_FinalTables_Vanilla = unimodaldata_modelcomparison_visualize(priors((end-2):end), noises((end-2):end), rescales((end-2):end), model_types((end-2):end), num_params((end-2):end), true, fontsize+1, model_path, data_path);
ax = gca; 
ax.FontSize = fontsize+1; 
% exportgraphics(gcf,figpath+'UJoint_ModelSelection_vanilla'+'.png','Resolution',png_dpi);
% exportgraphics(gcf,figpath+'UJoint_ModelSelection_vanilla'+'.pdf',"ContentType","vector");
save(analysis_path+'unimodaldata_modelcomparison_finaltables_vanilla','UnimodalData_ModelComparison_FinalTables_Vanilla');

%% All models on unimodal data 
figure('Position', [0 0 figsize(4)*0.9*5/4 figsize(4)*0.9]);
set(gcf, 'Color', 'w')
UnimodalData_ModelComparison_FinalTables = unimodaldata_modelcomparison_visualize(priors, noises, rescales, model_types, num_params, true, fontsize+1, model_path, data_path);
ax = gca; 
ax.FontSize = fontsize+1; 
% exportgraphics(gcf,figpath+'UJoint_ModelSelection'+'.png','Resolution',png_dpi);
% exportgraphics(gcf,figpath+'UJoint_ModelSelection'+'.pdf',"ContentType","vector");
save(analysis_path+'unimodaldata_modelcomparison_finaltables','UnimodalData_ModelComparison_FinalTables');

%%
function [UnimodalData_ModelComparison_FinalTables] = unimodaldata_modelcomparison_visualize(priors, noises, rescales, model_types, num_params, is_plot, fontsize, model_path, data_path)
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
    model_start_idx = 1;
    if(model_types(1)=="semi-parametric")
        load(model_path + "fittedparams_UJoint_Semiparam_rescalefree_lapseUniform.mat")
        [min_val, min_idx] = min(F_vals,[],1);
        [num_inits,num_params_nonparam]=size(theta_fitted);
        NLLs(1,:) = min_val;
        AICs(1,:) = 2.*min_val + 2.* num_params(1);
        BICs(1,:) = 2.*min_val + num_params(1).*log(n_data);
        model_start_idx = model_start_idx+1;
    end

    for model=model_start_idx:(num_models)
        load(model_path + "fittedparams_UJoint_"+noises(model)+"-"+priors(model)+"_rescale"+rescales(model)+"_lapseUniform.mat")
        [min_val, min_idx] = min(F_vals,[],2);
        NLLs(model,:) = min_val';
        AICs(model,:) = 2.*min_val' + 2.* num_params(model);
        BICs(model,:) = 2.*min_val' + num_params(model).*log(n_data);
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
    stat_names = ["$\Delta NLL$", "$\Delta AIC$", "$\Delta BIC$"];

    if(is_plot)
        tiledlayout(length(stat_names),1, 'TileSpacing', 'compact','Padding', 'none')
        for statistic=1:length(stat_names)
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
                bar(C,mean_stat,'FaceColor','k', 'FaceAlpha',0.2)
                errorbar(C,mean_stat,mean_stat-squeeze(bootstraps_errorbars(:, 1)),squeeze(bootstraps_errorbars(:, 2))-mean_stat, 'k.')
                %end
                if(statistic~=length(stat_names))
                    xticklabels([]);
                    set(gca,'xtick',[])
                end
                ylabel(stat_names(statistic), 'interpreter', 'latex')
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
    set(gca,'fontsize', fontsize)
end