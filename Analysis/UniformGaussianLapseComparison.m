model_path = "ModelFits\";
num_models = 2;
num_subjects=15;

NLLs = [];
load(model_path + "fittedparams_UJoint_exp-GaussianLaplaceBothFixedZero_rescalefree_lapseUniform.mat")
NLLs = [NLLs; min(F_vals,[],2)'];
load(model_path + "fittedparams_UJoint_exp-GaussianLaplaceBothFixedZero_rescalefree_lapseGaussian.mat")
NLLs = [NLLs; min(F_vals,[],2)'];
delta_NLLs = NLLs - NLLs(1,:);
NLL_allmodels_sumdiff = sum(delta_NLLs,2);

num_model_params = [14;15];
AICs = 2.*(NLLs + num_model_params);
delta_AICs = AICs - AICs(1,:);
AIC_allmodels_sumdiff = sum(delta_AICs,2);

load("data_stratified_UV.mat");
data_stratified_UV = data_stratified;
load("data_stratified_UA.mat");
data_stratified_UA = data_stratified;
data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
data_UA = data_stratified_to_data(data_stratified_UA, false, false);
n_data = zeros(1,15);
for subjidx=1:15
    n_data(subjidx) = length(data_UA{subjidx}) +length(data_UV{subjidx}); 
end
BICs = 2.*NLLs + num_model_params.*log(n_data);
delta_BICs = BICs - BICs(1,:);
BIC_allmodels_sumdiff = sum(delta_BICs,2);

% AIC/BIC bootstrapping
num_bootstrap_samps = 100000;
NLL_sum_bootstraps = zeros(num_models,num_bootstrap_samps);
AIC_sum_bootstraps = zeros(num_models,num_bootstrap_samps);
BIC_sum_bootstraps = zeros(num_models,num_bootstrap_samps);
rng('default')
rng(0)
for samp = 1:num_bootstrap_samps
    sampled_subj = datasample(1:num_subjects,num_subjects); 
    for model=1:num_models
        NLL_sum_bootstraps(model, samp) = sum(NLLs(model, sampled_subj));
        AIC_sum_bootstraps(model, samp) = sum(AICs(model, sampled_subj));
        BIC_sum_bootstraps(model, samp) = sum(BICs(model, sampled_subj));
    end
end
NLL_bootstraps_diff = NLL_sum_bootstraps - NLL_sum_bootstraps(1,:);
AIC_bootstraps_diff = AIC_sum_bootstraps - AIC_sum_bootstraps(1,:);
BIC_bootstraps_diff = BIC_sum_bootstraps - BIC_sum_bootstraps(1,:);
NLL_bootstraps_errorbars = prctile(NLL_bootstraps_diff,[2.5, 97.5], 2);
AIC_bootstraps_errorbars = prctile(AIC_bootstraps_diff,[2.5, 97.5], 2);
BIC_bootstraps_errorbars = prctile(BIC_bootstraps_diff,[2.5, 97,5], 2);

% Plot
bootstraps_errorbars_allstats = {NLL_bootstraps_errorbars, AIC_bootstraps_errorbars, BIC_bootstraps_errorbars};
allstats = {NLL_allmodels_sumdiff, AIC_allmodels_sumdiff, BIC_allmodels_sumdiff};
stat_names = ["$\Delta NLL$", "$\Delta AIC$", "$\Delta BIC$"];

if(false)
    tiledlayout(length(stat_names),1, 'TileSpacing', 'compact','Padding', 'none')
    for statistic=1:length(stat_names)
        nexttile
        bootstraps_errorbars = bootstraps_errorbars_allstats{statistic};
        mean_stat = allstats{statistic};
        hold on
        bar(1:(num_models),mean_stat,'FaceColor','k', 'FaceAlpha',0.2)
        errorbar(1:(num_models)',mean_stat,mean_stat-squeeze(bootstraps_errorbars(:, 1)),squeeze(bootstraps_errorbars(:, 2))-mean_stat, 'k.')
        if(statistic~=length(stat_names))
            xticklabels([])
        else
            xticks(1:(3*num_models));
            xticklabels(allmodel_xticklabels)
            xtickangle(20)
        end
        ylabel(stat_names(statistic), 'interpreter', 'latex')
        xlim([0, num_models+0.7])
    end
end

% Create a .mat file with the delta NLL, AIC, and BIC tables
AllData_ModelComparison_FinalTables = cell(1,3); 
for stat=1:3
    final_table = zeros(num_models,3);
    final_table(:,1) = bootstraps_errorbars_allstats{stat}(:,1); % 2.5% percentile
    final_table(:,2) = allstats{stat}; % Sum
    final_table(:,3) = bootstraps_errorbars_allstats{stat}(:,2); % 97.5% percentile
    AllData_ModelComparison_FinalTables{stat} = final_table;
end
