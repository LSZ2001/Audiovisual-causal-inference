clear all; close all;
cd('C:\Users\liu_s\Audiovisual-causal-inference')
fig_maxwidth_inches = 7.5;
fig_maxheight_inches = 8.75;
set(0,'units','inches');
Inch_SS = get(0,'screensize');
set(0,'units','pixels');
figsize = get(0, 'ScreenSize');
Res = figsize(3)./Inch_SS(3);
set(groot,'DefaultAxesFontName','Arial')

%figsize_RespDistr = [0,0,figsize(4)*4/3, figsize(4)];
figsize_RespDistr = [0,0,fig_maxwidth_inches, fig_maxheight_inches] .* Res; 

figformat = "svg";
figpath = "plots\"; %"newplots\"
fontsize=9; %9
png_dpi = 500;
plot_lapse = true;
lapse_type = "Uniform";
model_path = "modelfits\";
data_path = "data\";
analysis_path = "analysis\";
addpath(analysis_path,data_path,model_path,"utils\");

%% s_V, s_A generative distributions
linewidth = 1;

figure('Position', [0,0,figsize_RespDistr(3)./5,figsize_RespDistr(3)./8]); 
tiledlayout(1,1,'TileSpacing','none', 'Padding','none'); set(gca,'TickDir','out'); hold on;
sV_vals = [-25:1:-20,-20:1:20,20:1:25]; p_sV_vals = [zeros(size(-25:1:-20)),repmat(1/40,1,length(-20:1:20)),zeros(size(20:1:25))];
plot(sV_vals, p_sV_vals, "k-", 'LineWidth',linewidth);  area(sV_vals, p_sV_vals,'FaceColor','k', 'FaceAlpha',0.2); 
xlim([-25,25]);  ylim([0,0.03]); xticks(-20:20:20); yticks([0,1/40]); yticklabels(["0","1/40"])
xlabel("{\its}_V",'FontSize',fontsize+1); ylabel("p({\its}_V)",'FontSize',fontsize+1);
set(gca,'FontSize',fontsize)
exportgraphics(gcf,figpath+'sV_gen'+'.pdf',"ContentType","vector");

figure('Position', [0,0,figsize_RespDistr(3)./5,figsize_RespDistr(3)./8]); 
tiledlayout(1,1,'TileSpacing','none', 'Padding','tight'); set(gca,'TickDir','out'); hold on;
sA_vals = [-15:5:15]; p_sA_vals = repmat(1/length(sA_vals), 1,length(sA_vals)); 
h=stem(sA_vals, p_sA_vals,"k-", 'LineWidth',linewidth); set(h, 'Marker', 'none')
xlim([-20,20]); ylim([0,0.16]); xticks(-15:15:15); xtickangle(0); yticks([0,1/length(sA_vals)]); yticklabels(["0","1/7"])
xlabel("{\its}_A",'FontSize',fontsize); ylabel("p({\its}_A)",'FontSize',fontsize);
set(gca,'FontSize',fontsize)
exportgraphics(gcf,figpath+'sA_gen'+'.pdf',"ContentType","vector");

%% UAV data visualized, without model prediction ribbons.
prior = "NaN";
noise = "NaN";
aud_rescale = "NaN";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'UAV_dataonly'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'UAV_dataonly'+'.pdf',"ContentType","vector");


%% UJoint parametric model response distribution visualization
prior = "SingleGaussian";
noise = "constant";
aud_rescale = "1";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type, false, true);
exportgraphics(gcf,figpath+'Const-SingleGaussian_rescaleaud1'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Const-SingleGaussian_rescaleaud1'+'.pdf',"ContentType","vector");
%%
prior = "SingleGaussian";
noise = "constant";
aud_rescale = "free";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'Const-SingleGaussian'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Const-SingleGaussian'+'.pdf',"ContentType","vector");

%%
prior = "GaussianLaplaceBothFixedZero";
noise = "exp";
aud_rescale = "free";
%manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, four_by_three_figsize);
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'Exp-GaussianLaplace'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Exp-GaussianLaplace'+'.pdf',"ContentType","vector");

%% Individual-level
prior = "GaussianLaplaceBothFixedZero";
noise = "exp";
aud_rescale = "free";
save_name = "Exp-GaussianLaplace";
plot_individual = true;
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type, plot_individual);
figure(1)
saveas(gca, figpath+save_name+'_Individualmean.fig')
exportgraphics(gcf,figpath+save_name+'_Individualmean'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+save_name + '_Individualmean'+'.pdf',"ContentType","vector");
figure(2)
saveas(gca, figpath+save_name+'_IndividualSD.fig')
exportgraphics(gcf,figpath+save_name+'_IndividualSD'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+save_name + '_IndividualSD'+'.pdf',"ContentType","vector");

%% Exemplary subject
subjidx=7;
fitted_on_all_data = false;
allindvsubjplots_to_onesubjplot(save_name, subjidx, fitted_on_all_data, 10, [0 0 figsize_RespDistr(3) figsize_RespDistr(3)*0.5], figpath)
exportgraphics(gcf,figpath+'Exp-GaussianLaplace_Individual_example'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Exp-GaussianLaplace_Individual_example'+'.pdf',"ContentType","vector");

%% No Lapse 
prior = "GaussianLaplaceBothFixedZero";
noise = "exp";
aud_rescale = "free";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, false, lapse_type);
exportgraphics(gcf,figpath+'Exp-GaussianLaplace_nolapse'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Exp-GaussianLaplace_nolapse'+'.pdf',"ContentType","vector");

%% Gaussian
prior = "GaussianLaplaceBothFixedZero";
noise = "exp";
aud_rescale = "free";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, true, "Gaussian");
exportgraphics(gcf,figpath+'Exp-GaussianLaplace_Gaussianlapse'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Exp-GaussianLaplace_Gaussianlapse'+'.pdf',"ContentType","vector");

% Compute model comparison between uniform lapse model.
UnimodalData_ModelComparison_FinalTables_uniformgaussianlapse = unimodaldata_modelcomparison_visualize_uniformgaussianlapse(model_path, data_path);

%%
prior = "SingleGaussian";
noise = "exp";
aud_rescale = "free";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'Exp-SingleGaussian'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Exp-SingleGaussian'+'.pdf',"ContentType","vector");

prior = "GaussianLaplaceBothFixedZero";
noise = "constant";
aud_rescale = "free";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'Const-GaussianLaplace'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Const-GaussianLaplace'+'.pdf',"ContentType","vector");

prior = "TwoGaussiansBothFixedZero";
noise = "exp";
aud_rescale = "free";
manuscript_ujoint_respdistrvisualization(prior, noise, aud_rescale, fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'Exp-TwoGaussians'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Exp-TwoGaussians'+'.pdf',"ContentType","vector");

%% Unimodal semiparam model 
% Response distribution visualization
% manuscript_ujoint_respdistrvisualization_semiparam(fontsize, four_by_three_figsize);
manuscript_ujoint_respdistrvisualization_semiparam(fontsize, figsize_RespDistr, model_path, plot_lapse, lapse_type);
exportgraphics(gcf,figpath+'Semiparam_FittedRespDistr'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Semiparam_FittedRespDistr'+'.pdf',"ContentType","vector");

%%
% sigma(s), p(s) visualization
semiparam_sigmafun_prior_visualization(fontsize+1, figsize_RespDistr, model_path);
exportgraphics(gcf,figpath+'Semiparam_FittedParams'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'Semiparam_FittedParams'+'.pdf',"ContentType","vector");

%% Unimodal data ModelComparison
priors = ["","GaussianLaplaceBothFixedZero","GaussianLaplaceBothFixedZero","GaussianLaplaceBothFixedZero", "SingleGaussian", "GaussianLaplaceBothFixedZero", "TwoGaussiansBothFixedZero","SingleGaussian","SingleGaussian","SingleGaussian"];
noises = ["","exp", "exp", "exp", "exp", "constant", "exp", "constant","constant","constant"];
rescales = ["","free", "4over3", "1", "free","free","free","free","4over3","1"];
model_types = ["semiparametric","exp-GaussianLaplace", "exp-GaussianLaplace\_4/3","exp-GaussianLaplace\_1", "exp-SingleGaussian", "const-GaussianLaplace","exp-TwoGaussians","const-SingleGaussian","const-SingleGaussian\_4/3","const-SingleGaussian\_1"];
num_params = [40,14,13,13,12,10,14,8,7,7];

%% Vanila 3 models on unimodal data only
figure('Position', [0 0 5.2*Res figsize_RespDistr(4)*0.4]);
set(gcf, 'Color', 'w')
UnimodalData_ModelComparison_FinalTables_Vanilla = unimodaldata_modelcomparison_visualize(priors((end-2):end), noises((end-2):end), rescales((end-2):end), model_types((end-2):end), num_params((end-2):end), true, fontsize+1, model_path, data_path, false);
exportgraphics(gcf,figpath+'UJoint_ModelSelection_vanilla'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'UJoint_ModelSelection_vanilla'+'.pdf',"ContentType","vector");
save(analysis_path+'unimodaldata_modelcomparison_finaltables_vanilla','UnimodalData_ModelComparison_FinalTables_Vanilla');
        
% All models on unimodal data
keep_modelidx = [2,5,6,7,8,1];
figure('Position', [0 0 5.2*Res figsize_RespDistr(3)*0.6*0.5]);
set(gcf, 'Color', 'w')
UnimodalData_ModelComparison_FinalTables = unimodaldata_modelcomparison_visualize(priors(keep_modelidx), noises(keep_modelidx), rescales(keep_modelidx), model_types(keep_modelidx), num_params(keep_modelidx), true, fontsize+1, model_path, data_path, true);
exportgraphics(gcf,figpath+'UJoint_ModelSelection_BIC'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'UJoint_ModelSelection_BIC'+'.pdf',"ContentType","vector");

keep_modelidx = [2,3,4,5,6,7,8,9,10,1];
figure('Position', [0 0 5.2*Res figsize_RespDistr(3)*0.6*0.7]);
set(gcf, 'Color', 'w')
UnimodalData_ModelComparison_FinalTables = unimodaldata_modelcomparison_visualize(priors(keep_modelidx), noises(keep_modelidx), rescales(keep_modelidx), model_types(keep_modelidx), num_params(keep_modelidx), true, fontsize+1, model_path, data_path, true);
exportgraphics(gcf,figpath+'UJoint_ModelSelection_BIC_full'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'UJoint_ModelSelection_BIC_full'+'.pdf',"ContentType","vector");

figure('Position', [0 0 figsize_RespDistr(3) figsize_RespDistr(3)]);
set(gcf, 'Color', 'w')
UnimodalData_ModelComparison_FinalTables = unimodaldata_modelcomparison_visualize(priors, noises, rescales, model_types, num_params, true, fontsize+1, model_path, data_path, false);
exportgraphics(gcf,figpath+'UJoint_ModelSelection'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'UJoint_ModelSelection'+'.pdf',"ContentType","vector");
save(analysis_path+'unimodaldata_modelcomparison_finaltables','UnimodalData_ModelComparison_FinalTables');

%% AllData LiftedSemiparam Response distribution visualization
causal_inf_strategy = "ProbMatching";
save_name = "PM";
manuscript_allfits_respdistrvisualization_semiparaminsp_resc(causal_inf_strategy, fontsize, figsize_RespDistr, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type);

%%
causal_inf_strategy = "ModelSelection";
save_name = "MS";
manuscript_allfits_respdistrvisualization_semiparaminsp_resc(causal_inf_strategy, fontsize, figsize_RespDistr, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type);

causal_inf_strategy = "ModelAveraging";
save_name = "MA";
manuscript_allfits_respdistrvisualization_semiparaminsp_resc(causal_inf_strategy, fontsize, figsize_RespDistr, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type);

%% AllData parametric model response distributions
prior_type = "GaussianLaplaceBothFixedZero";
hetero_type = "exp";
causal_inf_strategy = "ProbMatching";
save_name = "exp-GaussianLaplace-PM";
manuscript_allfits_respdistrvisualization_resc(prior_type, hetero_type, causal_inf_strategy, fontsize, figsize_RespDistr, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type)

%% Individual-level plots for the above model
prior_type = "GaussianLaplaceBothFixedZero";
hetero_type = "exp";
causal_inf_strategy = "ProbMatching";
save_name = "exp-GaussianLaplace-PM";
manuscript_allfits_respdistrvisualization_resc(prior_type, hetero_type, causal_inf_strategy, fontsize, figsize_RespDistr, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type, true)

%% Exemplary subject
subjidx=7;
fitted_on_all_data = true;
allindvsubjplots_to_onesubjplot(save_name,subjidx, fitted_on_all_data, fontsize, figsize_RespDistr, figpath)
exportgraphics(gcf,figpath+'Exp-GaussianLaplace-PM_Individual_example'+'.pdf',"ContentType","vector");
exportgraphics(gcf,figpath+'Exp-GaussianLaplace-PM_Individual_example'+'.png','Resolution',png_dpi);

%%
prior_type = "GaussianLaplaceBothFixedZero";
hetero_type = "exp";
causal_inf_strategy = "ModelAveraging";
save_name = "exp-GaussianLaplace-MA";
manuscript_allfits_respdistrvisualization_resc(prior_type, hetero_type, causal_inf_strategy, fontsize, figsize_RespDistr, figpath, save_name, png_dpi, model_path, plot_lapse, lapse_type)

%% AllData ModelComparison
causal_inf_strategies = ["ModelSelection","ModelAveraging","ProbMatching"];
param_model_names = ["paramBest", "LiftedSemiparam","exp-GaussianLaplace","exp-SingleGaussian","const-GaussianLaplace","exp-TwoGaussians","const-SingleGaussian"];

figure('Position', [0 0 figsize_RespDistr(3) figsize_RespDistr(3)*0.6]);
set(gcf, 'Color', 'w')
alldata_modelcomparison_visualize(causal_inf_strategies, param_model_names, true, fontsize, model_path, data_path, true);
exportgraphics(gcf,figpath+'SemiparamIndv_ModelSelection_BIC'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'SemiparamIndv_ModelSelection_BIC'+'.pdf',"ContentType","vector");


figure('Position', [0 0 figsize_RespDistr(3) figsize_RespDistr(3)]);
set(gcf, 'Color', 'w')
AllData_ModelComparison_FinalTables = alldata_modelcomparison_visualize(causal_inf_strategies, param_model_names, true, fontsize, model_path, data_path, false);
exportgraphics(gcf,figpath+'SemiparamIndv_ModelSelection'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'SemiparamIndv_ModelSelection'+'.pdf',"ContentType","vector");
save(analysis_path+'alldata_modelcomparison_finaltables','AllData_ModelComparison_FinalTables');


%% sigma(s) and p(s) examples
figure('Position', [0 0 figsize_RespDistr(3) figsize_RespDistr(3)]);
set(gcf, 'Color', 'w')
sigmafun_prior_examples(fontsize);
exportgraphics(gcf,figpath+'SensoryNoisePriorParamFamilies'+'.png','Resolution',png_dpi);
exportgraphics(gcf,figpath+'SensoryNoisePriorParamFamilies'+'.pdf',"ContentType","vector");




























%% BELOW: Helper Functions 
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
                yticklabels(C)
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

%%
function [UnimodalData_ModelComparison_FinalTables] = unimodaldata_modelcomparison_visualize_uniformgaussianlapse(model_path, data_path) 
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

    load(data_path+"data_stratified_UV.mat");
    load(data_path+"data_stratified_UA.mat");
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
    stat_names = ["\DeltaNLL", "\DeltaAIC", "\DeltaBIC"];

    if(false)
        tiledlayout(length(stat_names),1, 'TileSpacing', 'tight','Padding', 'none')
        for statistic=1:length(stat_names)
            nexttile
            set(gca,'TickDir','out');
            bootstraps_errorbars = bootstraps_errorbars_allstats{statistic};
            mean_stat = allstats{statistic};
            hold on
            bar(1:(num_models),mean_stat,'FaceColor','k', 'FaceAlpha',0.2)
            errorbar(1:(num_models)',mean_stat,mean_stat-squeeze(bootstraps_errorbars(:, 1)),squeeze(bootstraps_errorbars(:, 2))-mean_stat, 'k.')
            if(statistic~=length(stat_names))
                xticklabels([])
                set(gca,'xtick',[])
            else
                xticks(1:(3*num_models));
                xticklabels(allmodel_xticklabels)
                xtickangle(20)
            end
            set(gca,'FontSize',9)
            ylabel(stat_names(statistic), 'FontSize',9)
            xlim([0, num_models+0.7])
        end
    end

    % Create a .mat file with the delta NLL, AIC, and BIC tables
    UnimodalData_ModelComparison_FinalTables = cell(1,3); 
    for stat=1:3
        final_table = zeros(num_models,3);
        final_table(:,1) = bootstraps_errorbars_allstats{stat}(:,1); % 2.5% percentile
        final_table(:,2) = allstats{stat}; % Sum
        final_table(:,3) = bootstraps_errorbars_allstats{stat}(:,2); % 97.5% percentile
        UnimodalData_ModelComparison_FinalTables{stat} = final_table;
    end

end

%%
function [AllData_ModelComparison_FinalTables] = alldata_modelcomparison_visualize(causal_inf_strategies, param_model_names, is_plot, fontsize, model_path, data_path, plot_BIC_only)
    num_strategies = length(causal_inf_strategies);
    num_models = length(param_model_names);
    num_subjects = 15;
    
    out_structs = cell(1,num_strategies);
    NLL_allmodels = zeros(num_models*num_strategies, 15);
    AIC_allmodels = zeros(num_models*num_strategies, 15);
    BIC_allmodels = zeros(num_models*num_strategies, 15);
    for causal_inf_strategy_idx=1:num_strategies
        causal_inf_strategy_idx
        out_struct = alldata_modelcomparison_visualize_helper(causal_inf_strategies(causal_inf_strategy_idx), model_path, data_path);
        out_structs{causal_inf_strategy_idx} = out_struct;
        NLL_allmodels(causal_inf_strategy_idx:3:end,:) = out_struct.NLLs;
        AIC_allmodels(causal_inf_strategy_idx:3:end,:) = out_struct.AICs;
        BIC_allmodels(causal_inf_strategy_idx:3:end,:) = out_struct.BICs;
    end

    % deltaNLL and deltaAIC across all 18 models
    NLL_allmodels_sum = sum(NLL_allmodels,2);
    NLL_allmodels_sumdiff = NLL_allmodels_sum - min(NLL_allmodels_sum);
    [~,NLL_min_model] = min(NLL_allmodels_sum);

    AIC_allmodels_sum = sum(AIC_allmodels,2);
    AIC_allmodels_sumdiff = AIC_allmodels_sum - min(AIC_allmodels_sum);
    [~,AIC_min_model] = min(AIC_allmodels_sum);

    BIC_allmodels_sum = sum(BIC_allmodels,2);
    BIC_allmodels_sumdiff = BIC_allmodels_sum - min(BIC_allmodels_sum);
    [~,BIC_min_model] = min(BIC_allmodels_sum);

    causal_inf_strategies_abbrev = "-"+["MS","MA","PM"];
    allmodel_xticklabels =  param_model_names' + causal_inf_strategies_abbrev;
    allmodel_xticklabels = allmodel_xticklabels';
    allmodel_xticklabels = allmodel_xticklabels(:);
    
    % AIC/BIC bootstrapping
    num_bootstrap_samps = 100000;
    NLL_sum_bootstraps = zeros(3*num_models,num_bootstrap_samps);
    AIC_sum_bootstraps = zeros(3*num_models,num_bootstrap_samps);
    BIC_sum_bootstraps = zeros(3*num_models,num_bootstrap_samps);
    rng('default')
    rng(0)
    for samp = 1:num_bootstrap_samps
        sampled_subj = datasample(1:num_subjects,num_subjects); 
        for model=1:(3*num_models)
            NLL_sum_bootstraps(model, samp) = sum(NLL_allmodels(model, sampled_subj));
            AIC_sum_bootstraps(model, samp) = sum(AIC_allmodels(model, sampled_subj));
            BIC_sum_bootstraps(model, samp) = sum(BIC_allmodels(model, sampled_subj));
        end
    end
    NLL_bootstraps_diff = NLL_sum_bootstraps - NLL_sum_bootstraps(NLL_min_model,:);
    AIC_bootstraps_diff = AIC_sum_bootstraps - AIC_sum_bootstraps(AIC_min_model,:);
    BIC_bootstraps_diff = BIC_sum_bootstraps - BIC_sum_bootstraps(BIC_min_model,:);
    NLL_bootstraps_errorbars = prctile(NLL_bootstraps_diff,[2.5, 97.5], 2);
    AIC_bootstraps_errorbars = prctile(AIC_bootstraps_diff,[2.5, 97.5], 2);
    BIC_bootstraps_errorbars = prctile(BIC_bootstraps_diff,[2.5, 97,5], 2);

    % Plot
    bootstraps_errorbars_allstats = {NLL_bootstraps_errorbars, AIC_bootstraps_errorbars, BIC_bootstraps_errorbars};
    allstats = {NLL_allmodels_sumdiff, AIC_allmodels_sumdiff, BIC_allmodels_sumdiff};
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
            bar(1:(3*num_models),mean_stat,'FaceColor','k', 'FaceAlpha',0.2)
            errorbar(1:(3*num_models)',mean_stat,mean_stat-squeeze(bootstraps_errorbars(:, 1)),squeeze(bootstraps_errorbars(:, 2))-mean_stat, 'k.')
            if(statistic~=length(stat_names))
                xticklabels([])
                set(gca,'xtick',[])
            else
                xticks(1:(3*num_models));
                xticklabels(allmodel_xticklabels)
                xtickangle(20)
            end
            set(gca,'FontSize',9)
            ylabel(stat_names(statistic), 'FontSize', 9)
            xlim([0, 3*num_models+0.7])
        end
    end

    % Create a .mat file with the delta NLL, AIC, and BIC tables
    AllData_ModelComparison_FinalTables = cell(1,3); 
    for stat=1:3
        final_table = zeros(num_models*3,3);
        final_table(:,1) = bootstraps_errorbars_allstats{stat}(:,1); % 2.5% percentile
        final_table(:,2) = allstats{stat}; % Sum
        final_table(:,3) = bootstraps_errorbars_allstats{stat}(:,2); % 97.5% percentile
        AllData_ModelComparison_FinalTables{stat} = final_table;
    end
end


%%
function [out_struct_allmodels] = alldata_modelcomparison_visualize_helper(causal_inf_strategy, model_path, data_path)
    if(nargin==0)
        causal_inf_strategy = "ModelSelection";
        is_plot=false;
    elseif(nargin==1)
        is_plot=false;
    end
    num_subjects = 15;
    
    % Get number of observations for BIC
    load(data_path+"BAV_data.mat");
    load(data_path+"BC_data.mat");
    load(data_path+"data_stratified_UV.mat");
    load(data_path+"data_stratified_UA.mat");
    data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
    data_UA = data_stratified_to_data(data_stratified_UA, false, false);
    UAV_data = cell(1,num_subjects);
    for i=1:num_subjects
        data_UA{i}(:,3) = 4;
        UAV_data{i} = [data_UV{i}; data_UA{i}];
    end
    n_data = zeros(1,15);
    n_data_bytasktype = zeros(4,15);
    for subjidx=1:num_subjects
        n_data(subjidx) = length(BAV_data{subjidx}) + length(BC_data{subjidx}) + length(UAV_data{subjidx}); 
        n_data_bytasktype(1,subjidx) = length(data_UV{subjidx});
        n_data_bytasktype(2,subjidx) = length(data_UA{subjidx});
        n_data_bytasktype(3,subjidx) = length(BC_data{subjidx});
        n_data_bytasktype(4,subjidx) = length(BAV_data{subjidx});
    end
    % save("NumTrials_allsubjects",'n_data_bytasktype','n_data')

    colors = brewermap(10,"Accent");
    colors = colors([1,2,3,5:10],:);
    
    % AllData parametric model fit NLLs for this causal_inf_strategy
    noises = ["exp","exp","constant","exp","constant"]
    priors = ["GaussianLaplaceBothFixedZero", "SingleGaussian", "GaussianLaplaceBothFixedZero", "TwoGaussiansBothFixedZero","SingleGaussian"];
    causal_inf = ["ModelSelection","ModelAveraging","ProbMatching"];
    num_models = length(priors);
    fitted_last3_params = zeros(num_models,15,3);
    NLLs_param = zeros(num_models,15);
    num_parametric_model_params = zeros(num_models,1);
    idx=0;
    for i=1:num_models 
            idx = idx+1;
            filename = "fittedparams_All_UBresc_"+noises(i)+"-"+priors(i)+"-"+causal_inf_strategy+"_rescalefree_lapseUniform.mat";
            load(model_path + filename)
            NLLs_param(idx,:) = min(F_vals,[],2);
            num_parametric_model_params(idx) = length(theta_fitted(1,:));
            if(noises(i)=="constant") % Remove place filler zeros in theta_fitted for k_vis, k_aud.
                num_parametric_model_params(idx) = num_parametric_model_params(idx)-2;
            end
    end
    [NLL_param_min, NLL_param_min_idx] = min(NLLs_param, [],1);

    % LiftedSemiparam model for this causal_inf_strategy
    filename_basis = "fittedparams_All_UBresc_SemiparamInspired_"+causal_inf_strategy+"_rescalefree_lapseUniform.mat";
    load(model_path + filename_basis)
    num_LiftedSemiparam_params = 9;
    [F_min_val, F_min_idx] = min(F_vals,[],1);

    % Plot best LiftedSemiparam fit NLL - best param fit NLL
    NLLs = [NLL_param_min;F_min_val;NLLs_param];

    % AIC
    AIC_LiftedSemiparam_best = 2.*(num_LiftedSemiparam_params+F_min_val);
    AIC_param_best = 2.*(num_parametric_model_params(NLL_param_min_idx)'+NLL_param_min);
    AICs = [AIC_param_best; AIC_LiftedSemiparam_best];
    for model=1:num_models
        AICs = [AICs; 2.*(num_parametric_model_params(repmat(model,1,num_subjects))'+NLLs_param(model,:))];
    end

    % BIC
    %load("NumTrials_allsubjects")
    BIC_LiftedSemiparam_best = 2.*F_min_val + num_LiftedSemiparam_params.*log(n_data);
    BIC_param_best = 2.*NLL_param_min + num_parametric_model_params(NLL_param_min_idx)'.*log(n_data);

    BICs = [BIC_param_best; BIC_LiftedSemiparam_best];
    for model=1:num_models
        BICs = [BICs; 2.*NLLs_param(model,:) + num_parametric_model_params(repmat(model,1,num_subjects))'.*log(n_data)];
    end

    out_struct_allmodels.NLLs = NLLs;
    out_struct_allmodels.AICs = AICs;
    out_struct_allmodels.BICs = BICs;
end


%% Plot nonparam fitted sigma(s) and p(s) shapes
function [] = semiparam_sigmafun_prior_visualization(fontsize, figspec, model_path)
    num_subjects = 15;
    num_params = 40;
    num_iters = 81*15;
    num_inits = 81;
    
    colors = brewermap(12,"Set3");
    colors2 = brewermap(8,"Set2");
    colors(2,:) = colors2(2,:);
    colors = [colors; colors2([4,7,8],:)];
    %colors = colors .* 0.9;

    filename_basis = "fittedparams_UJoint_Semiparam_rescalefree_lapseUniform.mat";
    load(model_path + filename_basis);
    theta_fitted_cmaes = theta_fitted;
    F_vals_cmaes = F_vals;

    ModelComponents.SPivot = [0,0.1,0.3,1,2,4,6,8,10,15,20,45];

    % Plot function shapes
    s_pivot = ModelComponents.SPivot;
    s_pivot_full = [-fliplr(s_pivot(2:end)), s_pivot];
    s_fine = linspace(0,15,2^7);
    s_fine_full = linspace(0,45,2^9);
    num_pivots = length(s_pivot);
    
    figure('Position', figspec);
    set(gcf, 'Color', 'w')
    T = tiledlayout(2,10,'TileSpacing','compact', "Padding","none");

    linewidth = 1;
    for fun_idx=1:3
        t = nexttile(T,[1,5]);
        set(t,'TickDir','out');
        hold(t,'on')
        pl = get(t, 'Position');
        switch fun_idx
            case 1
                h = axes('Parent', gcf, 'Position', [pl(1)+pl(3)*.61 pl(2)+pl(4)*0.72 pl(3)*0.35 pl(3)*0.35.*3/4]);
            case 2
                h = axes('Parent', gcf, 'Position', [pl(1)+pl(3)*.65 pl(2)+.07 pl(3)*0.33 pl(3)*0.35.*3/4]);
            case 3
                h = axes('Parent', gcf, 'Position', [pl(1)+pl(3)*.61 pl(2)+pl(4)*0.75 pl(3)*0.35 pl(3)*0.35.*3/4]);
        end
        %box(h,'on');
        hold(h,'on')
        set(h,'TickDir','out');
        for subj=1:num_subjects
            %[fun_idx, subj]

            theta=squeeze(theta_fitted_cmaes(subj,:));
            color = colors(subj,:);
            switch fun_idx
                case 1
                    sigma_fun_vis_rel_high_pivots = cumsum(theta(1:num_pivots));
                    sigma_fun_vis_rel_high_pivots = [fliplr(sigma_fun_vis_rel_high_pivots(2:end)), sigma_fun_vis_rel_high_pivots];
                    sigma_fun_vis = @(s)  min([repmat(45,length(s),1)' ; interp1(s_pivot_full, exp(sigma_fun_vis_rel_high_pivots), s, 'pchip')], [], 1);       
                    p=plot(t,s_fine, sigma_fun_vis(s_fine),'-','Color', color, 'LineWidth',linewidth);
                    p.Color(4)=0.9;
                    scatter1 = scatter(t,s_pivot_full, exp(sigma_fun_vis_rel_high_pivots),'o','MarkerFaceColor',color,'MarkerEdgeColor',color); 
                    %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
                    scatter1.SizeData = linewidth.*5;
                    
                    ylabel(t,"$\sigma_{\mathrm{V}}(s)$", 'interpreter','latex', 'FontSize', fontsize)
                    xlabel(t,"Visual stimulus location (\circ)", 'FontSize', fontsize)
                    xlim(t,[0,15])
                    ylim(t,[0,6.1])
                    xl = xlim(t); yl = ylim(t);
                    yticks(t,0:6)
                    t.XAxis.FontSize = 9;
                    t.YAxis.FontSize = 9;
                    
                    p=plot(h, s_fine_full, sigma_fun_vis(s_fine_full),'-', 'Color',color);
                    p.Color(4)=0.9;
                    scatter1 = scatter(h, s_pivot_full, min(45,exp(sigma_fun_vis_rel_high_pivots)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color); 
                    %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
                    scatter1.SizeData = 1;
                    rectangle('Position',[xl(1) yl(1) xl(2)-xl(1) yl(2)-yl(1)])
                    %ylabel(h,"$\sigma_{\mathrm{V}}(s)$", 'interpreter','latex', 'FontSize', fontsize)
                    xlim(h,[0,45])
                    ylim(h,[0,45])
                    xticks(h,0:15:45)
                    yticks(h,0:15:45)
                    h.XAxis.FontSize = 9;
                    h.YAxis.FontSize = 9;

                case 2
                    sigma_fun_aud_pivots = cumsum(theta((num_pivots+1):(2*num_pivots)));
                    sigma_fun_aud_pivots = [fliplr(sigma_fun_aud_pivots(2:end)), sigma_fun_aud_pivots];
                    sigma_fun_aud = @(s)  min([repmat(45,length(s),1)' ; interp1(s_pivot_full, exp(sigma_fun_aud_pivots), s, 'pchip')], [], 1);       
                    p=plot(t,s_fine, sigma_fun_aud(s_fine),'-', 'Color',color, 'LineWidth', linewidth);
                    p.Color(4)=0.9;
                    scatter1 = scatter(t,s_pivot_full, exp(sigma_fun_aud_pivots),'o','MarkerFaceColor',color,'MarkerEdgeColor',color); 
                    %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
                    scatter1.SizeData = linewidth.*5;
                    ylabel(t,"$\sigma_{\mathrm{A}}(s)$", 'interpreter','latex', 'FontSize', fontsize) 
                    xlabel(t,"Auditory stimulus location (\circ)", 'FontSize', fontsize)
                    xlim(t,[0,15])
                    ylim(t,[0,6.1])
                    xl = xlim(t); yl = ylim(t);
                    yticks(t,0:6)
                    t.XAxis.FontSize = 9;
                    t.YAxis.FontSize = 9;
                    
                    p=plot(h, s_fine_full, sigma_fun_aud(s_fine_full),'-', 'Color',color);
                    p.Color(4)=0.9;
                    scatter1 = scatter(h, s_pivot_full, min(45,exp(sigma_fun_aud_pivots)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color); 
                    %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
                    scatter1.SizeData = 1;
                    rectangle('Position',[xl(1) yl(1) xl(2)-xl(1) yl(2)-yl(1)])
                    %ylabel(h,"$\sigma_{\mathrm{A}}(s)$", 'interpreter','latex', 'FontSize', fontsize)
                    xlim(h,[0,45])
                    ylim(h,[0,45])
                    xticks(h,0:15:45)
                    yticks(h,0:15:45)
                    h.XAxis.FontSize = 9;
                    h.YAxis.FontSize = 9;
                case 3    
                    prior_pivots = cumsum([1,theta((2*num_pivots+1):(3*num_pivots-1))]);
                    prior_pivots = [fliplr(prior_pivots(2:end)), prior_pivots];
                    prior = @(s) exp(interp1(s_pivot_full, prior_pivots, s, 'pchip'));
                    
                    s_fine_full_width = s_fine_full(2) - s_fine_full(1);
                    normalization_constant = 1./(qtrapz(prior(s_fine_full).*s_fine_full_width));

                    p=plot(t,s_fine_full, (prior(s_fine_full).*normalization_constant),'-', 'Color',color, 'LineWidth', linewidth);
                    p.Color(4)=0.9;
                    %scatter(s_pivot_full, prior_pivots,"o", "MarkerFaceAlpha",0.1)
                    scatter1 = scatter(t,s_pivot_full, exp(prior_pivots).*normalization_constant,'o','MarkerFaceColor',color,'MarkerEdgeColor',color); 
                    %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
                    scatter1.SizeData = linewidth.*5;
                    ylabel(t,"$p(s)$", 'interpreter','latex', 'FontSize', fontsize) 
                    xlabel(t,"Visual/Auditory stimulus location (\circ)", 'FontSize', fontsize)
                    xlim(t,[0,3])
                    ylim(t,[0,5])
                    xl = xlim(t); yl = ylim(t);
                    xticks(t,0:1:3);
                    yticks(t, 0:1:5);
                    t.XAxis.FontSize = 9;
                    t.YAxis.FontSize = 9;
                    
                    p=plot(h, s_fine_full, log(prior(s_fine_full).*normalization_constant),'-', 'Color',color);
                    p.Color(4)=0.9;
                    scatter1 = scatter(h, s_pivot_full, log(exp(prior_pivots).*normalization_constant),'o','MarkerFaceColor',color,'MarkerEdgeColor',color); 
                    %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
                    scatter1.SizeData = 1;
                    %rectangle('Position',[xl(1) yl(1) xl(2)-xl(1) yl(2)-yl(1)])
                    xlim(h,[0,45])
                    ylim(h,[-20,3])
                    xticks(h,0:15:45)
                    yticks(h, -20:10:0);
                    h.XAxis.FontSize = 9;
                    h.YAxis.FontSize = 9;
                    ylabel(h,"$\log p(s)$", 'interpreter','latex', 'FontSize', fontsize) 
            end  
        end
        ax.XAxis.FontSize = fontsize;
        ax.YAxis.FontSize = fontsize;
    end
    
    x_labels_pos = (3*num_pivots):num_params;
    x_labels = "$"+{"\alpha_\mathrm{med}", "\alpha_\mathrm{low}", "\lambda","\sigma_\mathrm{motor}","\rho_\mathrm{A}"}+"$";

    t=nexttile(T,[1,3]);
    boxplot(theta_fitted_cmaes(:,x_labels_pos([1,2,5])), 'Color','k')
    hold on
    idx=0
    for param=[1,2,5]
        idx=idx+1
        for subj=1:num_subjects
            scatter1 = scatter(repmat(idx,1), theta_fitted_cmaes(subj,x_labels_pos(param)),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
            %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
            scatter1.SizeData = 5;
        end
    end
    xticks(1:3);
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    xticklabels(x_labels([1,2,5]));
    %xtickangle(45);
    ax=gca;
    ax.XAxis.FontSize = fontsize;
    xlim([1-0.5,idx+0.5])
    ylim([0, Inf])
    yticks(0:1:5)
    t.XAxis.FontSize = 9;
    t.YAxis.FontSize = 9;
    box off 
    set(gca,'TickDir','out');
        
    t=nexttile(T);
    boxplot(theta_fitted_cmaes(:,x_labels_pos(3)), 'Color','k')
    hold on
    for subj=1:num_subjects
        scatter1 = scatter(1, theta_fitted_cmaes(subj,x_labels_pos(3)),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
        %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
        scatter1.SizeData = 5;
    end
    % scatter1 = scatter(repmat(1,num_subjects,1), theta_fitted_cmaes(:,x_labels_pos(3)),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    % scatter1.SizeData = 5;
    xticks([1])
    xtickangle(0);
    yticks(0:0.005:0.015)
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    xticklabels(x_labels(3));
    ax=gca;
    ax.XAxis.FontSize = fontsize;
    xlim([0.5,1.5])
    ylim([0, 0.015])
    t.XAxis.FontSize = 9;
    t.YAxis.FontSize = 9;
    set(gca,'TickDir','out');
    box off 
    
    t=nexttile(T);
    boxplot(theta_fitted_cmaes(:,x_labels_pos(4)), 'Color','k')
    hold on
    for subj=1:num_subjects
        scatter1 = scatter(1, theta_fitted_cmaes(subj,x_labels_pos(4)),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
        %scatter1.MarkerFaceAlpha = .2; scatter1.MarkerEdgeAlpha = .2; 
        scatter1.SizeData = 5;
    end
    % scatter1 = scatter(repmat(1,num_subjects,1), theta_fitted_cmaes(:,x_labels_pos(4)),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    % scatter1.SizeData = 5;
    xticks([1])
    yticks(0:0.1:0.5);
    xtickangle(0);
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    xticklabels(x_labels(4));
    ax=gca;
    ax.XAxis.FontSize = fontsize;
    xlim([0.5,1.5])
    ylim([0, 0.5])
    t.XAxis.FontSize = 9;
    t.YAxis.FontSize = 9;
    set(gca,'TickDir','out');
    box off 
    
    %set(gca,'fontsize', fontsize)
end

%%
function [] = sigmafun_prior_examples(fontsize)
    s_grid = -45:0.1:45;
    sigma0_vals = [0.5,1,2,3];
    colors = brewermap(12,"Set1");

    hetero_type = "exp";
    sigma_fun_constant = heterotype_to_sigmafun("constant");
    sigma_fun_exp = heterotype_to_sigmafun("exp");

    subplot(2,3,1)
    set(gca,'TickDir','out');
    hold on
    for i=1:length(sigma0_vals)
        plot(s_grid, repmat(sigma0_vals(i), length(s_grid),1), "-", 'Color', colors(i,:));
    end
    lg = legend("$\sigma_0="+sigma0_vals+"$", 'Interpreter', 'latex', 'FontSize', fontsize);
    set(lg,'Box','off')
    ylim([0,6])
    yticks(0:1:6)
    xlabel("Visual/Auditory stimulus location (\circ)", 'FontSize', fontsize)
    ylabel("$\sigma(s)$", 'Interpreter', 'latex', 'FontSize', fontsize)
    title("Constant sensory noise", 'FontSize', fontsize+1)
    

    sigma0_vals = [0.5,1,1,1,3];
    k1_vals = [1,1,2,1,1];
    k2_vals = [0.1,0.1,0.1,0.5 0.5];
    subplot(2,3,2)
    set(gca,'TickDir','out');
    hold on
    for i=1:length(sigma0_vals)
        plot(s_grid, sigma_fun_exp(s_grid,sigma0_vals(i), [k1_vals(i),k2_vals(i)]), "-", 'Color', colors(i,:));
    end
    lg = legend("$\sigma_0="+sigma0_vals+", k_1="+k1_vals+", k_2="+k2_vals+"$", 'Interpreter', 'latex', 'FontSize', fontsize);
    set(lg,'Box','off')
    lg.Position(1:2) = [0.63,0.8];
    ylim([0,6])
    yticks(0:1:6)
    xlabel("Visual/Auditory stimulus location (\circ)", 'FontSize', fontsize)
    title("Exponential sensory noise", 'FontSize', fontsize+1)
    
    % Priors
    sigma_s_vals = [3, 5, 8, 10];
    subplot(2,3,4)
    set(gca,'TickDir','out');
    hold on
    for i=1:length(sigma_s_vals)
        plot(s_grid, normpdf(s_grid, 0, sigma_s_vals(i)), "-", 'Color', colors(i,:))
    end
    lg = legend("$\sigma_s="+sigma_s_vals+"$", 'Interpreter', 'latex', 'FontSize', fontsize);
    set(lg,'Box','off')
    lg.Position(1) = 0.25;
    xlabel("Visual/Auditory stimulus location (\circ)", 'FontSize', fontsize)
    ylabel("$p(s)$", 'Interpreter', 'latex', 'FontSize', fontsize)
    title("SingleGaussian prior", 'FontSize', fontsize+1)
    ylim([0,0.15])
    yticks(0:0.05:0.15)

    sigma_s_vals = [8, 8, 8, 15].*2;
    b_vals = [1,2,1,2].*2;
    w_vals = [0.3, 0.3, 0.5, 0.5];
    subplot(2,3,6)
    set(gca,'TickDir','out');
    hold on
    for i=1:length(sigma_s_vals)
        plot(s_grid, (1-w_vals(i)).*normpdf(s_grid, 0, sigma_s_vals(i)) + w_vals(i).*1./(2.*b_vals(i)).*exp(-abs(s_grid)./b_vals(i)), "-", 'Color', colors(i,:))
    end
    xlabel("Visual/Auditory stimulus location (\circ)", 'FontSize', fontsize)
    lg = legend("$\sigma_s="+sigma_s_vals+", b="+b_vals+", w="+w_vals+"$", 'Interpreter', 'latex', 'FontSize', fontsize);
    set(lg,'Box','off')
    lg.Position(1:2) = [0.69,0.48];
    title("GaussianLaplace prior", 'FontSize', fontsize+1)
    ylim([0,0.15])
    yticks(0:0.05:0.15)

    sigma_s_vals = [3,5,8,10];
    sigma_s2_vals = [8, 8, 8, 15].*2 - sigma_s_vals;
    w_vals = [0.3, 0.3, 0.7, 0.7];
    subplot(2,3,5)
    set(gca,'TickDir','out');
    hold on
    for i=1:length(sigma_s_vals)
        plot(s_grid, (1-w_vals(i)).*normpdf(s_grid, 0, sigma_s_vals(i)) + w_vals(i).*normpdf(s_grid, 0, sigma_s_vals(i)+sigma_s2_vals(i)), "-", 'Color', colors(i,:))
    end
    ylim([0,0.15])
    xlabel("Visual/Auditory stimulus location (\circ)", 'FontSize', fontsize)
    lg = legend("$\sigma_s="+sigma_s_vals+", \sigma_{\Delta}="+sigma_s2_vals+", w="+w_vals+"$", 'Interpreter', 'latex', 'FontSize', fontsize);
    set(lg,'Box','off')
    lg.Position(1:2) = [0.42, 0.348];
    title("TwoGaussians prior", 'FontSize', fontsize+1)
    ylim([0,0.15])
    yticks(0:0.05:0.15)
    
    %set(gca,'fontsize', fontsize)
end


%%
function [] = allindvsubjplots_to_onesubjplot(save_name, subjidx, fitted_on_all_data, fontsize, figspecs, figpath)
    % This function assumes that the individual-level plots have been saved
    % as .fig files.
    close all;
    
    if(~fitted_on_all_data)
        F1 = openfig(figpath + save_name + "_Individualmean.fig");
        t1 = nexttile(subjidx);
        ax1=gca;
        F2 = openfig(figpath + save_name + "_IndividualSD.fig");
        t2 = nexttile(subjidx);
        ax2=gca;
        
        figure('Position', figspecs);
        set(gcf, 'Color', 'w')
        T=tiledlayout(1,2,'Padding', 'tight', 'TileSpacing', 'tight');
        t1 = nexttile(1);
        set(gca,'TickDir','out');
        hold on
        plot([-20,20],[-20,20],"k--",'HandleVisibility','off');
        fig1 = get(ax1,'children');
        copyobj(fig1, t1);
        xlabel("Stimulus location (\circ)", 'FontSize', fontsize)
        ylabel("Bias (\circ)", 'FontSize', fontsize)
        ylim([-20,20])
        xticks(-20:10:20)
        set(gca,"FontSize",9)
        
        t2 = nexttile(2);
        fig2 = get(ax2,'children');
        set(gca,'TickDir','out');
        copyobj(fig2, t2);
        xlabel("Stimulus location (\circ)", 'FontSize', fontsize)
        ylabel("SD of location response (\circ)", 'FontSize', fontsize)
        h = findall(gca, 'LineStyle', '-');
        for i=1:4
            h(i).HandleVisibility="off";
        end
        lg = legend("Visual (high reliability)","Visual (med. reliability)", "Visual (low reliability)", "Auditory");
        set(lg,'Box','off')
        lg.FontSize = max(9,fontsize-1);
        lg.Location="northeast";
        lg.ItemTokenSize(1) = 10;
        xticks(-20:10:20)
        set(gca,"FontSize",9)
        
    else
        
        F1 = openfig(figpath + save_name + "-UAV_Individualmean.fig");
        t1 = nexttile(subjidx);
        ax1=gca;
        F2 = openfig(figpath + save_name + "-UAV_IndividualSD.fig");
        t2 = nexttile(subjidx);
        ax2=gca;
    
        F3 = openfig(figpath + save_name + "-BC_Individual.fig");
        ax3_center = F3.Children.Children((end-subjidx+1)).Children(2);
        ax3_periphery = F3.Children.Children((end-subjidx+1)).Children(1);
        
        F4 = openfig(figpath + save_name + "-BV_Individual.fig");
        ax4_right = F4.Children.Children((end-subjidx+1)).Children(1);
        ax4_center = F4.Children.Children((end-subjidx+1)).Children(2);
        ax4_left = F4.Children.Children((end-subjidx+1)).Children(3);
        
        F5 = openfig(figpath + save_name + "-BA_Individual.fig");
        ax5_right = F5.Children.Children((end-subjidx+1)).Children(1);
        ax5_center = F5.Children.Children((end-subjidx+1)).Children(2);
        ax5_left = F5.Children.Children((end-subjidx+1)).Children(3);
        
        %% Move to new plot
        figure('Position', figspecs);
        set(gcf, 'Color', 'w')
        T=tiledlayout(2,12,'Padding', 'tight', 'TileSpacing', 'tight');
        
        
        t12=tiledlayout(T,1,2, 'Padding','none','TileSpacing','tight');
        t12.Layout.Tile = 1;
        t12.Layout.TileSpan = [1 6];

        t1 = nexttile(t12);
        hold on
        set(gca,'TickDir','out');
        plot([-20,20],[-20,20],"k--",'HandleVisibility','off');
        fig1 = get(ax1,'children');
        copyobj(fig1, t1);
        xlabel("Stimulus location (\circ)", 'FontSize', fontsize)
        ylabel("Mean location response (\circ)", 'FontSize', fontsize)
        ylim([-20,20])
        ttl = title('(a)', "Fontsize", 10);
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = -0.3; % use negative values (ie, -0.1) to move further left
        ttl.HorizontalAlignment = 'left'; 
        xticks(-20:10:20)
        set(gca,"FontSize",9)
        
        t2 = nexttile(t12);
        set(gca,'TickDir','out');
        fig2 = get(ax2,'children');
        copyobj(fig2, t2);
        xlabel("Stimulus location (\circ)", 'FontSize', fontsize)
        ylabel("SD of location response (\circ)", 'FontSize', fontsize)
        ylim([0,9])
        h = findall(gca, 'LineStyle', '-');
        for i=1:4
            h(i).HandleVisibility="off";
        end
        lg = legend("Visual (high rel.)","Visual (med. rel.)", "Visual (low rel.)", "Auditory");
        set(lg,'Box','off')
        lg.FontSize = 9;
        lg.ItemTokenSize(1) = 10;
        xticks(-20:10:20)
        set(gca,"FontSize",9)
        
        % BC
        %t3 = nexttile([1,2]);
        t3=tiledlayout(T,1,2, 'Padding','none','TileSpacing','compact');
        t3.Layout.Tile = 7;
        t3.Layout.TileSpan = [1 6];
        xlabel(t3,"Stimulus location disparity, {\its}_A– {\its}_V (\circ)", 'FontSize',fontsize)
        ylabel(t3,{"{\rm \fontsize{9} {Proportion responding "+ '"'+'same'+ '"'+"}}"}, 'FontSize',fontsize);
        
        BC_strat_names = ["Center", "Periphery"];
        for strats=1:2
            tt = nexttile(t3);
            set(gca,'TickDir','out');
            hold on;
            if(strats==1)
                fig31 = get(ax3_center,'children');
                copyobj(fig31, tt);
                ttl = title('(b)', "Fontsize", 10);
                ttl.Units = 'Normalize'; 
                ttl.Position(1) = -0.3; % use negative values (ie, -0.1) to move further left
                ttl.HorizontalAlignment = 'left'; 
                subtitle(tt,"Center",'Fontsize', fontsize, 'FontWeight','bold')
                yticks(0:0.2:1)
            else
                fig32 = get(ax3_periphery,'children');
                copyobj(fig32, tt);
                subtitle(tt,"Periphery",'Fontsize', fontsize, 'FontWeight','bold')
                yticks([])
            end
            h = findall(gca, 'LineStyle', '-');
            for i=1:3
                h(i).HandleVisibility="off";
            end
            xlim([-30,30])
            ylim([0,1])
            xticks(-30:15:30)
            xtickangle(0)
            set(gca,"FontSize",9)
            
            lg = legend({"High vis. rel.","Med. vis. rel.","Low vis. rel."});
            set(lg,'Box','off')
            lg.FontSize = 9;
            lg.Location="south";
            lg.ItemTokenSize(1) = 10;
        end
        
        BAV_strat_names = ["Left","Center","Right"];
        t4=tiledlayout(T,1,3, 'Padding','none','TileSpacing','compact');
        t4.Layout.Tile = 13;
        t4.Layout.TileSpan = [1 6];
        xlabel(t4, "Stimulus location disparity, {\its}_A– {\its}_V (\circ)", 'FontSize',fontsize)
        ylabel(t4,"{\rm \fontsize{10} {Visual bias (\circ)}}");
        for strats=1:3
            tt = nexttile(t4);
            set(gca,'TickDir','out');
            hold on;
            if(strats==1)
                fig4 = get(ax4_left,'children');
                copyobj(fig4, tt);
                ttl = title('(c)', "Fontsize", 10);
                ttl.Units = 'Normalize'; 
                ttl.Position(1) = -0.4; % use negative values (ie, -0.1) to move further left
                ttl.HorizontalAlignment = 'left';
                subtitle(tt,"Left",'Fontsize', fontsize, 'FontWeight','bold')
            elseif(strats==2)
                yticks([])
                fig4 = get(ax4_center,'children');
                copyobj(fig4, tt);
                subtitle(tt,"Center",'Fontsize', fontsize, 'FontWeight','bold')
            else
                yticks([])
                fig4 = get(ax4_right,'children');
                copyobj(fig4, tt);
                subtitle(tt,"Right",'Fontsize', fontsize, 'FontWeight','bold')
            end
            
            h = findall(gca, 'LineStyle', '-');
            for i=1:3
                h(i).HandleVisibility="off";
            end
            xlim([-35,35])
            ylim([-15,15])
            xticks(-30:15:30)
            xtickangle(0)
            set(gca,"FontSize",9)
            lg = legend({"High vis. rel.","Med. vis. rel.","Low vis. rel."});
            set(lg,'Box','off')
            lg.FontSize = 9;
            lg.Location="north";
            lg.ItemTokenSize(1) = 10;
        end
        
        t5=tiledlayout(T,1,3, 'Padding','none','TileSpacing','compact');
        t5.Layout.Tile = 19;
        t5.Layout.TileSpan = [1 6];
        xlabel(t5, "Stimulus location disparity, {\its}_A– {\its}_V (\circ)", 'FontSize',fontsize)
        ylabel(t5,"{\rm \fontsize{10} {Auditory bias (\circ)}}");
        for strats=1:3
            tt = nexttile(t5);
            set(gca,'TickDir','out');
            hold on;
            if(strats==1)
                fig5 = get(ax5_left,'children');
                copyobj(fig5, tt);
                ttl = title('(d)', "Fontsize", 10);
                ttl.Units = 'Normalize'; 
                ttl.Position(1) = -0.4; % use negative values (ie, -0.1) to move further left
                ttl.HorizontalAlignment = 'left';
                subtitle(tt,"Left",'Fontsize', fontsize, 'FontWeight','bold')
            elseif(strats==2)
                yticks([])
                fig5 = get(ax5_center,'children');
                copyobj(fig5, tt);
                subtitle(tt,"Center",'Fontsize', fontsize, 'FontWeight','bold')
            else
                yticks([])
                fig5 = get(ax5_right,'children');
                copyobj(fig5, tt);
                subtitle(tt,"Right",'Fontsize', fontsize, 'FontWeight','bold')
            end
            h = findall(gca, 'LineStyle', '-');
            for i=1:3
                h(i).HandleVisibility="off";
            end
            set(gca,"FontSize",9)
            xlim([-35,35])
            ylim([-15,15])
            xticks(-30:15:30)
            xtickangle(0)
            
            lg = legend({"High vis. rel.","Med. vis. rel.","Low vis. rel."});
            set(lg,'Box','off')
            lg.FontSize = 9;
            lg.Location="north";
            lg.ItemTokenSize(1) = 10;
            switch strats
                case 2
                    lg.Position(1) = 0.46;
                    lg.Position(2) = 0.7725;
                case 3
                    lg.Position(1) = 0.74;
                    lg.Position(2) = 0.7725;
            end
        end
    end
end