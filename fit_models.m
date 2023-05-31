clear all; close all;
cd('C:\Users\liu_s\Audiovisual-causal-inference')
model_path = "\modelfits\";
model_path_temp = model_path+"temp\";
data_path = "\data\";
num_subjects = 15;

%% Parametric models on UA+UV data
prior_type = "SingleGaussian"; % Can also be "SingleGaussian", "TwoGaussiansBothFixedZero"
hetero_type = "constant"; % Can also be "exp";
lapse_type = "Uniform"; % Can also be "Gaussian";
rescale_aud = "1"; % Can also be "4/3" or "free";
num_inits_persubj = 10;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_ujointmodel_parametric(iter,prior_type,hetero_type, lapse_type, rescale_aud, num_inits_persubj, data_path, model_path_temp)
end
merge_parametric_ujointfits_files(prior_type,hetero_type,lapse_type,rescale_aud, num_inits_persubj, model_path, model_path_temp)

%% Semiparam models on UA+UV data
num_inits_persubj = 81;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_ujointmodel_semiparam(iter, num_inits_persubj, data_path, model_path_temp);
end
merge_semiparam_ujointfits_files(num_iters_persubj, model_path, model_path_temp)

%% SemiparamInsp models on all data
causal_inf_strategy = "ProbMatching"; % Can also either be "ModelSelection" or "ModelAveraging"
num_inits_persubj = 100;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_alldatamodel_semiparaminsp_resc(iter,causal_inf_strategy, num_inits_persubj, data_path, model_path, model_path_temp)
end
merge_semiparaminsp_alldatafits_files(causal_inf_strategy, num_iters_persubj, model_path, model_path_temp)

%% Parametric models on all data
prior_type = "GaussianLaplaceBothFixedZero"; % Can also be "SingleGaussian", "TwoGaussiansBothFixedZero"
hetero_type = "exp"; % Can also be "constant";
causal_inf_strategy = "ProbMatching"; % Can also either be "ModelSelection" or "ModelAveraging"
num_inits_persubj = 10;
for iter = 1:(num_inits_persubj * num_subjects)
    fit_alldatamodel_parametric_resc(iter,prior_type,hetero_type,causal_inf_strategy, num_inits_persubj, data_path, model_path_temp)
end
merge_parametric_alldatafits_files(prior_type,hetero_type,causal_inf_strategy, num_inits_persubj, model_path, model_path_temp)








%% Helper functions to merge different subjects/inits for the same model.
function [] = merge_parametric_ujointfits_files(prior_type,hetero_type,lapse_type,rescale_aud, num_inits_persubj, model_path, model_path_temp)
    filename = 'fittedparams_UJoint_'+hetero_type+"-"+prior_type+"_rescale"+rescale_aud+"_lapse"+lapse_type;
    load(model_path_temp+filename+"_1.mat");
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
        load(model_path_temp+filename+"_"+num2str(i));
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
    
    save(model_path+filename, "Theta_fitted", "T_Ends","F_vals", "theta_fitted");
end

%%
function [] = merge_semiparam_ujointfits_files(num_iters_persubj, model_path, model_path_temp)
    filename_basis = "fittedparams_UJoint_Semiparam_rescalefree_lapseUniform";
    num_subjects = 15;
    num_iters = num_iters_persubj*num_subjects;

    Theta_fitted = zeros(num_iters_persubj, num_subjects, num_params);
    F_vals = zeros(num_iters_persubj, num_subjects);
    T_Ends = zeros(num_iters_persubj, num_subjects);
    theta_inits = zeros(num_iters_persubj, num_subjects, num_params);

    for i=1:num_iters
        load(savepath+"1"+filename_basis+"_"+i+".mat")
        init = out_struct.init; subj = out_struct.subj;
        Theta_fitted(init, subj,:) = out_struct.Theta_fitted;
        theta_inits(init, subj,:) = out_struct.theta0;
        F_vals(init, subj) = out_struct.F_val;
        T_Ends(init, subj) = out_struct.T_end;
    end
    [F_min_val, min_idx] = min(F_vals,[],1);
    theta_fitted = zeros(num_subjects, num_params);
    for i=1:num_subjects
        theta_fitted(i,:) = Theta_fitted(min_idx(i),i,:);
    end
    save(filename_basis, 'Theta_fitted','F_vals','T_Ends', 'theta_fitted','theta_inits')
end
%%

function [] = merge_semiparaminsp_alldatafits_files(causal_inf_strategy, num_iters_persubj, model_path, model_path_temp)
    if(nargin==0)
        causal_inf_strategy = "ModelSelection";
    end
    num_subjects = 15;
    num_params = 9;
    colors = brewermap(10,"Accent");
    colors = colors([1,2,3,5],:);
    filename_basis = "fittedparams_All_UBresc_SemiparamInspired_"+causal_inf_strategy+"_rescalefree_lapseUniform";

    num_iters = num_subjects * num_iters_persubj;
    num_inits = ceil(num_iters/num_subjects);

    Theta_fitted = zeros(num_inits, num_subjects, num_params);
    F_vals = zeros(num_inits, num_subjects);
    T_Ends = zeros(num_inits, num_subjects);
    theta_inits = zeros(num_inits, num_subjects, num_params);
    SigmaFunsAllSubjs = cell(1,15);
    PriorUnnormalizedAllSubjs = zeros(201,15);
    for i=1:num_iters
        load(model_path_temp+filename_basis+"_"+i+".mat")
        init = out_struct.init; subj = out_struct.subjidx;
        Theta_fitted(init, subj,:) = out_struct.theta_fitted;
        theta_inits(init, subj,:) = out_struct.theta_init;
        F_vals(init, subj) = out_struct.F_vals;
        T_Ends(init, subj) = out_struct.T_Ends;
        if(mod(i,num_inits)==1)
            PriorUnnormalizedAllSubjs(:,subj) = ModelComponents.PriorUnnormalized;
            SigmaFunsAllSubjs{subj} = ModelComponents.SigmaFuns;
        end
    end
    fields = {'PriorUnnormalized','SigmaFuns'};
    ModelComponents = rmfield(ModelComponents,fields);
    ModelComponents.PriorUnnormalizedAllSubjs = PriorUnnormalizedAllSubjs;
    ModelComponents.SigmaFunsAllSubjs = SigmaFunsAllSubjs;

    [F_min_val, min_idx] = min(F_vals,[],1);
    theta_fitted = zeros(num_subjects, num_params);
    for i=1:num_subjects
        theta_fitted(i,:) = Theta_fitted(min_idx(i),i,:);
    end

    save(model_path + filename_basis, 'Theta_fitted','F_vals','T_Ends', 'theta_fitted','theta_inits', 'ModelComponents');
end

%%
function [] = merge_parametric_alldatafits_files(prior_type,hetero_type,causal_inf_strategy, num_inits_persubj, model_path, model_path_temp)
    filename = 'fittedparams_All_UBresc_'+hetero_type+"-"+prior_type+"-"+causal_inf_strategy+"_rescalefree_lapseUniform";
    load(model_path_temp+filename+"_1.mat");
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
        load(model_path_temp+filename+"_"+num2str(i));
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
    
    save(model_path+filename, "Theta_fitted", "T_Ends","F_vals", "theta_fitted");
end