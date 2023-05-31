clear all; close all;
%cd('C:\Users\liu_s\Audiovisual-causal-inference')
model_path = "fast_modelfit_example\";
model_path_temp = model_path+"temp\";
data_path = "data\";
num_subjects = 15;

%% Fit a model on UV+UA data quickly, using one random init per subject.
prior_type = "GaussianLaplaceBothFixedZero"; % Can also be "SingleGaussian", "TwoGaussiansBothFixedZero"
hetero_type = "exp"; % Can also be "constant";
lapse_type = "Uniform"; % Can also be "Gaussian";
rescale_aud = "free"; % Can also be "4/3" or "free";
num_inits_persubj = 1;
parfor iter = 1:(num_inits_persubj * num_subjects)
    fit_ujointmodel_parametric(iter,prior_type,hetero_type, lapse_type, rescale_aud, num_inits_persubj, data_path, model_path_temp)
end
merge_parametric_ujointfits_files_fast(prior_type,hetero_type,lapse_type,rescale_aud, num_inits_persubj, model_path, model_path_temp)

%% Visualize the fit.
manuscript_ujoint_respdistrvisualization(prior_type, hetero_type, rescale_aud, fontsize, [0 0 figsize(4)*4/3 figsize(4)], model_path, plot_lapse, lapse_type);


%% Helper function to merge different subjects/inits for the same model.
function [] = merge_parametric_ujointfits_files_fast(prior_type,hetero_type,lapse_type,rescale_aud, num_inits_persubj, model_path, model_path_temp)
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