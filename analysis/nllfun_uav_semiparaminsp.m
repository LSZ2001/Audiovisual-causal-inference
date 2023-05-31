% This is identical to
% midpoint_NLLfun_comprehensive_nonparam_mixture_indv2, but removes many
% function arguments due to constraints of the cmaes() optimization
% algorithm.
function [output] = nllfun_uav_semiparaminsp(theta, data, ModelComponents, return_predictive_samples, return_response_distr, plot_consider_lapse, lapse_type, Gaussian_lapse_SD)
    if(nargin==3) % default is to return NLL, not posterior preditive samples vectorized response distribition.
        return_predictive_samples = false; return_response_distr=false; plot_consider_lapse=true; lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==4) % if return posterior predictive samples and not vectorized response distribution, default is not include lapse mixture component.
        return_response_distr=false; plot_consider_lapse=true; lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif (nargin==5) % if return vectorized response distribution (need S to be one s-value and rel-level repmat, and R to be r_grid), default is not include lapse mixture component.
        plot_consider_lapse=true; lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==6)
        lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==7)
        Gaussian_lapse_SD = NaN;
    end
    
    s_pivot = [0,0.1,0.3,1,2,4,6,8,10,15,20,45];
    s_pivot_full = [-fliplr(s_pivot(2:end)), s_pivot];
    
    lapse_range = ModelComponents.LapseRange;
    has_rescale_vis = ModelComponents.RescaleVis=="free"; % fit k_vis
    has_rescale_aud = ModelComponents.RescaleAud=="free"; % fit k_vis
    num_vis_rels = 3; % Get number of visual reliability levels
    num_rels = num_vis_rels + 1; % plus the auditory -- for sigma0(cond);
    s_grid_integrate = ModelComponents.SGridIntegrate;
   
    % Dataset parsing
    s_matrix = data;
    S = s_matrix(:,1); R_raw = s_matrix(:,2); reliabilities = s_matrix(:,3);
    LL_idx_mat = zeros(length(S),1);
        
    %% Parameter extraction
    % [scale_rel_med, scale_rel_low, lapse, sig_motor, rescale_aud, p_same,
    % Bimod_rescale_vis, Bimod_rescale_aud, UB_shared_prior_power];
    rel_multiply_factors_uniq =[1, theta(1:2), 1];
    lapse = theta(3);  sigma_motor = theta(4);
    rescale_vis = 1; rescale_aud = theta(5);
    rescales = [repmat(rescale_vis,3,1); rescale_aud];
    UB_shared_prior_power = theta(end);
    if(~plot_consider_lapse) % Do not consider the lapse component.
        lapse=0;
    end
    
    %% For response distribution visualization -- need to change the sigma_fun and prior def. 
    if(return_predictive_samples) % If return predictive samples, cannot return response distribution despite input args.
        return_response_distr=false;
        num_samps = 100; % Number of samples for this trial. 
        
        % Assume input only contains one reliability level (1 to 4)!!
        rel_level = reliabilities(1);
        rescale = rescales(rel_level);
        rel_multiply_factor = rel_multiply_factors_uniq(rel_level);

        % Get normalized prior density.
        p_s_unnorm_base = ModelComponents.PriorUnnormalized;
        p_s_unnorm = p_s_unnorm_base .^ UB_shared_prior_power;
        p_s = p_s_unnorm / qtrapz(p_s_unnorm*(s_grid_integrate(2)-s_grid_integrate(1)));
        log_prior = log(p_s);

        sigma_fun = ModelComponents.SigmaFuns{rel_level};
        x_obs_vals = S + sigma_fun(S) .* rel_multiply_factor.*randn([length(S),num_samps]);
        
        x_obs_rescale = x_obs_vals(:)' .* rescale;            
        log_likelihood = -log(sigma_fun(s_grid_integrate).*rel_multiply_factor)-0.5.*log(2*pi) - ((x_obs_rescale-s_grid_integrate).^2)./(2.*((sigma_fun(s_grid_integrate).*rel_multiply_factor).^2));

        log_expr = log_likelihood + log_prior;
        log_expr_shifted = log_expr - max(log_expr);
        expr = exp(log_expr_shifted);

        s_hats_pm = qtrapz(s_grid_integrate.*expr,1) ./ qtrapz(expr,1);
        s_hats_PM = reshape(s_hats_pm', length(S), num_samps);

        % Add motor noise
        output = s_hats_PM + sigma_motor .* randn([length(S),num_samps]);    
        % Add in lapse trials
        lapse_idx = rand([length(S), num_samps]) < lapse;
        
        switch lapse_type
            case "Uniform"
                % Random responses (uniform in -45:1:45)
                lapse_val = min(lapse_range)+ (max(lapse_range)-min(lapse_range))*rand(sum(lapse_idx(:)),1);
            case "Gaussian"
                lapse_val = Gaussian_lapse_SD .* randn(sum(lapse_idx(:)),1);
        end
        output(find(lapse_idx~=0)) = lapse_val;
        return
    end
        
%% All below is NLL evaluation for model fitting
    for i=1:num_rels 
        % Extract trials with this reliablity/modality.
        rel_trial_idx = find(reliabilities==i);
        if(isempty(rel_trial_idx)) % if no trials in the data has this reliability, skip this iteration
            continue;
        end
        S_rel = S(rel_trial_idx); R_rel = R_raw(rel_trial_idx);

        % Get reliability and sigma(s) pivots for this modality.
        rescale = rescales(i);
        rel_multiply_factor = rel_multiply_factors_uniq(i);
        sigma_fun = ModelComponents.SigmaFuns{i};       

        % Bounds of integration for p(x|s)
        s_min = min(S_rel); s_max = max(S_rel);
        x_min = s_min(1) - 5*sigma_fun(s_min(1)).*rel_multiply_factor;
        x_max = s_max(1) + 5*sigma_fun(s_max(1)).*rel_multiply_factor; 

        %% This part borrowed from parametric NLL code to make x_grid identical.
        Nx = 2^10;
        x_grid = linspace(x_min, x_max, Nx); % Need to transpose due to matrix operation x-s.
        dx = x_grid(2) - x_grid(1);
        if(dx > 0.1) % Ensure that dx is much smaller than the LB (0.25) of sigma_motor.
            dx=0.1;
            x_grid = (x_min:dx:x_max);
            Nx = length(x_grid);
        end
        
        %%
        log_p_x_given_s = -log(sigma_fun(S_rel).*rel_multiply_factor) - 0.5.*log(2*pi) - 0.5.*((x_grid-S_rel)./(sigma_fun(S_rel).*rel_multiply_factor)).^2;
        log_offsets = max(log_p_x_given_s,2);
        log_p_x_given_s = log_p_x_given_s - log_offsets;
        p_x_given_s = exp(log_p_x_given_s) .* exp(log_offsets);

        % Get likelihood 
        x_grid_rescale = x_grid .* rescale;            
        log_likelihood = -log(sigma_fun(s_grid_integrate).*rel_multiply_factor)-0.5.*log(2*pi) - ((x_grid_rescale-s_grid_integrate).^2)./(2.*((sigma_fun(s_grid_integrate).*rel_multiply_factor).^2));

        % Get normalized prior density.
        p_s_unnorm_base = ModelComponents.PriorUnnormalized;
        p_s_unnorm = p_s_unnorm_base .^ UB_shared_prior_power;
        p_s = p_s_unnorm / qtrapz(p_s_unnorm*(s_grid_integrate(2)-s_grid_integrate(1)));
        log_prior = log(p_s);
        
        % Bayes rule 
        log_expr = log_likelihood + log_prior;
        log_expr_shifted = log_expr - max(log_expr);
        expr = exp(log_expr_shifted);

        %% TEST CODE
%             likelihood = exp(log_likelihood);
%             prior = p_s;
%             likelihood2 = 1./(sqrt(2.*pi).*sigma_fun(s_grid_integrate).*rel_multiply_factor) .* exp(-((x_grid_rescale-s_grid_integrate).^2)./(2.*((sigma_fun(s_grid_integrate).*rel_multiply_factor).^2)));
%             expr = prior .* likelihood2;
%             
         %% Param-nonparam comparison
%         load("fittedparams_UJoint_GaussianLaplaceBothFixedZeroPMexpdnoiseGaussian_rescalefree.mat");
%         param_subjidx = 5;
%         theta_param = theta_fitted(param_subjidx,:);
%         sigma_fun_param_basis = heterotype_to_sigmafun("exp");
%         sigma_fun_param = @(s) sigma_fun_param_basis(s,theta_param(1),theta_param(2:3));
%         prior_param = (1-theta_param(9)).*normpdf(s_grid_integrate, 0, theta_param(6)) + theta_param(9)./(2.*theta_param(8)).*exp(-abs(s_grid_integrate)./theta_param(8));       
%         
%         figure
%         subplot(1,2,1)
%         hold on
%         plot(s_grid_integrate, log(sigma_fun_param(s_grid_integrate).*rel_multiply_factor),"b-")
%         plot(s_grid_integrate, log(sigma_fun(s_grid_integrate).*rel_multiply_factor),"r--")
%         ylabel("log sigma(s)")
%         xlabel("s")
%         subplot(1,2,2)
%         hold on
%         plot(s_grid_integrate, log(prior_param), "b-")
%         plot(s_grid_integrate, log_prior,"r--")
%         ylabel("log p(s)")
%         xlabel("s")
%         sgtitle("UV: Subject "+param_subjidx)
%         
%         
%         likelihood_param = 1./(sqrt(2.*pi).*sigma_fun_param(s_grid_integrate).*rel_multiply_factor) .* exp(-((x_grid_rescale-s_grid_integrate).^2)./(2.*((sigma_fun_param(s_grid_integrate).*rel_multiply_factor).^2)));
%         expr_param = prior_param .* likelihood_param;
%         s_hats_PM_param = sum(s_grid_integrate.*expr_param,1) ./ sum(expr_param,1); % the binwidth term gets factored out.           
%         p_r_given_x_param = normpdf(R_rel, s_hats_PM_param, sigma_motor);        
%         prob_r_given_s_nolapse_param = qtrapz(p_r_given_x_param .* p_x_given_s .* dx,2);
%         ll_idx_param = lapse./(max(lapse_range)-min(lapse_range)) + (1-lapse).* prob_r_given_s_nolapse_param; % Account for lapse (Unif(lapse_range) part of the mixture)

        
        %%

        s_hats_PM = qtrapz(s_grid_integrate.*expr,1) ./ qtrapz(expr,1);
        p_r_given_x = normpdf(R_rel, s_hats_PM, sigma_motor);
        prob_r_given_s_nolapse = qtrapz(p_r_given_x .* p_x_given_s .* dx,2); %./ qtrapz(p_x_given_s,2);
        
        switch lapse_type
            case "Uniform"
                ll_idx = lapse./(max(lapse_range)-min(lapse_range)) + (1-lapse).* prob_r_given_s_nolapse; % Account for lapse (Unif(lapse_range) part of the mixture)
            case "Gaussian"
                ll_idx = lapse .* normpdf(R_rel, 0, Gaussian_lapse_SD) + (1-lapse).* prob_r_given_s_nolapse; 
        end
        LL_idx_mat(rel_trial_idx) = ll_idx;
        
        if(isnan(-sum(log(ll_idx))))
            ll_idx
        end
    end

if(return_response_distr)
    output=LL_idx_mat;
else
    NLL = -sum(log(LL_idx_mat));
    if(isnan(NLL))
        clc
        %p_x_given_s
        theta
        %LL_idx
    end
    output=NLL;
end  
return

         
