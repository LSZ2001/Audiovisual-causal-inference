% This function returns the negative log likelihood of observing some data 
% R under corresponding S, under the PM model with sigma(s) = sigma_fun.

% plot_consider_lapse=false would remove the lapse mixture component in the
% response_distr_vector (LL_trials) returned. In contrast, NLL always includes the
% lapse.

% Inputs: ModelComponents: struct with fields {
% str         PriorType: "SingleGaussian", "TwoGaussiansOneFixedZero", "TwoGaussiansBothFixedZero"
% str         SensoryNoise: "constant", "linear", "quad", "squarequad", "exp", "cos"
% float/str   Rescale: "free" means to be fitted. Also, can be 1 (no rescale) or 4/3.
% boolean     IsFixedPriorMean: (if true, fit a prior mean parameter called mu; if not, mu is fixed at 0).
% vector      LapseRange: (-45 to 45 for now)}
%
% theta, S, R, R_round, PMIntegrationParams=[a,b,num_bins], plot_consider_lapse.

% theta = [sigma0, k[...], sigma_s, mu?, lapse, rescale?, sigma_center?, w?] is a set of free params values raised by BADS.
%    for hetero_type="constant", still need to provide k (any value makes no difference, but not fitted.
%       the output result for the k column would be 0.
%    mu present only when ~IsFixedPriorMean (hence need ~IsFixedPriorMean for TwoGaussianOneFixedZero)
%    rescale only present when ~FixedRescale && IsFitRescale. 
%    (sigma_center, w) only present for PriorType="TwoGaussiansOneFixedZero".

% If want to return posterior predictive samples given trials s (R input argument not used), return_predictive_samples=true.
% elseif want to return p(r|s) for a grid of r got given s, return_predictive_samples=false, return_response_distr=true.
% elseif want to return NLL=-log(sum(p(r|s)) of all trials, return_predictive_samples=false, return_response_distr=false;


function [output] = nllfun_bc_ubresc_semiparaminsp(theta, data, ModelComponents, return_predictive_samples, return_response_distr, plot_consider_lapse)
    if(nargin==3) % default is to return NLL, not posterior preditive samples vectorized response distribition.
        return_predictive_samples = false; return_response_distr=false; plot_consider_lapse=true;
    elseif(nargin==4) % if return posterior predictive samples and not vectorized response distribution, default is not include lapse mixture component.
        return_response_distr=false; plot_consider_lapse=true;
    elseif (nargin==5) % if return vectorized response distribution (need S to be one s-value and rel-level repmat, and R to be r_grid), default is not include lapse mixture component.
        plot_consider_lapse=true;
    end
    causal_inf_strategy = ModelComponents.CausalInfStrategy;
    s_grid_integrate = ModelComponents.SGridIntegrate;
    ds = s_grid_integrate(2) - s_grid_integrate(1);
    num_bins = length(s_grid_integrate);
    num_rels_vis = 3; % 3 reliability levels for V.
    
    %% Parse dataset.
    S_V = data(:,[1,3]); % include reliability level info for each trial.
    S_A = data(:,2);
    R = data(:,4); % R is just C in this task type.
    [~,ns] = size(S_V);
    [~,nr] = size(R);
    if(ns>2)
        S_V=S_V';
        S_A=S_A';
    end
    if(nr>1)
        R=R';
    end
    num_trials = length(R);
    reliabilities = S_V(:,1);
    S_V = S_V(:,2);
    
    %% Parameter extraction
    % [scale_rel_med, scale_rel_low, lapse, sig_motor, rescale_aud, p_same,
    % Bimod_rescale_vis, Bimod_rescale_aud, UB_shared_prior_power];
    rel_multiply_factors_uniq =[1, theta(1:2), 1];
    lapse = theta(3); 
    rescale_vis = 1; rescale_aud = theta(5);
    p_same = theta(6);
    B_rescale_vis = theta(7); B_rescale_aud = theta(8);
    UB_shared_prior_power = theta(end);
    if(~plot_consider_lapse) % Do not consider the lapse component.
        lapse=0;
    end
    sigma_fun_vis_raw = ModelComponents.SigmaFuns{1};
    sigma_fun_aud_raw = ModelComponents.SigmaFuns{4};
    sigma_fun_vis = @(s) B_rescale_vis .* sigma_fun_vis_raw(s);
    sigma_fun_aud = @(s) B_rescale_aud .* sigma_fun_aud_raw(s);
    
    %% All below is NLL evaluation for model fitting
    LL_idx = zeros(num_trials,1);
    for i=1:num_rels_vis
        rel_multiply_factor = rel_multiply_factors_uniq(i);
        rel_trial_idx = find(reliabilities==i);
        if(isempty(rel_trial_idx)) % if no trials in the data has this reliability, skip this iteration
            continue;
        end
        S_V_rel = S_V(rel_trial_idx); S_A_rel = S_A(rel_trial_idx); R_rel = R(rel_trial_idx);

        % Bounds of V interation
        s_V_min = min(S_V_rel); s_V_max = max(S_V_rel);
        x_V_min = s_V_min(1) - 4*sigma_fun_vis(s_V_min(1))*rel_multiply_factor;
        x_V_max = s_V_max(1) + 4*sigma_fun_vis(s_V_max(1))*rel_multiply_factor;
        x_V_min = floor(x_V_min*2)/2; x_V_max = ceil(x_V_max*2)/2;

        x_V_grid = [x_V_min:0.5:-15.5, -15:0.25:-5.25, -5:0.1:5, 5.25:0.25:15 ,15.5:0.5:x_V_max];
        Nx_V = length(x_V_grid);
        dx_V_raw = x_V_grid(2:end) - x_V_grid(1:(end-1));
        dx_V_raw = [dx_V_raw, dx_V_raw(end)];
        dx_V = zeros(1,1,length(dx_V_raw));
        dx_V(1,1,:) = dx_V_raw;
        clear dx_V_raw

        % Bounds of A interation
        s_A_min = min(S_A_rel); s_A_max = max(S_A_rel);
        x_A_min = s_A_min(1) - 4*sigma_fun_aud(s_A_min(1)) - 5;
        x_A_max = s_A_max(1) + 4*sigma_fun_aud(s_A_max(1)) + 5;
        x_A_min = floor(x_A_min*2)/2; x_A_max = ceil(x_A_max*2)/2;

        x_A_grid = [x_A_min:0.5:-15.5, -15:0.25:-5.25, -5:0.1:5, 5.25:0.25:15 ,15.5:0.5:x_A_max];
        Nx_A = length(x_A_grid);
        dx_A = x_A_grid(2:end) - x_A_grid(1:(end-1));
        dx_A = [dx_A, dx_A(end)];

        % Change sigma_fun at the bounds to be larger than dx grid.
        grid_size_fun = @(s) 0.5.*(abs(s)>15) + 0.25.*((abs(s)<=15).*abs(s)>5) + 0.1.*(abs(s)<=5);
        sigma_fun_vis = @(s) max(sigma_fun_vis(s), grid_size_fun(s));
        sigma_fun_aud = @(s) max(sigma_fun_aud(s), grid_size_fun(s));


        % 3D tensors: (R_rel, xA, xV). 
        p_x_given_s_V_raw = normpdf(x_V_grid, S_V_rel, sigma_fun_vis(S_V_rel).*rel_multiply_factor);
        p_x_given_s_A_raw = normpdf(x_A_grid, S_A_rel, sigma_fun_aud(S_A_rel));

        p_x_given_s_V = reshape(p_x_given_s_V_raw, [length(S_V_rel), 1, Nx_V]);
        p_x_given_s_A = reshape(p_x_given_s_A_raw, [length(S_A_rel), Nx_A, 1]);
        clear p_x_given_s_V_raw p_x_given_s_A_raw

        % Get normalized prior density.
        p_s_unnorm_base = ModelComponents.PriorUnnormalized;
        p_s_unnorm = p_s_unnorm_base .^ UB_shared_prior_power;
        p_s = p_s_unnorm / qtrapz(p_s_unnorm*(s_grid_integrate(2)-s_grid_integrate(1)));
        log_prior = log(p_s);
        
        x_V_grid_rescale = x_V_grid .* rescale_vis;   
        x_A_grid_rescale = x_A_grid .* rescale_aud;  

        % These are non-Gaussian likelihood functions of s_grid_integrate.
        log_p_xV_given_s_integrate_raw = -log(sigma_fun_vis(s_grid_integrate).*rel_multiply_factor) - 0.5*log(2*pi) - 0.5.*(((x_V_grid_rescale-s_grid_integrate)./(sigma_fun_vis(s_grid_integrate).*rel_multiply_factor)).^2);
        log_p_xA_given_s_integrate_raw = -log(sigma_fun_aud(s_grid_integrate)) - 0.5*log(2*pi) - 0.5.*(((x_A_grid_rescale-s_grid_integrate)./(sigma_fun_aud(s_grid_integrate))).^2);

        log_p_xV_given_s_integrate = reshape(log_p_xV_given_s_integrate_raw , [num_bins, 1, Nx_V]);
        log_p_xA_given_s_integrate = reshape(log_p_xA_given_s_integrate_raw , [num_bins, Nx_A, 1]);
        clear log_p_xV_given_s_integrate_raw log_p_xA_given_s_integrate_raw

        posteriors_C = zeros(2, Nx_A, Nx_V);
        switch ModelComponents.PriorAssumption
            case "independent"               
                % Pr[C|x_V, x_A] for C={1,2}.
                if((p_same~=1) && (p_same~=0))
                    log_expr_C1 =  log_prior + log_p_xV_given_s_integrate + log_p_xA_given_s_integrate;
                    log_offset_C1 = max(log_expr_C1,[],1);
                    log_expr_C1 = log_expr_C1 - log_offset_C1;
                    expr_C1 = exp(log_expr_C1);
                    overall_log_likelihood_C1 = squeeze(log_offset_C1 + log(sum(expr_C1.*ds,1))); % This is now [x_A, x_V].           
                    %protoposteriors_C(1,:,:) = overall_likelihood_C1 .* p_same; %./ sum(overall_likelihood_C1(:));
                    clear log_expr_C1 log_offset_C1 expr_C1

                    log_expr_C2_V = log_prior + log_p_xV_given_s_integrate;
                    log_offset_C2_V = max(log_expr_C2_V,[],1);
                    log_expr_C2_V = log_expr_C2_V - log_offset_C2_V;
                    log_expr_C2_A = log_prior + log_p_xA_given_s_integrate;
                    log_offset_C2_A = max(log_expr_C2_A,[],1);
                    log_expr_C2_A = log_expr_C2_A - log_offset_C2_A;

                    expr_C2_V = exp(log_expr_C2_V);
                    expr_C2_A = exp(log_expr_C2_A);


                    expr_C2_V_term = log_offset_C2_V + log(sum(expr_C2_V.*ds,1));
                    expr_C2_A_term = log_offset_C2_A + log(sum(expr_C2_A.*ds,1));
                    overall_log_likelihood_C2 = squeeze(expr_C2_V_term + expr_C2_A_term);
                    %protoposteriors_C(2,:,:) = overall_likelihood_C2 .* (1-p_same);% ./ sum(overall_likelihood_C2(:));

                    clear  log_expr_C2_V  log_expr_C2_A log_offset_C2_V log_offset_C2_A expr_C2_V expr_C2_A expr_C2_V_term expr_C2_A_term

                    d = log(p_same/(1-p_same)) + overall_log_likelihood_C1 - overall_log_likelihood_C2;
                    posteriors_C(1,:,:) = 1 ./ (1 + exp(-d));
                    posteriors_C(2,:,:) = 1-posteriors_C(1,:,:);
                    clear d overall_log_likelihood_C1 overall_log_likelihood_C2
                else
                    posteriors_C(1,:,:)=p_same;
                    posteriors_C(2,:,:)=1-p_same;
                end

                % dimensions [trials, s_A, s_V]
                if(causal_inf_strategy=="ProbMatching")
                    expr = posteriors_C(R_rel,:,:) .* p_x_given_s_A .* p_x_given_s_V;
                else
                    expr = (posteriors_C(R_rel,:,:)>1/2) .* p_x_given_s_A .* p_x_given_s_V;
                end
                clear posteriors_C p_x_given_s_A p_x_given_s_V

                dx_A([1, end]) = dx_A([1, end])./2;
                dx_V([1, end]) = dx_V([1, end])./2;
                dx_A_changeidx = find(((x_A_grid==-15) + (x_A_grid==-5) + (x_A_grid==5) + (x_A_grid==15))>0);
                dx_A(dx_A_changeidx) = (dx_A(dx_A_changeidx+1) + dx_A(dx_A_changeidx-1)) ./2;
                dx_V_changeidx = find(((x_V_grid==-15) + (x_V_grid==-5) + (x_V_grid==5) + (x_V_grid==15))>0);
                dx_V(dx_V_changeidx) = (dx_V(dx_V_changeidx+1) + dx_V(dx_V_changeidx-1)) ./2;
                clear dx_V_changeidx dx_A_changeidx

                prob_r_given_s_nolapse = sum(sum(expr .* dx_A .* dx_V,2),3);
                clear expr dx_A dx_V

            case "empirical"
        end
        ll_idx = lapse./2 + (1-lapse).* prob_r_given_s_nolapse; % Account for lapse (Unif(lapse_range) part of the mixture)
        LL_idx(rel_trial_idx) = ll_idx;
        clear ll_idx
    end


    if(return_response_distr)
        output=LL_idx;
    else
        NLL = -sum(log(LL_idx));
        if(isnan(NLL))
            clc
            %p_x_given_s
            theta_raw
            %LL_idx
        end
        output=NLL;
    end  
    return

end