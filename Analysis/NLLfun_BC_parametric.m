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


function [output] = nllfun_bc_parametric(ModelComponents_V, ModelComponents_A, theta_raw, R, S_V, S_A, return_predictive_samples, return_response_distr, plot_consider_lapse, causal_inf_strategy)
    if(nargin<7) % default is to return NLL, not posterior preditive samples vectorized response distribition.
         return_predictive_samples = false; return_response_distr=false; causal_inf_strategy="ProbMatching";
    elseif(nargin<8) % if return posterior predictive samples and not vectorized response distribution, default is not include lapse mixture component.
        return_response_distr=false; plot_consider_lapse=true; causal_inf_strategy="ProbMatching";
    elseif (nargin<9) % if return vectorized response distribution (need S to be one s-value and rel-level repmat, and R to be r_grid), default is not include lapse mixture component.
        plot_consider_lapse=true; causal_inf_strategy="ProbMatching";
    elseif(nargin<10)
        causal_inf_strategy="ProbMatching";
    end
    PMIntegrationParams = ModelComponents_V.PMIntegrationParams;
    % Pop the parameter p_same from theta, so that V/A param distangling code from
    % unimodal trials can be applied directly on theta.
    UB_rescale_aud = theta_raw(end); % Rescale applied to bimodal cases, on sigma(s). 
    theta_raw = theta_raw(1:(end-1));
    UB_rescale_vis = theta_raw(end); % Rescale applied to bimodal cases, on sigma(s). 
    theta_raw = theta_raw(1:(end-1));
    
    p_same = theta_raw(end);
    theta = theta_raw(1:(end-1));
    % Get separate params for V and A.
    [LB_V, UB_V, PLB_V, PUB_V] = sigmafun_badsbounds_comprehensive(ModelComponents_V);
    num_V_params = length(LB_V);
    [LB_A, UB_A, PLB_A, PUB_A] = sigmafun_badsbounds_comprehensive(ModelComponents_A);
    [LB, ~, ~, ~, A_param_keep_idx] = merge_ujoint_badsbounds(LB_V,UB_V,PLB_V,PUB_V,LB_A,UB_A,PLB_A,PUB_A,ModelComponents_A);
    num_A_uniq_params = length(LB) - length(LB_A) - 1;
    theta_vis = theta(1:num_V_params);
    theta_aud = complete_thetaua_for_ujointfits(theta, A_param_keep_idx, ModelComponents_V.Rescale=="free");
      
    prior_type = ModelComponents_V.PriorType;
    hetero_type = ModelComponents_V.SensoryNoise;
    lapse_range = ModelComponents_V.LapseRange;
    is_fixed_prior_mean = ModelComponents_V.IsFixedPriorMean;
    num_rels = 3; % 3 reliability levels for V.
    
    rescale_vis = ModelComponents_V.Rescale; % visual rescale
    rescale_aud = ModelComponents_A.Rescale; % auditory rescale.
    
    motor_noise_type =  ModelComponents_V.MotorNoise;
    a=PMIntegrationParams(1); b=PMIntegrationParams(2); num_bins = PMIntegrationParams(3);
    sigma_fun_raw = heterotype_to_sigmafun(hetero_type);    
    sigma_fun_vis = @(s, sigma0, k) UB_rescale_vis .* sigma_fun_raw(s, sigma0, k);
    sigma_fun_aud = @(s, sigma0, k) UB_rescale_aud .* sigma_fun_raw(s, sigma0, k);
    
    % Ensure S and R inputs are column vectors.
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
    
    if(hetero_type=="exp") % there are 2 slope params, k(1), k(2).
        k_idx = [2,3];
    elseif(hetero_type=="constant")
        k_idx = 2;
    else
        k_idx = 2;
    end
    
    % Extract V parameters.
    sigma0_vis=theta_vis(1); k_vis=theta_vis(k_idx); sigma_s = theta_vis(k_idx(end)+(num_rels-1)+1);
    if(rescale_vis~="free") % not fit rescale param
        rescale_vis = str2num(rescale_vis);
        if(is_fixed_prior_mean) % PM model, mu is fixed at 0
            lapse=theta_vis(k_idx(end)+(num_rels-1)+2); mu=0;
        else % PMmufree model, mu is free param to be fitted
            mu=theta_vis(k_idx(end)+(num_rels-1)+2); lapse=theta_vis(k_idx(end)+(num_rels-1)+3);
        end
    else % fit rescale param; so override the provided, fixed rescale. Instead, use values in theta.
        if(is_fixed_prior_mean) % mu is fixed at 0
            lapse=theta_vis(k_idx(end)+(num_rels-1)+2); mu=0; rescale_vis=theta_vis(k_idx(end)+(num_rels-1)+3);
        else % mu is free param to be fitted
            mu=theta_vis(k_idx(end)+(num_rels-1)+2); lapse=theta_vis(k_idx(end)+(num_rels-1)+3); rescale_vis=theta_vis(k_idx(end)+(num_rels-1)+4);
        end
    end
    
    % Extract A unique parameters.
    sigma0_aud=theta_aud(1); k_aud=theta_aud(k_idx);
    if(rescale_aud~="free") % not fit rescale param
        rescale_aud = str2num(rescale_aud);
    else % fit rescale param; so override the provided, fixed rescale. Instead, use values in theta.
        if(is_fixed_prior_mean) % mu is fixed at 0
            rescale_aud=theta_aud(k_idx(end)+3);
        else % mu is free param to be fitted
            rescale_aud=theta_aud(k_idx(end)+4);
        end
    end
   
    % Get reliabilities and corresponding scales on sigma0_vis.
    reliabilities = S_V(:,1);
    S_V = S_V(:,2);
    rel_multiply_factors_uniq = [1,theta((k_idx(end)+1):(k_idx(end)+num_rels-1))]; % Reliability 1theta is most reliable, corresponding to sigma(s). 
    rel_multiply_factors = rel_multiply_factors_uniq(reliabilities)'; % sigma(s) multiply factor for each trial specifically

    
    if(hetero_type=="constant")
        k_vis=0; k_aud=0;
    end
    R_raw=R;

    if(motor_noise_type=="Gaussian")
        sigma_motor = theta_vis(end); % Last entry of theta_vis (shared with theta_aud) is the motor noise magnitude. 
        if(~plot_consider_lapse) % Do not consider the lapse component.
            lapse=0;
        end
      
        if(return_predictive_samples) % If return predictive samples, cannot return response distribution despite input args.
            return_response_distr=false;
            num_samps = 100; % Number of samples for this trial. 
            if(num_rels>1)
                rel_multiples_vis = rel_multiply_factors_uniq(reliabilities)';
            else
                rel_multiples_vis = rel_multiply_factors_uniq(reliabilities);
            end

            switch prior_type
                case "SingleGaussian"
                    % Use midpoint rule to get posterior mean
                    [posterior_C1,posterior_C2] = BC_getposteriors(sigma_fun, x_obs_vals, sigma0,k,sigma_s,mu, p_same, a, b, num_bins, "midpoint",rel_multiples_vis, prior_type);
                case "TwoGaussiansOneFixedZero"
                    sigma_center=theta(end-2); w=theta(end-1); % w is weight of second Gaussian.
                    [posterior_C1,posterior_C2] = BC_getposteriors(sigma_fun, x_obs_vals, sigma0,k,sigma_s,mu, p_same, a, b, num_bins, "midpoint",rel_multiples_vis, prior_type, sigma_center,w);
                    
                case "TwoGaussiansBothFixedZero"
                    delta_sigma=theta(end-2); w=theta(end-1);
                    sigma_s_larger = sigma_s + delta_sigma;
                    [posterior_C1,posterior_C2] = BC_getposteriors(sigma_fun, x_obs_vals, sigma0,k,sigma_s,mu, p_same, a, b, num_bins, "midpoint",rel_multiples_vis, prior_type, sigma_s_larger,w);
                    
                case "GaussianLaplaceBothFixedZero"
                    laplace_scale = theta(end-2); w=theta(end-1); % w is weight of Laplace.
                    [posterior_C1,posterior_C2] = BC_getposteriors(sigma_fun, x_obs_vals, sigma0,k,sigma_s,mu, p_same, a, b, num_bins, "midpoint",rel_multiples_vis, prior_type, laplace_scale,w);
            end

            % No motor noise for BC task; report category choice according
            % to causal inference strategy (that is, how to get from
            % posterior_C1 and posterior_C2 to a category choice)
            switch ModelComponents_V.CausalInfStrategy
                case "ModelAveraging"
                    output = (posterior_C1 < posterior_C2) + 1;
                case "ModelSelection"
                    output = (posterior_C1 < posterior_C2) + 1;
                case "ProbMatching"
                    output = (rand([length(S_V), num_samps]) > posterior_C1) + 1;
            end
                    
            
            % Add in lapse trials
            lapse_idx = rand([length(S_V), num_samps]) < lapse;
            lapse_idx_single = find(lapse_idx~=0);
            % Random responses (uniform in -45:1:45)
            output(lapse_idx_single) = randi(2,1,length(lapse_idx_single));
    
            return
        end
 
        
        LL_idx = zeros(num_trials,1);
        for i=1:num_rels
            rel_multiply_factor = rel_multiply_factors_uniq(i);
            rel_trial_idx = find(reliabilities==i);
            if(isempty(rel_trial_idx)) % if no trials in the data has this reliability, skip this iteration
                continue;
            end
            S_V_rel = S_V(rel_trial_idx); S_A_rel = S_A(rel_trial_idx); R_rel = R_raw(rel_trial_idx);
           
            % Bounds of V interation
            s_V_min = min(S_V_rel); s_V_max = max(S_V_rel);
            x_V_min = s_V_min(1) - 4*sigma_fun_vis(s_V_min(1), sigma0_vis,k_vis)*rel_multiply_factor;
            x_V_max = s_V_max(1) + 4*sigma_fun_vis(s_V_max(1),sigma0_vis,k_vis)*rel_multiply_factor;
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
            x_A_min = s_A_min(1) - 4*sigma_fun_aud(s_A_min(1), sigma0_aud,k_aud) - 5;
            x_A_max = s_A_max(1) + 4*sigma_fun_aud(s_A_max(1),sigma0_aud,k_aud) + 5;
            x_A_min = floor(x_A_min*2)/2; x_A_max = ceil(x_A_max*2)/2;
            
            x_A_grid = [x_A_min:0.5:-15.5, -15:0.25:-5.25, -5:0.1:5, 5.25:0.25:15 ,15.5:0.5:x_A_max];
            Nx_A = length(x_A_grid);
            dx_A = x_A_grid(2:end) - x_A_grid(1:(end-1));
            dx_A = [dx_A, dx_A(end)];
            
            % Change sigma_fun at the bounds to be larger than dx grid.
            grid_size_fun = @(s) 0.5.*(abs(s)>15) + 0.25.*((abs(s)<=15).*abs(s)>5) + 0.1.*(abs(s)<=5);
            sigma_fun_vis = @(s, sigma0, k) max(sigma_fun_vis(s, sigma0, k), grid_size_fun(s));
            sigma_fun_aud = @(s, sigma0, k) max(sigma_fun_aud(s, sigma0, k), grid_size_fun(s));

            
            % 3D tensors: (R_rel, xA, xV). 
            p_x_given_s_V_raw = normpdf(x_V_grid, S_V_rel, sigma_fun_vis(S_V_rel,sigma0_vis,k_vis).*rel_multiply_factor);
            p_x_given_s_A_raw = normpdf(x_A_grid, S_A_rel, sigma_fun_aud(S_A_rel,sigma0_aud,k_aud));

            p_x_given_s_V = reshape(p_x_given_s_V_raw, [length(S_V_rel), 1, Nx_V]);
            p_x_given_s_A = reshape(p_x_given_s_A_raw, [length(S_A_rel), Nx_A, 1]);
            clear p_x_given_s_V_raw p_x_given_s_A_raw
            
            % Trapezoid Rule to get s_hat_PMs along x_grid.
            binedges = linspace(a,b,num_bins);
            binwidth =  binedges(2)-binedges(1);
            s_grid_integrate = (binedges(1:end-1) + binwidth/2)';
            ds = s_grid_integrate(2) - s_grid_integrate(1);
            switch prior_type
                case "SingleGaussian"
                    %prior = normpdf(s_grid_integrate, mu, sigma_s);
                    log_prior = -log(sigma_s) - 0.5*log(2*pi) - 0.5.*(((s_grid_integrate-mu)/sigma_s).^2);
                    %log_offset_prior = max(log_prior);
                    %log_prior_new = log_prior - log_offset_prior;
                    %prior = exp(log_prior_new) .* exp(log_offset_prior);
                 case "TwoGaussiansOneFixedZero"
                    % The free moving gaussian is N(mu, sigma_s), both defined
                    % above. (Hence prior_type="TwoGaussiansOneFixedZero" must be paired with IsFixedPriorMean=false)
                    % sigma_center is the center of the gaussian fixed at mu_fixed=0.
                    % w is the weight of the fix gaussian as a mixture component for the prior.
                    sigma_center=theta_vis(end-2); w=theta_vis(end-1);
                    % Since 1/sqrt(2pi*sigma_s) and 1/sqrt(2pi*sigma_center) are
                    % constants wrt to s, hence absorbed into proportionality sign.
                    prior = (1-w).*normpdf(s_grid_integrate, mu, sigma_s) + w.* normpdf(s_grid_integrate, 0, sigma_center);
                case "TwoGaussiansBothFixedZero"
                    delta_sigma=theta_vis(end-2); w=theta_vis(end-1);
                    sigma_s_larger = sigma_s + delta_sigma;
                    % Since 1/sqrt(2pi*sigma_s) and 1/sqrt(2pi*sigma_center) are
                    % constants wrt to s, hence absorbed into proportionality sign.
                    log_prior = log((1-w).*normpdf(s_grid_integrate, 0, sigma_s) + w.* normpdf(s_grid_integrate, 0, sigma_s_larger));
                case "GaussianLaplaceBothFixedZero"
                    % simga_s is width of Gaussian component, with mixture
                    % weight (1-w).
                    laplace_scale = theta_vis(end-2); w=theta_vis(end-1);
                    log_prior = log((1-w).*normpdf(s_grid_integrate, 0, sigma_s) + w./(2.*laplace_scale).*exp(-abs(s_grid_integrate-0)./laplace_scale));       
            end
            %prior(prior==0) = 10^(-20);
            
            
            x_V_grid_rescale = x_V_grid .* rescale_vis;   
            x_A_grid_rescale = x_A_grid .* rescale_aud;  
            
            % These are non-Gaussian likelihood functions of s_grid_integrate.
            log_p_xV_given_s_integrate_raw = -log(sigma_fun_vis(s_grid_integrate,sigma0_vis,k_vis).*rel_multiply_factor) - 0.5*log(2*pi) - 0.5.*(((x_V_grid_rescale-s_grid_integrate)./(sigma_fun_vis(s_grid_integrate,sigma0_vis,k_vis).*rel_multiply_factor)).^2);
            log_p_xA_given_s_integrate_raw = -log(sigma_fun_aud(s_grid_integrate,sigma0_aud,k_aud)) - 0.5*log(2*pi) - 0.5.*(((x_A_grid_rescale-s_grid_integrate)./(sigma_fun_aud(s_grid_integrate,sigma0_aud,k_aud))).^2);

            log_p_xV_given_s_integrate = reshape(log_p_xV_given_s_integrate_raw , [num_bins-1, 1, Nx_V]);
            log_p_xA_given_s_integrate = reshape(log_p_xA_given_s_integrate_raw , [num_bins-1, Nx_A, 1]);
            clear log_p_xV_given_s_integrate_raw log_p_xA_given_s_integrate_raw
            
            posteriors_C = zeros(2, Nx_A, Nx_V);
            switch ModelComponents_V.PriorAssumption
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
         
end