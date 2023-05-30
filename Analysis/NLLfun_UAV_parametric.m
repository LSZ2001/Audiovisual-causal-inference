function [output] = NLLfun_UAV_parametric(ModelComponents, theta, R, S, return_predictive_samples, return_response_distr, plot_consider_lapse, lapse_type, Gaussian_lapse_SD)
    
    if(nargin==4) % default is to return NLL, not posterior preditive samples vectorized response distribition.
         return_predictive_samples = false; return_response_distr=false; lapse_type = "Uniform";
    elseif(nargin==5) % if return posterior predictive samples and not vectorized response distribution, default is not include lapse mixture component.
        return_response_distr=false; plot_consider_lapse=true; lapse_type = "Uniform";
    elseif (nargin==6) % if return vectorized response distribution (need S to be one s-value and rel-level repmat, and R to be r_grid), default is not include lapse mixture component.
        plot_consider_lapse=true; lapse_type = "Uniform";
    elseif(nargin==7)
        lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==8)
        Gaussian_lapse_SD = NaN;
    end
    
    if(lapse_type=="Gaussian")
%         sigma_lapse = 10; % For plotting models that do not have
%         sigma_lapse, and want to visualize some other fixed Gaussian
%         lapse distribution with some predefined SD value.
        sigma_lapse = theta(end);
        theta = theta(1:(end-1));
    end
    
    if(ModelComponents.IsMLModel)
        print("ALERT: Intput is ML model. Do not use this function, which is for models involving posteriors only!!!")
        output=NaN;
        return
    end    
    PMIntegrationParams = ModelComponents.PMIntegrationParams;
    prior_type = ModelComponents.PriorType;
    hetero_type = ModelComponents.SensoryNoise;
    lapse_range = ModelComponents.LapseRange;
    is_fixed_prior_mean = ModelComponents.IsFixedPriorMean;
    if(~isfield(ModelComponents, "NumReliabilityLevels")) % For auditory info, no reliability levels by default.
        num_rels = 1;
    else
        num_rels = ModelComponents.NumReliabilityLevels;
    end
    
    rescale = ModelComponents.Rescale;
    motor_noise_type =  ModelComponents.MotorNoise;
    a=PMIntegrationParams(1); b=PMIntegrationParams(2); num_bins = PMIntegrationParams(3);
    sigma_fun = heterotype_to_sigmafun(hetero_type);
    
    % Ensure S and R inputs are column vectors.
    [~,ns] = size(S);
    [~,nr] = size(R);
    if(ns>2)
        S=S';
    end
    if(nr>1)
        R=R';
    end
    num_trials = length(R);

    % Ensure S and R are col vectors here.
    
    if(hetero_type=="exp") % there are 2 slope params, k(1), k(2).
        k_idx = [2,3];
    elseif(hetero_type=="constant")
        k_idx = 2;
    else
        k_idx = 2;
    end
    
    if(rescale~="free") % not fit rescale param
        rescale = str2num(rescale);
        if(is_fixed_prior_mean) % PM model, mu is fixed at 0
            sigma0=theta(1); k=theta(k_idx); sigma_s = theta(k_idx(end)+(num_rels-1)+1); lapse=theta(k_idx(end)+(num_rels-1)+2); mu=0;
        else % PMmufree model, mu is free param to be fitted
            sigma0=theta(1); k=theta(k_idx); sigma_s = theta(k_idx(end)+(num_rels-1)+1); mu=theta(k_idx(end)+(num_rels-1)+2); lapse=theta(k_idx(end)+(num_rels-1)+3);
        end
    else % fit rescale param; so override the provided, fixed rescale. Instead, use values in theta.
        if(is_fixed_prior_mean) % mu is fixed at 0
            sigma0=theta(1); k=theta(k_idx); sigma_s = theta(k_idx(end)+(num_rels-1)+1); lapse=theta(k_idx(end)+(num_rels-1)+2); mu=0; rescale=theta(k_idx(end)+(num_rels-1)+3);
        else % mu is free param to be fitted
            sigma0=theta(1); k=theta(k_idx); sigma_s = theta(k_idx(end)+(num_rels-1)+1); mu=theta(k_idx(end)+(num_rels-1)+2); lapse=theta(k_idx(end)+(num_rels-1)+3); rescale=theta(k_idx(end)+(num_rels-1)+4);
        end
    end
   
    reliabilities = repmat(1,num_trials,1);
    if(num_rels>1) % If visual data, extract multiply factors on sigma of higher reliabilities.
        reliabilities = S(:,2);
        S = S(:,1);
        rel_multiply_factors_uniq = [1,theta((k_idx(end)+1):(k_idx(end)+num_rels-1))]; % Reliability 1theta is most reliable, corresponding to sigma(s). 
        rel_multiply_factors = rel_multiply_factors_uniq(reliabilities)'; % sigma(s) multiply factor for each trial specifically
    else % If auditory data, sigma(s) no need scaling.
        rel_multiply_factors_uniq = 1;
        rel_multiply_factors = reliabilities;
    end
    if(return_predictive_samples && num_rels==1) % For visualizing UA response distributions using post pred samples
        reliabilities = S(:,2);
        S = S(:,1);
    end
    
    if(hetero_type=="constant")
        k=0;
    end
    R_raw=R;

    if(motor_noise_type=="Gaussian")
        sigma_motor = theta(end); % Last entry of theta is the motor noise magnitude. 
        if(~plot_consider_lapse) % Do not consider the lapse component.
            lapse=0;
        end
        
        
        
        
        if(return_predictive_samples) % If return predictive samples, cannot return response distribution despite input args.
            return_response_distr=false;
            num_samps = 100; % Number of samples for this trial. 
            if(num_rels>1)
                rel_multiples = rel_multiply_factors_uniq(reliabilities)';
            else % Audition -- no reliability levels.
                rel_multiples = 1;
            end
            
            sigma = sigma_fun(S, sigma0, k) .* rel_multiples;
            
            
            x_obs_vals = S + sigma.*randn([length(S),num_samps]);
            
            % STILL NEEDS CHANGING! For the TwoGaussian and GaussianLaplace
            % mixture priors, cannot draw x|s ~ SingleGaussian!
            switch prior_type
                case "SingleGaussian"
                    % Use midpoint rule to get posterior mean
                    [s_hats_PM,~] = midpoint_trepazoid_getpostmean(sigma_fun, x_obs_vals .* rescale, sigma0,k,sigma_s,mu, a, b, num_bins, "midpoint",rel_multiples, prior_type);
                case "TwoGaussiansOneFixedZero"
                    sigma_center=theta(end-2); w=theta(end-1); % w is weight of second Gaussian.
                    [s_hats_PM,~] = midpoint_trepazoid_getpostmean(sigma_fun, x_obs_vals .* rescale, sigma0,k,sigma_s,mu, a, b, num_bins, "midpoint",rel_multiples, prior_type, sigma_center,w);
                    
                case "TwoGaussiansBothFixedZero"
                    delta_sigma=theta(end-2); w=theta(end-1);
                    sigma_s_larger = sigma_s + delta_sigma;
                    [s_hats_PM,~] = midpoint_trepazoid_getpostmean(sigma_fun, x_obs_vals .* rescale, sigma0,k,sigma_s,mu, a, b, num_bins, "midpoint",rel_multiples, prior_type, sigma_s_larger,w);
                    
                case "GaussianLaplaceBothFixedZero"
                    laplace_scale = theta(end-2); w=theta(end-1); % w is weight of Laplace.
                    [s_hats_PM,~] = midpoint_trepazoid_getpostmean(sigma_fun, x_obs_vals .* rescale, sigma0,k,sigma_s,mu, a, b, num_bins, "midpoint",rel_multiples, prior_type, laplace_scale,w);
            end
            
            % Add motor noise
            output = s_hats_PM + sigma_motor .* randn([length(S),num_samps]);
            
            % Add in lapse trials
            lapse_idx = rand([length(S), num_samps]) < lapse;
            switch lapse_type
                case "Uniform"
                    % Random responses (uniform in -45:1:45)
                    lapse_val = min(lapse_range)+ (max(lapse_range)-min(lapse_range))*rand(sum(lapse_idx(:)),1);
                case "Gaussian"
                    %lapse_val = Gaussian_lapse_SD .* randn(sum(lapse_idx(:)),1);
                    low_trunc = repmat(min(lapse_range)/sigma_lapse, sum(lapse_idx(:)), 1);
                    high_trunc = repmat(max(lapse_range)/sigma_lapse, sum(lapse_idx(:)), 1);
                    lapse_val = sigma_lapse .* trandn(low_trunc, high_trunc);
                    %lapse_val = sigma_lapse .* randn(sum(lapse_idx(:)),1);
            end
            output(find(lapse_idx~=0)) = lapse_val;
    
            return
        end
        
        
          
        
        LL_idx = zeros(num_trials,1);
        for i=1:num_rels
            rel_multiply_factor = rel_multiply_factors_uniq(i);
            rel_trial_idx = find(reliabilities==i);
            if(isempty(rel_trial_idx)) % if no trials in the data has this reliability, skip this iteration
                continue;
            end
            S_rel = S(rel_trial_idx); R_rel = R_raw(rel_trial_idx);
            
            % Bounds of interation for p(x|s)
            s_min = min(S_rel); s_max = max(S_rel);
            x_min = s_min(1) - 5*sigma_fun(s_min(1), sigma0,k)*rel_multiply_factor;
            x_max = s_max(1) + 5*sigma_fun(s_max(1),sigma0,k)*rel_multiply_factor;
            Nx = 2^10;
            x_grid = linspace(x_min, x_max, Nx); % Need to transpose due to matrix operation x-s.
            dx = x_grid(2) - x_grid(1);
            if(dx > 0.1) % Ensure that dx is much smaller than the LB (0.25) of sigma_motor.
                dx=0.1;
                x_grid = (x_min:dx:x_max);
                Nx = length(x_grid);
            end
            
            p_x_given_s = normpdf(x_grid, S_rel, sigma_fun(S_rel,sigma0,k).*rel_multiply_factor);
            %p_x_given_s = p_x_given_s ./ sum(p_x_given_s,2);
            
            % Midpoint Rule to get s_hat_PMs along x_grid.
            binedges = linspace(a,b,num_bins+1);
            binwidth =  binedges(2)-binedges(1);
            s_grid_integrate = (binedges(1:end-1) + binwidth/2)';
            switch prior_type
                case "SingleGaussian"
                    log_prior = -log(sigma_s) - 0.5*log(2*pi) - 0.5.*(((s_grid_integrate-mu)/sigma_s).^2);
                case "TwoGaussiansOneFixedZero"
                    % The free moving gaussian is N(mu, sigma_s), both defined
                    % above. (Hence prior_type="TwoGaussiansOneFixedZero" must be paired with IsFixedPriorMean=false)
                    % sigma_center is the center of the gaussian fixed at mu_fixed=0.
                    % w is the weight of the fix gaussian as a mixture component for the prior.
                    sigma_center=theta(end-2); w=theta(end-1);
                    % Since 1/sqrt(2pi*sigma_s) and 1/sqrt(2pi*sigma_center) are
                    % constants wrt to s, hence absorbed into proportionality sign.
                    prior = (1-w).*normpdf(s_grid_integrate, mu, sigma_s) + w.* normpdf(s_grid_integrate, 0, sigma_center);
                case "TwoGaussiansBothFixedZero"
                    delta_sigma=theta(end-2); w=theta(end-1);
                    sigma_s_larger = sigma_s + delta_sigma;
                    % Since 1/sqrt(2pi*sigma_s) and 1/sqrt(2pi*sigma_center) are
                    % constants wrt to s, hence absorbed into proportionality sign.
                    log_prior = log((1-w).*normpdf(s_grid_integrate, 0, sigma_s) + w.* normpdf(s_grid_integrate, 0, sigma_s_larger));
                case "GaussianLaplaceBothFixedZero"
                    % simga_s is width of Gaussian component, with mixture
                    % weight (1-w).
                    laplace_scale = theta(end-2); w=theta(end-1);
                    log_prior = log((1-w).*normpdf(s_grid_integrate, 0, sigma_s) + w./(2.*laplace_scale).*exp(-abs(s_grid_integrate-0)./laplace_scale));       

            end
            x_grid_rescale = x_grid .* rescale;            
            log_likelihood = -log(sigma_fun(s_grid_integrate, sigma0,k).*rel_multiply_factor)-0.5.*log(2*pi) - ((x_grid_rescale-s_grid_integrate).^2)./(2.*((sigma_fun(s_grid_integrate, sigma0,k).*rel_multiply_factor).^2));
            log_expr = log_likelihood + log_prior;
            log_expr_shifted = log_expr - max(log_expr);
            expr = exp(log_expr_shifted);
        %%
        
            s_hats_PM = sum(s_grid_integrate.*expr,1) ./ sum(expr,1); % the binwidth term gets factored out.           
            p_r_given_x = normpdf(R_rel, s_hats_PM, sigma_motor);        
            prob_r_given_s_nolapse = qtrapz(p_r_given_x .* p_x_given_s .* dx,2);
            switch lapse_type
                case "Uniform"
                    ll_idx = lapse./(max(lapse_range)-min(lapse_range)) + (1-lapse).* prob_r_given_s_nolapse; % Account for lapse (Unif(lapse_range) part of the mixture)
                case "Gaussian"
                    truncated_Gaussian_area_perc = normcdf(max(lapse_range), 0, sigma_lapse) - normcdf(min(lapse_range), 0, sigma_lapse);
                    ll_idx = lapse .* normpdf(R_rel, 0, sigma_lapse)./truncated_Gaussian_area_perc + (1-lapse).* prob_r_given_s_nolapse; 
                    %ll_idx = lapse .* normpdf(R_rel, 0, Gaussian_lapse_SD) + (1-lapse).* prob_r_given_s_nolapse; 
            end
            LL_idx(rel_trial_idx) = ll_idx;
        end
        
        
        if(return_response_distr)
            output=LL_idx;
        else
            NLL = -sum(log(LL_idx));
            if(isnan(NLL))
                clc
                %p_x_given_s
                theta
                %LL_idx
            end
            output=NLL;
        end  
        return
    end
         
end