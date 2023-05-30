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

% Takes in either BA or BV data (but not both) and returns NLL. (fittype
% either "BA" or "BV"
function [output] = NLLfun_BAV_UBresc_semiparamInsp(theta, data, ModelComponents, return_predictive_samples, return_response_distr, plot_consider_lapse, lapse_type, Gaussian_lapse_SD)
    if(nargin==3) % default is to return NLL, not posterior preditive samples vectorized response distribition.
        return_predictive_samples = false; return_response_distr=false; plot_consider_lapse=true;
        lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==4) % if return posterior predictive samples and not vectorized response distribution, default is not include lapse mixture component.
        return_response_distr=false; plot_consider_lapse=true;
        lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif (nargin==5) % if return vectorized response distribution (need S to be one s-value and rel-level repmat, and R to be r_grid), default is not include lapse mixture component.
        plot_consider_lapse=true;
        lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==6)
        lapse_type = "Uniform";  Gaussian_lapse_SD = NaN;
    elseif(nargin==7)
        Gaussian_lapse_SD = NaN;
    end
    lapse_range = ModelComponents.LapseRange;
    causal_inf_strategy = ModelComponents.CausalInfStrategy;
    s_grid_integrate = ModelComponents.SGridIntegrate;
    ds = s_grid_integrate(2) - s_grid_integrate(1);
    num_bins = length(s_grid_integrate);
    num_rels = 4; % 3 reliability levels for V.
    
    
    %% Parse dataset.
    S_V = data(:,[1,3]); % include reliability level info for each trial.
    S_A = data(:,2);
    R = data(:,[4,5]); % R is just C in this task type.
    [~,ns] = size(S_V);
    [~,nr] = size(R);
    if(ns>2)
        S_V=S_V';
        S_A=S_A';
    end
    if(nr>2)
        R=R';
    end
    response_types = R(:,2);
    R = R(:,1);
    reliabilities = S_V(:,1);
    S_V = S_V(:,2);
    num_trials = length(S_A);
    
    
    %% Parameter extraction
    % [scale_rel_med, scale_rel_low, lapse, sig_motor, rescale_aud, p_same,
    % Bimod_rescale_vis, Bimod_rescale_aud, UB_shared_prior_power];
    rel_multiply_factors_uniq =[1, theta(1:2), 1];
    lapse = theta(3); sigma_motor = theta(4);
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
    
    % Get normalized prior density.
    p_s_unnorm_base = ModelComponents.PriorUnnormalized;
    p_s_unnorm = p_s_unnorm_base .^ UB_shared_prior_power;
    prior = p_s_unnorm / qtrapz(p_s_unnorm*(s_grid_integrate(2)-s_grid_integrate(1)));
    log_prior = log(prior);
        

    %% For response distribution visualization -- need to change the sigma_fun and prior def. 
    if(return_predictive_samples) % If return predictive samples, cannot return response distribution despite input args.
        return_response_distr=false;
        num_samps = 1; % Number of samples for this trial. 

        task_types = ["BV", "BA"];
        % Assume all inputs have the same task type as Trial 1.
        tasktype = task_types(response_types(1));
        % Assume all inputs have the same reliability level as Trial 1.
        rel_multiply_factor = rel_multiply_factors_uniq(reliabilities(1));
        sigma_vis = sigma_fun_vis(S_V) .* rel_multiply_factor;
        sigma_aud = sigma_fun_aud(S_A);

        x_V = S_V + sigma_vis.*randn([num_trials,num_samps]);
        x_A = S_A + sigma_aud.*randn([num_trials,num_samps]);

        ds = s_grid_integrate(2) - s_grid_integrate(1);

        x_V_rescale = x_V';   
        x_A_rescale = x_A' .* rescale_aud;  

        % These are nonGaussian likelihood functions of s_grid_integrate.
        p_xV_given_s_integrate = 1./(sqrt(2.*pi).*sigma_fun_vis(s_grid_integrate).*rel_multiply_factor) .* exp(-((x_V_rescale-s_grid_integrate).^2)./(2.*((sigma_fun_vis(s_grid_integrate).*rel_multiply_factor).^2)));
        p_xA_given_s_integrate = 1./(sqrt(2.*pi).*sigma_fun_aud(s_grid_integrate)) .* exp(-((x_A_rescale-s_grid_integrate).^2)./(2.*((sigma_fun_aud(s_grid_integrate)).^2)));

        protoposteriors_C = zeros(2, num_trials);
        expr_C1 =  prior .* p_xV_given_s_integrate .* p_xA_given_s_integrate;
        overall_likelihood_C1 = squeeze(sum(expr_C1 .* ds,1)); % This is now [x_A, x_V].           
        protoposteriors_C(1,:) = overall_likelihood_C1 .* p_same; %./ sum(overall_likelihood_C1(:));

        clear expr_C1


        expr_C2_V = sum(prior .* p_xV_given_s_integrate .* ds,1);
        expr_C2_A = sum(prior .* p_xA_given_s_integrate .* ds,1);
        overall_likelihood_C2 = squeeze(expr_C2_V .* expr_C2_A);
        protoposteriors_C(2,:) = overall_likelihood_C2 .* (1-p_same); %./ sum(overall_likelihood_C2(:));

        % Pr[C|x_V, x_A] for C={1,2}.
        protoposteriors_C =  protoposteriors_C';
        posteriors_C = protoposteriors_C ./ sum(protoposteriors_C,2);

        PMintegrand_C1 =  prior .* p_xV_given_s_integrate .* p_xA_given_s_integrate;
        %E[s|C=1, xV, xA]
        sPM_C1 = squeeze((sum(PMintegrand_C1.*s_grid_integrate.*ds,1)))'; % This is now [x_A, x_V].                                     
        sPM_C1 = sPM_C1 ./ squeeze((sum(PMintegrand_C1.*ds,1)))';
        clear PMintegrand_C1

        switch tasktype
            case "BV"
                PMintegrand_C2 =  prior .* p_xV_given_s_integrate;
            case "BA"
                PMintegrand_C2 =  prior .* p_xA_given_s_integrate;
        end
        %E[s (either V or A)|C=2, xV, xA]
        sPM_C2 = squeeze((sum(PMintegrand_C2.*s_grid_integrate.*ds,1)))'; % This is now [x_A, x_V].           
        sPM_C2 = sPM_C2 ./ squeeze((sum(PMintegrand_C2.*ds,1)))';
        clear PMintegrand_C2

        switch causal_inf_strategy
            case "ProbMatching"
                believe_C1 = rand(num_trials,1) < squeeze(posteriors_C(:,1));
                s_PM = believe_C1.* sPM_C1 + (1-believe_C1).* sPM_C2;
            case "ModelSelection"
                believe_C1 = (squeeze(posteriors_C(:,1))>=0.5);
                s_PM = believe_C1.* sPM_C1 + (1-believe_C1).* sPM_C2;
            case "ModelAveraging"
                % 1*115*132*1 1*115*132*1
                s_PM = squeeze(posteriors_C(:,1)) .* sPM_C1 + squeeze(posteriors_C(:,2)) .* sPM_C2;                           
        end
        lapse_idx = rand(size(s_PM)) < lapse;
        % Random responses (uniform in -45:1:45)
        switch lapse_type
            case "Uniform"
                % Random responses (uniform in -45:1:45)
                lapse_val = min(lapse_range) + (max(lapse_range)-min(lapse_range))*rand(sum(lapse_idx),1);
            case "Gaussian"
                lapse_val = Gaussian_lapse_SD .* randn(sum(lapse_idx),1);
        end
        s_PM(lapse_idx) = lapse_val;
        % I now store the x_obs that lead to s\hat.
        output = s_PM;
        return
    end

    
    
    %% All below is NLL evaluation for model fitting
    LL_idx = zeros(num_trials,1);
    for resp_type = 1:2 % BV, then BA
    for i=1:num_rels
        rel_multiply_factor = rel_multiply_factors_uniq(i);
        rel_trial_idx = find(((reliabilities==i).* (response_types==resp_type))==1);
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

%             Nx_V = 2^6;
%             x_V_grid = linspace(x_V_min, x_V_max, Nx_V); % Need to transpose due to matrix operation x-s.
%             dx_V = x_V_grid(2) - x_V_grid(1);
%             if(dx_V > dx_max) % Ensure that dx is much smaller than the LB (0.25) of sigma_motor.
%                 dx_V=dx_max;
%                 x_V_grid = (x_V_min:dx_V:x_V_max);
%                 Nx_V = length(x_V_grid);
%             end

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
        sigma_fun_vis = @(s, sigma0, k) max(sigma_fun_vis(s), grid_size_fun(s));
        sigma_fun_aud = @(s, sigma0, k) max(sigma_fun_aud(s), grid_size_fun(s));



%             Nx_A = 2^6;
%             x_A_grid = linspace(x_A_min, x_A_max, Nx_A); % Need to transpose due to matrix operation x-s.
%             dx_A = x_A_grid(2) - x_A_grid(1);
%             if(dx_A > dx_max) % Ensure that dx is much smaller than the LB (0.25) of sigma_motor.
%                 dx_A=dx_max;
%                 x_A_grid = (x_A_min:dx_A:x_A_max);
%                 Nx_A = length(x_A_grid);
%             end

        % 3D tensors: (R_rel, xA, xV). 
        p_x_given_s_V_raw = normpdf(x_V_grid, S_V_rel, sigma_fun_vis(S_V_rel).*rel_multiply_factor);
        p_x_given_s_A_raw = normpdf(x_A_grid, S_A_rel, sigma_fun_aud(S_A_rel));

        p_x_given_s_V = reshape(p_x_given_s_V_raw, [length(S_V_rel), 1, Nx_V]);
        p_x_given_s_A = reshape(p_x_given_s_A_raw, [length(S_A_rel), Nx_A, 1]);

        clear p_x_given_s_V_raw p_x_given_s_A_raw

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

                    %log p(xV,xA|C=1)
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

                    d = log(p_same/(1-p_same)) + overall_log_likelihood_C1 - overall_log_likelihood_C2;
                    posteriors_C(1,:,:) = 1 ./ (1 + exp(-d));
                    posteriors_C(2,:,:) = 1-posteriors_C(1,:,:);

                    clear  d overall_log_likelihood_C1 log_expr_C2_A log_offset_C2_A log_expr_C2_V log_offset_C2_V expr_C2_V_term expr_C2_A expr_C2_V expr_C2_A_term overall_log_likelihood_C2
                else
                    posteriors_C(1,:,:)=p_same;
                    posteriors_C(2,:,:)=1-p_same;
                end



                % Below, dimensions [s_grid_integrate, x_A_grid, x_V_grid]
                log_PMintegrand_C1 =  log_prior + log_p_xV_given_s_integrate + log_p_xA_given_s_integrate;
                log_offset_PMintegrand_C1 = max(log_PMintegrand_C1,[],1);
                log_PMintegrand_C1 = log_PMintegrand_C1 - log_offset_PMintegrand_C1;
                PMintegrand_C1 = exp(log_PMintegrand_C1);

                clear log_PMintegrand_C1 log_offset_PMintegrand_C1

                %E[s|C=1, xV, xA]
                % Both numerator and denom have factor exp(log_offset_PMintegrand_C1), hence it is factored out.
                sPM_C1 = squeeze((qtrapz(PMintegrand_C1.*s_grid_integrate.*ds,1))); % This is now [x_A, x_V].                                     
                sPM_C1 = sPM_C1 ./ squeeze((qtrapz(PMintegrand_C1.*ds,1))); % This is now [x_A, x_V].                                    
                clear log_PMintegrand_C1 log_offset_PMintegrand_C1 PMintegrand_C1

                switch resp_type
                    case 1
                        log_PMintegrand_C2 =  log_prior + log_p_xV_given_s_integrate;
                    case 2
                        log_PMintegrand_C2 =  log_prior + log_p_xA_given_s_integrate;
                end
                log_offset_PMintegrand_C2 = max(log_PMintegrand_C2,[],1);
                log_PMintegrand_C2 = log_PMintegrand_C2 - log_offset_PMintegrand_C2;
                PMintegrand_C2 = exp(log_PMintegrand_C2);

                clear log_prior log_PMintegrand_C2 log_offset_PMintegrand_C2 log_p_xV_given_s_integrate log_p_xA_given_s_integrate

                %E[s (either V or A)|C=2, xV, xA]
                sPM_C2 = squeeze((qtrapz(PMintegrand_C2.*s_grid_integrate.*ds,1)))'; % This is now [x_A, x_V].           
                sPM_C2 = sPM_C2 ./ squeeze((qtrapz(PMintegrand_C2.*ds,1)))'; % This is now [x_A, x_V].                                     
                clear PMintegrand_C2

                switch causal_inf_strategy
                    case "ProbMatching"
                        sPM_C1_mat = zeros(1, Nx_A, Nx_V); 
                        sPM_C1_mat(1,:,:) = sPM_C1; 
                        if(length(sPM_C2)==Nx_V)
                            sPM_C2_mat = zeros(1, 1, Nx_V);
                            sPM_C2_mat(1,1,:) = sPM_C2;
                        else
                            sPM_C2_mat = zeros(1, Nx_A, 1);
                            sPM_C2_mat(1,:,1) = sPM_C2;
                        end
                        p_r_given_x = posteriors_C(1,:,:) .* normpdf(R_rel, sPM_C1_mat, max(grid_size_fun(sPM_C1_mat),sigma_motor)) + posteriors_C(2,:,:).*normpdf(R_rel, sPM_C2_mat,  max(grid_size_fun(sPM_C2_mat),sigma_motor));

                        clear sPM_C1_mat sPM_C2_mat

                    case "ModelSelection"
                        s_hats_PM = zeros(1, Nx_A, Nx_V);
                        believe_C1 = (squeeze(posteriors_C(1,:,:))>=0.5);
                        s_hats_PM(1,:,:) = believe_C1.* sPM_C1 + (1-believe_C1).* sPM_C2;
                        %180*132*115         180*1*1  1*132*115    1*1*1 
                        p_r_given_x = normpdf(R_rel, s_hats_PM, max(grid_size_fun(s_hats_PM),sigma_motor));
                    case "ModelAveraging"
                        s_hats_PM = zeros(1, Nx_A, Nx_V);
                        % 1*115*132*1 1*115*132*1
                        s_hats_PM(1,:,:) = squeeze(posteriors_C(1,:,:)) .* sPM_C1 + squeeze(posteriors_C(2,:,:)) .* sPM_C2;                           
                        %180*132*115         180*1*1  1*132*115    1*1*1 
                        p_r_given_x = normpdf(R_rel, s_hats_PM, max(grid_size_fun(s_hats_PM),sigma_motor));
                end

                clear sPM_C1 sPM_C2 posteriors_C

                %p_r_given_x = 1/(sqrt(2*pi)*sigma_motor) .* exp(-0.5.*((R_rel-s_hats_PM)./sigma_motor).^2);
                %180*132*115            180*115*1        180*1*132
                Expr = p_r_given_x .* p_x_given_s_A .* p_x_given_s_V;
                clear p_r_given_x p_x_given_s_A p_x_given_s_V

                dx_A([1, end]) = dx_A([1, end])./2;
                dx_V([1, end]) = dx_V([1, end])./2;

                dx_A_changeidx = find(((x_A_grid==-15) + (x_A_grid==-5) + (x_A_grid==5) + (x_A_grid==15))>0);
                dx_A(dx_A_changeidx) = (dx_A(dx_A_changeidx+1) + dx_A(dx_A_changeidx-1)) ./2;
                dx_V_changeidx = find(((x_V_grid==-15) + (x_V_grid==-5) + (x_V_grid==5) + (x_V_grid==15))>0);
                dx_V(dx_V_changeidx) = (dx_V(dx_V_changeidx+1) + dx_V(dx_V_changeidx-1)) ./2;

                clear dx_V_changeidx dx_A_changeidx
                %180*1*1
                prob_r_given_s_nolapse = sum(sum(Expr .* dx_A .* dx_V,2),3);

            case "empirical"
        end

        clear Expr
        switch lapse_type
            case "Uniform"
                ll_idx = lapse./(max(lapse_range)-min(lapse_range)) + (1-lapse).* prob_r_given_s_nolapse; % Account for lapse (Unif(lapse_range) part of the mixture)
            case "Gaussian"
                ll_idx = lapse .* normpdf(R_rel, 0, Gaussian_lapse_SD) + (1-lapse).* prob_r_given_s_nolapse; 
        end
        LL_idx(rel_trial_idx) = ll_idx;
    end
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