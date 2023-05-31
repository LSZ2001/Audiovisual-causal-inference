% This code assumes a free audivisual rescaling parameter for auditory
% trials.

function [out_struct] = fit_ujointmodel_semiparam(iter, num_inits, data_path, model_path_temp)
    if(nargin==0)
        iter=1; num_inits = 81; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==1)
        num_inits = 81; data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==2)
        data_path = "data\"; model_path_temp = "modelfits\temp\";
    elseif(nargin==3)
        model_path_temp = "modelfits\temp\";
    end
    lapse_type = "Uniform";
    num_params = 40;
    
    %%
    s_a_range = -15:5:15;
    s_v_range = -20:1:20;

    model_type="PM"; % two gaussians components of the prior both centered at 0. 
    PMIntegrationParams = [-45,45,201]; % PM midpoint Rule bounds and numbins.
    consider_lapse=true; % fit a lapse parameter.
    randomize_trials=false;
    mu = 0; % prior distribution mean

    % UV data model components
    ModelComponents.SPivot = [0,0.1,0.3,1,2,4,6,8,10,15,20,45];
    ModelComponents.LapseRange = [-45,45];
    ModelComponents.RescaleVis="1"; % fit k_vis
    ModelComponents.RescaleAud="free"; % fit k_aud
    ModelComponents.NumReliabilityLevels = 3; % Get number of visual reliability levels
    ModelComponents.MotorNoise = "Gaussian";

    num_pivots = length(ModelComponents.SPivot);

    load(data_path+"data_stratified_uv.mat");
    load(data_path+"data_stratified_ua.mat");
    [~, num_subjects] = size(data_stratified_UA);
    % Unstratify data. 
    data_UV = data_stratified_to_data(data_stratified_UV, false, true); % last argument is is_visual.
    data_UA = data_stratified_to_data(data_stratified_UA, false, false);
    UAV_data = cell(1,num_subjects);
    Gaussian_lapse_SDs = zeros(1,num_subjects); % If assume Gaussian Lapse, then its SD is just the SD across all UAV trials for that subject. 
    for i=1:num_subjects
        data_UA{i}(:,3) = 4;
        UAV_data{i} = [data_UV{i}; data_UA{i}];
        Gaussian_lapse_SDs(i) = std(UAV_data{i}(:,2));
    end

    %% Get inits with flat prior and/or flat sigma(s) vis and/or flat sigma(s) aud. 
        
    % I am creating theta_inits for all the iters per iteration, but in
    % fact only one of them is used. Then in the next iteration, everything
    % is recreated. This is redundant, but not very time-consuming.
    theta_inits = zeros(num_inits, num_subjects, num_params);
    
    % Flat prior
    theta_inits(:,:,(2*num_pivots+1):(3*num_pivots-1)) = 0;
    
    % Flat sigma(s) vis, sigma(s) aud
    sigmafun_const_vals = [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]; 
    num_sigmafun_vals = length(sigmafun_const_vals);
    for i=1:num_sigmafun_vals
        theta_inits((num_sigmafun_vals *(i-1)+1):(num_sigmafun_vals *i),:,1) = sigmafun_const_vals(i);
    end
    theta_inits(:,:,2:num_pivots) = 0;

    for i=1:num_inits
        i_mod = mod(i,num_sigmafun_vals);
        if(i_mod==0)
            i_mod = num_sigmafun_vals;
        end
        theta_inits(i,:,(1+num_pivots)) = sigmafun_const_vals(i_mod);
    end
    theta_inits(:,:,(2+num_pivots):(2*num_pivots)) = 0;
    
    
    % 5 remaining params -- uniform between PLB and PUB
    plb = [1,1,0,0.25,1]'; pub = [5,10,0.05,1,2]';
    PLB = permute(repmat(plb,1,num_inits, num_subjects), [2,3,1]);
    PUB = permute(repmat(pub,1,num_inits, num_subjects), [2,3,1]);
    theta_inits(:,:,(end-4):end) = rand(size(PLB)).*(PUB-PLB) + PLB;

    % insigma
    insigma = zeros(num_params,1);
    insigma(1:(3*num_pivots-1)) = 0.5;
    insigma((3*num_pivots):end) = (pub-plb)./sqrt(12);

    %% Plot the inits
% s_pivot = ModelComponents.SPivot;
% s_pivot_full = [-fliplr(s_pivot(2:end)), s_pivot];
% s_fine = linspace(0,15,2^6);
% num_pivots = length(s_pivot);
% 
% for fun_idx=1:3
%     figure
%     tiledlayout(3,5,'TileSpacing','compact');
%     for subj=1:num_subjects
%         [fun_idx, subj]
%         nexttile;
%         hold on
%         title("Subject "+subj)
%         xlabel("$s$", 'interpreter','latex') 
%         for init=1:num_inits_new
%             theta=squeeze(theta_inits(init,subj,:))';
%             col = "c:";
%             if(sum(theta(2:num_pivots)==0)>=11)
%                 col = "r:";
%             elseif(sum(theta((2+num_pivots):(2*num_pivots))==0)>=11)
%                 col = "r:";
%             elseif(sum(theta((1+2*num_pivots):(3*num_pivots-1))==0)>=11)
%                 col = "r:";
%                 
%             else
%                 continue;
%             end
%             
%             switch fun_idx
%                 case 1
%                     sigma_fun_vis_rel_high_pivots = cumsum(theta(1:num_pivots));
%                     sigma_fun_vis_rel_high_pivots = [fliplr(sigma_fun_vis_rel_high_pivots(2:end)), sigma_fun_vis_rel_high_pivots];
%                     sigma_fun_vis = @(s)  min([repmat(45,length(s),1)' ; interp1(s_pivot_full, exp(sigma_fun_vis_rel_high_pivots), s, 'pchip')], [], 1);       
%                     h(init)=plot(s_fine, log(sigma_fun_vis(s_fine)),col);
%                     ylabel("$\log \sigma_{vis}(s)$", 'interpreter','latex') 
%                 case 2
%                     sigma_fun_aud_pivots = cumsum(theta((num_pivots+1):(2*num_pivots)));
%                     sigma_fun_aud_pivots = [fliplr(sigma_fun_aud_pivots(2:end)), sigma_fun_aud_pivots];
%                     sigma_fun_aud = @(s)  min([repmat(45,length(s),1)' ; interp1(s_pivot_full, exp(sigma_fun_aud_pivots), s, 'pchip')], [], 1);       
%                     h(init)=plot(s_fine, log(sigma_fun_aud(s_fine)),col);
%                     ylabel("$\log \sigma_{aud}(s)$", 'interpreter','latex') 
%                 case 3    
%                     prior_pivots = cumsum([1,theta((2*num_pivots+1):(3*num_pivots-1))]);
%                     prior_pivots = [fliplr(prior_pivots(2:end)), prior_pivots];
%                     prior = @(s) exp(interp1(s_pivot_full, prior_pivots, s, 'pchip'));
%                     h(init)=plot(s_fine, log(prior(s_fine)),col);
%                     ylabel("$\log p(s)$", 'interpreter','latex') 
%             end
% 
%         end
%     end
% end
% 
% x_labels_pos = (3*num_pivots):num_params;
% x_labels = {"scale med", "scale low", "lapse","sigma motor","rescale aud"};
% 
% figure
% tiledlayout(3,5,'TileSpacing','compact');
% for subj=1:num_subjects
%     nexttile;
%     hold on
%     
%     plot(repmat(x_labels_pos,1,1), squeeze(theta_inits(:,subj,x_labels_pos)), "c.");
%     xticks(x_labels_pos);
%     xticklabels(x_labels);
%     xtickangle(45);
%     ax=gca;
%     ax.XAxis.FontSize = 7;
%     title("Subject " + subj)
% 
% end

    
    %% LB and UB of params
    % [sigma(s)_vis pivots, sigma(s)_aud pivots, prior_pivots, scale_vis_med,
    % scale_vis_low, lapse, sigma_motor, rescale_aud]. 
    lb = [-5,repmat(-eps,1,num_pivots-1),-5,repmat(-eps,1,num_pivots-1), repmat(-15,1,num_pivots-2),-70,0.7,0.7,0,0.249,1/3];
    ub = [4,repmat(5,1,num_pivots-1),4,repmat(5,1,num_pivots-1), repmat(eps,1,num_pivots-1), 20,20,1,2,3];


    theta_inits = min(max(theta_inits, permute(repmat(lb'+2*eps,1,num_inits, num_subjects),[2,3,1])),permute(repmat(ub'-2*eps,1,num_inits, num_subjects),[2,3,1]));
    
    %% Check that theta_inits are within LB, UB.
%     for i=1:(num_subjects*num_inits_new)
%         init = mod(i, num_inits_new);
%         if(init==0)
%             init=num_inits_new;
%         end
%         subj = ceil(i/num_inits_new);
%         theta0 = squeeze(theta_inits(init,subj,:));
%         indx = [subj,init];
%         if(prod(theta0>lb)~=1)
%             indx
%             q = "LB: "
%             theta0>lb
%         elseif(prod(theta0<ub)~=1)
%             indx
%             q = "UB: "
%             theta0<ub
%             theta0
%         end
%     end
%     
   

    %% Fit nonparam models starting from nonparam inits.
    fun_name = convertStringsToChars("nllfun_uav_semiparam");

    init = mod(iter, num_inits);
    if(init==0)
        init=num_inits;
    end
    subj = ceil(iter/num_inits);

    % Setup options for CMA-ES optimization
    cmaes_opts = cmaes_modded('defaults');
    cmaes_opts.MaxFunEvals = 10*10^4;
    cmaes_opts.TolX = '1e-11*max(insigma)';
    cmaes_opts.TolFun = 1e-12;
    cmaes_opts.TolHistFun = 1e-13;
    cmaes_opts.EvalParallel = 'no';
    cmaes_opts.DispFinal = 'off';
    cmaes_opts.SaveVariables = 'off';
    cmaes_opts.DispModulo = Inf;
    cmaes_opts.LogModulo = 0;
    cmaes_opts.CMA.active = 1;      % Use Active CMA (generally better)

    cmaes_opts.DispFinal = 'on';
    cmaes_opts.DispModulo = 10;
    cmaes_opts.LBounds = squeeze(lb)';
    cmaes_opts.UBounds = squeeze(ub)';
    cmaes_opts.Seed = subj*100+init;

    indx=[subj, init]   
    %nllfun = @(theta) midpoint_NLLfun_comprehensive_nonparam_mixture_indv_CMAES(theta, {data{j}});
    theta0 = squeeze(theta_inits(init,subj,:)); %+ randn(1,num_params); 
    tStart = tic;
    [theta_fitted_iter, f_val,COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes(fun_name,...
                                                theta0, ...    % objective variables initial point, determines N
                                                insigma, ...   % initial coordinate wise standard deviation(s)
                                                cmaes_opts, ...    % options struct, see defopts below
                                                {UAV_data{subj}}, false, false, true, lapse_type, Gaussian_lapse_SDs(subj) ...   % arguments passed to objective function 
                                                )
    t_end = toc(tStart);

    out_struct.iter = iter;
    out_struct.subj = subj;
    out_struct.init = init;
    out_struct.Theta_fitted = theta_fitted_iter;
    out_struct.F_val = f_val;
    out_struct.T_end = t_end;
    out_struct.theta0 = theta0;
    
    save(model_path_temp+"fittedparams_UJoint_Semiparam_rescalefree_lapseUniform_"+iter, 'out_struct', 'cmaes_opts');
    
    Q = "Code finished."

end
