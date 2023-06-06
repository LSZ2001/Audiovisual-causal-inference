function manuscript_unimodalfits_visualization(data_stratified, fitted_params_PM, ModelComponents, is_visual, colors, s_range, model_family, plot_lapse, is_pred_samples, fontsize, figspecs, lapse_type, Gaussian_lapse_SDs)
    %cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Analysis')
    return_predictive_samples = true;
    return_response_distribution=false;
    
    if(nargin<12)
        lapse_type = "Uniform"; Gaussian_lapse_SDs = NaN(1,15);
    elseif(nargin<13)
        Gaussian_lapse_SDs = NaN(1,15);
    end
    
    if(~isnan(figspecs))
        figure('Position', figspecs);
        set(gcf, 'Color', 'w')
        tcl=tiledlayout(4,4,'Padding', 'compact', 'TileSpacing', 'compact');
    end

    if(model_family=="parametric")
        hetero_type = ModelComponents.SensoryNoise;
    end
    
    num_rels = ModelComponents.NumReliabilityLevels;
    if(~is_visual)
        num_rels=4;
    end
    reliability_titles=["High","Med","Low"];

    
    mu=0; 
    [num_subjects, num_params] = size(fitted_params_PM);
    num_s_bins = length(data_stratified{1,1});
    R_round_decimals = 2; % R round to nearest integer for binning.
    dr = 10^(-R_round_decimals);
    R_grid = -45:dr:45;
    
    if(~is_visual)
        for i=1:num_subjects
            for j =1:num_s_bins
                data_stratified{i}{j}(:,3)=4;
            end
        end
    end
    
    subj_idx = cell(1,num_subjects);
    for i=1:num_subjects
        subj_idx{i} = "subj "+num2str(i);
    end
    
    hist_range = -60:3:60;
    binedges = [-Inf, hist_range(1:(end-1)), +Inf] + 1.5;
    
    PMR_distr_params = zeros(num_subjects, num_s_bins, 3);
    s_cond_s_hat_means = zeros(num_subjects, num_s_bins);
    s_cond_s_hat_stds = zeros(num_subjects, num_s_bins);
    s_cond_s_hat_iqrs = zeros(num_subjects, num_s_bins);
    raw_data_histcounts = zeros(num_subjects, num_s_bins, length(binedges)-1);
    
    PMR_vals = zeros(num_subjects, num_s_bins, length(binedges)-1);
    
    %%
    if(~is_visual) % Do not visualize UA model ribbons 3 times!
        rel_levels=4;
    else
        rel_levels = 1:num_rels;
    end
    for l=rel_levels  % If want to plot UA model ribbons, use "for l=4". 
        for j=1:num_s_bins
            if(is_pred_samples || is_visual) % If want to plot UA model ribbons, change argument to true.
                return_predictive_samples = true;
                return_response_distribution=false;

                dataset_size_bysubj = zeros(1,num_subjects);
                data = cell(1,num_subjects);
                for i=1:num_subjects
                    dataset_size_bysubj(i) = length(data_stratified{i}{j});
                    data_allrels = data_stratified{i}{j};
                    data{i} = data_allrels(data_allrels(:,3)==l,:);
                end
                dataset_startidx_bysubj = [0,cumsum(dataset_size_bysubj)];
                
                for i=1:num_subjects
                    Gaussian_lapse_SD = Gaussian_lapse_SDs(i);
                    switch model_family                  
                        case "parametric"
                            S_vals = data{i}(:,[1,3]);
                            samples = nllfun_uav_parametric(ModelComponents,fitted_params_PM(i,:),R_grid,S_vals, return_predictive_samples, return_response_distribution, plot_lapse, lapse_type, Gaussian_lapse_SD);
                        case "semiparam"
                            samples = nllfun_uav_semiparam(fitted_params_PM(i,:)', {data{i}}, return_predictive_samples, return_response_distribution, plot_lapse, lapse_type, Gaussian_lapse_SD);
                            samples = samples{1};
                        case "semiparamInsp"
                            ModelComponents.PriorUnnormalized = ModelComponents.PriorUnnormalizedAllSubjs(:,i);
                            ModelComponents.SigmaFuns = ModelComponents.SigmaFunsAllSubjs{i};
                            samples = nllfun_uav_semiparaminsp(fitted_params_PM(i,:),data{i},ModelComponents,return_predictive_samples, return_response_distribution,plot_lapse, lapse_type, Gaussian_lapse_SD);
                    end

                    PMR_distr_params(i,j,:) = [mean(samples(:)),mean(std(samples,[],1)), mean(iqr(samples,1))];
                    PMRs = [];
                    for k=1:(length(binedges)-1)
                        [l,i,j,k]
                        sample_counts = ((int8(samples(:)>=binedges(k)) .* int8(samples(:)<binedges(k+1)))~=0);
                        PMRs = [PMRs, sum(sample_counts)];
                    end
                    PMR_vals(i,j,:) = PMRs ./ sum(PMRs);
                end

            else % Auditory
                return_predictive_samples = false;
                return_response_distribution=true;

                dataset_size_bysubj = zeros(1,num_subjects);
                data = cell(1,num_subjects);
                for i=1:num_subjects
                    s = data_stratified{i}{j}(1,1);
                    dataset_size_bysubj(i) = length(R_grid);
                    data{i} = [repmat(s, length(R_grid),1), R_grid', repmat(4, length(R_grid),1)];
                end
                dataset_startidx_bysubj = [0,cumsum(dataset_size_bysubj)];
                
                for i=1:num_subjects
                    switch model_family
                        case "parametric"
                            S_vals = data{i}(:,1);
                            Gaussian_lapse_SD = Gaussian_lapse_SDs(i);
                            %PMRs_vector = midpoint_postmean_NLLfun_comprehensive(ModelComponents,fitted_params_PM(i,:),R_grid,S_vals, NaN,[-45,45,201], false, true,plot_lapse, lapse_type, Gaussian_lapse_SD);
                            PMRs_vector = nllfun_uav_parametric(ModelComponents,fitted_params_PM(i,:),R_grid,S_vals, false, true,plot_lapse, lapse_type, Gaussian_lapse_SD);
                        case "semiparam"
                            PMRs_vector = nllfun_uav_semiparam(fitted_params_PM(i,:)', {data{i}}, return_predictive_samples,return_response_distribution);
                        case "semiparamInsp"
                            ModelComponents.PriorUnnormalized = ModelComponents.PriorUnnormalizedAllSubjs(:,i);
                            ModelComponents.SigmaFuns = ModelComponents.SigmaFunsAllSubjs{i};
                            PMRs_vector = nllfun_uav_semiparaminsp(fitted_params_PM(i,:),data{i},ModelComponents,return_predictive_samples, return_response_distribution,plot_lapse);
                    end

                    Q1 = find(cumsum(PMRs_vector) <= 0.25*sum(PMRs_vector)); Q3 = find(cumsum(PMRs_vector) > 0.75*sum(PMRs_vector));
                    Q1 = R_grid(Q1(end));
                    Q3 = R_grid(Q3(1));
                    PMR_distr_params(i,j,:) = [sum(R_grid'.*PMRs_vector.*dr),sqrt(sum(R_grid'.^2.*PMRs_vector.*dr) - sum(R_grid'.*PMRs_vector.*dr).^2), Q3-Q1];
                    PMRs = [];
                    for k=1:(length(binedges)-1)
                        [l,i,j,k]
                        relevant_vector_entries = find((int8(R_grid>=binedges(k)) .* int8(R_grid<binedges(k+1)))~=0);
                        PMRs = [PMRs, sum(PMRs_vector(relevant_vector_entries))];
                    end
                    PMR_vals(i,j,:) = PMRs ./ sum(PMRs);
                end

            end
        end

        for i=1:num_subjects
            for j=1:num_s_bins
                s_matrix = data_stratified{i}{j};

                if(is_visual) % if is visual data; then use center of bins that we grouped S into.
                    reliabilities_alltrials = s_matrix(:,3);
                    rel_trial_idx = find(reliabilities_alltrials==l);
                    S_vals = [s_matrix(rel_trial_idx,1), s_matrix(rel_trial_idx,3)];
                    r_relevant = s_matrix(rel_trial_idx,2);
                else % Auditory, no reliability levels.
                    s = s_matrix(1,1);
                    S_vals = repmat(s, 1,length(R_grid));
                    r_relevant = s_matrix(:,2);
                end
                
                [N,edges] = histcounts(r_relevant, binedges, 'Normalization','pdf');
                raw_data_histcounts(i,j,:) = N ./ sum(N);
                s_cond_s_hat_means(i,j) = mean(r_relevant);
                s_cond_s_hat_stds(i,j) = std(r_relevant);
                s_cond_s_hat_iqrs(i,j) = iqr(r_relevant);
            end

        end

        % Loop again to plot things
        raw_data_histcounts_meanoversubj = squeeze(mean(raw_data_histcounts, 1));
        raw_data_histcounts_SDoversubj = squeeze(std(raw_data_histcounts));
        raw_data_histcounts_SEMoversubj = raw_data_histcounts_SDoversubj ./ sqrt(num_subjects);

        model_vals = PMR_vals; 
        model_params = PMR_distr_params;
        
        model_vals_meanoversubj = squeeze(mean(model_vals, 1));
        model_vals_SDoversubj = squeeze(std(model_vals));
        model_vals_SEMoversubj = model_vals_SDoversubj ./ sqrt(num_subjects);
        Curves1 = model_vals_meanoversubj + model_vals_SEMoversubj;
        Curves2 = model_vals_meanoversubj - model_vals_SEMoversubj;

        nexttile([1,2]);
        for j=1:num_s_bins
            hold on
            model_mean_vals = model_vals_meanoversubj(j,:);
            curve1 = squeeze(Curves1(j,:)); curve2 = squeeze(Curves2(j,:));
            x2 = [hist_range, fliplr(hist_range)];
            inBetween = [curve1, fliplr(curve2)];
            p=fill(x2, inBetween, colors(j,:),'FaceAlpha',0.25, 'EdgeColor', 'none');
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';

            keep_x_idx = find(raw_data_histcounts_meanoversubj(j,:) ~=0);

            errorbar(squeeze(hist_range(keep_x_idx)), raw_data_histcounts_meanoversubj(j,keep_x_idx), raw_data_histcounts_SEMoversubj(j,keep_x_idx), ".", 'Color', colors(j,:));

        end
%         if(is_visual) % visual data
%             xlim([-35,35])
%             ylim([0,0.7])
%         end
        xlim([-35,35]) %[-32,42]
        ylim([0,0.73])
        if(is_visual)
            xlabel("Estimated visual stimulus location (\circ)", 'FontSize', fontsize)
        else
            xlabel("Estimated auditory stimulus location (\circ)", 'FontSize', fontsize)
        end
        %ylabel({"\fontsize{10} \bf{Visual}","\rm{Proportion}"}, 'FontSize', fontsize)
        
        if(is_visual)
            ylabel({"{\fontsize{10} \bf{Visual ("+reliability_titles(l)+" Reliability)}}","{\rm \fontsize{9} {Proportion}}"})

            %title("Visual (Reliability " +reliability_titles(l)+")", 'FontSize', fontsize+1)
        else
            ylabel({"{\fontsize{10} \bf{Auditory}}","{\rm \fontsize{9} {Proportion}}"})
            %title("Auditory", 'FontSize', fontsize+1)
        end
        
        if(is_visual)
            %lg = legend("$s \in "+{"[-20, -100/7]","[-100/7,-60/7]","[-60/7,-20/7]","[-20/7,+20/7]","[+20/7,+60/7]","[+60/7,+100/7]","[+100/7,+20]"}+"$",'interpreter', 'latex', 'Position', [0.8,0.05,0.15,0.8], 'color','none');
            lg = legend("$s \in "+{"[-20, -14.3]","[-14.3,-8.6]","[-8.6,-2.9]","[-2.9,+2.9]","[+2.9,+8.6]","[+8.6,+14.3]","[+14.3,+20]"}+"$",'interpreter', 'latex', 'Position', [0.4053    0.8044    0.1307    0.1418], 'color','none');

        else
            lg = legend("$s = "+{"-15","-10","-5","0","+5","+10","+15"}+"$",'interpreter', 'latex', 'Location', 'northeast', 'color','none');
        end
        % set(lg, 'Position', [0.4053    0.8044    0.1307    0.1418])
        %lg.Position = [0.4053    0.8044    0.1307    0.1418];
        lg.ItemTokenSize(1) = 6;
        set(lg,'Box','off')
        lg.FontSize = 6;
        
        if(l==1)
            ttl = title('(a)', "Fontsize", 10);
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = -0.05; % use negative values (ie, -0.1) to move further left
            ttl.HorizontalAlignment = 'left';  
        end
        
%         nexttile;
%         xticklabels([]); yticklabels([]);
%         axis off

        % Means of histograms
        nexttile;
        hold on
        model_params_mu = model_params(:,:,1);
        model_params_mu_meanoversubj = squeeze(mean(model_params_mu,1));
        model_params_mu_SEMsoversubj = squeeze(std(model_params_mu))./ sqrt(num_subjects);
        s_cond_s_hat_means_meanoversubj = squeeze(mean(s_cond_s_hat_means,1));
        s_cond_s_hat_means_SEMsoversubj = squeeze(std(s_cond_s_hat_means)) ./ sqrt(num_subjects);
        curve1 = model_params_mu_meanoversubj + model_params_mu_SEMsoversubj;
        curve2 = model_params_mu_meanoversubj - model_params_mu_SEMsoversubj;
        x2 = [s_range, fliplr(s_range)];
        inBetween = [curve1, fliplr(curve2)];
        p=fill(x2, inBetween, "k",'FaceAlpha',0.2, 'EdgeColor', 'none');
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h2 = plot([s_range(1),s_range(end)], [s_range(1),s_range(end)], "k--");
        h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        for j=1:num_s_bins
            errorbar(squeeze(s_range), s_cond_s_hat_means_meanoversubj, s_cond_s_hat_means_SEMsoversubj, ".", 'Color', 'k');
        end
        xlim([-20,20])
        ylabel("Mean loc estimate (\circ)", 'FontSize', fontsize)
        if(is_visual)
            xlabel("Visual stimulus location (\circ)", 'FontSize', fontsize)
        else
            xlabel("Auditory stimulus location (\circ)", 'FontSize', fontsize)
        end
        if(l==1)
            ttl = title('(b)', "Fontsize", 10);
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = -0.29; % use negative values (ie, -0.1) to move further left
            ttl.HorizontalAlignment = 'left';  
        end

        % SDs of histograms
        nexttile;
        hold on
        model_params_sigma = model_params(:,:,2);
        model_params_sigma_meanoversubj = squeeze(mean(model_params_sigma,1));
        model_params_sigma_SEMsoversubj = squeeze(std(model_params_sigma))./ sqrt(num_subjects);
        s_cond_s_hat_stds_meanoversubj = squeeze(mean(s_cond_s_hat_stds,1));
        s_cond_s_hat_stds_SEMsoversubj = squeeze(std(s_cond_s_hat_stds)) ./ sqrt(num_subjects);
        curve1 = model_params_sigma_meanoversubj + model_params_sigma_SEMsoversubj;
        curve2 = model_params_sigma_meanoversubj - model_params_sigma_SEMsoversubj;
        x2 = [s_range, fliplr(s_range)];
        inBetween = [curve1, fliplr(curve2)];
        p=fill(x2, inBetween,"k",'FaceAlpha',0.2, 'EdgeColor', 'none');
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        for j=1:num_s_bins
            errorbar(squeeze(s_range), s_cond_s_hat_stds_meanoversubj, s_cond_s_hat_stds_SEMsoversubj, ".", 'Color', 'k');
        end
        xlim([-20,20])
        ylim([0,6])
        ylabel("SD of loc estimate (\circ)", 'FontSize', fontsize)
        if(is_visual)
            xlabel("Visual stimulus location (\circ)", 'FontSize', fontsize)
        else
            xlabel("Auditory stimulus location (\circ)", 'FontSize', fontsize)
        end    
        if(l==1)
            ttl = title('(c)', "Fontsize", 10);
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = -0.23; % use negative values (ie, -0.1) to move further left
            ttl.HorizontalAlignment = 'left';  
        end
        
%         % IQRs of histograms
%         nexttile;
%         hold on
%         model_params_iqr = model_params(:,:,3);
%         model_params_iqr_meanoversubj = squeeze(mean(model_params_iqr,1));
%         model_params_iqr_SEMsoversubj = squeeze(std(model_params_iqr))./ sqrt(num_subjects);
%         s_cond_s_hat_iqrs_meanoversubj = squeeze(mean(s_cond_s_hat_iqrs,1));
%         s_cond_s_hat_iqrs_SEMsoversubj = squeeze(std(s_cond_s_hat_iqrs)) ./ sqrt(num_subjects);
%         curve1 = model_params_iqr_meanoversubj + model_params_iqr_SEMsoversubj;
%         curve2 = model_params_iqr_meanoversubj - model_params_iqr_SEMsoversubj;
%         x2 = [s_range, fliplr(s_range)];
%         inBetween = [curve1, fliplr(curve2)];
%         p=fill(x2, inBetween,"k",'FaceAlpha',0.1, 'EdgeColor', 'none');
%         p.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         for j=1:num_s_bins
%             errorbar(squeeze(s_range), s_cond_s_hat_iqrs_meanoversubj, s_cond_s_hat_iqrs_SEMsoversubj, ".", 'Color', 'k');
%         end
%         xlim([-20,20])
%         ylim([0,8])
%         ylabel("IQR of Loc Estimate (\circ)", 'FontSize', fontsize)
%         if(is_visual)
%             xlabel("Visual Stimulus Location (\circ)", 'FontSize', fontsize)
%         else
%             xlabel("Auditory Stimulus Location (\circ)", 'FontSize', fontsize)
%         end
        %cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Analysis')
    end

end