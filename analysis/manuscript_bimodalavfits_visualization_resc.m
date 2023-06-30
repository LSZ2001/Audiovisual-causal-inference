function manuscript_bimodalavfits_visualization_resc(BAV_data, fitted_params_PM, ModelComponents_V, ModelComponents_A, PMIntegrationParams, dx_max, model_family, plot_lapse, causal_inf_strategy, unique_bins_subj, fontsize, figspecs, figpath, save_name, png_dpi, lapse_type, Gaussian_lapse_SDs, plot_individual)
    % Note: if model_family=="semiparamInsp", then ModelComponents_A,
    % PMIntegrationParams, dx_max, causal_inf_strategy arguments are useless.
    if(nargin<15)
        lapse_type = "Uniform"; Gaussian_lapse_SDs = NaN(1,15); plot_individual=false;
    elseif(nargin<16)
        Gaussian_lapse_SDs = NaN(1,15); plot_individual=false;
    elseif(nargin<17)
        plot_individual=false;
    end
    
    return_predictive_samples = true;
    return_response_distr = false;
    num_rels = 3;
    reliability_titles=["High","Med","Low"];
    switch model_family
        case "parametric"
            hetero_type = ModelComponents_V.SensoryNoise;
        case "semiparamInsp"
            % ModelComponents is shared between V and A.
            ModelComponents = ModelComponents_V;
            dx_max = NaN; % This param not needed.
            % ModelComponents already contains causal_inf_strategy
            causal_inf_strategy = ModelComponents.CausalInfStrategy;
    end
    
    num_subjects = length(BAV_data); % length(BAV_data)
    R_round_decimals = 2; % R round to nearest integer for binning.
    dr = 10^(-R_round_decimals);
    R_grid = -45:dr:45;
    
    subj_idx = cell(1,num_subjects);
    s_diff_alltrials = [];
    s_sum_alltrials = [];
    BAV_data_alltrials = [];
    for i=1:num_subjects
        subj_idx{i} = "subj "+num2str(i);
        BAV_data_alltrials = [BAV_data_alltrials; BAV_data{i}];
        s_diff_alltrials = [s_diff_alltrials; BAV_data{i}(:,2) - BAV_data{i}(:,3)];
        s_sum_alltrials = [s_sum_alltrials; BAV_data{i}(:,2)+ BAV_data{i}(:,3)];
    end
%     binedges = quantile(s_diff_alltrials(abs(s_diff_alltrials)>0.01), linspace(0,1,7));
%     bincenters = binedges(1:(end-1)) + (binedges(2:end) - binedges(1:(end-1))) ./2;    
    num_binedges = 8;
  
    stratify_labels = ["({\its}_A+{\its}_V) on left", "({\its}_A+{\its}_V) in center", "({\its}_A+{\its}_V) on right"];
    stratify_thres = quantile(s_sum_alltrials, [0.25, 0.75]);
    for response_type=1:2
    bincenters_stratrels = zeros(3,3,num_subjects,num_binedges-1);
    Mean_Biases = zeros(3,3, num_subjects, num_binedges-1);
    Mean_Biases_modelfit = Mean_Biases;
    
    for l=1:3 %3, Need changing back to 1:3
    for strats=1:3 %3, Need changing back to 1:3
    %figure
        
        for i=1:num_subjects % 4, Need changing back to 1:num_subjects
            
            % Get s\hat bins for this subject
%            s_diff_alltrials_subj = BAV_data{i}(:,2)- BAV_data{i}(:,3);
%             binedges = quantile(s_diff_alltrials_subj(abs(s_diff_alltrials_subj)>0.01), linspace(0,1,num_binedges));
%             bincenters = binedges(1:(end-1)) + (binedges(2:end) - binedges(1:(end-1))) ./2;

%           subplot(3,5,i)
            hold on
            cond0 = ((BAV_data{i}(:,5) == response_type)==1);
            COND0 = ((BAV_data_alltrials(:,5) == response_type)==1);
            cond1 = ((BAV_data{i}(:,1) == l)==1);
            COND1 = ((BAV_data_alltrials(:,1) == l)==1);
            if(strats==1)
                cond2 = ((BAV_data{i}(:,2)+BAV_data{i}(:,3))<stratify_thres(1));
                COND2 = ((BAV_data_alltrials(:,2)+BAV_data_alltrials(:,3))<stratify_thres(1));
                cond3 = 1;
                COND3 = 1;
                trials_condsubj_idx = find(cond0 .* cond1 .* cond2.*cond3~=0);
            elseif(strats==2)
                cond2 = ((BAV_data{i}(:,2)+BAV_data{i}(:,3))>=stratify_thres(1));
                COND2 = ((BAV_data_alltrials(:,2)+BAV_data_alltrials(:,3))>=stratify_thres(1));
                cond3 = ((BAV_data{i}(:,2)+BAV_data{i}(:,3))<stratify_thres(2));
                COND3 = ((BAV_data_alltrials(:,2)+BAV_data_alltrials(:,3))<stratify_thres(2));
                trials_condsubj_idx = find(cond0 .* cond1 .* cond2 .* cond3 ~=0);
                %trials_condsubj_idx = find(((BAV_data{i}(:,1) == l)==1 .* ((BAV_data{i}(:,2)+BAV_data{i}(:,3))>=stratify_thres(1)) .* ((BAV_data{i}(:,2)+BAV_data{i}(:,3))<stratify_thres(2)))~=0);
            else  
                cond2 = ((BAV_data{i}(:,2)+BAV_data{i}(:,3))>stratify_thres(2));
                COND2 = ((BAV_data_alltrials(:,2)+BAV_data_alltrials(:,3))>stratify_thres(2));
                cond3 = 1;
                COND3=1;
                trials_condsubj_idx = find(cond0 .* cond1 .* cond2.*cond3 ~=0);
                %trials_condsubj_idx = find(((BAV_data{i}(:,1) == l)==1 .* ((BAV_data{i}(:,2)+BAV_data{i}(:,3))>=stratify_thres(2)))~=0);
            end
            trials_consubj = BAV_data{i}(trials_condsubj_idx,:);
            
            
            
% ---------------------------------------------------------------------
            if(unique_bins_subj)
                % Option 1: Each subject has its own binedges.
                trials_binedges_idx = find(cond0.* cond1 .* cond2 .* cond3 ~=0);
                trials_binedges = BAV_data{i}(trials_binedges_idx,:);
            else
                % Option 2: All subjects share the same binedges. 
                TRIALS_binedges_idx = find(COND0.* COND1 .* COND2 .* COND3 ~=0);
                trials_binedges = BAV_data_alltrials(TRIALS_binedges_idx,:);
            end
% ---------------------------------------------------------------------
            
            s_diff_alltrials_subj = trials_binedges(:,2)- trials_binedges(:,3);
            binedges = quantile(s_diff_alltrials_subj(abs(s_diff_alltrials_subj)>0.01), linspace(0,1,num_binedges));
            % Prevent percentiles from exactly landing on one of the s_diff
            % vals, causing difficulty in binning (double count, never
            % count)
            for b = 1:num_binedges
                if(~isempty(find(s_diff_alltrials_subj==binedges(b))==1))
                    if(b~=1)
                         binedges(b) = binedges(b)+0.001;
                    else
                        binedges(b) = binedges(b)-0.001;
                    end
                end
            end
                    
            
            bincenters = binedges(1:(end-1)) + (binedges(2:end) - binedges(1:(end-1))) ./2;
 
            % For model ribbon visualization
            S_A_condsubj = trials_consubj(:,2);
            S_V_condsubj = trials_consubj(:,[1,3]);
            R_condsubj = trials_consubj(:,[4,5]); 
            
            
% --------------For model fit visualization -- omitted for now!-----------
            num_samps = 100;
            shats_given_s_modelfit = zeros(length(trials_consubj(:,4)),num_samps);
            for samp=1:num_samps
                switch model_family
                    case "parametric"
                        shats_given_s_modelfit(:, samp) = nllfun_bav_parametric(ModelComponents_V, ModelComponents_A, fitted_params_PM(i,:), R_condsubj, S_V_condsubj, S_A_condsubj, PMIntegrationParams, return_predictive_samples, return_response_distr, plot_lapse, dx_max, causal_inf_strategy, lapse_type, Gaussian_lapse_SDs(i));
                    case "semiparamInsp"	
                        ModelComponents.PriorUnnormalized = ModelComponents.PriorUnnormalizedAllSubjs(:,i);
                        ModelComponents.SigmaFuns = ModelComponents.SigmaFunsAllSubjs{i};
                        shats_given_s_modelfit(:, samp) = nllfun_bav_ubresc_semiparaminsp(fitted_params_PM(i,:), trials_consubj, ModelComponents, return_predictive_samples, return_response_distr, plot_lapse, lapse_type, Gaussian_lapse_SDs(i));
                end
            end

            mean_biases_modelfit = zeros(1,length(binedges)-1);
            
            % For human data visualization
            s_diffs = trials_consubj(:,2)- trials_consubj(:,3);
            mean_biases = zeros(1,length(binedges)-1);
            
            for k=1:(num_binedges-1) % 5; Need change back to 1:(num_binedges-1)
                [response_type,strats, l,i,k]
                % Find relevant trials with (s_A-s_V) disparity inside this bin.
                trials_condsubjbin_idx = find((int8(s_diffs(:)>=binedges(k)) .* int8(s_diffs(:)<binedges(k+1)))~=0);
                trials_condsubjbin = trials_consubj(trials_condsubjbin_idx,:);
                
                % Bias = response - stim (visual stim for visual resp, vice versa);
                stim_response_bias = trials_condsubjbin(:,4) - trials_condsubjbin(:,3-(response_type==2));
                
                % For model fits     
                shats_given_s_modelfit_rel = shats_given_s_modelfit(trials_condsubjbin_idx,:);
                stim_response_bias_modelfits = shats_given_s_modelfit_rel- trials_condsubjbin(:,3-(response_type==2));
                                
                mean_biases(k) = mean(stim_response_bias);
                mean_biases_modelfit(k) = mean(stim_response_bias_modelfits(:));
            end
            bincenters_stratrels(l,strats,i,:) = bincenters;
            %Prob_C1(isnan(Prob_C1))=0;
            Mean_Biases(l,strats,i,:) = mean_biases;
            Mean_Biases_modelfit(l,strats,i,:) = mean_biases_modelfit;
    
        end
    end 
    end
        
    %%
        if(plot_individual)
            figure('Position', figspecs);
            % figsize = get(0, 'ScreenSize');
            % figure('Position', figsize);
            
            set(gcf, 'Color', 'w')
            T=tiledlayout(4,4,'Padding', 'none', 'TileSpacing', 'tight');
            rel_level_alphas = [1, 0.5, 0.25];
            colors_AV = brewermap(10,"Set1");
            colors_AV = colors_AV([5,4,7],:); 
            for i=1:num_subjects
                t=tiledlayout(T,1,3, 'Padding', 'none', 'TileSpacing', 'tight');
                t.Layout.Tile = i;
                t.Layout.TileSpan = [1 1];
                for strats=1:3
                    color_AV = colors_AV(strats,:);
                    tt=nexttile(t);
                    set(gca,'TickDir','out');
                    hold on
                    title(t,"Subject " + i, 'FontSize',fontsize, 'FontWeight', 'bold')
                    for l=1:3
                        alpha_level = rel_level_alphas(l);
                        scatter(squeeze(bincenters_stratrels(l,strats,i,:)), squeeze(Mean_Biases(l,strats,i,:)),10,'MarkerFaceColor',color_AV,'MarkerEdgeColor','none','MarkerFaceAlpha',alpha_level,'MarkerEdgeAlpha',alpha_level);
                        patchline(squeeze(bincenters_stratrels(l,strats,i,:)), squeeze(Mean_Biases_modelfit(l,strats,i,:)), 'edgecolor',color_AV,'linewidth',1,'edgealpha',alpha_level);
%                         p = plot(squeeze(bincenters_stratrels(l,strats,i,:)), squeeze(Mean_Biases_modelfit(l,strats,i,:)), 'LineStyle', '-', 'Color',color_AV,'LineWidth',1);
%                         p.Color(4) = alpha_level;
                    end
                    xlim([-35,35])
                    ylim([-15,15])
                    if((mod(i,4)==1) && (strats==1))
                        yticks(-15:5:15)
                    elseif(strats==1)
                        set(tt,'YTickMode', 'manual', 'YTick',-15:5:15,'Yticklabel',[]);
                    else
                        yticks([])
                    end
                    xticks([-20,0,20]);
                    tt.TickLength(1)=0.03;
                    xtickangle(0)
                    tt.XAxis.FontSize = 9;
                    tt.YAxis.FontSize = 9;
                end
            end
            
            % Create legend
            idx=num_subjects+1;
            i = min(idx, num_subjects);
            nexttile(T);
            set(gca,'TickDir','out');
            hold on
            for strats=1:3
                color_AV = colors_AV(strats,:);
                for l=1:3
                    alpha_level = rel_level_alphas(l);
                    scatter(squeeze(bincenters_stratrels(l,strats,i,:)), squeeze(Mean_Biases(l,strats,i,:)),10,'MarkerFaceColor',color_AV,'MarkerEdgeColor',color_AV,'MarkerFaceAlpha',alpha_level,'MarkerEdgeAlpha',alpha_level);
                end
            end
            axis off
            fill([-40,-40,40,40],[-40,40,40,-40],"w", 'EdgeColor','none','HandleVisibility','off')
            lg = legend("Left, High visual rel","Left, Med visual rel","Left, Low visual rel","Center, High visual rel","Center, Med visual rel","Center, Low visual rel", "Right, High visual rel", "Right, Med visual rel", "Right, Low visual rel");
            set(lg,'Box','off')
            lg.FontSize = 9;
            lg.ItemTokenSize(1) = 10;
            lg.Position(2) = 0.05;
            xlim([-35,35])
            ylim([-15,15])
            yticks(-15:5:15)
            
            
            xlabel(T, "Stimulus location disparity, {\its}_A– {\its}_V (\circ)", 'FontSize',fontsize)
            switch response_type
                case 1
                    ylabel(T, "{\rm \fontsize{9} {Visual bias (\circ)}}");
                case 2
                    ylabel(T, "{\rm \fontsize{9} {Auditory bias (\circ)}}");
            end
                

        else % Default plots -- mean+-SEM across subjects
            figure('Position', figspecs);
            set(gcf, 'Color', 'w')
            % bottom_tick = 0.02;
            % top_tick = 0.98;
            % bracket_y_coord = 0.023;
            % annotation("line",[bracket_y_coord,bracket_y_coord],[bottom_tick,top_tick], 'LineWidth',1)
            % annotation("line",[bracket_y_coord,0.047],[bottom_tick,bottom_tick], 'LineWidth',1)
            % annotation("line",[bracket_y_coord,0.047],[top_tick,top_tick], 'LineWidth',1)

            tcl=tiledlayout(3,3,'Padding', 'none', 'TileSpacing', 'compact');
            for l=1:3
                for strats=1:3
                    nexttile
                    set(gca,'TickDir','out');
                    % Plot average plots
                    bincenters_stratrels_avg = squeeze(mean(bincenters_stratrels, 3));
                    bincenters_stratrels_std = squeeze(std(bincenters_stratrels, [], 3));

                    Mean_Biases_avg = squeeze(nanmean(Mean_Biases, 3));
                    Mean_Biases_sem = squeeze(nanstd(Mean_Biases, [], 3)) ./ sqrt(squeeze(sum(~isnan(Mean_Biases), 3)));   
                    Mean_Biases_modelfit_avg = squeeze(nanmean(Mean_Biases_modelfit, 3));
                    Mean_Biases_modelfit_sem = squeeze(nanstd(Mean_Biases_modelfit, [], 3)) ./ sqrt(squeeze(sum(~isnan(Mean_Biases_modelfit), 3)));   

                    hold on
                    % Plot human data
                    errorbar(squeeze(bincenters_stratrels_avg(l,strats,:)), squeeze(Mean_Biases_avg(l,strats,:)), squeeze(Mean_Biases_sem(l,strats,:)), "k.", 'CapSize', 3)

                    % Plot model ribbons
                    curve1 = Mean_Biases_modelfit_avg(l,strats,:) + Mean_Biases_modelfit_sem(l,strats,:);
                    curve2 = Mean_Biases_modelfit_avg(l,strats,:) - Mean_Biases_modelfit_sem(l,strats,:);
                    curve1(isnan(curve1))=0; curve2(isnan(curve2))=0;
                    x2 = [squeeze(bincenters_stratrels_avg(l,strats,:))', fliplr(squeeze(bincenters_stratrels_avg(l,strats,:))')];
                    inBetween = [squeeze(curve1)', fliplr(squeeze(curve2)')];
                    p=fill(x2, inBetween,"k",'FaceAlpha',0.2, 'EdgeColor', 'none');
                    p.Annotation.LegendInformation.IconDisplayStyle = 'off';

                    if(l==1)
                        title(stratify_labels(strats), 'FontSize', 10)
                    elseif(l==3)
                        xlabel("Stimulus location disparity, {\its}_A"+char(8722)+"{\its}_V (\circ)", 'FontSize',fontsize)
                    end
                    sensory_bias_type = ["Visual","Auditory"];
                    if(strats==1)
                        ylabel({"{\fontsize{10} \bf{"+reliability_titles(l)+" Visual Reliability}}","{\rm \fontsize{9} {"+sensory_bias_type(response_type)+" bias (\circ)}}"})
                        % if(l==2)
                        %     ylabel({"{\fontsize{10} \bf{Visual}}", "{\fontsize{9} \bf{"+reliability_titles(l)+" Reliability}}","{\rm \fontsize{9} {Mean "+sensory_bias_type(response_type)+" bias}}"})
                        % else
                        %     ylabel({"","{\fontsize{9} \bf{"+reliability_titles(l)+" Reliability}}","{\rm \fontsize{9} {Mean "+sensory_bias_type(response_type)+" bias}}"})
                        % end
                    end
                    xlim([-35,35])
                    ylim([-15,15])
                    xticks(-30:10:30)
                end
            end

        end

    %%
    
    switch response_type
    case 1
        if(plot_individual)
            saveas(gca, figpath+save_name+'-BV_Individual.fig')
            exportgraphics(gcf,figpath+save_name+'-BV_Individual.png','Resolution',png_dpi);
            exportgraphics(gcf,figpath+save_name+'-BV_Individual.pdf',"ContentType","vector");
        else
            exportgraphics(gcf,figpath+save_name+'-BV.png','Resolution',png_dpi);
            exportgraphics(gcf,figpath+save_name+'-BV.pdf',"ContentType","vector");
        end
    case 2
        if(plot_individual)
            saveas(gca, figpath+save_name+'-BA_Individual.fig')
            exportgraphics(gcf,figpath+save_name+'-BA_Individual.png','Resolution',png_dpi);
            exportgraphics(gcf,figpath+save_name+'-BA_Individual.pdf',"ContentType","vector");
        else
            exportgraphics(gcf,figpath+save_name+'-BA.png','Resolution',png_dpi);
            exportgraphics(gcf,figpath+save_name+'-BA.pdf',"ContentType","vector");
        end
    end
    
    end
end
