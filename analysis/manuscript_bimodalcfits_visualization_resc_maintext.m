function manuscript_bimodalcfits_visualization_resc_maintext(T,BC_data, fitted_params_PM, ModelComponents_V, ModelComponents_A, model_family, plot_lapse, causal_inf_strategy, fontsize, figspecs, plot_individual)
    % Note: if model_family=="semiparamInsp", then ModelComponents_A,
    % PMIntegrationParams, dx_max, causal_inf_strategy arguments are useless.

    return_predictive_samples = false;
    return_response_distr = true;
    num_rels = 3;
    reliability_titles=["High","Med.","Low"];
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
    
    num_subjects = length(BC_data);
    R_round_decimals = 2; % R round to nearest integer for binning.
    dr = 10^(-R_round_decimals);
    R_grid = -45:dr:45;
    
    subj_idx = cell(1,num_subjects);
    s_diff_alltrials = [];
    s_sum_alltrials = [];
    for i=1:num_subjects
        subj_idx{i} = "subj "+num2str(i);
        s_diff_alltrials = [s_diff_alltrials; BC_data{i}(:,2)- BC_data{i}(:,3)];
        s_sum_alltrials = [s_sum_alltrials; BC_data{i}(:,2)+ BC_data{i}(:,3)];
    end
    binedges = quantile(s_diff_alltrials(abs(s_diff_alltrials)>0.01), linspace(0,1,10));
    bincenters = binedges(1:(end-1)) + (binedges(2:end) - binedges(1:(end-1))) ./2;
    Probs_C1s = zeros(3, 2, num_subjects, length(binedges)-1);
    Probs_C1s_modelfit = Probs_C1s;
  
    stratify_labels = ["|{\its}_A+{\its}_V| center", "|{\its}_A+{\its}_V| periphery"];
    stratify_thres = quantile(abs(s_sum_alltrials), 0.5);
    for l=1:3
        for strats=1:2
            for i=1:num_subjects

                if(strats==1)
                    trials_condsubj_idx = find((BC_data{i}(:,1) == l)==1 .* (abs(BC_data{i}(:,2)+BC_data{i}(:,3))<stratify_thres));
                else
                    trials_condsubj_idx = find((BC_data{i}(:,1) == l)==1 .* (abs(BC_data{i}(:,2)+BC_data{i}(:,3))>=stratify_thres));
                end
                trials_consubj = BC_data{i}(trials_condsubj_idx,:);

                % For model ribbon visualization
                S_A_condsubj = trials_consubj(:,2);
                S_V_condsubj = trials_consubj(:,[1,3]);
                R_condsubj = repmat(1,length(S_A_condsubj),1); 

                switch model_family   
                    case "parametric"
                        Prob_C1_given_s_modelfit = nllfun_bc_parametric(ModelComponents_V, ModelComponents_A, fitted_params_PM(i,:), R_condsubj, S_V_condsubj, S_A_condsubj, return_predictive_samples, return_response_distr, plot_lapse, causal_inf_strategy);
                    case "semiparamInsp"
                        ModelComponents.PriorUnnormalized = ModelComponents.PriorUnnormalizedAllSubjs(:,i);
                        ModelComponents.SigmaFuns = ModelComponents.SigmaFunsAllSubjs{i};
                        data = trials_consubj;
                        data(:,4) = R_condsubj;
                        Prob_C1_given_s_modelfit = nllfun_bc_ubresc_semiparaminsp(fitted_params_PM(i,:), data, ModelComponents, return_predictive_samples, return_response_distr, plot_lapse);
                end

                % For human data visualization
                s_diffs = trials_consubj(:,2)- trials_consubj(:,3);

                Prob_C1 = zeros(1,length(binedges)-1);
                Prob_C1_modelfit = zeros(1,length(binedges)-1);

                for k=1:(length(binedges)-1)
                    [l,strats,i,k]
                    % Find relevant trials with (s_A-s_V) disparity inside this bin.
                    trials_condsubjbin_idx = find((int8(s_diffs(:)>=binedges(k)) .* int8(s_diffs(:)<binedges(k+1)))~=0);
                    trials_condsubjbin = trials_consubj(trials_condsubjbin_idx,:);

                    % for model fits        
                    p_resp_modelfits_condsubjbin = Prob_C1_given_s_modelfit(trials_condsubjbin_idx);
                    %model_predictions = 1+binornd(1,1-relevant_modelpreds,[length(sample_counts),1]);

                    Prob_C1(k) = mean(trials_condsubjbin(:,4)==1);
                    Prob_C1_modelfit(k) = mean(p_resp_modelfits_condsubjbin);
                end

                Probs_C1s(l,strats,i,:) = Prob_C1;
                Probs_C1s_modelfit(l,strats,i,:) = Prob_C1_modelfit;
            end
        end
    end


        % Plot average plots
        rel_colors = brewermap(9,"Set1");
        rel_colors = rel_colors([1,2,3],:);
        Probs_C1s_avg = squeeze(mean(Probs_C1s, 3));
        Probs_C1s_sem = squeeze(std(Probs_C1s, [], 3)) ./ sqrt(num_subjects);

        Probs_C1s_modelfit_avg = squeeze(mean(Probs_C1s_modelfit, 3));
        Probs_C1s_modelfit_sem = squeeze(std(Probs_C1s_modelfit, [], 3)) ./ sqrt(num_subjects);

        
        % Default plots -- mean+-SEM across subjects
        for l=1:3
            t=tiledlayout(T,1,2, 'Padding', 'none', 'TileSpacing', 'tight');
            t.Layout.Tile = 1+(l-1)*10;
            t.Layout.TileSpan = [1 4];
            for strats=1:2
                %nexttile((strats*2-1)+((l-1)*11), [1,2]);
                tt=nexttile(t, strats);

                set(gca,'TickDir','out');
                hold on
                % Plot human data
                errorbar(bincenters, squeeze(Probs_C1s_avg(l,strats,:)), squeeze(Probs_C1s_sem(l,strats,:)), "k.", 'CapSize', 3)

                % Plot model ribbons
                curve1 = Probs_C1s_modelfit_avg(l,strats,:) + Probs_C1s_modelfit_sem(l,strats,:);
                curve2 = Probs_C1s_modelfit_avg(l,strats,:) - Probs_C1s_modelfit_sem(l,strats,:);
                curve1(isnan(curve1))=0; curve2(isnan(curve2))=0;
                x2 = [bincenters, fliplr(bincenters)];
                inBetween = [squeeze(curve1)', fliplr(squeeze(curve2)')];
                p=fill(x2, inBetween,"k",'FaceAlpha',0.2, 'EdgeColor', 'none');
                p.Annotation.LegendInformation.IconDisplayStyle = 'off';

                if(l==1)
                    title(stratify_labels(strats), 'FontSize',10)
                    if(strats==1)
                        text(-0.3, 1.15, '(a)', 'FontWeight', 'bold', 'FontSize', fontsize+1, 'HorizontalAlignment', 'left', 'Units', 'normalized');
                    end
                end
                yticks(0:0.2:1)
                if(strats==1)
                    ylabel({"{\fontsize{10} \bf{"+reliability_titles(l)+" Visual Reliability}}","{\rm \fontsize{9} {Proportion responding "+ '"'+'same'+ '"'+"}}"})
                else
                    yticklabels([])
                end
                xlim([-30,30])
                ylim([0,1])
            end
        end

end