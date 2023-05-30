%% Extrace Bimondal location estimate (BAV) task data 
% load("alldata.mat")
% num_subjects = length(data);
% BAV_data = cell(1,num_subjects);
% 
% % Cols: [Vis_Rel, sV, sA, s\hat, s\hat_is_visual]
% for i=1:num_subjects
%     BV_subj = [data{i}.dataV, repmat(1,length(data{i}.dataV),1)];
%     BA_subj = [data{i}.dataA, repmat(2,length(data{i}.dataA),1)];
%     BAV_data{i} = [BV_subj; BA_subj];
% end
% save('BAV_data.mat', 'BAV_data');
% 
% % Col 1: visual stimulus reliability (most reliable if value is 1)
% % Col 2: auditory stimulus true location
% % Col 3: visual stimulus true location
% % Col 4: subject's stimulus location estimate
% % Col 5: whether Col is visual (1) or auditory (2) location estimate.

%%
close all; clear all;
load("BAV_data.mat")
num_subjects = length(BAV_data);

s_a_range = -15:5:15;
s_v_range_binedges = linspace(-20,20,8);
s_v_range = (s_v_range_binedges(2:end) - s_v_range_binedges(1:(end-1)))./2 + s_v_range_binedges(1:(end-1));
num_sA_bins = length(s_a_range);
num_sV_bins = length(s_v_range);

s_hat_binedges = -45:5:45;
s_hat_bincenters = (s_hat_binedges(2:end) - s_hat_binedges(1:(end-1)))./2 + s_hat_binedges(1:(end-1));
num_s_hat_bins = length(s_hat_bincenters);

colors = brewermap(length(s_v_range)+1,"Set1");
colors = colors([1,2,3,4,5,7,8],:);
response_types_names = ["BV", "BA"];
%%
responses_bysubj_stats = zeros(2, num_sV_bins, num_sA_bins, 3, num_subjects, num_s_hat_bins);
for response_type=1:2 % Separate figures for BV/BA trials
    figure;
    tiledlayout(num_sV_bins, num_sA_bins, 'TileSpacing', 'compact')
    
    for v_binidx = 1:num_sV_bins % Subplot rows are sV bins
        for a_binidx = 1:num_sA_bins % Subplot cols are sA bins
            
            nexttile
            [response_type, v_binidx, a_binidx]
            hold on
            
            
            for l=1:3
                for subjidx=1:num_subjects
                    subj_trials = BAV_data{subjidx};
                
                    relevant_subjrelbin_trials_idx = (subj_trials(:,5)==response_type) .* (subj_trials(:,1)==l) .* (abs(subj_trials(:,2)-s_a_range(a_binidx))<=0.1);
                    relevant_subjrelbin_trials_idx = relevant_subjrelbin_trials_idx .* (subj_trials(:,3)>=s_v_range_binedges(v_binidx)) .* (subj_trials(:,3)<s_v_range_binedges(v_binidx+1));
                    relevant_subjrelbin_trials_idx = find(relevant_subjrelbin_trials_idx==1);
                    
                    relevant_trials = subj_trials(relevant_subjrelbin_trials_idx,:);
                    
                    [s_hat_histcounts,~] = histcounts(relevant_trials(:,4),s_hat_binedges);

                    if(sum(s_hat_histcounts)~=0)
                        subj_vals = s_hat_histcounts ./ sum(s_hat_histcounts);
                    else
                        subj_vals = s_hat_histcounts;
                    end
                    responses_bysubj_stats(response_type, v_binidx, a_binidx, l, subjidx,:) = subj_vals;
                    %subj_vals(subj_vals==0)=NaN;
                    plot(s_hat_bincenters, subj_vals, ".-", 'Color', [colors(l,:),0.1])
                end
            end
            xlim([-35,35])
            if(v_binidx==7)
                xlabel("Location Estimate")
            end
            ylim([0,1])
            if(a_binidx==1)
                ylabel("Rel Freq")
            end
            title("sV binidx "+ num2str(v_binidx)+ ", sA="+ num2str(s_a_range(a_binidx)))
        end
    end
    sgtitle("Human subjects: "+response_types_names(response_type)+" trials")  
    
end

%% Plot raw data
mean_responses_oversubj_stats = mean(responses_bysubj_stats, 5);
sem_responses_oversubj_stats = std(responses_bysubj_stats,[],5) ./ sqrt(num_subjects);
     
for response_type=1:2 % Separate figures for BV/BA trials
    figure;
    tiledlayout(num_sV_bins, num_sA_bins, 'TileSpacing', 'compact')
    for v_binidx = 1:num_sV_bins % Subplot rows are sV bins
        for a_binidx = 1:num_sA_bins % Subplot cols are sA bins
            [response_type, v_binidx, a_binidx]
            nexttile
            hold on
            for l=1:3
                % Vectors over s_hat_grid. 
                means = squeeze(mean_responses_oversubj_stats(response_type, v_binidx, a_binidx, l, :, :))';
                sems = squeeze(sem_responses_oversubj_stats(response_type, v_binidx, a_binidx, l, :, :))';
                
                % Remove all errorbars corresponding to no data.
                nanidx = find((means~=0)+(sems~=0)==0);
%                 if(~isempty(nanidx))
%                     means(nanidx)=NaN;
%                 end
                errorbar(s_hat_bincenters, means, sems, ".-", 'Color', colors(l,:));
            end
            xlim([-35,35])
            if(v_binidx==7)
                xlabel("Location Estimate")
            end
            ylim([0,1])
            if(a_binidx==1)
                ylabel("Rel Freq")
            end
            title("$s_V \in ["+ num2str(round(s_v_range_binedges(v_binidx),1))+", "+ num2str(round(s_v_range_binedges(v_binidx+1),1))+ "), s_A="+ num2str(s_a_range(a_binidx)) + "$", 'interpreter', 'latex')
        end
    end
    sgtitle("Human subjects agg: "+response_types_names(response_type)+" trials")
        
end
       
                
                
                
