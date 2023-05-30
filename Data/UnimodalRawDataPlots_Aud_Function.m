% Shuze Liu
% This code visualizes the audio-visual effect and heterskedasticity. 
close all; clear all;
cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Data')
load("alldata.mat");
cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Plots')

%% Calls function
%names = {'aaron' 'allison' 'avauna' 'clare' 'connie' 'daisy' 'david' 'gabriela' 'jeff' 'jobyeon' 'julia' 'katherine' 'priyanka' 'stephanie' 'wendy'} ;
num_subjects = length(data);
subj_idx = cell(1,15);
for i=1:num_subjects
    subj_idx{i} = "subj "+num2str(i);
end

colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840 ];

s_a_range = -15:5:15;
s_v_range = -20:1:20;

is_rescale = false; rescale = 4/3;
is_bimodal=false; is_save=false; 

% No longer need this -- negative pixel value issue solved!
discard_outliers=false; discard_thres = 50;

[data_stratified, s_cond_s_hat_means, s_cond_s_hat_stds] = Aud_RawDataPlots(data, colors, subj_idx, s_a_range, is_bimodal, is_save, discard_outliers, discard_thres, is_rescale, rescale);
cd('C:\Users\liu_s\OneDrive\桌面\MATLAB\AudioVisual\Analysis')
save('data_stratified_UA', 'data_stratified')
%%
function [data_stratified, s_cond_s_hat_means, s_cond_s_hat_stds] = Aud_RawDataPlots(data, colors, subj_idx, s_range, is_bimodal, is_save, discard_outliers, discard_thres, is_rescale, rescale)
    % For bimodal, this does not separate according to 3 visual coherence
    % levels!
    s_idx=2;
    fns = fieldnames(data{1});
    num_subjects = length(subj_idx);
    c = [0 0.4470 0.7410];
    
    % Create cell to store data for each subject -> each s_level. 
    % First level: 15 subjects.
    % Second level: trials with different s_values. 
    % If is_outlier/is_rescale, will resemble data after these changes. 
    % Discarding outliers always occur before rescaling!
    data_stratified = cell(1,15);
    
    % Specify which of the 5 tasks to consider. Get relevant column of s_A (Col 2)
    save_suffix = "_";
    if(is_bimodal)
        save_suffix = save_suffix + "BA";  
        fns_idx=2;
    else
        save_suffix = save_suffix + "UA";
        fns_idx=4;
    end
    if(is_rescale)
        save_suffix = save_suffix + "rescale";  
    end
    if(discard_outliers)
        save_suffix = save_suffix + "-noOutliers";
    end
    
    %% Make scatterplot (Fig 1) and response histograms (Fig 2)
    s_range = -15:5:15;
    s_cond_s_hat_means = zeros(num_subjects,length(s_range));
    s_cond_s_hat_stds = zeros(num_subjects,length(s_range));
    s_cond_s_hat_iqrs = zeros(num_subjects,length(s_range));
    
    hist_range = -30:3:30;
    binedges = [-Inf, hist_range(1:(end-1)), +Inf] + 0.5;

    Fig1=figure('Position', get(0, 'Screensize'));
    Fig2=figure('Position', get(0, 'Screensize'));
    for i=1:num_subjects
        % Common to Fig 1 and 2
        dataA_subject = data{i}.(fns{fns_idx});
        s_vals = dataA_subject(:,s_idx);
        s_hat_vals = dataA_subject(:,end)
        
        data_stratified{i} = cell(1, length(s_range));
        if(discard_outliers) % Discard outlier data when needed.
            is_not_outlier = abs(s_hat_vals) < discard_thres;
            s_vals = s_vals(is_not_outlier);
            s_hat_vals = s_hat_vals(is_not_outlier);
        end
        if(is_rescale)
            s_hat_vals = s_hat_vals .* rescale;
            rescale
        end
            
        
        % Fig 1 stuff in this iteration
        figure(1)
        subplot(3,5,i)
        hold on
        plot(s_vals, s_hat_vals, ".", 'Color',c)
        plot([-40,40], [-40,40], "k--")
        xlabel("s")
        ylabel("s hat")
        title(subj_idx{i})
        
        % Fig 2 stuff in this iteration
        figure(2)
        subplot(3,5,i)
        hold on
        for j=1:length(s_range)
            relevant_trials = s_vals==s_range(j);
            s_hats_relevant = s_hat_vals(relevant_trials);
            data_stratified{i}{j} = [s_vals(relevant_trials), s_hats_relevant];
            [N,edges] = histcounts(s_hats_relevant, binedges, 'Normalization','pdf');
            N = N ./ sum(N);
            %edges
            %hist_range = edges(2:end) - (edges(2)-edges(1))/2;
            plot(hist_range, N, ".-", 'Color', colors(j,:));
            s_cond_s_hat_means(i,j) = mean(s_hats_relevant);
            s_cond_s_hat_stds(i,j) = std(s_hats_relevant);
            s_cond_s_hat_iqrs(i,j) = iqr(s_hats_relevant);
        end
        plot([s_range(1),s_range(1)],[0,0.2], "k--")
        plot([s_range(end),s_range(end)],[0,0.2], "k--")
        plot([-20, -20],[0,0.2], "k--")
        plot([20,20],[0,0.2], "k--")
        xlim([-40,40])
        xlabel("s hat")
        ylabel("Relative Density")
        if(i==num_subjects)
            legend({"s=-15", "s=-10","s=-5","s=0","s=5","s=10", "s=15"})
        end
        title(subj_idx{i})
    end

    %% E[s_hat | s] for all subjects.
    figure;
    for i=1:length(s_range)
        subplot(2,4,i);
        histogram(s_cond_s_hat_means(:,i),5,'Normalization','probability')
        xlabel("E[s hat | s]")
        ylabel("Rel Freq")
        title({"s = ", num2str(s_range(i))})
    end

    % Pr[s_hat | s] means for each subject.
    Fig3=figure('Position', get(0, 'Screensize'));
    colormap(parula(num_subjects))
    for i=1:num_subjects
        plot(s_range, s_cond_s_hat_means(i,:),".-")
        hold on
    end
    xlabel("s")
    ylabel("E[s hat | s]")
    plot(s_range, mean(s_cond_s_hat_means,1), "k:", 'LineWidth',5)
    legend(subj_idx)
    h1=plot([-20,20], [-20,20], "k--");
    h2=plot([-20,20], [s_range(1),s_range(1)], "k--");
    h3=plot([-20,20], [s_range(end),s_range(end)], "k--");
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    %% SD[s_hat | s] for all subjects.
    Fig4=figure('Position', get(0, 'Screensize'));
    for i=1:num_subjects
        subplot(1,2,1)
        plot(s_range, s_cond_s_hat_stds(i,:),".-")
        hold on
        subplot(1,2,2)
        plot(s_range, s_cond_s_hat_iqrs(i,:),".-")
        hold on
    end
    subplot(1,2,1)
    colormap(jet(num_subjects))
    plot(s_range, mean(s_cond_s_hat_stds,1), "k:", 'LineWidth',5)
    legend(subj_idx)
    xlabel("s")
    ylabel("SD[s hat | s]")
    title("SD[s hat | s]")
    xlim([-20,20])
    ylim([0,10])
    
    subplot(1,2,2)
    colormap(jet(num_subjects))
    plot(s_range, mean(s_cond_s_hat_iqrs,1), "k:", 'LineWidth',5)
    legend(subj_idx)
    xlabel("s")
    ylabel("IQR[s hat | s]")
    title("IQR[s hat | s]")
    xlim([-20,20])
    %ylim([0,10])
    
    % Save when needed.
    if(is_save)
        saveas(Fig1,'s_hat_VS_s_scatter'+save_suffix+'.png')
        saveas(Fig2,'ResponseDistrs_VS_s_hat'+save_suffix+'.png')
        saveas(Fig3,'E_ResponseDistr_VS_s'+save_suffix+'.png')
        saveas(Fig4,'SD_ResponseDistr_VS_s'+save_suffix+'.png')
    end
end





