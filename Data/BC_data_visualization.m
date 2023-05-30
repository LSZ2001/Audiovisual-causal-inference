%% Extrace Bimondal category (BC) choice task data 
% load("alldata.mat")
% num_subjects = length(data);
% BC_data = cell(1,num_subjects);
% for i=1:num_subjects
%     BC_data{i} = data{i}.dataC;
% end
%save('BC_data.mat', 'BC_data');

% Col 1: visual stimulus reliability (most reliable if value is 1)
% Col 2: auditory stimulus true location
% Col 3: visual stimulus true location
% Col 4: subject's unity judgment -- 1 is same source, 2 is different sources.

%% Visualize results
close all; clear all;
load("BC_data.mat")
num_subjects = length(BC_data);

reliability_names = ["High", "Med", "Low"];
rel_colors = brewermap(11,"BrBG");
rel_colors = rel_colors([1,3,4],:);

%% Bin the data and make similar plots
binedges = linspace(-20-15, 20+15, 31); % binedges for grouping s_diffs.
bincenters = binedges(1:(end-1)) + (binedges(2:end) - binedges(1:(end-1))) ./2;
PMR_vals = zeros(3, num_subjects, 2, length(binedges)-1); % Rels -> subjects -> C=1 or C=2 -> bins.
Probs_C1s = zeros(3, num_subjects, length(binedges)-1);
for l=1:3
    figure
    for i=1:num_subjects
        subplot(3,5,i)
        hold on
        relevant_trial_idx = find((BC_data{i}(:,1) == l)==1);
        relevant_trials = BC_data{i}(relevant_trial_idx,:);
        
        C1_relevant_trials = relevant_trials(relevant_trials(:,4)==1,:);
        C2_relevant_trials = relevant_trials(relevant_trials(:,4)==2,:);
        s_diffs_C1 = C1_relevant_trials(:,2)-C1_relevant_trials(:,3);
        s_diffs_C2 = C2_relevant_trials(:,2)-C2_relevant_trials(:,3);
        
        PMRs_C1 = [];
        PMRs_C2 = [];
        for k=1:(length(binedges)-1)
            [l,i,k]
            sample_counts_C1 = find((int8(s_diffs_C1(:)>=binedges(k)) .* int8(s_diffs_C1(:)<binedges(k+1)))~=0);
            sample_counts_C2 = find((int8(s_diffs_C2(:)>=binedges(k)) .* int8(s_diffs_C2(:)<binedges(k+1)))~=0);
            PMRs_C1 = [PMRs_C1, sum(sample_counts_C1)];
            PMRs_C2 = [PMRs_C2, sum(sample_counts_C2)];
        end
        PMR_vals(l,i,1,:) = PMRs_C1 ./ sum(PMRs_C1);
        PMR_vals(l,i,2,:) = PMRs_C2 ./ sum(PMRs_C2);
        Prob_C1s(l,i,:) = squeeze(PMR_vals(l,i,1,:) ./ (PMR_vals(l,i,1,:)+PMR_vals(l,i,2,:)));
      
%         plot(bincenters, squeeze(PMR_vals(l,i,1,:)), ".-")
%         plot(bincenters, squeeze(PMR_vals(l,i,2,:)), ".-")
        plot(bincenters, squeeze(Prob_C1s(l,i,:)), "k.")

        xlabel("s_A - s_V")
        ylabel("Pr[C hat = 1]")
        title("Subject "+num2str(i))
    end
    sgtitle("True subjects BC trials: Response histograms of Visual Reliability "+reliability_names(l))
end


%% Average over subject plots -> errorbars
PMR_vals_avg = squeeze(mean(PMR_vals, 2)); % Average over subjects
PMR_vals_sem = squeeze(std(PMR_vals, [], 2)) ./ sqrt(num_subjects); % Average over subjects, take SEM

%Prob_C1s = PMR_vals(:,:,1,:) ./ (PMR_vals(:,:,1,:) + PMR_vals(:,:,2,:));
Prob_C1s_avg = squeeze(mean(Prob_C1s, 2));
Prob_C1s_sem = squeeze(std(Prob_C1s, [], 2)) ./ sqrt(num_subjects);


figure;
for l=1:3
    subplot(1,3,l)
    hold on
    errorbar(bincenters, squeeze(PMR_vals_avg(l,1,:)), squeeze(PMR_vals_sem(l,1,:)), ".")
    errorbar(bincenters, squeeze(PMR_vals_avg(l,2,:)), squeeze(PMR_vals_sem(l,2,:)), ".")   
    xlabel("s_A - s_V")
    ylim([0,0.45])
    ylabel("Relative Frequency")
    title("Visual Reliability " + reliability_names(l));
    sgtitle("True subjects BC trials: Response histograms")

end

figure
hold on
for l=1:3   
    errorbar(bincenters, Prob_C1s_avg(l,:), Prob_C1s_sem(l,:), ".", 'Color', rel_colors(l,:), 'LineWidth', 2)
    xlabel("s_A - s_V")
    ylabel("Pr[C hat = 1]")
    title("True subjects BC trials: Probability of choosing unity, by Reliability")
end
legend(reliability_names)

%% Heatmaps: Pr[C_hat = 1 | s_V, s_A]


binedges_sV = linspace(-20, 20, 41);
binedges_sA = -17.5:5:17.5;
bincenters_sV = binedges_sV(1:(end-1)) + (binedges_sV(2:end) - binedges_sV(1:(end-1))) ./2;
bincenters_sA = binedges_sA(1:(end-1)) + (binedges_sA(2:end) - binedges_sA(1:(end-1))) ./2;

PMR_vals_s = zeros(3, num_subjects, 2, length(binedges_sV)-1, length(binedges_sA)-1); % Rels -> subjects -> C=1 or C=2 -> bins.
Prob_C1s = zeros(3, num_subjects, length(binedges_sV)-1, length(binedges_sA)-1);
for l=1:3
    figure
    for i=1:num_subjects
        subplot(3,5,i)
        hold on
        relevant_trial_idx = find((BC_data{i}(:,1) == l)==1);
        relevant_trials = BC_data{i}(relevant_trial_idx,:);
        
        C1_relevant_trials = relevant_trials(relevant_trials(:,4)==1,:);
        C2_relevant_trials = relevant_trials(relevant_trials(:,4)==2,:);
        sA_C1 = C1_relevant_trials(:,2);
        sV_C1 = C1_relevant_trials(:,3);
        sA_C2 = C2_relevant_trials(:,2);
        sV_C2 = C2_relevant_trials(:,3);
        
        PMRs_C1 = [];
        PMRs_C2 = [];
        for kv=1:(length(binedges_sV)-1)
            for ka=1:(length(binedges_sA)-1)
                [l,i,kv,ka]
                PMR_vals_s(l,i,1,kv,ka) = sum(find((int8(sV_C1(:)>=binedges_sV(kv)) .* int8(sV_C1(:)<binedges_sV(kv+1)) .* int8(sA_C1(:)>=binedges_sA(ka)) .* int8(sA_C1(:)<binedges_sA(ka+1)))~=0));
                PMR_vals_s(l,i,2,kv,ka) = sum(find((int8(sV_C2(:)>=binedges_sV(kv)) .* int8(sV_C2(:)<binedges_sV(kv+1)) .* int8(sA_C2(:)>=binedges_sA(ka)) .* int8(sA_C2(:)<binedges_sA(ka+1)))~=0));
            end
        end
        PMR_vals_s(l,i,1,:,:) = PMR_vals_s(l,i,1,:,:) ./ sum(PMR_vals_s(l,i,1,:,:));
        PMR_vals_s(l,i,2,:,:) = PMR_vals_s(l,i,2,:,:) ./ sum(PMR_vals_s(l,i,2,:,:));
        PMR_vals_s(isnan(PMR_vals_s)) = 0;
        
        Prob_C1s(l,i,:,:) = squeeze(PMR_vals_s(l,i,1,:,:) ./ (PMR_vals_s(l,i,1,:,:)+PMR_vals_s(l,i,2,:,:)));
        Prob_C1s(isnan(Prob_C1s)) = 0;
        
        [sA_mesh, sV_mesh] = meshgrid(bincenters_sA, bincenters_sV); 
        contour(sA_mesh, sV_mesh,squeeze(Prob_C1s(l,i,:,:)));
        %imagesc(bincenters_sA, bincenters_sV, 1-Prob_C1s(l,i,:,:));
        newmap = contrast(Prob_C1s(l,i,:,:));
        colormap(newmap)
        xlabel("s_A")
        ylabel("s_V")
        title("Subject "+num2str(i))

    end
    sgtitle("True subjects BC trials: $Pr[\hat{C}=1 | s_V, s_A]$  of Visual Reliability "+reliability_names(l), 'interpreter', 'latex')
end

%%
figure
for l=1:3
    Prob_C1s_mean = squeeze(mean(Prob_C1s(l,:,:,:),2));
    subplot(1,3,l)
    imagesc(bincenters_sA, (bincenters_sV), (Prob_C1s_mean));
    caxis([0 1])
    newmap = contrast((Prob_C1s_mean));
    colormap(flipud(newmap))
    colorbar
    hold on
    plot([-15,15],[-15,15], "r--")
    xlabel("s_A")
    ylabel("s_V")
    title("Visual Reliability " + reliability_names(l));
end
sgtitle("True subjects BC trials: $Pr[\hat{C}=1 | s_V, s_A]$", 'interpreter', 'latex')
