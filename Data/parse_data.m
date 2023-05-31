clear all;
% cd('C:\Users\liu_s\Audiovisual-causal-inference\Analysis')
data_path = "..\Data\";
load(data_path+"alldata.mat");

s_a_range = -15:5:15;
s_v_binedges = linspace(-20,20,8);

%% BC data
num_subjects = length(data);
BC_data = cell(1,num_subjects);
for i=1:num_subjects
    BC_data{i} = data{i}.dataC;
end
save(data_path + 'BC_data.mat', 'BC_data');

%% BAV data
BAV_data = cell(1,num_subjects);
for i=1:num_subjects
    BV_subj = [data{i}.dataV, repmat(1,length(data{i}.dataV),1)];
    BA_subj = [data{i}.dataA, repmat(2,length(data{i}.dataA),1)];
    BAV_data{i} = [BV_subj; BA_subj];
end
save(data_path + 'BAV_data.mat', 'BAV_data');

%% UA data
data_stratified = cell(1,15);
for i=1:num_subjects
    data_stratified{i} = cell(1, length(s_a_range));
    dataA_subject = data{i}.dataA_uni;
    s_vals = dataA_subject(:,2); % save visual stimulus location, and reliability level of each trial.
    s_hat_vals = dataA_subject(:,end);
    
    for j=1:length(s_a_range)
        relevant_trials = s_vals==s_a_range(j);
        s_hats_relevant = s_hat_vals(relevant_trials);
        data_stratified{i}{j} = [s_vals(relevant_trials), s_hats_relevant];
    end
end
data_stratified_UA = data_stratified;
save(data_path+'data_stratified_UA', 'data_stratified_UA')

%% UV data
s_v_bin_centers = (s_v_binedges(1:end-1) + s_v_binedges(2:end)) ./ 2;
data_stratified = cell(1,15);
for i=1:num_subjects
    data_stratified{i} = cell(1, length(s_v_binedges)-1);
    dataV_subject = data{i}.dataV_uni;
    s_vals = dataV_subject(:,3); % save visual stimulus location, and reliability level of each trial.
    reliability_vals = dataV_subject(:,1);
    s_hat_vals = dataV_subject(:,end);
    
    for j=1:(length(s_v_binedges)-1)
        relevant_trials = (s_vals>=s_v_binedges(j)).*(s_vals<s_v_binedges(j+1));
        relevant_trials = find(relevant_trials ~=0);
        s_hats_relevant = s_hat_vals(relevant_trials);
        data_stratified{i}{j} = [s_vals(relevant_trials), s_hats_relevant, reliability_vals(relevant_trials)];
    end
end
data_stratified_UV = data_stratified;
save(data_path+'data_stratified_UV', 'data_stratified_UV')
