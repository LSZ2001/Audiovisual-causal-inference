% This function converts stratified data (1x15 cell -> each cell is 1x7 cell) 
% to data (only 1x15 cell -> not addtionally stratified by 7 true s_values)
function [data] = data_stratified_to_data(data_stratified, randomize_trials, is_vision)
    if(nargin==3) % By default, assume input is auditory data.
        is_vision=false;
    end
    
    num_s_val_strats = length(data_stratified{1});
    
    [~,num_subjects] = size(data_stratified);
    if(~ is_vision) % Auditory trials
        data = cell(1,num_subjects);
        for i=1:num_subjects
            data_vals = [];
            for j=1:num_s_val_strats
                data_vals = [data_vals; data_stratified{i}{j}];
            end
            if(randomize_trials)
                data_vals = data_vals(randperm(size(data_vals, 1)), :);
            end
            data{i} = data_vals;
        end
        
    else % Vision trials
        data = cell(1,num_subjects);
        for i=1:num_subjects
            data_vals = [];
            for j=1:num_s_val_strats
                data_vals = [data_vals; data_stratified{i}{j}];
            end
            if(randomize_trials)
                data_vals = data_vals(randperm(size(data_vals, 1)), :);
            end
            data{i} = data_vals;
        end
    end
end

