% Shuze Liu:
% Col 1,2,11: always 0.
% Col 3: visual coherence level idx in {1,2,3} (1 is smallest sigma). 0 if audio stimulus only.
% Col 4: the demanded task response: 1 for visual, 2 for auditory, 3 for same/diff judgment.
% Col 5: auditory source location in pixels -> converted to deg. Can only take on 7 values.
% Col 6: 0 when auditory only -> position of visual stimulus in pixels -> converted to deg. 
% Col 7: 1 for single stimulus present, 2 for both audio/visual stimulus.
% Col 8?: -45 when Col4 ~= 2. Otherwise, looks random. -> Audio stim position estimate of subject.
% Col 9?: -45 when Col4 ~= 1. Otherwise, looks random. -> Visual stim position estimate of subject.
% Col 10: report for same/diff judgment, either 1 (same) or 2 (diff). 0 for trials with other demands.
close all; clear all;

addpath ../matlab
names = {'aaron' 'allison' 'avauna' 'clare' 'connie' 'daisy' 'david' 'gabriela' 'jeff' 'jobyeon' 'julia' 'katherine' 'priyanka' 'stephanie' 'wendy'} ;

mup = 507; % center of the screen -- where the white fixation cross lies.
pxtodeg = 0.0892857142; % pixels to degrees conversion scale

% Col 1,2 are never used...
for subjidx = 1:length(names)
    subjname = names{subjidx}
    fname = [subjname 'userdata' '.mat'];
    load(fname);
    D = userdata;
    
    % Remove negative pixels trials!
    D = D(D(:,8)>=0, :);
    D = D(D(:,9)>=0, :);
    
    % Col 3 is visual reliability levels in [1,2,3]. But why can it take on
    % values of 0? -> because these trials are auditory only?!
    
    %On each trial, the std of the Gaussian spot was randomly chosen from
    %three possible values: 21, 130, and 250 pixels.
    D(D(:,3)==21,3)=1;
    D(D(:,3)==130,3)=2;
    D(D(:,3)==250,3)=3;
    
    % Col 5,6,8,9 must be some coordinates -> need to transfer to deg. 
    D(:,[5 6 8 9]) = (D(:,[5 6 8 9])-mup)*pxtodeg;
    
    % Col 5 need rounding -> maybe because of 7 discrete audio speakers?
    D(:,5)= round(D(:,5));
    
    % For visual trials, replace 0's in aud response Col with NaN. 
    % For auditory trials, replace 0's in vis response Col with NaN. 
    
    % Col 7 can only take on values in {1,2} -> binomial trials always 2 ->
    % must be the # of stimuli in the task!
    % Col 4 seems to indicate the demanded task reponse.
    data{subjidx}.dataC = D(D(:,4)==3 & D(:,7)==2,[3 5 6 10]); % bimodal trials, unity judgment
    data{subjidx}.dataA = D(D(:,4)==2 & D(:,7)==2,[3 5 6 8]); % bimodal trials, auditory localization
    data{subjidx}.dataV = D(D(:,4)==1 & D(:,7)==2,[3 5 6 9]); % bimodal trials, visual localization
    
    data{subjidx}.dataA_uni = D(D(:,4)==2 & D(:,7)==1,[3 5 6 8]); % unimodal trials, auditory localization
    data{subjidx}.dataV_uni = D(D(:,4)==1 & D(:,7)==1,[3 5 6 9]); % unimodal trials, visual localization
    
    data{subjidx}.dataA_uni(:,[1,3]) = NaN; % Remove irrelevant visual stim info 
    data{subjidx}.dataV_uni(:,2) = NaN; % Remove irrelevant auditory stim info 
end

save alldata data