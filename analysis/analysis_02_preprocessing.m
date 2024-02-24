%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% Preprocess the Drum Trainer EEG
% 
% Other m-files required: 
% EEGLAB toolbox
% ccn_prep.m
% ccn_check.m
% find_artifacts.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% change according to needs
gitDir = '/Users/rh/Documents/GitHub/drumtrainer';
MATLABfunctionsDir = '/Users/rh/Documents/MATLAB';
eegLABDir = '/Users/rh/Documents/eeglab';
bvaDir = '/Users/rh/Documents/GitHub/bva-io';

% Set raw and preprocessed folders - set as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152/';
end

% cd gitDir;
addpath(genpath(gitDir)); 
addpath(genpath(bvaDir));
addpath(genpath(MATLABfunctionsDir));
addpath(genpath(eegLABDir));
addpath(dataFolder);

eeglab;

% Participant strings
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% ICA settings - which triggers? Window size?
icaTriggers = {'S 21','S 22','S 23','S 24','S 25','S 26','S 27','S 28','S 29','S 30','S 31','S 32'};
icaWindow = [0  30]; % In seconds

% Artifact check settings - which triggers? Window size? Baseline?
checkTriggers = {'S 81','S 82','S 83','S 84','S 85','S 86','S 87','S 88','S 89','S 90','S 91','S 92'};
checkWindow = [-0.2,0.6]; % In seconds
baseline = [-200 0]; % In milliseconds

% High-pass filter setting - use 0.1 for FRN and 0.01 for temporal scaling
highPass = 0.1;  
appendString = '';

toRemove = {};
allBadChannels = {};

artifactSettings.maxMin = 100;
artifactSettings.level = 100;
artifactSettings.step = 40;
artifactSettings.lowest = 0.1;

% Loop through participants
for p = 1:length(ps)

    rng(2021); % Set for consistency

    subName = ['sub-' ps{p}];
    taskName = 'drumtrainer';
    
    filters = [0.1 20];
    reference = {'TP9','TP10'};
    
    % Specify bad channels, overwrite as needed
    badChannels = {};
    if strcmp(ps{p},'21')
        badChannels = {'Fp1','Fp2'};
    end

    ccn_prep(dataFolder,subName,taskName,icaTriggers,icaWindow,filters,reference,badChannels,appendString);
    allBadChannels{p} = ccn_check(dataFolder,subName,taskName,checkTriggers,checkWindow,baseline,artifactSettings);
end