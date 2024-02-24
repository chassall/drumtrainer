%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% "First response" triggers did not record properly beyond first block
% Load the preprocessed data, add the triggers, and resave
%
% Other m-files required: 
% EEGLAB toolbox

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Setup
close all; clear all;

ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder
dataFolder = '/Users/rh/Documents/ds004152/';

% Trigger settings (triggers are base + condition code)

% Base triggersx
readyBaseTrigger = 20;
firstResponseTrigger = 40;
earlyBaseTrigger = 60;
onTimeBaseTrigger = 80;
lateBaseTrigger = 100;
redXBaseTrigger = 120;
baseTriggers = [40 60 80 100 120];

% Condition triggers (1-12)
conditionBaseTriggers = 1:12;
fast = 1:4; % Fast tempo
medium = 5:8; % Medium tempo
slow = 9:12; % Slow tempo
p1 = [1,2,5,6,9,10]; % Pattern 1
p2 = [3,4,7,8,11,12]; % Pattern 2

% Convert triggers to BV format
cueTriggers = num2bv(1:12);
readyTriggers = num2bv(1:12,readyBaseTrigger);
p1triggers = num2bv(p1, baseTriggers);
p2triggers = num2bv(p2, baseTriggers);
readyP1Triggers = num2bv(p1,readyBaseTrigger);
readyP2Triggers = num2bv(p2,readyBaseTrigger);
respTriggers = num2bv(conditionBaseTriggers,baseTriggers);
firstRespTriggers = num2bv(conditionBaseTriggers,firstResponseTrigger);
earlyTriggers = num2bv(conditionBaseTriggers,earlyBaseTrigger);
onTimeTriggers = num2bv(conditionBaseTriggers,onTimeBaseTrigger);
lateTriggers = num2bv(conditionBaseTriggers,lateBaseTrigger);
redXTriggers = num2bv(conditionBaseTriggers,redXBaseTrigger);

% Loop through participants
for p = 1:length(ps)
    subName = ['sub-' ps{p}];

    % Load beh data
    rawFile = [subName '_task-drumtrainer_beh.tsv'];
    thisData = readtable(fullfile(dataFolder,subName,'beh',rawFile),'FileType','text');

    % Find missing trials using behavioural record
    missingTrials = thisData.trialLoop_thisRepN == 71;
    missingRTs = thisData.trialResp_rt(missingTrials);
    missingTypes = thisData.trigger(missingTrials);
    
    % Load preprocessed EEG
    prepFileName = [subName '_task-drumtrainer_eegprep.mat'];
    EEG = pop_loadset(fullfile(dataFolder,'derivatives','eegprep',subName,prepFileName));
        
    % Add in the first response trigger
    for i = 1:length(EEG.event)-1
        if ismember(EEG.event(i).type,readyTriggers)
            thisCondition = str2num(EEG.event(i).type(2:end)) - readyBaseTrigger;
            thisTrigger = num2bv(thisCondition,firstResponseTrigger);
            EEG.event(i+1).type = thisTrigger{:};
        end
    end
    
    % Add in the missing last trigger
    newTriggerLatencies = [];
    newTriggerTypes = [];
    missingTriggerCount = 0;
    for i = 1:length(EEG.event)
        if ismember(EEG.event(i).type,firstRespTriggers) 
            missingTriggerCount = missingTriggerCount + 1;
            newTriggerLatencies = [newTriggerLatencies EEG.event(i+71).latency + missingRTs(missingTriggerCount) * EEG.srate];
        end
    end
    disp(missingTriggerCount);
   
    for t = 1:length(newTriggerLatencies)
        numEvents = length(EEG.event)
        thisTrigger = num2bv(missingTypes(t));
        EEG.event(numEvents+1).type = thisTrigger{:};
        EEG.event(numEvents+1).latency = newTriggerLatencies(t);
        EEG.event(numEvents+1).urevent = numEvents + 1;
    end
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    % Save over old prep file
    save(fullfile(dataFolder,'derivatives','eegprep',subName,prepFileName),'EEG');
end