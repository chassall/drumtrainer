%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% This script adds files and information to a BIDS folder that has been 
% created by Brain Visions's BV2BIDS command line tool
% 
% Other m-files required: 
% bids-matlab toolbox (https://github.com/bids-standard/bids-matlab)

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% specify BIDS folder
bidsFolder =  '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';

% specify source folder (contain the raw behavioural files
% that will be bidsified below)
sourceFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data_source';

%% read dataset_description.json, created by BV2BIDS, and modify
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
datasetFile = fullfile(bidsFolder, 'dataset_description.json');
dataset = bids.util.jsondecode(datasetFile);
dataset.DatasetType = 'raw';
% dataset.ReferencesAndLinks = {''};
options.indent = '  '; % Adds white space, easier to read
bids.util.jsonwrite(datasetFile, dataset, options);

%% write participant tsv files, not created by BV2BIDS
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html

% load participant file
participantFile = fullfile(bidsFolder, 'participants.tsv');
participants.participant_id = {'sub-01'; 'sub-02'; 'sub-03'; 'sub-04'; 'sub-05'; 'sub-06'; 'sub-07'; 'sub-08'; 'sub-09'; 'sub-10'; 'sub-11'; 'sub-13'; 'sub-14'; 'sub-15'; 'sub-16'; 'sub-17'; 'sub-18'; 'sub-19'; 'sub-20'; 'sub-21'};
participants.participant = {'01'; '02';'03'; '04';'05';'06';'07';'08';'09'; '10'; '11'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'};
participants.datetime = {'19-May-2021 11:56:37'; '21-May-2021 10:24:48'; '26-May-2021 09:36:23';... 
'27-May-2021 10:34:13'; '28-May-2021 11:00:31'; '03-Jun-2021 11:12:24';...
'10-Jun-2021 10:43:34'; '14-Jun-2021 12:27:05'; '16-Jun-2021 09:56:57';...
'16-Jun-2021 12:57:28'; '10-Nov-2021 15:30:09'; ...
'19-Nov-2021 14:31:18'; '26-Nov-2021 14:34:36'; '03-Dec-2021 13:41:07';...
'08-Dec-2021 14:35:15'; '06-Jan-2022 15:11:54'; '02-Feb-2022 10:23:30';...
'10-Feb-2022 11:20:46'; '16-Feb-2022 10:44:16'; '25-Feb-2022 13:03:09'};

participants.age = [23; 22; 23; 25; 27; 29; 26; 41; 33; 25; 24; 24; 22; 28; 27; 21; 23; 24; 26; 24];
participants.sex = {'F';'M'; 'M'; 'F'; 'F'; 'F'; 'F'; 'M'; 'F'; 'M'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'F'; 'M'; 'F'};

participants.handedness = {'R'; 'L'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'L'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'};

bids.util.tsvwrite(participantFile, participants);

%% write participant json file, not created by BV2BIDS
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
participantJSONFile = fullfile(bidsFolder, 'participants.json');
pInfoDesc.participant.Description = 'participant number as input by the tester';
pInfoDesc.datetime.Description = 'date and time at start of task';
pInfoDesc.age.Description = 'self-reported age of participant';
pInfoDesc.age.Units = 'years';
pInfoDesc.sex.Description = 'self-reported sex of participant';
pInfoDesc.sex.Levels.M = 'male';
pInfoDesc.sex.Levels.F = 'female';
pInfoDesc.handedness.Description = 'self-reported handedness of participant';
pInfoDesc.handedness.Levels.L = 'left-handed';
pInfoDesc.handedness.Levels.R = 'right-handed';
pInfoDesc.handedness.Levels.LR = 'ambidextrous';
options.indent = '  '; % Adds white space, easier to read
bids.util.jsonwrite(participantJSONFile,pInfoDesc,options);

%% read eeg json for each participant, created by BV2BIDS, and modify
% https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#electroencephalography
whichPs = [1:11 13:21];
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    eegJSONFile = fullfile(bidsFolder, pString,'eeg',[pString '_task-drumtrainer_eeg.json']);
    eeg = bids.util.jsondecode(eegJSONFile);
    eeg.InstitutionName = 'University of Oxford';
    eeg.InstitutionAddress = 'Warneford Hospital, Oxford, OX3 7JX';
    eeg.InstitutionalDepartmentName = 'Department of Psychiatry';
    eeg.ManufacturersModelName = 'actiCHamp Plus';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(eegJSONFile, eeg, options);
end

%% write events json for each participant, not created by BV2BIDS
whichPs = [1:11 13:21];
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    eventsJSONFile = fullfile(bidsFolder, pString,'eeg',[pString '_task-drumtrainer_events.json']);

    eInfoDesc.trial_type.Description = 'BrainVision event number and type (string)';
    eInfoDesc.value.Description = 'Event value (string)';
    eInfoDesc.value.Levels.S1 = 'Condition 1 metronome beat';
    eInfoDesc.value.Levels.S2 = 'Condition 2 metronome beat';
    eInfoDesc.value.Levels.S3 = 'Condition 3 metronome beat';
    eInfoDesc.value.Levels.S4 = 'Condition 4 metronome beat';
    eInfoDesc.value.Levels.S5 = 'Condition 5 metronome beat';
    eInfoDesc.value.Levels.S6 = 'Condition 6 metronome beat';
    eInfoDesc.value.Levels.S7 = 'Condition 7 metronome beat';
    eInfoDesc.value.Levels.S8 = 'Condition 8 metronome beat';
    eInfoDesc.value.Levels.S9 = 'Condition 9 metronome beat';
    eInfoDesc.value.Levels.S10 = 'Condition 10 metronome beat';
    eInfoDesc.value.Levels.S11 = 'Condition 11 metronome beat';
    eInfoDesc.value.Levels.S12 = 'Condition 12 metronome beat';
    eInfoDesc.value.Levels.S21 = 'Condition 1 ready screen';
    eInfoDesc.value.Levels.S22 = 'Condition 2 ready screen';
    eInfoDesc.value.Levels.S23 = 'Condition 3 ready screen';
    eInfoDesc.value.Levels.S24 = 'Condition 4 ready screen';
    eInfoDesc.value.Levels.S25 = 'Condition 5 ready screen';
    eInfoDesc.value.Levels.S26 = 'Condition 6 ready screen';
    eInfoDesc.value.Levels.S27 = 'Condition 7 ready screen';
    eInfoDesc.value.Levels.S28 = 'Condition 8 ready screen';
    eInfoDesc.value.Levels.S29 = 'Condition 9 ready screen';
    eInfoDesc.value.Levels.S30 = 'Condition 10 ready screen';
    eInfoDesc.value.Levels.S31 = 'Condition 11 ready screen';
    eInfoDesc.value.Levels.S32 = 'Condition 12 ready screen';
    eInfoDesc.value.Levels.S41 = 'Condition 1 first response';
    eInfoDesc.value.Levels.S42 = 'Condition 2 first response';
    eInfoDesc.value.Levels.S43 = 'Condition 3 first response';
    eInfoDesc.value.Levels.S44 = 'Condition 4 first response';
    eInfoDesc.value.Levels.S45 = 'Condition 5 first response';
    eInfoDesc.value.Levels.S46 = 'Condition 6 first response';
    eInfoDesc.value.Levels.S47 = 'Condition 7 first response';
    eInfoDesc.value.Levels.S48 = 'Condition 8 first response';
    eInfoDesc.value.Levels.S49 = 'Condition 9 first response';
    eInfoDesc.value.Levels.S50 = 'Condition 10 first response';
    eInfoDesc.value.Levels.S51 = 'Condition 11 first response';
    eInfoDesc.value.Levels.S52 = 'Condition 12 first response';
    eInfoDesc.value.Levels.S61 = 'Condition 1 early response';
    eInfoDesc.value.Levels.S62 = 'Condition 2 early response';
    eInfoDesc.value.Levels.S63 = 'Condition 3 early response';
    eInfoDesc.value.Levels.S64 = 'Condition 4 early response';
    eInfoDesc.value.Levels.S65 = 'Condition 5 early response';
    eInfoDesc.value.Levels.S66 = 'Condition 6 early response';
    eInfoDesc.value.Levels.S67 = 'Condition 7 early response';
    eInfoDesc.value.Levels.S68 = 'Condition 8 early response';
    eInfoDesc.value.Levels.S69 = 'Condition 9 early response';
    eInfoDesc.value.Levels.S70 = 'Condition 10 early response';
    eInfoDesc.value.Levels.S71 = 'Condition 11 early response';
    eInfoDesc.value.Levels.S72 = 'Condition 12 early response';
    eInfoDesc.value.Levels.S81 = 'Condition 1 on time response';
    eInfoDesc.value.Levels.S82 = 'Condition 2 on time response';
    eInfoDesc.value.Levels.S83 = 'Condition 3 on time response';
    eInfoDesc.value.Levels.S84 = 'Condition 4 on time response';
    eInfoDesc.value.Levels.S85 = 'Condition 5 on time response';
    eInfoDesc.value.Levels.S86 = 'Condition 6 on time response';
    eInfoDesc.value.Levels.S87 = 'Condition 7 on time response';
    eInfoDesc.value.Levels.S88 = 'Condition 8 on time response';
    eInfoDesc.value.Levels.S89 = 'Condition 9 on time response';
    eInfoDesc.value.Levels.S90 = 'Condition 10 on time response';
    eInfoDesc.value.Levels.S91 = 'Condition 11 on time response';
    eInfoDesc.value.Levels.S92 = 'Condition 12 on time response';
    eInfoDesc.value.Levels.S101 = 'Condition 1 late response';
    eInfoDesc.value.Levels.S102 = 'Condition 2 late response';
    eInfoDesc.value.Levels.S103 = 'Condition 3 late response';
    eInfoDesc.value.Levels.S104 = 'Condition 4 late response';
    eInfoDesc.value.Levels.S105 = 'Condition 5 late response';
    eInfoDesc.value.Levels.S106 = 'Condition 6 late response';
    eInfoDesc.value.Levels.S107 = 'Condition 7 late response';
    eInfoDesc.value.Levels.S108 = 'Condition 8 late response';
    eInfoDesc.value.Levels.S109 = 'Condition 9 late response';
    eInfoDesc.value.Levels.S110 = 'Condition 10 late response';
    eInfoDesc.value.Levels.S111 = 'Condition 11 late response';
    eInfoDesc.value.Levels.S112 = 'Condition 12 late response';
    eInfoDesc.value.Levels.S121 = 'Condition 1 red X (incorrect response)';
    eInfoDesc.value.Levels.S122 = 'Condition 2 red X (incorrect response)';
    eInfoDesc.value.Levels.S123 = 'Condition 3 red X (incorrect response)';
    eInfoDesc.value.Levels.S124 = 'Condition 4 red X (incorrect response)';
    eInfoDesc.value.Levels.S125 = 'Condition 5 red X (incorrect response)';
    eInfoDesc.value.Levels.S126 = 'Condition 6 red X (incorrect response)';
    eInfoDesc.value.Levels.S127 = 'Condition 7 red X (incorrect response)';
    eInfoDesc.value.Levels.S128 = 'Condition 8 red X (incorrect response)';
    eInfoDesc.value.Levels.S129 = 'Condition 9 red X (incorrect response)';
    eInfoDesc.value.Levels.S130 = 'Condition 10 red X (incorrect response)';
    eInfoDesc.value.Levels.S131 = 'Condition 11 red X (incorrect response)';
    eInfoDesc.value.Levels.S132 = 'Condition 12 red X (incorrect response)';
    eInfoDesc.onset.Description = 'Event onset';
    eInfoDesc.onset.Units = 'milisecond';
    eInfoDesc.duration.Description = 'Event duration';
    eInfoDesc.duration.Units = 'milisecond';
    eInfoDesc.channel.Description = 'Channel number';
    eInfoDesc.channel.Levels.x0 = 'All channels';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(eventsJSONFile,eInfoDesc,options);
end

%% load behavioural data and save as tsv files
whichPs = [1:11 13:21];
for p = 1:length(whichPs)
    pString = num2str(whichPs(p),'%0.02i');
    pBIDSString = ['sub-' pString];
    thisDir = dir(fullfile(sourceFolder,pBIDSString,'beh','*.csv')); % This file has a header
    thisData = readtable(fullfile(sourceFolder,pBIDSString,'beh',thisDir.name));

    % Remove metronome data (since no behaviour)
    thisData(isnan(thisData.trialLoop_thisRepN),:) = [];

    thisData.blockType = thisData.blockType + 1; % Switch to index from 1

%     thisData.firstKey(strcmp(thisData.firstKey,'')) = {'n/a'};
%     toChange = thisData.readyKey_corrv
%     thisData.readyKey_corr = num2cell(thisData.readyKey_corr);
%     thisData.readyKey_corr(thisData.readyKey_corr) = {'n/a'};

    % Remove some unused columns that were inserted by PsychoPy
    thisData.instLoop_thisRepN = [];
    thisData.instLoop_thisTrialN = [];
    thisData.instLoop_thisN = [];
    thisData.instLoop_thisIndex = [];
    thisData.blockLoop_thisTrialN = [];
    thisData.blockLoop_thisIndex = [];
    thisData.blockLoop_thisN = [];
    thisData.cueLoop_thisRepN = [];
    thisData.cueLoop_thisTrialN = [];
    thisData.cueLoop_thisIndex = [];
    thisData.cueLoop_thisN = [];
    thisData.cuePort_started = [];
    thisData.cuePort_stopped = [];
    thisData.cueSnare_started = [];
    thisData.cueSnare_stopped = [];
    thisData.cueBass_started = [];
    thisData.cueBass_stopped = [];
    thisData.trialLoop_thisTrialN = [];
    thisData.trialLoop_thisIndex = [];
    thisData.trialLoop_thisN = [];
    thisData.leftSpace_started = [];
    thisData.leftSpace_stopped = [];
    thisData.leftSquare_started = [];
    thisData.leftSquare_stopped = [];
    thisData.rightSquare_started = [];
    thisData.rightSpace_started = [];
    thisData.rightSpace_stopped = [];
    thisData.rightSquare_stopped = [];
    thisData.readyKey_keys = [];
    readyKey_started.readyKey_started = [];
    thisData.readyPort_started = [];
    thisData.readyPort_stopped = [];
    thisData.readyKey_started = [];
    thisData.readyKey_stopped = [];
    thisData.trialPort_started = [];	
    thisData.trialPort_stopped = [];
    thisData.trialFeedback_stopped = [];
    thisData.trialFixation_stopped = [];
    thisData.trialResp_stopped = [];	
    thisData.trialSnare_started = [];	
    thisData.trialSnare_stopped = [];	
    thisData.trialHihat_started = [];	
    thisData.trialHihat_stopped = [];	
    thisData.trialBass_started = [];	
    thisData.trialBass_stopped = [];
    thisData.colourCode = [];
    thisData.redX_started = [];	
    thisData.redX_stopped = [];
    thisData.trialFeedback_started = [];
    thisData.trialFixation_started = [];
    thisData.trialResp_started = [];
    thisData.participant = [];	
    thisData.date = [];	
    thisData.expName = [];
    thisData.psychopyVersion = [];
    thisData.frameRate = [];
    

    % Make a struct out of the behavioural data
    beh.blockLoop_thisRepN = thisData.blockLoop_thisRepN;
    beh.trialLoop_thisRepN = thisData.trialLoop_thisRepN;
    beh.trigger = thisData.trigger;
    beh.firstKey = thisData.firstKey;
    beh.readyKey_corr = thisData.readyKey_corr;
    beh.readyKey_rt = thisData.readyKey_rt;
    beh.targetKey = thisData.targetKey;
    beh.firstKey = thisData.firstKey;
    beh.readyKey_corr = thisData.readyKey_corr;
    beh.readyKey_rt = thisData.readyKey_rt;
    beh.targetKey = thisData.targetKey;
    beh.blockType = thisData.blockType;
    beh.margin = thisData.margin;
    beh.outcome = thisData.outcome;
    beh.trialResp_keys = thisData.trialResp_keys;
    beh.trialResp_corr = thisData.trialResp_corr;
    beh.trialResp_rt = thisData.trialResp_rt;
    
    behFolder = fullfile(bidsFolder,pBIDSString,'beh');
    behFile = [pBIDSString '_task-drumtrainer_beh.tsv'];
    bids.util.tsvwrite(fullfile(behFolder,behFile), beh);

%     behFolder = fullfile(bidsFolder,pBIDSString,'beh');
%     behFile = [pBIDSString '_task-drumtrainer_beh.tsv'];
%     if ~exist(behFolder, 'dir')
%         mkdir(behFolder);
%     end
%     writetable(thisData,fullfile(behFolder,behFile),'Delimiter','\t','FileType','text');
end

%% write beh json for each participant
whichPs = [1:11 13:21];
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    behJSONFile = fullfile(bidsFolder, pString,'beh',[pString '_task-drumtrainer_beh.json']);

    bInfoDesc.blockLoop_thisRepN.Description = 'Block number (integer)';	
	bInfoDesc.trialLoop_thisRepN.Description = 'Trial number (integer)';
    bInfoDesc.trigger.Description	 = 'EEG trigger for this response (integer)';
    bInfoDesc.firstKey.Description = 'First key pressed this block (character)';
    bInfoDesc.firstKey.Levels.f = 'Left';
    bInfoDesc.firstKey.Levels.j = 'Right';
    bInfoDesc.readyKey_corr.Description = 'First key was correct (boolean)';
    bInfoDesc.readyKey_corr.Levels.x0 = 'First key press this block was incorrect';
    bInfoDesc.readyKey_corr.Levels.x1 = 'First key press this block was correct';
    bInfoDesc.readyKey_rt.Description = 'First key RT (seconds)';
    bInfoDesc.targetKey.Description = 'Target key for this round (character)';
    bInfoDesc.targetKey.Levels.f = 'Left';
    bInfoDesc.targetKey.Levels.j = 'Right';
    bInfoDesc.blockType.Description = 'Condition code (integer)';
    bInfoDesc.blockType.Levels.x1 = 'Fast, pattern 1, left-hand start';
    bInfoDesc.blockType.Levels.x2 = 'Fast, pattern 1, right-hand start';
    bInfoDesc.blockType.Levels.x3 = 'Fast, pattern 2, left-hand start';
    bInfoDesc.blockType.Levels.x4 = 'Fast, pattern 2, right-hand start';
    bInfoDesc.blockType.Levels.x5 = 'Medium, pattern 1, left-hand start';
    bInfoDesc.blockType.Levels.x6 = 'Medium, pattern 1, right-hand start';
    bInfoDesc.blockType.Levels.x7 = 'Medium, pattern 2, left-hand start';
    bInfoDesc.blockType.Levels.x8 = 'Medium, pattern 2, right-hand start';
    bInfoDesc.blockType.Levels.x9 = 'Slow, pattern 1, left-hand start';
    bInfoDesc.blockType.Levels.x10 = 'Slow, pattern 1, right-hand start';
    bInfoDesc.blockType.Levels.x11 = 'Slow, pattern 2, left-hand start';
    bInfoDesc.blockType.Levels.x12 = 'Slow, pattern 2, right-hand start';
    bInfoDesc.margin.Description	= 'Margin around target time (seconds)';
    bInfoDesc.outcome.Description	= 'Trial outcome (integer)';
    bInfoDesc.outcome.x0 = 'Early';
    bInfoDesc.outcome.x1 = 'On time';
    bInfoDesc.outcome.x2 = 'Late';
    bInfoDesc.trialResp_keys.Description = 'Response (character)';
    bInfoDesc.trialResp_keys.Levels.f = 'Left';
    bInfoDesc.trialResp_keys.Levels.j = 'Right';
    bInfoDesc.trialResp_corr.Description = 'Correct response flag (boolean)';
    bInfoDesc.trialResp_corr.Levels.x0 = 'Incorrect key was pressed';
    bInfoDesc.trialResp_corr.Levels.x1 = 'Correct key was pressed';
    bInfoDesc.trialResp_rt.Description = 'Response time (seconds)';


    options.indent = '  '; % Adds white space, easier to read
    delete(behJSONFile);
    bids.util.jsonwrite(behJSONFile,bInfoDesc,options);
end