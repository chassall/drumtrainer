
close all; clear all; clc;

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';
end

% Conditions
% 1, fast, pattern 1, LRRR
% 2, fast, pattern 1, RLLL
% 3, fast, pattern 2, LRRRRR
% 4, fast, pattern 2, RLLLLL
% 5, medium, pattern 1, LRRR
% 6, medium, pattern 1, RLLL
% 7, medium, pattern 2, LRRRRR
% 8, medium, pattern 2, RLLLLL
% 9, slow, pattern 1, LRRR
% 10, slow, pattern 1, RLLL
% 11, slow, pattern 2, LRRRRR
% 12, slow, pattern 2, RLLLLL

% Loop through participants
fitPs = [];
for p = 1:length(ps)

    oldFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data_old/derivatives/eegtbt';
    % Load trial-by-trial scores 
    subName = ['sub-' ps{p}];
    tbtFolder = [dataFolder '/derivatives/erptbt/' subName];
    tbtFile = [subName '_task-drumtrainer_erptbt.mat'];
    load(fullfile(tbtFolder,tbtFile),'erpRewp','erpArtifact','rERPRewp','rERPArtifact');
    %load(fullfile(oldFolder,subName,[subName '_task-dt_eegtbt_0_0.mat']));

    % Load behavioural data
    rawFile = [subName '_task-drumtrainer_beh.tsv'];
    beh = readtable(fullfile(dataFolder,subName,'beh',rawFile),'FileType','text');
    
    rt = beh.trialResp_rt;
    isBadRT = rt < 0.1 | rt > 2;
    isBadPrevRT = [0; isBadRT(1:end-1)];
    rtAdjust = [NaN; beh.trialResp_rt(2:end) - beh.trialResp_rt(1:end-1)];
    isFirstTrial = beh.trialLoop_thisRepN == 0;
    
    isFastTempo = ismember(beh.blockType,1:4);
    isMediumTempo = ismember(beh.blockType,5:8);
    isSlowTempo = ismember(beh.blockType,9:12);
    
    tooFastI = beh.outcome == 0 & ~isFirstTrial & ~isBadRT & ~isBadPrevRT;
    onTimeI = beh.outcome == 1 & ~isFirstTrial & ~isBadRT & ~isBadPrevRT;
    tooSlowI = beh.outcome == 2 & ~isFirstTrial & ~isBadRT & ~isBadPrevRT;

    subplot(1,2,1);
    plot(rt(~isBadRT));
    subplot(1,2,2);
    bar([mean(rtAdjust(tooFastI)) mean(rtAdjust(onTimeI)) mean(rtAdjust(tooSlowI))]);
    pause();

end