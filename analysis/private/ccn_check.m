function badChannels = ccn_check(dataFolder,subName,taskName,triggers,window,baseline,artifactSettings)
%CCN_CHECK Summary of this function goes here
%   Detailed explanation goes here

preprocessedFolder = [dataFolder '/derivatives/eegprep/' subName];
eyeCorrFile = [subName '_task-' taskName '_eegprep.mat'];

load(fullfile(preprocessedFolder,eyeCorrFile),'EEG');


[~, allArtifacts,~] = make_erp(EEG,triggers,window,baseline,artifactSettings);

allLabels = {EEG.chanlocs.labels};
isBad = allArtifacts > 0.10;
%badChannels = {};
%theseBadChannelIs = find(isBad);
badChannels= allLabels(isBad);

if isempty(badChannels)
    disp('all channels good');
else
    disp(['bad channels: ' badChannels]);
end

