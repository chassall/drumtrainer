% number of ocular artifacts
clear all; clc;
% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152';
end

NOcular = [];

for p = 1:length(ps)
    % Load the preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');
    NOcular = [NOcular EEG.numOcular];
end

%% Load and combine the beta weights

close all; clear all; % Variables will be loaded below

if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data'
    outputFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\analysis\output'
else
    dataFolder = '/Users/rh/Documents/ds004152/';
    outputFolder = '/Users/rh/Documents/ds004152/R/';
end

fontSize = 7;
lineWidth = 1;
condColours = brewermap(5,'Dark2');

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152/';
end

sElectrode = 'FCz'; % electrode of interest
% tWindow = [0.228, 0.304]; % time window of Rewp, in s 50% of peak
tWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)

% Load participant data and get beta values
allBetasDC = [];
allBetasNoDC = [];
for p = 1:length(ps)

    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');

    betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
    betaFile = [subName '_task-drumtrainer_eegbeta.mat'];
    load(fullfile(betaFolder,betaFile),'EEG');
    
    allBetasDC(p,:,:,:) = EEG.unfold.beta_dc;
    allBetasNoDC(p,:,:,:) = EEG.unfold.beta_nodc;

end

% Reminder about order of betas
for i = 1:length(EEG.unfold.eventtypes)
    disp(EEG.unfold.eventtypes{i});
end

meanBetasDC = squeeze(mean(allBetasDC,1));
meanBetasNoDC = squeeze(mean(allBetasNoDC,1));

baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');
% iWindow = dsearchn(EEG.unfold.times',tWindow');

meanBetasDC = meanBetasDC - mean(meanBetasDC(:,baselineI(1):baselineI(2),:),2);
meanBetasNoDC = meanBetasNoDC - mean(meanBetasNoDC(:,baselineI(1):baselineI(2),:),2);

%% plot the ERN vs. FRN for all tempos collapsed
outputFolder = '/Users/rh/Documents/ds004152/R/';

plotColours = cbrewer('qual','Set1',3); % requires cbrewer
sElectrode = 'FCz'; % electrode of interest
iElectrodeFRN = eeg_chaninds(EEG,sElectrode);
sElectrode = 'Fz'; % electrode of interest
iElectrodeERN = eeg_chaninds(EEG,sElectrode);
% alltimes = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
times = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
tWindow1 = [0.240 0.340]; % Sambrook and Goslin (2015)
iWindowFRN =  dsearchn(times',tWindow1');
tWindow2 = [0.070 0.110]; %ERN window around 100ms
iWindowERN =  dsearchn(times',tWindow2');

ontime = squeeze(mean(meanBetasDC(:,:,[2 5 8]),3));
notOntime = squeeze(mean(meanBetasDC(:,:,[1 3 4 6 7 9]),3));
redX = meanBetasDC(:,:,10);
meanDiffERN = redX-ontime;
meanTopoERN = mean(meanDiffERN(:,iWindowERN(1):iWindowERN(2)),2);
meanDiffFRN = notOntime-ontime;
meanTopoFRN = mean(meanDiffFRN(:,iWindowFRN(1):iWindowFRN(2)),2);

makefigure(26,8);
subplot(1,6,1:2);
plot(times,ontime(iElectrodeFRN,:),'Color',plotColours(2,:),'lineWidth',2); hold on;
plot(times,notOntime(iElectrodeFRN,:),'Color',plotColours(3,:),'lineWidth',2); hold on;
plot(times,redX(iElectrodeERN,:),'--','Color',plotColours(1,:),'lineWidth',2);
ax = gca;
area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
xlim([-0.2,0.6]);
legend('Average on-time wave','Average not on-time wave','RedX','Box','off','Location','NorthWest');

plotColours = cbrewer('qual','Dark2',2); % requires cbrewer
subplot(1,6,3:4);
plot(times,meanDiffERN(iElectrodeFRN,:),'Color',plotColours(1,:),'lineWidth',2); hold on;
plot(times,meanDiffFRN(iElectrodeERN,:),'Color',plotColours(2,:),'lineWidth',2); 
ax = gca;
% area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
xlim([-0.2,0.6]);
legend('ERN','FRN','Box','off','Location','NorthWest');
title('Difference wave');

subplot(1,6,5);
tp = topoplot(meanTopoERN,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('ERN');

subplot(1,6,6);
tp = topoplot(meanTopoFRN,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('FRN');

print(fullfile(outputFolder,'rnr_fig1.tiff'),'-dtiff','-r600');

%% make ERN for each tempo
% Setup
close all; clear all; init_unfold(); rng(2022);

% High-level flags
collapseFB = 0; % 0 = all fb separate, 1 = collapse all, 2 = collapse incorrect only
collapseT = 0; % 0 = separate tempos, 1 = collapse across tempos

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152/';
end

% Variables
allBetasDC = [];
allERPs = [];
allRERPs = [];
allERPArtifacts = [];
allBetasNoDC = [];

% Base triggers
% 20 ready screen
% 40 first response
% 60 early
% 80 on time
% 100 late
% 120 red x
readyBaseTrigger = 20;
firstResponseTrigger = 40;
earlyBaseTrigger = 60;
onTimeBaseTrigger = 80;
lateBaseTrigger = 100;
feedbackTriggers = [earlyBaseTrigger onTimeBaseTrigger lateBaseTrigger];
redXBaseTrigger = 120;
baseTriggers = [40 60 80 100 120];

% Condition triggers (1-12)
conditionBaseTriggers = 1:12;
fast = 1:4; % Fast tempo
medium = 5:8; % Medium tempo
slow = 9:12; % Slow tempo
% p1 = [1,2,5,6,9,10]; % Pattern 1
% p2 = [3,4,7,8,11,12]; % Pattern 2

% Convert triggers to BV format
% p1triggers = num2bv(p1, baseTriggers);
% p2triggers = num2bv(p2, baseTriggers);
% readyP1Triggers = num2bv(p1,readyBaseTrigger);
% readyP2Triggers = num2bv(p2,readyBaseTrigger);

respTriggers = num2bv(conditionBaseTriggers,baseTriggers);
firstRespTriggers = num2bv(conditionBaseTriggers,firstResponseTrigger);
earlyTriggers = num2bv(conditionBaseTriggers,earlyBaseTrigger);
onTimeTriggers = num2bv(conditionBaseTriggers,onTimeBaseTrigger);
lateTriggers = num2bv(conditionBaseTriggers,lateBaseTrigger);
redXTriggers = num2bv(conditionBaseTriggers,redXBaseTrigger);

% Fast
respTriggersF = num2bv(fast,baseTriggers);
feedbackTriggersF = num2bv(fast,feedbackTriggers);
firstRespTriggersF = num2bv(fast,firstResponseTrigger);
earlyTriggersF = num2bv(fast,earlyBaseTrigger);
onTimeTriggersF = num2bv(fast,onTimeBaseTrigger);
lateTriggersF = num2bv(fast,lateBaseTrigger);
redXTriggersF = num2bv(fast,redXBaseTrigger);

% Medium
respTriggersM = num2bv(medium,baseTriggers);
feedbackTriggersM = num2bv(medium,feedbackTriggers);
firstRespTriggersM = num2bv(medium,firstResponseTrigger);
earlyTriggersM = num2bv(medium,earlyBaseTrigger);
onTimeTriggersM = num2bv(medium,onTimeBaseTrigger);
lateTriggersM = num2bv(medium,lateBaseTrigger);
redXTriggersM = num2bv(medium,redXBaseTrigger);

% Slow
respTriggersS = num2bv(slow,baseTriggers);
feedbackTriggersS = num2bv(slow,feedbackTriggers);
firstRespTriggersS = num2bv(slow,firstResponseTrigger);
earlyTriggersS = num2bv(slow,earlyBaseTrigger);
onTimeTriggersS = num2bv(slow,onTimeBaseTrigger);
lateTriggersS = num2bv(slow,lateBaseTrigger);
redXTriggersS = num2bv(slow,redXBaseTrigger);

% GLM
% Compute indices of pattern start/end RELATIVE TO POSITION 0
% p1I = 3:4:72; % pattern 1
% p2I = 4:6:72; % pattern 2
ufTime = [-1.5 1.5]; % time window for GLM, in seconds
% ufTime = [-0.9 0.6];

%% Loop through participants to get ERPs BY TEMPO
% Variables
allBetasDC = [];
allERPs = [];
allRERPs = [];
allERPArtifacts = [];
allBetasNoDC = [];

% Base triggers
% 20 ready screen
% 40 first response
% 60 early
% 80 on time
% 100 late
% 120 red x
readyBaseTrigger = 20;
firstResponseTrigger = 40;
earlyBaseTrigger = 60;
onTimeBaseTrigger = 80;
lateBaseTrigger = 100;
feedbackTriggers = [earlyBaseTrigger onTimeBaseTrigger lateBaseTrigger];
redXBaseTrigger = 120;
baseTriggers = [40 60 80 100 120];

% Condition triggers (1-12)
conditionBaseTriggers = 1:12;
fast = 1:4; % Fast tempo
medium = 5:8; % Medium tempo
slow = 9:12; % Slow tempo
% p1 = [1,2,5,6,9,10]; % Pattern 1
% p2 = [3,4,7,8,11,12]; % Pattern 2

% Convert triggers to BV format
% p1triggers = num2bv(p1, baseTriggers);
% p2triggers = num2bv(p2, baseTriggers);
% readyP1Triggers = num2bv(p1,readyBaseTrigger);
% readyP2Triggers = num2bv(p2,readyBaseTrigger);

respTriggers = num2bv(conditionBaseTriggers,baseTriggers);
firstRespTriggers = num2bv(conditionBaseTriggers,firstResponseTrigger);
earlyTriggers = num2bv(conditionBaseTriggers,earlyBaseTrigger);
onTimeTriggers = num2bv(conditionBaseTriggers,onTimeBaseTrigger);
lateTriggers = num2bv(conditionBaseTriggers,lateBaseTrigger);
redXTriggers = num2bv(conditionBaseTriggers,redXBaseTrigger);

% Fast
respTriggersF = num2bv(fast,baseTriggers);
feedbackTriggersF = num2bv(fast,feedbackTriggers);
firstRespTriggersF = num2bv(fast,firstResponseTrigger);
earlyTriggersF = num2bv(fast,earlyBaseTrigger);
onTimeTriggersF = num2bv(fast,onTimeBaseTrigger);
lateTriggersF = num2bv(fast,lateBaseTrigger);
redXTriggersF = num2bv(fast,redXBaseTrigger);

% Medium
respTriggersM = num2bv(medium,baseTriggers);
feedbackTriggersM = num2bv(medium,feedbackTriggers);
firstRespTriggersM = num2bv(medium,firstResponseTrigger);
earlyTriggersM = num2bv(medium,earlyBaseTrigger);
onTimeTriggersM = num2bv(medium,onTimeBaseTrigger);
lateTriggersM = num2bv(medium,lateBaseTrigger);
redXTriggersM = num2bv(medium,redXBaseTrigger);

% Slow
respTriggersS = num2bv(slow,baseTriggers);
feedbackTriggersS = num2bv(slow,feedbackTriggers);
firstRespTriggersS = num2bv(slow,firstResponseTrigger);
earlyTriggersS = num2bv(slow,earlyBaseTrigger);
onTimeTriggersS = num2bv(slow,onTimeBaseTrigger);
lateTriggersS = num2bv(slow,lateBaseTrigger);
redXTriggersS = num2bv(slow,redXBaseTrigger);

% GLM
ufTime = [-1.5 1.5]; % time window for GLM, in seconds


for p = 1:length(ps)
    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');
    
    % Relabel event types for unfold (all early, on time, late, and red X)
    intervalTimes = [];
    patternTimes = [];
    trialCount = 0;
    trialTypes = []; % 0 = early, 1 = on time, 2 = late, 3 = red x, 4 = first resp
    for i = 1:length(EEG.event)
     % Keep a record of trial types
     if ismember(EEG.event(i).type,[earlyTriggers])
         trialTypes = [trialTypes; 0];
     elseif ismember(EEG.event(i).type,[onTimeTriggers])
         trialTypes = [trialTypes; 1];
     elseif ismember(EEG.event(i).type,[lateTriggers])
         trialTypes = [trialTypes; 2];
     elseif ismember(EEG.event(i).type,[redXTriggers])
         trialTypes = [trialTypes; 3];
     end
            % Relabel
         if ismember(EEG.event(i).type,[earlyTriggersF])
                    EEG.event(i).type = 'fastEarly'; trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[onTimeTriggersF])
                    EEG.event(i).type = 'fastOnTime';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[lateTriggersF])
                    EEG.event(i).type = 'fastLate';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[earlyTriggersM])
                    EEG.event(i).type = 'medEarly';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[onTimeTriggersM])
                    EEG.event(i).type = 'medOnTime';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[lateTriggersM])
                    EEG.event(i).type = 'medLate';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[earlyTriggersS])
                    EEG.event(i).type = 'slowEarly';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[onTimeTriggersS])
                    EEG.event(i).type = 'slowOnTime';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[lateTriggersS])
                    EEG.event(i).type = 'slowLate';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[redXTriggersF])
                    EEG.event(i).type = 'fastRedX';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[redXTriggersM])
                    EEG.event(i).type = 'medRedX';trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[redXTriggersS])
                    EEG.event(i).type = 'slowRedX';trialCount = trialCount + 1;
         end
    end
    
    EEG = uf_designmat(EEG,'eventtypes',{'fastEarly','fastOnTime','fastLate','fastRedX','medEarly','medOnTime','medLate','medRedX',...
        'slowEarly','slowOnTime','slowLate','slowRedX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1',...
            'y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1', 'y ~ 1'});

    % Make design matrix for GLM using unfold
    EEG = uf_timeexpandDesignmat(EEG,'timelimits',ufTime);
    
    winrej = uf_continuousArtifactDetect(EEG, 'amplitudeThreshold', 250);
    
    % Make copy of EEG to eventually hold the residual EEG
    rEEG = EEG;
    
    % Solve GLM with Unfold
    EEG= uf_glmfit(EEG);

    % Epoch with Unfold (this step turned the EEG.data dimension from 30 * datapoints into 30 * 750 *
    % events
    EEGall = uf_epoch(EEG,struct('winrej',[],'timelimits',ufTime));
    EEG = uf_epoch(EEG,struct('winrej',winrej,'timelimits',ufTime));


    % Make traditional ERPs with Unfold
    EEG = uf_glmfit_nodc(EEG);

    % Save this participant's EEG, which is now epoched and contains Unfold
    % data
    ufFolder = [dataFolder '/derivatives/eegbeta/' subName];
    if ~exist(ufFolder,'dir')
        mkdir(ufFolder);
    end
    ufFile = [subName '_task-drumtrainer_eegbeta_sepERN.mat'];
    save(fullfile(ufFolder,ufFile),'EEG');
end

return;

%% Loop through participants to get ERPs (Collapsed vs. redX) for locating ERN
% Variables
allBetasDC = [];
allERPs = [];
allRERPs = [];
allERPArtifacts = [];
allBetasNoDC = [];

% Base triggers
% 20 ready screen
% 40 first response
% 60 early
% 80 on time
% 100 late
% 120 red x
readyBaseTrigger = 20;
firstResponseTrigger = 40;
earlyBaseTrigger = 60;
onTimeBaseTrigger = 80;
lateBaseTrigger = 100;
feedbackTriggers = [earlyBaseTrigger onTimeBaseTrigger lateBaseTrigger];
redXBaseTrigger = 120;
baseTriggers = [40 60 80 100 120];

% Condition triggers (1-12)
conditionBaseTriggers = 1:12;
fast = 1:4; % Fast tempo
medium = 5:8; % Medium tempo
slow = 9:12; % Slow tempo
% p1 = [1,2,5,6,9,10]; % Pattern 1
% p2 = [3,4,7,8,11,12]; % Pattern 2

% Convert triggers to BV format
% p1triggers = num2bv(p1, baseTriggers);
% p2triggers = num2bv(p2, baseTriggers);
% readyP1Triggers = num2bv(p1,readyBaseTrigger);
% readyP2Triggers = num2bv(p2,readyBaseTrigger);

respTriggers = num2bv(conditionBaseTriggers,baseTriggers);
firstRespTriggers = num2bv(conditionBaseTriggers,firstResponseTrigger);
earlyTriggers = num2bv(conditionBaseTriggers,earlyBaseTrigger);
onTimeTriggers = num2bv(conditionBaseTriggers,onTimeBaseTrigger);
lateTriggers = num2bv(conditionBaseTriggers,lateBaseTrigger);
redXTriggers = num2bv(conditionBaseTriggers,redXBaseTrigger);

for p = 1:length(ps)
    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');
    
    % Relabel event types for unfold (all early, on time, late, and red X)
    intervalTimes = [];
    patternTimes = [];
    trialCount = 0;
    trialTypes = []; % 0 = early, 1 = on time, 2 = late, 3 = red x, 4 = first resp
    for i = 1:length(EEG.event)
            % Relabel
         if ismember(EEG.event(i).type,[earlyTriggers onTimeTriggers lateTriggers])
                    EEG.event(i).type = 'noRedX'; trialCount = trialCount + 1;
         elseif ismember(EEG.event(i).type,[redXTriggers])
                    EEG.event(i).type = 'RedX';trialCount = trialCount + 1;
         end
    end
    
    EEG = uf_designmat(EEG,'eventtypes',{'noRedX','RedX'},'formula',{'y ~ 1','y ~ 1'});

    % Make design matrix for GLM using unfold
    EEG = uf_timeexpandDesignmat(EEG,'timelimits',ufTime);
    
    winrej = uf_continuousArtifactDetect(EEG, 'amplitudeThreshold', 250);
    
    % Make copy of EEG to eventually hold the residual EEG
    rEEG = EEG;
    
    % Solve GLM with Unfold
    EEG= uf_glmfit(EEG);

    % Epoch with Unfold (this step turned the EEG.data dimension from 30 * datapoints into 30 * 750 *
    % events
    EEGall = uf_epoch(EEG,struct('winrej',[],'timelimits',ufTime));
    EEG = uf_epoch(EEG,struct('winrej',winrej,'timelimits',ufTime));


    % Make traditional ERPs with Unfold
    EEG = uf_glmfit_nodc(EEG);

    % Save this participant's EEG, which is now epoched and contains Unfold
    % data
    ufFolder = [dataFolder '/derivatives/eegbeta/' subName];
    if ~exist(ufFolder,'dir')
        mkdir(ufFolder);
    end
    ufFile = [subName '_task-drumtrainer_eegbeta_collapsedERN.mat'];
    save(fullfile(ufFolder,ufFile),'EEG');
end

return;


%% read in collapsed ERN for locating an ERN
close all; clear all; % Variables will be loaded below

fontSize = 7;
lineWidth = 1;
condColours = brewermap(5,'Dark2');

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152/';
end

% Load participant data and get beta values
allBetasDC = [];
allBetasNoDC = [];
for p = 1:length(ps)

    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');

    betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
    betaFile = [subName '_task-drumtrainer_eegbeta_collapsedERN.mat'];
    load(fullfile(betaFolder,betaFile),'EEG');
    
    allBetasDC(p,:,:,:) = EEG.unfold.beta_dc;
    allBetasNoDC(p,:,:,:) = EEG.unfold.beta_nodc;

end

meanBetasDC = squeeze(mean(allBetasDC,1));
meanBetasNoDC = squeeze(mean(allBetasNoDC,1));

baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');

meanBetasDC = meanBetasDC - mean(meanBetasDC(:,baselineI(1):baselineI(2),:),2);
meanBetasNoDC = meanBetasNoDC - mean(meanBetasNoDC(:,baselineI(1):baselineI(2),:),2);

%% plot the collapsed ERN
sElectrode = 'FCz'; % electrode of interest
iElectrodeERN = eeg_chaninds(EEG,sElectrode);
% alltimes = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
times = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
tWindow = [0.080 0.120]; % ERN window around 100ms
iWindowERN =  dsearchn(times',tWindow');

plotColours = cbrewer('qual','Set1',2); % requires cbrewer

makefigure(18,10);
subplot(1,3,1:2);
plot(times,meanBetasDC(iElectrodeERN,:,2),'Color',plotColours(1,:),'lineWidth',2); hold on;
plot(times,meanBetasDC(iElectrodeERN,:,1),'Color',plotColours(2,:),'lineWidth',2); hold on;
xlim([-0.2,0.6]);
ylim([-1.5 13]);

ax = gca;
area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

legend('red X','all other trials collapsed','Box','off','Location','NorthWest');
title('Average RedX wave');
xlabel('time (s)');
ylabel('Regression coef');

%meanTopo = mean(meanBetasDC(:,iWindowERN(1):iWindowERN(2),2) - meanBetasDC(:,iWindowERN(1):iWindowERN(2),1), 2);
meanTopo = mean(meanBetasDC(:,iWindowERN(1):iWindowERN(2),2), 2);
subplot(1,3,3);
tp = topoplot(meanTopo,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('redX condition');
c = colorbar;
c.Label.String = 'regression coef.';

outputFolder = '/Users/rh/Documents/ds004152/R/';
print(fullfile(outputFolder,'rnr_fig4.tiff'),'-dtiff','-r600');
%% tally the frequency of RedX
close all; clear all; % Variables will be loaded below

fontSize = 7;
lineWidth = 1;
condColours = brewermap(5,'Dark2');

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152/';
end

redXfreq = [];
names = {'fastRedX','medRedX','slowRedX'};

ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

for p = 1:length(ps)

    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');

    betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
    betaFile = [subName '_task-drumtrainer_eegbeta_sepERN.mat'];
    load(fullfile(betaFolder,betaFile),'EEG');
    events = struct2table(EEG.urevent);
    
    for t = 1:3
        redXfreq(p,t) = sum(strcmp(events.type, string(names(t))));
    end
    
    clear EEG events
end

% Load participant data and get beta values
allBetasDC = [];
allBetasNoDC = [];
for p = 1:length(ps)

    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');

    betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
    betaFile = [subName '_task-drumtrainer_eegbeta_sepERN.mat'];
    load(fullfile(betaFolder,betaFile),'EEG');
    
    allBetasDC(p,:,:,:) = EEG.unfold.beta_dc;
    allBetasNoDC(p,:,:,:) = EEG.unfold.beta_nodc;

end

% Reminder about order of betas
for i = 1:length(EEG.unfold.eventtypes)
    disp(EEG.unfold.eventtypes{i});
end

n_censor = 8;
allBetasDC(redXfreq(:,1) <= n_censor,:,:,4) = NaN;
allBetasDC(redXfreq(:,2) <= n_censor,:,:,8) = NaN;
allBetasDC(redXfreq(:,3) <= n_censor,:,:,12) = NaN;

meanBetasDC = squeeze(mean(allBetasDC,1, 'omitnan'));
meanBetasNoDC = squeeze(mean(allBetasNoDC,1, 'omitnan'));

baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');

meanBetasDC = meanBetasDC - mean(meanBetasDC(:,baselineI(1):baselineI(2),:),2);
meanBetasNoDC = meanBetasNoDC - mean(meanBetasNoDC(:,baselineI(1):baselineI(2),:),2);

allBetasDC = allBetasDC - mean(allBetasDC(:,:,baselineI(1):baselineI(2),:),3, 'omitnan');
allBetasNoDC = allBetasNoDC - mean(allBetasNoDC(:,:,baselineI(1):baselineI(2),:),3, 'omitnan');


% Raw signal for topo - RedX

outputFolder = '/Users/rh/Documents/ds004152/R/';
plotColours = [205,183, 158; 139,125,107; 0, 0, 0;]/255; % requires cbrewer
plotColours2 = cbrewer('qual','Dark2',2); % requires cbrewer

sElectrode = 'FCz'; % electrode of interest
iElectrodeERN = eeg_chaninds(EEG,sElectrode);
% alltimes = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
times = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
tWindow1 = [0.240 0.340]; % Sambrook and Goslin (2015)
iWindowFRN =  dsearchn(times',tWindow1');
tWindow2 = [0.080 0.120]; % ERN window around 100ms
iWindowERN =  dsearchn(times',tWindow2');

earlyF = squeeze(meanBetasDC(:,:,1));
ontimeF = squeeze(meanBetasDC(:,:,2));
lateF = squeeze(meanBetasDC(:,:,3));
redXF = squeeze(meanBetasDC(:,:,4));
earlyM = squeeze(meanBetasDC(:,:,5));
ontimeM = squeeze(meanBetasDC(:,:,6));
lateM = squeeze(meanBetasDC(:,:,7));
redXM = squeeze(meanBetasDC(:,:,8));
earlyS = squeeze(meanBetasDC(:,:,9));
ontimeS = squeeze(meanBetasDC(:,:,10));
lateS = squeeze(meanBetasDC(:,:,11));
redXS = squeeze(meanBetasDC(:,:,12));

Ontime = squeeze(mean(meanBetasDC(:,:,[2 6 10]),3));
notOntime = squeeze(mean(meanBetasDC(:,:,[1 3 5 7 9 11]),3));
redX = squeeze(mean(meanBetasDC(:,:,[4 8 12]),3));

meanTopoERNF = mean(redXF(:,iWindowERN(1):iWindowERN(2)),2);

meanTopoERNM = mean(redXM(:,iWindowERN(1):iWindowERN(2)),2);

meanTopoERNS = mean(redXS(:,iWindowERN(1):iWindowERN(2)),2);

letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
%%

makefigure(26,18);

subplot(4,4,[1 2 5 6]);
plot(times,redXF(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2); hold on;
plot(times,redXM(iElectrodeERN,:),'Color',plotColours(2,:),'lineWidth',2); hold on;
plot(times,redXS(iElectrodeERN,:),'Color',plotColours(3,:),'lineWidth',2); hold on;
ax = gca;
area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
xlim([-0.2,0.6]);
ylim([-1.5 20]);
legend('Fast','Medium','Slow','Box','off','Location','NorthWest');
title('Average RedX wave');
xlabel('time (s)');
ylabel('Regression coef');

letter = letters(1);
xlims = get(ax, 'XLim');
ylims = get(ax, 'YLim');
text(ax.XLim(1)-0.2, ax.YLim(2)+0.3, letter, 'FontSize', 14);

min1 = min([meanTopoERNF meanTopoERNM meanTopoERNS],[],'all');
max1 = max([meanTopoERNF meanTopoERNM meanTopoERNS],[],'all');

subplot(4,4,9);
tp = topoplot(meanTopoERNF,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill','maplimits',[min1,max1]);
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('Fast');
letter = letters(3);
text(-1.2, 0.6, letter, 'FontSize', 14);

subplot(4,4,10);
tp = topoplot(meanTopoERNM,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill','maplimits',[min1,max1]);
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('Medium');


subplot(4,4,11);
tp = topoplot(meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill','maplimits',[min1,max1]);
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('Slow');
c = colorbar;
c.Label.String = 'regression coef.';

subplot(4,4,[3 4 7 8]);
allRawF = mean(squeeze(allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),4)),2);
allRawF(redXfreq(:,1)<=n_censor) = NaN;
allRawM = mean(squeeze(allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),8)),2);
allRawM(redXfreq(:,2)<=n_censor) = NaN;
allRawS = mean(squeeze(allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),12)),2);
allRawS(redXfreq(:,3)<=n_censor) = NaN;

nbp = notBoxPlot([allRawF allRawM allRawS]); 
formatNBP(nbp); 
set(gca,'XTickLabel',{'Fast','Medium','Slow'});
xlabel('Tempo');
ylabel('ERN amplitude');
letter = letters(2);
xlims = get(ax, 'XLim');
ylims = get(ax, 'YLim');
text(-0.5, 14.5, letter, 'FontSize', 14);

min2 = min([meanTopoERNF-meanTopoERNS,meanTopoERNM-meanTopoERNS],[],'all');
max2 = max([meanTopoERNF-meanTopoERNS,meanTopoERNM-meanTopoERNS],[],'all');
subplot(4,4,13);
tp = topoplot(meanTopoERNF-meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill','maplimits',[min2,max2]);
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('Fast - Slow');
letter = letters(4);
text(-1.2, 0.6, letter, 'FontSize', 14);


subplot(4,4,14);
tp = topoplot(meanTopoERNM-meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill','maplimits',[min2,max2]);
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
colormap(topoMap);
title('Medium - Slow');
c = colorbar;
c.Label.String = 'regression coef.';
% c = colorbar;
% c.Label.String = 'amplitude';

print(fullfile(outputFolder,'rnr_fig3.tiff'),'-dtiff','-r600');

%%
find(meanTopoERNM-meanTopoERNS == max(meanTopoERNM-meanTopoERNS))

find(meanTopoERNF-meanTopoERNS == max(meanTopoERNF-meanTopoERNS))
%%
[H,P,CI,STATS] = ttest(allRawF,allRawM)
mean(allRawF,'omitNaN') - mean(allRawM,'omitNaN')
[H,P,CI,STATS] = ttest(allRawF,allRawS)
mean(allRawF,'omitNaN') - mean(allRawS,'omitNaN')
%%
% %% plot the ERN vs. FRN for all tempos SEPARATELY
% outputFolder = '/Users/rh/Documents/ds004152/R/';
% 
% purpleGold = [85, 26, 139; 255, 215, 0; 205, 150, 205;]/255;
% plotColours = cbrewer('qual','Set1',1); % requires cbrewer
% plotColours2 = cbrewer('qual','Dark2',2); % requires cbrewer
% 
% sElectrode = 'FCz'; % electrode of interest
% iElectrodeFRN = eeg_chaninds(EEG,sElectrode);
% sElectrode = 'FCz'; % electrode of interest
% iElectrodeERN = eeg_chaninds(EEG,sElectrode);
% % alltimes = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
% times = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
% tWindow1 = [0.240 0.340]; % Sambrook and Goslin (2015)
% iWindowFRN =  dsearchn(times',tWindow1');
% tWindow2 = [0.070 0.110]; % ERN window around 100ms
% iWindowERN =  dsearchn(times',tWindow2');
% 
% earlyF = squeeze(meanBetasDC(:,:,1));
% ontimeF = squeeze(meanBetasDC(:,:,2));
% lateF = squeeze(meanBetasDC(:,:,3));
% redXF = squeeze(meanBetasDC(:,:,4));
% earlyM = squeeze(meanBetasDC(:,:,5));
% ontimeM = squeeze(meanBetasDC(:,:,6));
% lateM = squeeze(meanBetasDC(:,:,7));
% redXM = squeeze(meanBetasDC(:,:,8));
% earlyS = squeeze(meanBetasDC(:,:,9));
% ontimeS = squeeze(meanBetasDC(:,:,10));
% lateS = squeeze(meanBetasDC(:,:,11));
% redXS = squeeze(meanBetasDC(:,:,12));
% 
% Ontime = squeeze(mean(meanBetasDC(:,:,[2 6 10]),3));
% notOntime = squeeze(mean(meanBetasDC(:,:,[1 3 5 7 9 11]),3));
% redX = squeeze(mean(meanBetasDC(:,:,[4 8 12]),3));
% 
% meanDiffERNF = redXF-ontimeF;
% meanTopoERNF = mean(meanDiffERNF(:,iWindowERN(1):iWindowERN(2)),2);
% meanDiffFRNF = ((earlyF-ontimeF) + (lateF-ontimeF))/2;
% meanTopoFRNF = mean(meanDiffFRNF(:,iWindowFRN(1):iWindowFRN(2)),2);
% 
% meanDiffERNM = redXM-ontimeM;
% meanTopoERNM = mean(meanDiffERNM(:,iWindowERN(1):iWindowERN(2)),2);
% meanDiffFRNM = ((earlyM-ontimeM) + (lateM-ontimeM))/2;
% meanTopoFRNM = mean(meanDiffFRNM(:,iWindowFRN(1):iWindowFRN(2)),2);
% 
% meanDiffERNS = redXS-ontimeS;
% meanTopoERNS = mean(meanDiffERNS(:,iWindowERN(1):iWindowERN(2)),2);
% meanDiffFRNS = ((earlyS-ontimeS) + (lateS-ontimeS))/2;
% meanTopoFRNS = mean(meanDiffFRNS(:,iWindowFRN(1):iWindowFRN(2)),2);
% 
% makefigure(26,18);
% subplot(3,6,1:2);
% plot(times,earlyF(iElectrodeFRN,:),'Color',purpleGold(1,:),'lineWidth',2); hold on;
% plot(times,ontimeF(iElectrodeFRN,:),'Color',purpleGold(2,:),'lineWidth',2); hold on;
% plot(times,lateF(iElectrodeFRN,:),'Color',purpleGold(3,:),'lineWidth',2); hold on;
% plot(times,redXF(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2);
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 20]);
% legend('early','on-time','late','redX','Box','off','Location','NorthWest');
% title('Fast tempo');
% 
% subplot(3,6,3:4);
% plot(times,meanDiffERNF(iElectrodeFRN,:),'Color',plotColours2(1,:),'lineWidth',2); hold on;
% plot(times,meanDiffFRNF(iElectrodeERN,:),'Color',plotColours2(2,:),'lineWidth',2); 
% ax = gca;
% % area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% % area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-5 10]);
% legend('ERN (redX - on-time)','FRN (not on time - on-time)','Box','off','Location','NorthWest');
% title('Difference wave');
% 
% subplot(3,6,5);
% tp = topoplot(meanTopoERNF,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('ERN');
% 
% subplot(3,6,6);
% tp = topoplot(meanTopoFRNF,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('FRN');
% 
% subplot(3,6,7:8);
% plot(times,earlyM(iElectrodeFRN,:),'Color',purpleGold(1,:),'lineWidth',2); hold on;
% plot(times,ontimeM(iElectrodeFRN,:),'Color',purpleGold(2,:),'lineWidth',2); hold on;
% plot(times,lateM(iElectrodeFRN,:),'Color',purpleGold(3,:),'lineWidth',2); hold on;
% plot(times,redXM(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2);
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 20]);
% legend('early','on-time','late','redX','Box','off','Location','NorthWest');
% title('Medium tempo');
% 
% subplot(3,6,9:10);
% plot(times,meanDiffERNM(iElectrodeFRN,:),'Color',plotColours2(1,:),'lineWidth',2); hold on;
% plot(times,meanDiffFRNM(iElectrodeERN,:),'Color',plotColours2(2,:),'lineWidth',2); 
% ax = gca;
% % area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% % area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-5 10]);
% legend('ERN (redX - on-time)','FRN (not on time - on-time)','Box','off','Location','NorthWest');
% title('Difference wave');
% 
% subplot(3,6,11);
% tp = topoplot(meanTopoERNM,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('ERN');
% 
% subplot(3,6,12);
% tp = topoplot(meanTopoFRNM,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('FRN');
% 
% subplot(3,6,13:14);
% plot(times,earlyS(iElectrodeFRN,:),'Color',purpleGold(1,:),'lineWidth',2); hold on;
% plot(times,ontimeS(iElectrodeFRN,:),'Color',purpleGold(2,:),'lineWidth',2); hold on;
% plot(times,lateS(iElectrodeFRN,:),'Color',purpleGold(3,:),'lineWidth',2); hold on;
% plot(times,redXS(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2);
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 20]);
% legend('early','on-time','late','redX','Box','off','Location','NorthWest');
% title('Slow tempo');
% 
% plotColours = cbrewer('qual','Dark2',2); % requires cbrewer
% subplot(3,6,15:16);
% plot(times,meanDiffERNS(iElectrodeFRN,:),'Color',plotColours2(1,:),'lineWidth',2); hold on;
% plot(times,meanDiffFRNS(iElectrodeERN,:),'Color',plotColours2(2,:),'lineWidth',2); 
% ax = gca;
% % area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% % area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-5 10]);
% legend('ERN (redX - on-time)','FRN (not on time - on-time)','Box','off','Location','NorthWest');
% title('Difference wave');
% 
% subplot(3,6,17);
% tp = topoplot(meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('ERN');
% 
% subplot(3,6,18);
% tp = topoplot(meanTopoFRNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('FRN');
% 
% print(fullfile(outputFolder,'rnr_fig2.tiff'),'-dtiff','-r600');
% 
% %% Different time window, same difference wave at FCz
% outputFolder = '/Users/rh/Documents/ds004152/R/';
% 
% purpleGold = [85, 26, 139; 255, 215, 0; 205, 150, 205;]/255;
% plotColours = cbrewer('qual','Set1',1); % requires cbrewer
% plotColours2 = cbrewer('qual','Dark2',2); % requires cbrewer
% 
% sElectrode = 'FCz'; % electrode of interest
% iElectrodeFRN = eeg_chaninds(EEG,sElectrode);
% sElectrode = 'FCz'; % electrode of interest
% iElectrodeERN = eeg_chaninds(EEG,sElectrode);
% % alltimes = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
% times = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
% tWindow1 = [0.240 0.340]; % Sambrook and Goslin (2015)
% iWindowFRN =  dsearchn(times',tWindow1');
% tWindow2 = [0.050 0.100]; % ERN window around 100ms
% iWindowERN =  dsearchn(times',tWindow2');
% 
% earlyF = squeeze(meanBetasDC(:,:,1));
% ontimeF = squeeze(meanBetasDC(:,:,2));
% lateF = squeeze(meanBetasDC(:,:,3));
% redXF = squeeze(meanBetasDC(:,:,4));
% earlyM = squeeze(meanBetasDC(:,:,5));
% ontimeM = squeeze(meanBetasDC(:,:,6));
% lateM = squeeze(meanBetasDC(:,:,7));
% redXM = squeeze(meanBetasDC(:,:,8));
% earlyS = squeeze(meanBetasDC(:,:,9));
% ontimeS = squeeze(meanBetasDC(:,:,10));
% lateS = squeeze(meanBetasDC(:,:,11));
% redXS = squeeze(meanBetasDC(:,:,12));
% 
% Ontime = squeeze(mean(meanBetasDC(:,:,[2 6 10]),3));
% notOntime = squeeze(mean(meanBetasDC(:,:,[1 3 5 7 9 11]),3));
% redX = squeeze(mean(meanBetasDC(:,:,[4 8 12]),3));
% 
% meanDiffFRNF = ((earlyF-ontimeF) + (lateF-ontimeF))/2;
% meanTopoFRNF = mean(meanDiffFRNF(:,iWindowFRN(1):iWindowFRN(2)),2);
% meanTopoERNF = mean(meanDiffFRNF(:,iWindowERN(1):iWindowERN(2)),2);
% 
% meanDiffFRNM = ((earlyM-ontimeM) + (lateM-ontimeM))/2;
% meanTopoFRNM = mean(meanDiffFRNM(:,iWindowFRN(1):iWindowFRN(2)),2);
% meanTopoERNM = mean(meanDiffFRNM(:,iWindowERN(1):iWindowERN(2)),2);
% 
% meanDiffFRNS = ((earlyS-ontimeS) + (lateS-ontimeS))/2;
% meanTopoFRNS = mean(meanDiffFRNS(:,iWindowFRN(1):iWindowFRN(2)),2);
% meanTopoERNS = mean(meanDiffFRNS(:,iWindowERN(1):iWindowERN(2)),2);
% 
% makefigure(26,18);
% subplot(3,6,1:2);
% plot(times,earlyF(iElectrodeFRN,:),'Color',purpleGold(1,:),'lineWidth',2); hold on;
% plot(times,ontimeF(iElectrodeFRN,:),'Color',purpleGold(2,:),'lineWidth',2); hold on;
% plot(times,lateF(iElectrodeFRN,:),'Color',purpleGold(3,:),'lineWidth',2); hold on;
% plot(times,redXF(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2);
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 20]);
% legend('early','on-time','late','redX','Box','off','Location','NorthWest');
% title('Fast tempo');
% 
% subplot(3,6,3:4);
% plot(times,meanDiffFRNF(iElectrodeERN,:),'Color',plotColours2(2,:),'lineWidth',2); hold on;
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-5 5]);
% legend('not on time - on-time','Box','off','Location','NorthWest');
% title('Difference wave');
% 
% subplot(3,6,5);
% tp = topoplot(meanTopoERNF,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('ERN');
% 
% subplot(3,6,6);
% tp = topoplot(meanTopoFRNF,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('FRN');
% 
% subplot(3,6,7:8);
% plot(times,earlyM(iElectrodeFRN,:),'Color',purpleGold(1,:),'lineWidth',2); hold on;
% plot(times,ontimeM(iElectrodeFRN,:),'Color',purpleGold(2,:),'lineWidth',2); hold on;
% plot(times,lateM(iElectrodeFRN,:),'Color',purpleGold(3,:),'lineWidth',2); hold on;
% plot(times,redXM(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2);
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 20]);
% legend('early','on-time','late','redX','Box','off','Location','NorthWest');
% title('Medium tempo');
% 
% subplot(3,6,9:10);
% plot(times,meanDiffFRNM(iElectrodeERN,:),'Color',plotColours2(2,:),'lineWidth',2); hold on;
% ax = gca;
% % area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% % area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-5 5]);
% legend('not on time - on-time','Box','off','Location','NorthWest');
% title('Difference wave');
% 
% subplot(3,6,11);
% tp = topoplot(meanTopoERNM,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('ERN');
% 
% subplot(3,6,12);
% tp = topoplot(meanTopoFRNM,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('FRN');
% 
% subplot(3,6,13:14);
% plot(times,earlyS(iElectrodeFRN,:),'Color',purpleGold(1,:),'lineWidth',2); hold on;
% plot(times,ontimeS(iElectrodeFRN,:),'Color',purpleGold(2,:),'lineWidth',2); hold on;
% plot(times,lateS(iElectrodeFRN,:),'Color',purpleGold(3,:),'lineWidth',2); hold on;
% plot(times,redXS(iElectrodeERN,:),'Color',plotColours(1,:),'lineWidth',2);
% ax = gca;
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 20]);
% legend('early','on-time','late','redX','Box','off','Location','NorthWest');
% title('Slow tempo');
% 
% plotColours = cbrewer('qual','Dark2',2); % requires cbrewer
% subplot(3,6,15:16);
% plot(times,meanDiffFRNS(iElectrodeFRN,:),'Color',plotColours2(2,:),'lineWidth',2); hold on;
% ax = gca;
% % area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% % area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow1, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-5 5]);
% legend('not on time - on-time','Box','off','Location','NorthWest');
% title('Difference wave');
% 
% subplot(3,6,17);
% tp = topoplot(meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('ERN');
% 
% subplot(3,6,18);
% tp = topoplot(meanTopoFRNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('FRN');
% 
% %% Raw signal for topo - incorrect
% 
% outputFolder = '/Users/rh/Documents/ds004152/R/';
% plotColours = [205,183, 158; 139,125,107; 0, 0, 0;]/255; % requires cbrewer
% plotColours2 = cbrewer('qual','Dark2',2); % requires cbrewer
% 
% sElectrode = 'FCz'; % electrode of interest
% iElectrodeFRN = eeg_chaninds(EEG,sElectrode);
% sElectrode = 'FCz'; % electrode of interest
% iElectrodeERN = eeg_chaninds(EEG,sElectrode);
% % alltimes = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
% times = -1.5:(1/EEG.srate):(1.5-1/EEG.srate);
% tWindow1 = [0.240 0.340]; % Sambrook and Goslin (2015)
% iWindowFRN =  dsearchn(times',tWindow1');
% tWindow2 = [0.130 0.170]; % ERN window around 100ms
% iWindowERN =  dsearchn(times',tWindow2');
% 
% earlyF = squeeze(meanBetasDC(:,:,1));
% ontimeF = squeeze(meanBetasDC(:,:,2));
% lateF = squeeze(meanBetasDC(:,:,3));
% redXF = squeeze(meanBetasDC(:,:,4));
% earlyM = squeeze(meanBetasDC(:,:,5));
% ontimeM = squeeze(meanBetasDC(:,:,6));
% lateM = squeeze(meanBetasDC(:,:,7));
% redXM = squeeze(meanBetasDC(:,:,8));
% earlyS = squeeze(meanBetasDC(:,:,9));
% ontimeS = squeeze(meanBetasDC(:,:,10));
% lateS = squeeze(meanBetasDC(:,:,11));
% redXS = squeeze(meanBetasDC(:,:,12));
% 
% Ontime = squeeze(mean(meanBetasDC(:,:,[2 6 10]),3));
% notOntime = squeeze(mean(meanBetasDC(:,:,[1 3 5 7 9 11]),3));
% redX = squeeze(mean(meanBetasDC(:,:,[4 8 12]),3));
% 
% meanDiffFRNF = ((earlyF-ontimeF) + (lateF-ontimeF))/2;
% meanRawF = (earlyF + lateF)/2;
% meanTopoFRNF = mean(meanDiffFRNF(:,iWindowFRN(1):iWindowFRN(2)),2);
% meanTopoERNF = mean(meanRawF(:,iWindowERN(1):iWindowERN(2)),2);
% 
% meanDiffFRNM = ((earlyM-ontimeM) + (lateM-ontimeM))/2;
% meanRawM = (earlyM + lateM)/2;
% meanTopoFRNM = mean(meanDiffFRNM(:,iWindowFRN(1):iWindowFRN(2)),2);
% meanTopoERNM = mean(meanRawM(:,iWindowERN(1):iWindowERN(2)),2);
% 
% meanDiffFRNS = ((earlyS-ontimeS) + (lateS-ontimeS))/2;
% meanRawS = (earlyS + lateS)/2;
% meanTopoFRNS = mean(meanDiffFRNS(:,iWindowFRN(1):iWindowFRN(2)),2);
% meanTopoERNS = mean(meanRawS(:,iWindowERN(1):iWindowERN(2)),2);
% 
% makefigure(26,18);
% subplot(4,4,[1 2 5 6]);
% plot(times,meanRawF(iElectrodeFRN,:),'Color',plotColours(1,:),'lineWidth',2); hold on;
% plot(times,meanRawM(iElectrodeFRN,:),'Color',plotColours(2,:),'lineWidth',2); hold on;
% plot(times,meanRawS(iElectrodeFRN,:),'Color',plotColours(3,:),'lineWidth',2); hold on;
% ax = gca;
% area(tWindow2, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% area(tWindow2, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
% xlim([-0.2,0.6]);
% ylim([-1.5 15]);
% legend('Fast','Medium','Slow','Box','off','Location','NorthWest');
% title('Average incorrect wave');
% 
% 
% subplot(4,4,9);
% tp = topoplot(meanTopoERNF,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('Fast');
% 
% subplot(4,4,10);
% tp = topoplot(meanTopoERNM,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('Medium');
% 
% 
% 
% subplot(4,4,11);
% tp = topoplot(meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('Slow');
% 
% subplot(4,4,[3 4 7 8]);
% allRawF = mean(squeeze((allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),1) + allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),3))/2),2);
% allRawM = mean(squeeze((allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),5) + allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),7))/2),2);
% allRawS = mean(squeeze((allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),9) + allBetasDC(:,iElectrodeFRN,iWindowERN(1):iWindowERN(2),11))/2),2);
% 
% nbp = notBoxPlot([allRawF allRawM allRawS]); formatNBP(nbp); %,['Fast','Medium','Slow']
% set(gca,'XTickLabel',{'Fast','Medium','Slow'});
% 
% subplot(4,4,13);
% tp = topoplot(meanTopoERNF-meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('Fast - Slow');
% 
% 
% subplot(4,4,14);
% tp = topoplot(meanTopoERNM-meanTopoERNS,EEG.chanlocs,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
% tp.Parent.XLim = [-0.6 0.6];
% tp.Parent.YLim = [-0.6 0.6];
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% colormap(topoMap);
% title('Medium - Slow');
% 
% print(fullfile(outputFolder,'rnr_fig3.tiff'),'-dtiff','-r600');


