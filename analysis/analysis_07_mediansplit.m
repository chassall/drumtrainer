%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% Make regression-ERPs based on a median split of RT
%
% Other m-files required: 
% EEGLAB toolbox
% Unfold toolbox
% num2bv.m
% make_erp();

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Setup
close all; clear all; init_unfold(); rng(2022);

% High-level flags
collapseFB = 0; % 0 = all fb separate, 1 = collapse all, 2 = collapse incorrect only
collapseT = 0; % 0 = separate tempos, 1 = collapse across tempos
vartype = 'rt_cat'; %rt_cat rt_adjust_cat
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

% Load RT info
rtInfo = readtable([dataFolder 'R/rt_vector.csv']);

% Loop through participants
for p = 1:length(ps)

    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');
    
    % Get this participant's RT info
    thisRTInfo = rtInfo(rtInfo.participant_ID == str2num(ps{p}),:);

    % Relabel event types for unfold (all early, on time, late, and red X)
    intervalTimes = [];
    patternTimes = [];
    trialCount = 1;
    trialTypes = []; % 0 = early, 1 = on time, 2 = late, 3 = red x, 4 = first resp
    for i = 1:length(EEG.event)
        if ismember(EEG.event(i).type,[earlyTriggers onTimeTriggers lateTriggers])
            % Keep a record of trial types
            if ismember(EEG.event(i).type,[earlyTriggers  ])
                trialTypes = [trialTypes; 0];
            elseif ismember(EEG.event(i).type,[onTimeTriggers])
                trialTypes = [trialTypes; 1];
            elseif ismember(EEG.event(i).type,[lateTriggers])
                trialTypes = [trialTypes; 2];
            end
            % Relabel
            if collapseFB && collapseT
                switch collapseFB
                    case 1
                        EEG.event(i).type = 'feedback';
                    case 2
                        if ismember(EEG.event(i).type,[earlyTriggers lateTriggers])
                            EEG.event(i).type = 'incorrect';
                        elseif ismember(EEG.event(i).type,[ onTimeTriggers ])
                            EEG.event(i).type = 'correct';
                        end
                end
            elseif collapseFB && ~collapseT
                switch collapseFB
                    case 1
                        if ismember(EEG.event(i).type,[feedbackTriggersF])
                            EEG.event(i).type = 'fast';
                        elseif ismember(EEG.event(i).type,[ feedbackTriggersM ])
                            EEG.event(i).type = 'medium';
                        elseif ismember(EEG.event(i).type,[  feedbackTriggersS])
                            EEG.event(i).type = 'slow';
                        end
                    case 2
                        error('flag combo not supported');
                end
            elseif ~collapseFB && collapseT
                if ismember(EEG.event(i).type,[earlyTriggers])
                    EEG.event(i).type = 'early';
                elseif ismember(EEG.event(i).type,[ onTimeTriggers ])
                    EEG.event(i).type = 'onTime';
                elseif ismember(EEG.event(i).type,[  lateTriggers])
                    EEG.event(i).type = 'late';
                end

            elseif ~collapseFB && ~collapseT && trialCount <= size(thisRTInfo,1)

                thisRTAdjustCat = num2str(thisRTInfo{trialCount,vartype});
%                 thisRTAdjustCat = num2str(thisRTInfo{trialCount,'rt_cat'});

                if ismember(EEG.event(i).type,[earlyTriggersF])
                    EEG.event(i).type = ['fastEarly' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[onTimeTriggersF])
                    EEG.event(i).type = ['fastOnTime' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[lateTriggersF])
                    EEG.event(i).type = ['fastLate' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[earlyTriggersM])
                    EEG.event(i).type = ['medEarly' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[onTimeTriggersM])
                    EEG.event(i).type = ['medOnTime' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[lateTriggersM])
                    EEG.event(i).type = ['medLate' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[earlyTriggersS])
                    EEG.event(i).type = ['slowEarly' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[onTimeTriggersS])
                    EEG.event(i).type = ['slowOnTime' thisRTAdjustCat]; trialCount = trialCount + 1;
                elseif ismember(EEG.event(i).type,[lateTriggersS])
                    EEG.event(i).type = ['slowLate' thisRTAdjustCat]; trialCount = trialCount + 1;
                end

            end
        elseif ismember(EEG.event(i).type,redXTriggers)
            trialTypes = [trialTypes; 3];
            EEG.event(i).type = 'redX'; trialCount = trialCount + 1;
        end
    end

    % Make design matrix for GLM using unfold
    if collapseFB && collapseT
        switch collapseFB
            case 1
                EEG = uf_designmat(EEG,'eventtypes',{'feedback','redX'},'formula',{'y ~ 1','y ~ 1'});
            case 2
                EEG = uf_designmat(EEG,'eventtypes',{'incorrect','correct','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1'});
        end
    elseif collapseFB && ~collapseT
        switch collapseFB
            case 1
                EEG = uf_designmat(EEG,'eventtypes',{'fast','medium','slow','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
            case 2
                error('flag combo not supported');
        end
    elseif ~collapseFB && collapseT
        EEG = uf_designmat(EEG,'eventtypes',{'early','onTime','late','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
    elseif ~collapseFB && ~collapseT
        
%         EEG = uf_designmat(EEG,'eventtypes',{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly',...
%             'slowOnTime','slowLate','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
        EEG = uf_designmat(EEG,'eventtypes',{'fastEarly0','fastOnTime0','fastLate0','medEarly0','medOnTime0','medLate0','slowEarly0',...
            'slowOnTime0','slowLate0','fastEarly1','fastOnTime1','fastLate1','medEarly1','medOnTime1','medLate1','slowEarly1',...
            'slowOnTime1','slowLate1','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
    end
    EEG = uf_timeexpandDesignmat(EEG,'timelimits',ufTime);
   
    
    winrej = uf_continuousArtifactDetect(EEG, 'amplitudeThreshold', 250);
    
    % Make copy of EEG to eventually hold the residual EEG
    rEEG = EEG;
    
    % Solve GLM with Unfold
    EEG= uf_glmfit(EEG);

    % Epoch with Unfold
    EEG = uf_epoch(EEG,struct('winrej',winrej,'timelimits',ufTime));

    % Make traditional ERPs with Unfold
    EEG = uf_glmfit_nodc(EEG);

    % Save this participant's EEG, which is now epoched and contains Unfold
    % data
    ufFolder = [dataFolder '/derivatives/eegbetartms/' subName];
    if ~exist(ufFolder,'dir')
        mkdir(ufFolder);
    end
    ufFile = [subName '_task-drumtrainer_eegbetartms.mat'];
    save(fullfile(ufFolder,ufFile),'EEG');

end

return;

%% Load and combine the beta weights

close all; clear all; % Variables will be loaded below

if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data'
    outputFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\analysis\output'
else
    dataFolder = '/Users/rh/Documents/ds004152/';
    outputFolder = '/Users/rh/Documents/ds004152/output/';
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

% Load participant data and get beta values
allBetasDC = [];
allBetasNoDC = [];
for p = 1:length(ps)

    % Load preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG');
    
    % Load beta values
    betaFolder = [dataFolder '/derivatives/eegbetartms/' subName];
    betaFile = [subName '_task-drumtrainer_eegbetartms.mat'];
    load(fullfile(betaFolder,betaFile),'EEG');
    allBetasDC(p,:,:,:) = EEG.unfold.beta_dc;
    allBetasNoDC(p,:,:,:) = EEG.unfold.beta_nodc;
end

% Reminder about order of betas
for i = 1:length(EEG.unfold.eventtypes)
    disp(EEG.unfold.eventtypes{i});
end

% Baseline correction
baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');
allBetasDC = allBetasDC - mean(allBetasDC(:,:,baselineI(1):baselineI(2),:),3);
allBetasNoDC = allBetasNoDC - mean(allBetasNoDC(:,:,baselineI(1):baselineI(2),:),3);

meanBetasDC = squeeze(mean(allBetasDC,1));
meanBetasNoDC = squeeze(mean(allBetasNoDC,1));

%% Make the figure 
% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.2 0.08], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

sElectrode = 'FCz'; % electrode of interest
iElectrode = eeg_chaninds(EEG,sElectrode);
% tWindow = [0.228, 0.304]; % time window of Rewp, in s 50% of peak
tWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)
iWindow = dsearchn(EEG.unfold.times',tWindow');

close all; 
makefigure(18,10);
%all conditions:
whichBs = [1 3 10 12; 4 6 13 15; 7 9 16 18];
%only slow condition:


axs = {};
fbColours = [];
for i = 1:3
    axs{i} = subplot(2,3,i);
    phs = [];
    for j = 1:2
        phs(j) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth); hold on;
        %plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
    end
    for j = 3:4
        phs(j) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth,'LineStyle',':'); hold on;
        %plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
    end

    xlim([-0.2 0.6]);
    ylim([-1 15]);
    area(tWindow, [axs{i}.YLim(1) axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{i}.YLim(2) axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend(phs,{'Early-Fast','Late-Fast','Early-Slow','Late-Slow'},'Box','off','Location','NorthWest');
    
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ax = gca;
    
    colororder(ax,condColours([1 2],:));

    text(axs{i}.XLim(1),axs{i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
end

for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end

% print(fullfile(outputFolder,'erps_mediansplitonrtadjust_inc.tiff'),'-dtiff','-r600');


% Make the figure (Correct Only)

% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.2 0.08], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

sElectrode = 'FCz'; % electrode of interest
iElectrode = eeg_chaninds(EEG,sElectrode);
tWindow = [0.228, 0.304]; % time window of Rewp, in s 50% of peak
tWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)

iWindow = dsearchn(EEG.unfold.times',tWindow');

whichBs = [2 11; 5 14; 8 17]; 
fbColours = [];
for i = 1:3
    axs{3+i} = subplot(2,3,3+i);
    phs = [];
    phs(1) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,1)),'LineWidth',lineWidth,'Color',condColours(3,:)); hold on;
    phs(2) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,2)),'LineWidth',lineWidth,'LineStyle',':','Color',condColours(3,:));
    
    area(tWindow, [axs{3+i}.YLim(1) axs{3+i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{3+i}.YLim(2) axs{3+i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

    legend(phs,{'Fast','Slow'},'Box','off','Location','NorthWest');
    xlim([-0.2 0.6]);
    ylim([-1 15]);
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ax = gca;
   
    text(axs{3+i}.XLim(1),axs{3+i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
end

for i = 1:length(axs)
    axs{3+i}.Box = 'off';
    axs{3+i}.FontSize = fontSize;
end

vartype = 'rt_cat'; %rt_cat rt_adjust_cat

if vartype ==  "rt_adjust_cat"
    print(fullfile(outputFolder,'erps_mediansplitonrtadjust.tiff'),'-dtiff','-r600');
elseif vartype ==  "rt_cat"
     print(fullfile(outputFolder,'erps_mediansplitonrt.tiff'),'-dtiff','-r600');
end


%% plot difference wave for fast minus slow
for i = 1:length(EEG.unfold.eventtypes)
    disp(EEG.unfold.eventtypes{i});
end

whichBs = [1 3 10 12; 4 6 13 15; 7 9 16 18]; %3conditions * 4 feedbacks
EEG.unfold.eventtypes{whichBs(1,:)}

axs = {};
fbColours = [];
for i = 1:3
    axs{i} = subplot(2,3,i);
    phs = [];
    phs(1) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,3))-meanBetasDC(iElectrode,:,whichBs(i,1)),'LineWidth',lineWidth); hold on;
        %plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
    phs(2) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,2))-meanBetasDC(iElectrode,:,whichBs(i,4)),'LineWidth',lineWidth); hold on;
        %plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);

    xlim([-0.2 0.6]);
    ylim([-5 10]);
    area(tWindow, [axs{i}.YLim(1) axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{i}.YLim(2) axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend(phs,{'earlySlow - earlyFast','lateFast - lateSlow'},'Box','off','Location','NorthWest');
    title("correct minus incorrect trials");
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ax = gca;
    
    colororder(ax,condColours([1 2],:));

    text(axs{i}.XLim(1),axs{i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
end


whichBs = [2 11; 5 14; 8 17]; 
EEG.unfold.eventtypes{whichBs(1,:)}
fbColours = [];
for i = 1:3
    axs{3+i} = subplot(2,3,3+i);
    phs = [];
    phs(1) = plot(EEG.unfold.times,(meanBetasDC(iElectrode,:,whichBs(i,2))-meanBetasDC(iElectrode,:,whichBs(i,1))),'LineWidth',lineWidth,'Color',condColours(3,:)); hold on;
    
    area(tWindow, [axs{3+i}.YLim(1) axs{3+i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{3+i}.YLim(2) axs{3+i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

    legend(phs,{'fast-slow'},'Box','off','Location','NorthWest');
    xlim([-0.2 0.6]);
    ylim([-5 10]);
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ax = gca;
   
    text(axs{3+i}.XLim(1),axs{3+i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
end
print(fullfile(outputFolder,'erps_mediansplitonrtadjust_diff.tiff'),'-dtiff','-r600');