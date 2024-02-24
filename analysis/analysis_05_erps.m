%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% Make ERPs and regression-ERPs for the Drum Trainer project
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
collapseFB = 1; % 0 = all fb separate, 1 = collapse all, 2 = collapse incorrect only
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

%% Loop through participants to get ERPs
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
                        if ismember(EEG.event(i).type,[earlyTriggersF lateTriggersF])
                            EEG.event(i).type = 'fastIncorrect';
                        elseif ismember(EEG.event(i).type,[onTimeTriggersF])
                            EEG.event(i).type = 'fastCorrect';
                        elseif ismember(EEG.event(i).type,[earlyTriggersM lateTriggersM])
                            EEG.event(i).type = 'medIncorrect';
                        elseif ismember(EEG.event(i).type,[onTimeTriggersM])
                            EEG.event(i).type = 'medCorrect';
                        elseif ismember(EEG.event(i).type,[earlyTriggersS lateTriggersS])
                            EEG.event(i).type = 'slowIncorrect';
                        elseif ismember(EEG.event(i).type,[onTimeTriggersS])
                            EEG.event(i).type = 'slowCorrect';
                        end
                    case 2
                        error('flag combo not supported');
                end
            elseif ~collapseFB && collapseT
                if ismember(EEG.event(i).type,[earlyTriggers])
                    EEG.event(i).type = 'early';
                elseif ismember(EEG.event(i).type,[ onTimeTriggers ])
                    EEG.event(i).type = 'onTime';
                elseif ismember(EEG.event(i).type,[ lateTriggers])
                    EEG.event(i).type = 'late';
                end

            elseif ~collapseFB && ~collapseT
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
                EEG = uf_designmat(EEG,'eventtypes',{'fastIncorrect','fastCorrect','medIncorrect','medCorrect','slowIncorrect','slowCorrect','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
            case 2
                error('flag combo not supported');
        end
    elseif ~collapseFB && collapseT
        EEG = uf_designmat(EEG,'eventtypes',{'early','onTime','late','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
    elseif ~collapseFB && ~collapseT
        EEG = uf_designmat(EEG,'eventtypes',{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly',...
            'slowOnTime','slowLate','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
    end
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
    ufFile = [subName '_' num2str(collapseFB) num2str(collapseT) '_task-drumtrainer_eegbeta.mat'];
    save(fullfile(ufFolder,ufFile),'EEG');
    
    % Compute and save the residual
        %rEEG
    resEEG = rEEG;
    beta_rEEG = reshape(EEG.unfold.beta_dc,size(EEG.unfold.beta_dc,1),size(EEG.unfold.beta_dc,2)*size(EEG.unfold.beta_dc,3));
    resEEG.data = rEEG.data - (EEG.unfold.Xdc*beta_rEEG')';
    resFolder = [dataFolder '/derivatives/eegresidual/' subName];
    if ~exist(resFolder,'dir')
        mkdir(resFolder);
    end
    resFile = [subName '_' num2str(collapseFB) num2str(collapseT) '_task-drumtrainer_eegresidual_rerp.mat'];
    save(fullfile(resFolder,resFile),'resEEG');
    
        %traditional ERP
    resEEG = rEEG;
    %EEG.unfold.X is a one hot matrix denoting trial type. repeat this for
    %each electrode and event.
    sizeX = size(EEGall.unfold.X);
    sizeB = size(EEG.unfold.beta_nodc);
    beta_traditional = reshape(EEG.unfold.beta_nodc, [sizeB(1)*sizeB(2), sizeB(3)]);

    % Perform matrix multiplication
    common = reshape(beta_traditional * EEGall.unfold.X', [sizeB(1), sizeB(2), sizeX(1)]);
    resEEG_traditional = EEGall.data - common;
    resEEG.data = resEEG_traditional;
    %plot(1:750,resEEG_traditional(1,:,1)); hold on; plot(1:750,EEG.data(1,:,1));

    resFile = [subName '_' num2str(collapseFB) num2str(collapseT) '_task-drumtrainer_eegresidual_traditional.mat'];
    save(fullfile(resFolder,resFile),'resEEG');

%     figure();
%     whichChannel = 6;
%     plot(squeeze(EEG.unfold.beta_dc(whichChannel,:,:)));
    
%     ufresult = uf_condense(EEG);
%     ufresultE = uf_condense(eEEG);
%     cfg = [];
%     cfg.channel = 6;
%     ax = uf_plotParam(ufresult,cfg);

%     figure();
%     %first plot the deconvoluted betas
%     g = uf_plotParam(ufresult,'channel',6,'deconv',1,'baseline',[-0.2 0]);
% 
%     g = uf_plotParam(ufresultE,'channel',6,'deconv',0,'baseline',[-0.2 0],'gramm',g);

%     % Solve GLM via lsmr
%     % Requires lsmr
%     % Set bad samples to 0 in the design matrix
%     for i = 1:size(winrej,1)
%         X(winrej(i,1):winrej(i,2),:) = 0;
%         data(:,winrej(i,1):winrej(i,2)) = 0;
%     end
%     whichMethod = 'lsmr';
%     lsmriterations = 400;
%     beta = [];
%     for e = 1:EEG.nbchan
%         [beta(:,e),ISTOP,ITN] = lsmr(X,double(data(e,:)'),[],10^-8,10^-8,[],lsmriterations);
%     end
    
%     figure();
%     n = length(EEG.unfold.times);
%     for i = 1:4
%         theseSamples = (1+(i-1)*n):i*n;
%         subplot(1,4,i);
%         plot(EEG.unfold.times,beta(theseSamples,whichChannel)); hold on;
%         plot(EEG.unfold.times,EEG.unfold.beta_dc(whichChannel,:,i));
%         legend('lsmr','unfold');
%     end
%     % Compute the residual
%     rEEG = EEG;
%     rEEG.data = EEG.data - (X*beta)';

%     pause();

end

return;

%% Load and combine the beta weights

close all; clear all; % Variables will be loaded below

collapseFB = 1; % 0 = all fb separate, 1 = collapse all, 2 = collapse incorrect only
collapseT = 0; % 0 = separate tempos, 1 = collapse across tempos

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
allBetasDCcollapsedFB = [];
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
    
    betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
    betaFile = [subName '_' num2str(collapseFB) num2str(collapseT) '_task-drumtrainer_eegbeta.mat'];
    load(fullfile(betaFolder,betaFile),'EEG');
    allBetasDCcollapsedFB(p,:,:,:) = EEG.unfold.beta_dc;

end

% Reminder about order of betas
for i = 1:length(EEG.unfold.eventtypes)
    disp(EEG.unfold.eventtypes{i});
end

whichBs = [1 2 3; 4 5 6; 7 8 9];

meanBetasDC = squeeze(mean(allBetasDC,1));
meanBetasNoDC = squeeze(mean(allBetasNoDC,1));
meanBetasDCcollapsedFB = squeeze(mean(allBetasDCcollapsedFB,1));

baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');
iWindow = dsearchn(EEG.unfold.times',tWindow');

meanBetasDC = meanBetasDC - mean(meanBetasDC(:,baselineI(1):baselineI(2),:),2);
meanBetasNoDC = meanBetasNoDC - mean(meanBetasNoDC(:,baselineI(1):baselineI(2),:),2);
meanBetasDCcollapsedFB = meanBetasDCcollapsedFB - mean(meanBetasDCcollapsedFB(:,baselineI(1):baselineI(2),:),2);

diffBetasDCearly = meanBetasDC(:,:,whichBs(:,2))-meanBetasDC(:,:,whichBs(:,1));
diffBetasDClate = meanBetasDC(:,:,whichBs(:,2))-meanBetasDC(:,:,whichBs(:,3));
% compute SE for plotting the shadings
BetaDC_se = 2.093 * squeeze(std(allBetasDC,0,1,'omitnan'))./sqrt(20);

diffBetasnoDCearly = meanBetasNoDC(:,:,whichBs(:,2))-meanBetasNoDC(:,:,whichBs(:,1));
diffBetasnoDClate = meanBetasNoDC(:,:,whichBs(:,2))-meanBetasNoDC(:,:,whichBs(:,3));

diffBetasDCcollapsedFast = meanBetasDCcollapsedFB(:,:,2) - meanBetasDCcollapsedFB(:,:,1);
diffBetasDCcollapsedMedium = meanBetasDCcollapsedFB(:,:,4) - meanBetasDCcollapsedFB(:,:,3);
diffBetasDCcollapsedSlow = meanBetasDCcollapsedFB(:,:,6) - meanBetasDCcollapsedFB(:,:,5);

diffBetasDCcollapsedFBfast_se = squeeze(std(allBetasDCcollapsedFB(:,:,:,2) - allBetasDCcollapsedFB(:,:,:,1),0,1,'omitnan'))./sqrt(20);
diffBetasDCcollapsedFBmedium_se = squeeze(std(allBetasDCcollapsedFB(:,:,:,4) - allBetasDCcollapsedFB(:,:,:,3),0,1,'omitnan'))./sqrt(20);
diffBetasDCcollapsedFBslow_se = squeeze(std(allBetasDCcollapsedFB(:,:,:,6) - allBetasDCcollapsedFB(:,:,:,5),0,1,'omitnan'))./sqrt(20);

diffBetasDCearly_se = squeeze(std(allBetasDC(:,:,:,whichBs(:,2)) - allBetasDC(:,:,:,whichBs(:,1)),0,1,'omitnan'))./sqrt(20);
diffBetasDClate_se = squeeze(std(allBetasDC(:,:,:,whichBs(:,2)) - allBetasDC(:,:,:,whichBs(:,3)),0,1,'omitnan'))./sqrt(20);

BetanoDC_se = 2.093 * squeeze(std(allBetasNoDC,0,1,'omitnan'))./sqrt(20);

diffBetasnoDCearly_se = squeeze(std(allBetasNoDC(:,:,:,whichBs(:,2)) - allBetasNoDC(:,:,:,whichBs(:,1)),0,1,'omitnan'))./sqrt(20);
diffBetasnoDClate_se = squeeze(std(allBetasNoDC(:,:,:,whichBs(:,2)) - allBetasNoDC(:,:,:,whichBs(:,3)),0,1,'omitnan'))./sqrt(20);

%% Make the figure (DC)

% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.05 0.08], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

% purple4 #551A8B plum3 #CD96CD gold #FFD700

purpleGold = [85, 26, 139; 255, 215, 0; 205, 150, 205;]/255;

close all; 
makefigure(20,16);
whichBs = [1 2 3; 4 5 6; 7 8 9];
axs = {};
iElectrode = eeg_chaninds(EEG,sElectrode);

% Define the letters to use as annotations
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
tempos_labs = ["Fast","Medium","Slow"];

for i = 1:3
    axs{i} = subplot(3,3,i);
    phs = [];
    for j = 1:3
        phs(j) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth,...
        'Color',purpleGold(j,:)); hold on;
%         ciplot(squeeze(meanBetasDC(iElectrode,:,whichBs(i,j))-BetaDC_se(iElectrode,:,whichBs(i,j))),...
%             meanBetasDC(iElectrode,:,whichBs(i,j))+squeeze(BetaDC_se(iElectrode,:,whichBs(i,j))),...
%             EEG.unfold.times,purpleGold(j,:),0.15);
%         plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
    end
    xlim([-0.2 0.6]);
    ylim([-1 15]);
    area(tWindow, [axs{i}.YLim(1) axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{i}.YLim(2) axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend(phs,{'Early','On time','Late'},'Box','off','Location','NorthWest');
    xlabel('Time (s)');
    ylabel('Regression Coef');
    title(tempos_labs(i));
    ax = gca;
%     colororder(ax,purpleGold);
    letter = letters(i);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    text(axs{i}.XLim(1)-0.2, axs{i}.YLim(2)-0.2, letter, 'FontSize', 14);
end

for i = 1:3
    axs{3+i} = subplot(3,3,3+i); 
    plot(EEG.unfold.times,diffBetasDCearly(iElectrode,:,i),'LineWidth',lineWidth, 'Color',purpleGold(1,:)); hold on;
    plot(EEG.unfold.times,diffBetasDClate(iElectrode,:,i),'LineWidth',lineWidth,'Color',purpleGold(3,:)); hold on;
    ciplot(diffBetasDCearly(iElectrode,:,i)-diffBetasDCearly_se(iElectrode,i),...
            diffBetasDCearly(iElectrode,:,i)+diffBetasDCearly_se(iElectrode,i),...
            EEG.unfold.times,purpleGold(1,:),0.2);
            ciplot(diffBetasDClate(iElectrode,:,i)-diffBetasDClate_se(iElectrode,i),...
            diffBetasDClate(iElectrode,:,i)+diffBetasDClate_se(iElectrode,i),...
            EEG.unfold.times,purpleGold(3,:),0.2);
    xlim([-0.2 0.6]);
    ylim([-4 6]);
    area(tWindow, [axs{3+i}.YLim(1) axs{3+i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{3+i}.YLim(2) axs{3+i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend('On time - Early','On time - Late','Box','off','Location','NorthWest');
    xlabel('Time (s)');
    title(tempos_labs(i));
    ylabel('Regression Coef');
    ax = gca;
%     colororder(ax,purpleGold([1 3],:)); 
    letter = letters(i+3);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    text(axs{3+i}.XLim(1)-0.2, axs{3+i}.YLim(2)-0.2, letter, 'FontSize', 14);
end


topos = squeeze(mean(meanBetasDC(:,iWindow(1):iWindow(2),:),2));
topoMap = brewermap(128,'PuOr');
topoMap = flip(topoMap);
for i = 1:3
    axs{6+i} = subplot(3,3,6+i);

    meanIncTopo = mean(topos(:,whichBs(i,[1 3])),2);
    tp = topoplot(topos(:,whichBs(i,2))- meanIncTopo,EEG.chanlocs,'maplimits',[-3 3],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
    tp.Parent.XLim = [-0.6 0.6];
    tp.Parent.YLim = [-0.6 0.6];
    colormap(topoMap);
    letter = letters(i+6);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    title(tempos_labs(i));
    text(axs{6+i}.XLim(1)-0.2, axs{6+i}.YLim(2), letter, 'FontSize', 14);
    c = colorbar;
    c.Label.String = 'amplitude (regression coef.)';
end

for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end

print(fullfile(outputFolder,'erps.tiff'),'-dtiff','-r600');


%% Make the figure (DC, with shading from permutation)
M = readtable([outputFolder,'cluster.csv']);
M.Start = M.Start./1000;
M.End = M.End./1000;
% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.05 0.08], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

% purple4 #551A8B plum3 #CD96CD gold #FFD700

close all; 
makefigure(20,6);
whichBs = [1 2 3; 4 5 6; 7 8 9];
axs = {};
iElectrode = eeg_chaninds(EEG,sElectrode);

% Define the letters to use as annotations
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
tempos_labs = ["Fast","Medium","Slow"];

for i = 1:3
    axs{i} = subplot(1,3,i);
    phs = [];
    for k = 1:3
        phs(k) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,k)),'LineWidth',lineWidth,'Color',purpleGold(k,:)); hold on;
%         ciplot(squeeze(meanBetasDC(iElectrode,:,whichBs(i,k))-BetaDC_se(iElectrode,:,whichBs(i,k))),...
%             meanBetasDC(iElectrode,:,whichBs(i,k))+squeeze(BetaDC_se(iElectrode,:,whichBs(i,k))),...
%             EEG.unfold.times,purpleGold(k,:),0.15);
    end
    
    if i == 3
        for j = 1:3
            time_id = (EEG.unfold.times >= M.Start(j) & EEG.unfold.times <= M.End(j));
%             area([M.Start(j) M.End(j)], [axs{i}.YLim(1) axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
%             area([M.Start(j) M.End(j)], [axs{i}.YLim(2) axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
            plot(EEG.unfold.times(time_id),meanBetasDC(iElectrode,time_id,whichBs(3,2)),...
                'lineWidth',4,'Color',purpleGold(2,:));
            plot(EEG.unfold.times(time_id),zeros(length(EEG.unfold.times(time_id)),1),...
                'lineWidth',1,'Color',purpleGold(2,:));
            time_1id = find(time_id== 1);
            text(EEG.unfold.times(time_1id(1)),0.5, ['p = ' num2str(round(M.cluster_p_values(j),2))],'FontSize',6);
        end
    end
    xlim([-0.2 0.6]);
    ylim([-1 15]);
    area(tWindow, [axs{i}.YLim(1) axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{i}.YLim(2) axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend(phs,{'Early','On time','Late'},'Box','off','Location','NorthWest');
    xlabel('Time (s)');
    ylabel('Regression Coef');
    title(tempos_labs(i)); 
    ax = gca;
    letter = letters(i);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    text(axs{i}.XLim(1)-0.2, axs{i}.YLim(2)-0.2, letter, 'FontSize', 14);
end

for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end

print(fullfile(outputFolder,'erps_with_shading.tiff'),'-dtiff','-r600');

%% Make the figure (noDC)

% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.05 0.08], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

% purple4 #551A8B plum3 #CD96CD gold #FFD700

purpleGold = [85, 26, 139; 255, 215, 0; 205, 150, 205;]/255;

close all; 
makefigure(20,16);
whichBs = [1 2 3; 4 5 6; 7 8 9];
axs = {};
iElectrode = eeg_chaninds(EEG,sElectrode);

% Define the letters to use as annotations
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

for i = 1:3
    axs{i} = subplot(3,3,i);
    phs = [];
    for j = 1:3
        phs(j) = plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth); hold on;
%         plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
    end
    xlim([-0.2 0.6]);
    ylim([-3 10]);
    area(tWindow, [axs{i}.YLim(1) axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{i}.YLim(2) axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend(phs,{'Early','On time','Late'},'Box','off','Location','NorthWest');
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    title(tempos_labs(i));
    ax = gca;
    colororder(ax,purpleGold);
    letter = letters(i);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    text(axs{i}.XLim(1)-0.2, axs{i}.YLim(2)-0.2, letter, 'FontSize', 14);
end

for i = 1:3
    axs{3+i} = subplot(3,3,3+i); 
    plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,2))-meanBetasNoDC(iElectrode,:,whichBs(i,1)),...
        'LineWidth',lineWidth,'Color',purpleGold(1,:)); hold on;
    plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,2))-meanBetasNoDC(iElectrode,:,whichBs(i,3)),...
        'LineWidth',lineWidth,'Color',purpleGold(3,:)); hold on;
    ciplot(diffBetasnoDCearly(iElectrode,:,i)-diffBetasnoDCearly_se(iElectrode,i),...
            diffBetasnoDCearly(iElectrode,:,i)+diffBetasnoDCearly_se(iElectrode,i),...
            EEG.unfold.times,purpleGold(1,:),0.2);
            ciplot(diffBetasnoDClate(iElectrode,:,i)-diffBetasnoDClate_se(iElectrode,i),...
            diffBetasnoDClate(iElectrode,:,i)+diffBetasnoDClate_se(iElectrode,i),...
            EEG.unfold.times,purpleGold(3,:),0.2);
    xlim([-0.2 0.6]);
    ylim([-4 6]);
    area(tWindow, [axs{3+i}.YLim(1) axs{3+i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{3+i}.YLim(2) axs{3+i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend('On time - Early','On time - Late','Box','off','Location','NorthWest');
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    title(tempos_labs(i));
    ax = gca;
    letter = letters(i+3);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    text(axs{3+i}.XLim(1)-0.2, axs{3+i}.YLim(2)-0.2, letter, 'FontSize', 14);
end


topos = squeeze(mean(meanBetasNoDC(:,iWindow(1):iWindow(2),:),2));
topoMap = brewermap(128,'PuOr');
topoMap = flip(topoMap);
for i = 1:3
    axs{6+i} = subplot(3,3,6+i);

    meanIncTopo = mean(topos(:,whichBs(i,[1 3])),2);
    tp = topoplot(topos(:,whichBs(i,2))- meanIncTopo,EEG.chanlocs,'maplimits',[-3 3],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
    tp.Parent.XLim = [-0.6 0.6];
    tp.Parent.YLim = [-0.6 0.6];
    colormap(topoMap);
    letter = letters(i+6);
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    title(tempos_labs(i));
    text(axs{6+i}.XLim(1)-0.2, axs{6+i}.YLim(2), letter, 'FontSize', 14);
    c = colorbar;
    c.Label.String = 'Voltage (\muV)';
end

for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end

print(fullfile(outputFolder,'erps_noDC.tiff'),'-dtiff','-r600');
%% difference wave for all individuals (by feedback)

close all; 
makefigure(20,14);

baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');
iWindow = dsearchn(EEG.unfold.times',tWindow');
allBetasDC = allBetasDC - mean(allBetasDC(:,:,baselineI(1):baselineI(2),:),3);

allDiffWave = [];
meanDiffWave = [];

for i = 1:3
    allDiffWave(i,1,:,:) = squeeze(allBetasDC(:,iElectrode,:,whichBs(i,2))- allBetasDC(:,iElectrode,:,whichBs(i,1)));
    allDiffWave(i,2,:,:) = squeeze(allBetasDC(:,iElectrode,:,whichBs(i,2))- allBetasDC(:,iElectrode,:,whichBs(i,3)));
    meanDiffWave(i,1,:) = squeeze(meanBetasDC(iElectrode,:,whichBs(i,2))- meanBetasDC(iElectrode,:,whichBs(i,1)));
    meanDiffWave(i,2,:) = squeeze(meanBetasDC(iElectrode,:,whichBs(i,2))- meanBetasDC(iElectrode,:,whichBs(i,3)));
end


colors = ['r','b'];
for k = 1:2 %early,late
    colr = colors(k);
    for i = 1:3
        for j = 1:20
            axs{(k-1)*3 + i} = subplot(2,3,(k-1)*3 + i); 
            p1 = plot(EEG.unfold.times,squeeze(allDiffWave(i,k,j,:)),'LineWidth',lineWidth,'Color',colr); hold on;
            p1.Color(4) = 0.2;%alpha
        end
        plot(EEG.unfold.times,squeeze(meanDiffWave(i,k,:)),'LineWidth',2,'Color',colr); hold on;
        xlim([-0.2 0.6]);
        ylim([-4 10]);
        area(tWindow, [ axs{(k-1)*3 + i}.YLim(1)  axs{(k-1)*3 + i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
        area(tWindow, [ axs{(k-1)*3 + i}.YLim(2)  axs{(k-1)*3 + i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
        if k == 1
            xlabel('Time (s)');
        end
        ylabel('Voltage (\muV)');
        ax = gca;
        colororder(ax,condColours(4:5,:)); 
    end
end


print(fullfile(outputFolder,'erps_ind.tiff'),'-dtiff','-r600');

%% difference wave for all individuals (collapsed FB)
close all; 
makefigure(20,6);

baselineT = [-0.04 0];
baselineI = dsearchn(EEG.unfold.times',baselineT');
iWindow = dsearchn(EEG.unfold.times',tWindow');


for i = 1:3
        for j = 1:20
            axs{i} = subplot(1,3,i); 
            p1 = plot(EEG.unfold.times,mean(squeeze(allDiffWave(i,:,j,:)),1),'LineWidth',lineWidth,'Color','k'); hold on;
            p1.Color(4) = 0.4;%alpha
        end
        plot(EEG.unfold.times,squeeze(mean(squeeze(allDiffWave(i,:,:,:)),[1 2])),'LineWidth',2,'Color','r'); hold on;
        xlim([-0.2 0.6]);
        ylim([-4 8]);
        area(tWindow, [ axs{i}.YLim(1)  axs{i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
        area(tWindow, [ axs{i}.YLim(2)  axs{i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
        xlabel('Time (s)');
        ylabel('Voltage (\muV)');
        ax = gca;
        colororder(ax,condColours(4:5,:)); 
end

print(fullfile(outputFolder,'erps_ind_collapsed.tiff'),'-dtiff','-r600');

%% make the manuscript figure

% % Overwrite subplot with subtightplot for more control
% subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.05 0.08], [0.08 0.05]);
% %     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
% %              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
% %              relatively large axis. 
% %     marg_h  margins in height in normalized units (0...1)
% %              or [lower uppper] for different lower and upper margins 
% %     marg_w  margins in width in normalized units (0...1)
% %              or [left right] for different left and right margins 
% 
% close all; 
% makefigure(16,24);
% whichBs = [1 2 3; 4 5 6; 7 8 9];
% axs = {};
% iElectrode = eeg_chaninds(EEG,sElectrode);
% fbColours = [];
% for i = 1:3
%     axs{i} = subplot(5,3,i);
%     phs = [];
%     for j = 1:3
%         phs(j) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth); hold on;
%         plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
%     end
%     legend(phs,{'Early','On time','Late'},'Box','off','Location','NorthWest');
%     xlim([-0.2 0.6]);
%     ylim([-1 15]);
%     xlabel('Time (s)');
%     ylabel('Voltage (\muV)');
%     ax = gca;
%     colororder(ax,condColours(1:3,:));
% 
%     text(axs{i}.XLim(1),axs{i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
% end
% 
% for i = 1:3
%     axs{3+i} = subplot(5,3,3+i); 
%     plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,2))-meanBetasDC(iElectrode,:,whichBs(i,1)),'LineWidth',lineWidth); hold on;
%     plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,2))-meanBetasDC(iElectrode,:,whichBs(i,3)),'LineWidth',lineWidth); hold on;
%     xlim([-0.2 0.6]);
%     ylim([-4 8]);
%     area(tWindow, [axs{3+i}.YLim(1) axs{3+i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
%     area(tWindow, [axs{3+i}.YLim(2) axs{3+i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
%     legend('Early','Late','Box','off','Location','NorthWest');
%     xlabel('Time (s)');
%     ylabel('Voltage (\muV)');
%     ax = gca;
%     colororder(ax,condColours(4:5,:)); 
%     text(axs{3+i}.XLim(1),axs{3+i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
% end
% 
% 
% for i = 1:length(axs)
%     axs{i}.Box = 'off';
%     axs{i}.FontSize = fontSize;
% end
% 
% % difference wave for all individuals
% baselineT = [-0.04 0];
% baselineI = dsearchn(EEG.unfold.times',baselineT');
% iWindow = dsearchn(EEG.unfold.times',tWindow');
% allBetasDC = allBetasDC - mean(allBetasDC(:,:,baselineI(1):baselineI(2),:),3);
% 
% allDiffWave = [];
% meanDiffWave = [];
% 
% for i = 1:3
%     allDiffWave(i,1,:,:) = squeeze(allBetasDC(:,iElectrode,:,whichBs(i,2))- allBetasDC(:,iElectrode,:,whichBs(i,1)));
%     allDiffWave(i,2,:,:) = squeeze(allBetasDC(:,iElectrode,:,whichBs(i,2))- allBetasDC(:,iElectrode,:,whichBs(i,3)));
%     meanDiffWave(i,1,:) = squeeze(meanBetasDC(iElectrode,:,whichBs(i,2))- meanBetasDC(iElectrode,:,whichBs(i,1)));
%     meanDiffWave(i,2,:) = squeeze(meanBetasDC(iElectrode,:,whichBs(i,2))- meanBetasDC(iElectrode,:,whichBs(i,3)));
% end
% 
% 
% 
% topos = squeeze(mean(meanBetasDC(:,iWindow(1):iWindow(2),:),2));
% topoMap = brewermap(128,'RdBu');
% topoMap = flip(topoMap);
% for i = 1:3
%     axs{6+i} = subplot(5,3,6+i);
% 
%     meanIncTopo = mean(topos(:,whichBs(i,[1 3])),2);
%     tp = topoplot(topos(:,whichBs(i,2))- meanIncTopo,EEG.chanlocs,'maplimits',[-3 3],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
%     tp.Parent.XLim = [-0.6 0.6];
%     tp.Parent.YLim = [-0.6 0.6];
%     colormap(topoMap);
% end
% 
% colors = ['r','b'];
% 
% for k = 1:2 %early,late
%     colr = colors(k);
%     for i = 1:3
%         for j = 1:20
%             axs{9+(k-1)*3 +i} = subplot(5,3,9+(k-1)*3 +i); 
%             p1 = plot(EEG.unfold.times,squeeze(allDiffWave(i,k,j,:)),'LineWidth',lineWidth,'Color',colr); hold on;
%             p1.Color(4) = 0.2;%alpha
%         end
%         plot(EEG.unfold.times,squeeze(meanDiffWave(i,k,:)),'LineWidth',2,'Color',colr); hold on;
%         xlim([-0.2 0.6]);
%         ylim([-4 10]);
%         area(tWindow, [ axs{9+(k-1)*3 +i}.YLim(1)  axs{9+(k-1)*3 +i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
%         area(tWindow, [ axs{9+(k-1)*3 +i}.YLim(2)  axs{9+(k-1)*3 +i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
%         if k == 1
%             xlabel('Time (s)');
%         end
%         ylabel('Voltage (\muV)');
%         ax = gca;
%         colororder(ax,condColours(4:5,:)); 
%     end
% end
% 
% print(fullfile(outputFolder,'erps_manuscript.tiff'),'-dtiff','-r600');

%% topograph of difference wave
topos = squeeze(mean(meanBetasDC(:,iWindow(1):iWindow(2),:),2));
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
for i = 1:3
    axs{6+i} = subplot(5,3,6+i);

    meanIncTopo = mean(topos(:,whichBs(i,[1 3])),2);
    tp = topoplot(topos(:,whichBs(i,2))- meanIncTopo,EEG.chanlocs,'maplimits',[-3 3],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
    tp.Parent.XLim = [-0.6 0.6];
    tp.Parent.YLim = [-0.6 0.6];
    colormap(topoMap);
end
%% output csv file for analysis
rERP_file = reshape(mean(allDiffWave(:,:,:,iWindow(1):iWindow(2)),4),[6 20]);
rERP_file =  array2table(rERP_file');
rERP_file.Properties.VariableNames = {'fast_early','medium_early','slow_early','fast_late','medium_late','slow_late'};
rERP_file = [array2table(str2double(ps)'),rERP_file];
rERP_file.Properties.VariableNames(1) = {'participant_ID'};

writetable(rERP_file,[dataFolder '/R/rERP_meta.csv']);

%% plot individual traditional ERPs
allBetasNoDC = allBetasNoDC - mean(allBetasNoDC(:,:,baselineI(1):baselineI(2),:),3);

allDiffWave = [];
meanDiffWave = [];

for i = 1:3
    allDiffWave(i,1,:,:) = squeeze(allBetasNoDC(:,iElectrode,:,whichBs(i,2))- allBetasNoDC(:,iElectrode,:,whichBs(i,1)));
    allDiffWave(i,2,:,:) = squeeze(allBetasNoDC(:,iElectrode,:,whichBs(i,2))- allBetasNoDC(:,iElectrode,:,whichBs(i,3)));
    meanDiffWave(i,1,:) = squeeze(meanBetasNoDC(iElectrode,:,whichBs(i,2))- meanBetasNoDC(iElectrode,:,whichBs(i,1)));
    meanDiffWave(i,2,:) = squeeze(meanBetasNoDC(iElectrode,:,whichBs(i,2))- meanBetasNoDC(iElectrode,:,whichBs(i,3)));
end

colors = ['r','b'];
for k = 1:2 %early,late
    colr = colors(k);
    for i = 1:3
        for j = 1:20
            axs{(k-1)*3 + i} = subplot(2,3,(k-1)*3 + i); 
            p1 = plot(EEG.unfold.times,squeeze(allDiffWave(i,k,j,:)),'LineWidth',lineWidth,'Color',colr); hold on;
            p1.Color(4) = 0.1;%alpha
        end
        plot(EEG.unfold.times,squeeze(meanDiffWave(i,k,:)),'LineWidth',2,'Color',colr); hold on;
        xlim([-0.2 0.6]);
        ylim([-4 8]);
        area(tWindow, [ axs{(k-1)*3 + i}.YLim(1)  axs{(k-1)*3 + i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
        area(tWindow, [ axs{(k-1)*3 + i}.YLim(2)  axs{(k-1)*3 + i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
        xlabel('Time (s)');
        ylabel('Voltage (\muV)');
        title("traditional ERP");
        ax = gca;
        colororder(ax,condColours(4:5,:)); 
    end
end

ERP_file = reshape(mean(allDiffWave(:,:,:,iWindow(1):iWindow(2)),4),[6 20]);
ERP_file =  array2table(ERP_file');
ERP_file.Properties.VariableNames = {'fast_early','medium_early','slow_early','fast_late','medium_late','slow_late'};
ERP_file = [array2table(str2double(ps)'),ERP_file];
ERP_file.Properties.VariableNames(1) = {'participant_ID'};

writetable(ERP_file,[dataFolder '/R/traditional_ERP_meta.csv']);

%% New Figure (DC)
makefigure(15,10);
plotColours = [205,183, 158; 139,125,107; 0, 0, 0;]/255; % requires cbrewer
plot(EEG.unfold.times,diffBetasDCcollapsedFast(iElectrode,:),'LineWidth',lineWidth, 'Color',plotColours(1,:),'lineWidth',2); hold on;
plot(EEG.unfold.times,diffBetasDCcollapsedMedium(iElectrode,:),'LineWidth',lineWidth, 'Color',plotColours(2,:),'lineWidth',2); hold on;
plot(EEG.unfold.times,diffBetasDCcollapsedSlow(iElectrode,:),'LineWidth',lineWidth, 'Color',plotColours(3,:),'lineWidth',2); hold on;
xlim([-0.2 0.6]);
ylim([-4 5]);
xlabel('Time (s)');
ylabel('Voltage (\muV)');
title("Difference wave (correct - incorrect)");
ax=gca;
area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

ciplot(diffBetasDCcollapsedFast(iElectrode,:)-diffBetasDCcollapsedFBfast_se(iElectrode,:),...
    diffBetasDCcollapsedFast(iElectrode,:)+diffBetasDCcollapsedFBfast_se(iElectrode,:),...
            EEG.unfold.times,plotColours(1,:),0.2);
ciplot(diffBetasDCcollapsedMedium(iElectrode,:)-diffBetasDCcollapsedFBmedium_se(iElectrode,:),...
    diffBetasDCcollapsedMedium(iElectrode,:)+diffBetasDCcollapsedFBmedium_se(iElectrode,:),...
            EEG.unfold.times,plotColours(2,:),0.2);
ciplot(diffBetasDCcollapsedSlow(iElectrode,:)-diffBetasDCcollapsedFBslow_se(iElectrode,:),...
    diffBetasDCcollapsedSlow(iElectrode,:)+diffBetasDCcollapsedFBslow_se(iElectrode,:),...
            EEG.unfold.times,plotColours(3,:),0.2);
legend({'Fast','Medium','Slow'},'Box','off','Location','NorthWest');

print(fullfile(outputFolder,'erps_dich.tiff'),'-dtiff','-r600');
