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
collapseFB = 0; % 0 = all fb separate, 1 = collapse all, 2 = collapse incorrect only
collapseT = 0; % 0 = separate tempos, 1 = collapse across tempos

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';
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

% Loop through participants
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
                EEG = uf_designmat(EEG,'eventtypes',{'fast','medium','slow','redX'},'formula',{'y ~ 1','y ~ 1','y ~ 1','y ~ 1'});
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

    % Epoch with Unfold
    EEG = uf_epoch(EEG,struct('winrej',winrej,'timelimits',ufTime));

    % Make traditional ERPs with Unfold
    EEG = uf_glmfit_nodc(EEG);

    % Save this participant's EEG, which is now epoched and contains Unfold
    % data
    ufFolder = [dataFolder '/derivatives/eegbeta/' subName];
    if ~exist(ufFolder,'dir')
        mkdir(ufFolder);
    end
    ufFile = [subName '_task-drumtrainer_eegbeta.mat'];
    save(fullfile(ufFolder,ufFile),'EEG');
    
    % Compute and save the residual
    beta = reshape(EEG.unfold.beta_dc,size(EEG.unfold.beta_dc,1),size(EEG.unfold.beta_dc,2)*size(EEG.unfold.beta_dc,3));
    rEEG.data = rEEG.data - (EEG.unfold.Xdc*beta')';
    resFolder = [dataFolder '/derivatives/eegresidual/' subName];
    if ~exist(resFolder,'dir')
        mkdir(resFolder);
    end
    resFile = [subName '_task-drumtrainer_eegresidual.mat'];
    save(fullfile(resFolder,resFile),'rEEG');


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

if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data'
    outputFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\analysis\output'
else
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';
    outputFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/analysis/output';
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
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';
end

sElectrode = 'FCz'; % electrode of interest
tWindow = [0.228, 0.304]; % time window of Rewp, in s 50% of peak
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
iWindow = dsearchn(EEG.unfold.times',tWindow');
meanBetasDC = meanBetasDC - mean(meanBetasDC(:,baselineI(1):baselineI(2),:),2);
meanBetasNoDC = meanBetasNoDC - mean(meanBetasNoDC(:,baselineI(1):baselineI(2),:),2);

%% Make the figure

% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.05 0.08], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

close all; 
makefigure(18,12);
whichBs = [1 2 3; 4 5 6; 7 8 9];
axs = {};
iElectrode = eeg_chaninds(EEG,sElectrode);
fbColours = [];
for i = 1:3
    axs{i} = subplot(3,3,i);
    phs = [];
    for j = 1:3
        phs(j) = plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth); hold on;
        plot(EEG.unfold.times,meanBetasNoDC(iElectrode,:,whichBs(i,j)),'LineWidth',lineWidth/2,'Color',[0.8 0.8 0.8]);
    end
    legend(phs,{'Early','On time','Late'},'Box','off','Location','NorthWest');
    xlim([-0.2 0.6]);
    ylim([-1 15]);
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ax = gca;
    colororder(ax,condColours(1:3,:));

    text(axs{i}.XLim(1),axs{i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
end

for i = 1:3
    axs{3+i} = subplot(3,3,3+i); 
    plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,2))-meanBetasDC(iElectrode,:,whichBs(i,1)),'LineWidth',lineWidth); hold on;
    plot(EEG.unfold.times,meanBetasDC(iElectrode,:,whichBs(i,2))-meanBetasDC(iElectrode,:,whichBs(i,3)),'LineWidth',lineWidth); hold on;
    xlim([-0.2 0.6]);
    ylim([-4 8]);
    area(tWindow, [axs{3+i}.YLim(1) axs{3+i}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    area(tWindow, [axs{3+i}.YLim(2) axs{3+i}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
    legend('Early','Late','Box','off','Location','NorthWest');
    xlabel('Time (s)');
    ylabel('Voltage (\muV)');
    ax = gca;
    colororder(ax,condColours(4:5,:)); 
    text(axs{3+i}.XLim(1),axs{3+i}.YLim(2),['  ' sElectrode],'FontSize',fontSize);
end


topos = squeeze(mean(meanBetasDC(:,iWindow(1):iWindow(2),:),2));
topoMap = brewermap(128,'RdBu');
topoMap = flip(topoMap);
for i = 1:3
    axs{6+i} = subplot(3,3,6+i);

    meanIncTopo = mean(topos(:,whichBs(i,[1 3])),2);
    tp = topoplot(topos(:,whichBs(i,2))- meanIncTopo,EEG.chanlocs,'maplimits',[-3 3],'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
    tp.Parent.XLim = [-0.6 0.6];
    tp.Parent.YLim = [-0.6 0.6];
    colormap(topoMap);
end

for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end

print(fullfile(outputFolder,'erps.tiff'),'-dtiff','-r600');