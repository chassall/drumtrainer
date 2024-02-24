%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% Load data for sample participant (sub-01)

close all; clear all; clc;

% Set data folder - change as needed
dataFolder = '/Users/rh/Documents/ds004152/';

% Load the EEG
subName = 'sub-01';
eegFolder = [dataFolder '/derivatives/eegprep/' subName];
eegFile = [subName '_task-drumtrainer_eegprep.mat'];
load(fullfile(eegFolder,eegFile),'EEG');

% Load the residual EEG;
resFolder = [dataFolder '/derivatives/eegresidual/' subName];
resFile = [subName '_task-drumtrainer_eegresidual.mat'];
load(fullfile(resFolder,resFile),'rEEG');

% Load the betas
betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
betaFile = [subName '_task-drumtrainer_eegbeta.mat'];
S = load(fullfile(betaFolder,betaFile),'EEG');
betaEEG = S.EEG;
beta = reshape(betaEEG.unfold.beta_dc,size(betaEEG.unfold.beta_dc,1),size(betaEEG.unfold.beta_dc,2)*size(betaEEG.unfold.beta_dc,3));
pEEG = (betaEEG.unfold.Xdc*beta')'; % Predicted EEG

% Analysis settings
sElectrode = 'FCz'; % electrode of interest
iElectrode = eeg_chaninds(rEEG,sElectrode);
rewpWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)
interval = [-0.2, 0.6]; % epoch time, in s
baseline = [-40 0];

% Get response times
rTimes = [];
for i = 1:length(rEEG.event)
    if any(strcmp(rEEG.event(i).type,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'}))
        rTimes = [rTimes; rEEG.event(i).latency];
    end
end
rTimes = rTimes / EEG.srate; % Convert from samples to time in seconds
%% Display a few seconds of raw EEG, predicted EEG, and residual EEG

% Plot settings
lineColours = cbrewer('qual','Dark2',4);
lineWidth = 1;
fontSize = 7;
analysisWindow = [0.250 0.340];

outputFolder = '/Users/rh/Documents/ds004152/output/';

makefigure(8,2);
times = (1:EEG.pnts)/EEG.srate;
p1 = plot(times,movmean(EEG.data(iElectrode,:),20),'LineWidth',lineWidth,'Color',lineColours(1,:)); hold on;
p2 = plot(times,movmean(pEEG(iElectrode,:),20),'LineWidth',lineWidth,'Color',lineColours(3,:));
p3 = plot(times,movmean(rEEG.data(iElectrode,:),20),'LineWidth',lineWidth,'Color',lineColours(2,:),'LineStyle',':');
xlabel('Time (s)');
ylabel('Voltage (\muV');

ax = gca;
ax.Box = 'off';
ax.XLim =  [498  502.5]; % In seconds

% Responses
whichTimes = rTimes(rTimes >= ax.XLim(1) & rTimes <= ax.XLim(2));
% xlines = xline(whichTimes);

% Shade in analysis areas
for i = 1:length(whichTimes)
    thisTimeWindow = whichTimes(i) + analysisWindow;
    area(thisTimeWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none','ShowBaseLine','off');
    area(thisTimeWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none','ShowBaseLine','off');
end

legend('Raw','Prediction','Residual','Box','off','Location','EastOutside');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.FontSize = fontSize;
print(fullfile(outputFolder,'eegwithresidual.tiff'),'-dtiff','-r600');