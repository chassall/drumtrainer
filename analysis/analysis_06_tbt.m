% Load residual EEG and compute trial-by-trial RewP scores for the Drum
% Trainer task

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';
end

% Loop through participants
for p = 1:length(ps)

    % Load the residual EEG
    subName = ['sub-' ps{p}];
    resFolder = [dataFolder '/derivatives/eegresidual/' subName];
    resFile = [subName '_task-drumtrainer_eegresidual.mat'];
    load(fullfile(resFolder,resFile),'rEEG');

    % Analysis settings
    sElectrode = 'FCz'; % electrode of interest
    iElectrode = eeg_chaninds(rEEG,sElectrode);
    rewpWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)
    interval = [-0.2, 0.6]; % epoch time, in s
    baseline = [-40 0];

    % Make ERPs with the residual EEG
    [~,~,times,rerpdata, isRERPArtifact] = make_erp(rEEG,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    iRewpWindow =  dsearchn(times',rewpWindow');
    rERPArtifact = isRERPArtifact(iElectrode,:)';

    % Define a RewP score
    rERPRewp = squeeze(mean(rerpdata(iElectrode,iRewpWindow(1):iRewpWindow(2),:),2));

    % * * * Now do the same for the regular EEG * * *

     % Load the preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG'); 

    % Make ERPs
    EEG.event = rEEG.event; % Swap event names to readable version
    EEG.urevent = rEEG.urevent;
    interval = [-0.2, 0.6]; % epoch time, in s
    baseline = [-40 0]; % baseline, in ms
    [~,~,times,erpdata, isERPArtifact] = make_erp(EEG,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    erpArtifact = isERPArtifact(iElectrode,:)';

    % Define a RewP score
    erpRewp = squeeze(mean(erpdata(iElectrode,iRewpWindow(1):iRewpWindow(2),:),2));

    % Save 
    tbtFolder = [dataFolder '/derivatives/erptbt/' subName];
    tbtFile = [subName '_task-drumtrainer_erptbt.mat'];
    if ~exist(tbtFolder,'dir')
        mkdir(tbtFolder);
    end
    save(fullfile(tbtFolder,tbtFile),'erpRewp','erpArtifact','rERPRewp','rERPArtifact');
end

return;

%% Load trial-by-trial data and format for further analysis in R

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';

% Loop through participants
for p = 1:length(ps)

    % Load trial-by-trial scores 
    subName = ['sub-' ps{p}];
    tbtFolder = [dataFolder '/derivatives/erptbt/' subName];
    tbtFile = [subName '_task-drumtrainer_erptbt.mat'];
    load(fullfile(tbtFolder,tbtFile),'erpRewp','erpArtifact','rERPRewp','rERPArtifact');

    % Load behavioural data
    rawFile = [subName '_task-drumtrainer_beh.tsv'];
    beh = readtable(fullfile(dataFolder,subName,'beh',rawFile),'FileType','text');
    
    % * * * FORMAT DATA HERE * * *
end