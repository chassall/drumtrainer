%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% Load residual EEG and compute trial-by-trial RewP scores for the Drum
% Trainer task

% Author: Cameron Hassall and Ryan Yan, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk; ryany98@stanford.edu
% Website: http://www.cameronhassall.com

clear all; clc;
% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = 'E:\OneDrive - Nexus365\Projects\2021_EEG_DrumTrainer_Hassall\data';
else
    dataFolder = '/Users/rh/Documents/ds004152';
end

% Loop through participants
for p = 1:length(ps)

    % Load the residual EEG
    subName = ['sub-' ps{p}];
    resFolder = [dataFolder '/derivatives/eegresidual/' subName];
    resFile = [subName '_task-drumtrainer_eegresidual_rerp.mat'];
    resEEGr = load(fullfile(resFolder,resFile)).resEEG;
    resFile = [subName '_task-drumtrainer_eegresidual_traditional.mat'];
    resEEGt = load(fullfile(resFolder,resFile)).resEEG;

    % Analysis settings
    sElectrode = 'FCz'; % electrode of interest
    iElectrode = eeg_chaninds(resEEGr,sElectrode);
    rewpWindow = [240, 340]; % Sambrook and Goslin (2015)
%     rewpWindow = [228, 304]; % time window of Rewp, in s 50% of peak
    interval = [-0.2, 0.6]; % epoch time, in s
    baseline = [-40 0];
    tWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)
    

    % =========================================================
    % ============== make ERPs with residual rEEG =============
    % =========================================================
%     [~,~,times,rerpdata, isRERPArtifact] = make_erp(rEEG,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    [~,~,times,rerpdata, isRERPArtifact] = make_erp(resEEGr,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    iRewpWindow =  dsearchn(times',rewpWindow');
    rERPArtifact = isRERPArtifact(iElectrode,:)';

    % Define a RewP score
    rERPRewp = squeeze(mean(rerpdata(iElectrode,iRewpWindow(1):iRewpWindow(2),:),2));
    
    % =========================================================
    % ========= Now do the same for the regular EEG ===========
    % =========================================================

     % Load the preprocessed EEG
    subName = ['sub-' ps{p}];
    prepFolder = [dataFolder '/derivatives/eegprep/' subName];
    prepFile = [subName '_task-drumtrainer_eegprep.mat'];
    load(fullfile(prepFolder,prepFile),'EEG'); 
    EEG.event = resEEGr.event; % Swap event names to readable version
    EEG.urevent = resEEGr.urevent;
%     betaFolder = [dataFolder '/derivatives/eegbeta/' subName];
%     betaFile = [subName '_task-drumtrainer_eegbeta.mat'];
%     load(fullfile(betaFolder,betaFile),'EEG');
    
    % take grand mean waveform
    [~,~,times,erpdata, isERPArtifact] = make_erp(EEG,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    erpArtifact = isERPArtifact(iElectrode,:)';
    tWindow = [0.240, 0.340]; % Sambrook and Goslin (2015)
    iRewpWindow =  dsearchn(resEEGt.unfold.times',tWindow');

    % Define a RewP score
    %erpRewp = squeeze(mean(erpdata(iElectrode,iRewpWindow(1):iRewpWindow(2),:),2));
    
    %get residual EEG
    erpRewp = squeeze(mean(resEEGt.data(iElectrode,iRewpWindow(1):iRewpWindow(2),:),2));

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
dataFolder = '/Users/rh/Documents/ds004152';
tbt_file=[];
% Loop through participants
for p = 1:length(ps)

    % Load trial-by-trial scores 
    subName = ['sub-' ps{p}];
    tbtFolder = [dataFolder '/derivatives/erptbt/' subName];
    tbtFile = [subName '_task-drumtrainer_erptbt.mat'];
    thisData = load(fullfile(tbtFolder,tbtFile),'erpRewp','erpArtifact','rERPRewp','rERPArtifact');

    % Load behavioural data
    rawFile = [subName '_task-drumtrainer_beh.tsv'];
    beh = readtable(fullfile(dataFolder,subName,'beh',rawFile),'FileType','text');
    
    % * * * FORMAT DATA HERE * * *
    
    % assign impossible value to artifacts to exclude later
    thisData.erpRewp(thisData.erpArtifact) = -99; 
    thisData.rERPRewp(thisData.rERPArtifact) = -99;
   
    if size(beh,1) > 1728 %participant 1
        thisData.erpRewp = thisData.erpRewp(1:1728);
        thisData.rERPRewp = thisData.rERPRewp(1:1728);
        beh = beh(1:1728,:);
    end
    
    tbt_file = [tbt_file; repmat(str2double(ps{p}),1728,1),beh.blockLoop_thisRepN,beh.trialLoop_thisRepN,thisData.erpRewp,thisData.rERPRewp,beh.blockType,beh.outcome,...
        beh.trialResp_corr,beh.trialResp_rt];
    

end

tbt_file =  array2table(tbt_file);
tbt_file.Properties.VariableNames = {'participant_ID','blockLoop_thisRepN','trialLoop_thisRepN','erpRewp','rERPRewp','blockType','outcome',...
    'trialResp_corr','trialResp_rt'};
writetable(tbt_file,[dataFolder '/R/tbt_meta.csv']);

%% gather data for data-driven analysis (this may take a while)
data_driven_mat = zeros(length(ps),size(rerpdata,3),10,length(times)); %20 1728 10 200
    
for p = 1:length(ps)

    % Load the residual EEG
    subName = ['sub-' ps{p}];
    resFolder = [dataFolder '/derivatives/eegresidual/' subName];
    resFile = [subName '_task-drumtrainer_eegresidual_rerp.mat'];
    load(fullfile(resFolder,resFile),'resEEG');

    % Analysis settings
    sElectrode = 'FCz'; % electrode of interest
    iElectrode = eeg_chaninds(resEEG,sElectrode);
     rewpWindow = [240, 340]; % Sambrook and Goslin (2015)
%     rewpWindow = [228, 304]; % time window of Rewp, in s 50% of peak
    interval = [-0.2, 0.6]; % epoch time, in s
    baseline = [-40 0];

    % rERP
%     [~,~,times,rerpdata, isRERPArtifact] = make_erp(rEEG,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    [~,~,times,rerpdata, isRERPArtifact] = make_erp(resEEG,{'fastEarly','fastOnTime','fastLate','medEarly','medOnTime','medLate','slowEarly','slowOnTime','slowLate','redX'},interval,baseline);
    rerpdata(isRERPArtifact) = -nan;
    if size(rerpdata,3) > 1728 %participant 1 has twice as much trials
        rerpdata = rerpdata(:,:,1:1728);
    end
    
    % Load behavioural data
    rawFile = [subName '_task-drumtrainer_beh.tsv'];
    beh = readtable(fullfile(dataFolder,subName,'beh',rawFile),'FileType','text');
    
    if size(beh,1) > 1728 %participant 1
        thisData.erpRewp = thisData.erpRewp(1:1728);
        thisData.rERPRewp = thisData.rERPRewp(1:1728);
        beh = beh(1:1728,:);
    end
    

    beh_temp = [repmat(str2double(ps{p}),1728,1),beh.blockLoop_thisRepN,beh.trialLoop_thisRepN,thisData.erpRewp,thisData.rERPRewp,beh.blockType,beh.outcome,...
        beh.trialResp_corr,beh.trialResp_rt]; %1728 9
    %rt adjustment
    for t = 0:(length(beh_temp)-1)
         if t == 0
            beh_temp(t+1,10) = nan;
         elseif beh_temp(t,2) == beh_temp(t+1,2) %if same block, calc
            beh_temp(t+1,10) = beh_temp(t+1,9)-beh_temp(t,9);
         else
             beh_temp(t+1,10) = nan;
         end
    end
    
    beh_mat_temp = zeros(1,size(rerpdata,3),10,length(times)); %1 1728 9 200
    rerp_mat_temp = zeros(1,size(rerpdata,3),1,length(times)); %1 1728 1 200
    for i = 1:length(times)
        beh_mat_temp(1,:,:,i) = beh_temp;
    end
    
    for i = 1:length(times)
        rerp_mat_temp(1,:,1,i) = rerpdata(iElectrode,i,:);
    end
    
        data_driven_mat(p,:,1:10,:) = beh_mat_temp;
        data_driven_mat(p,:,11,:) = rerp_mat_temp;
end
        

%% for each time point, generate a df for r
data_driven_out = [];

for t = 1:length(times)
    for p = 1:length(ps)
        data_driven_out = [data_driven_out;[squeeze(data_driven_mat(p,:,:,t)),repmat(times(t),1728,1)]];
    end
    
end

data_driven_out =  array2table(data_driven_out);
data_driven_out.Properties.VariableNames = {'participant_ID','blockLoop_thisRepN','trialLoop_thisRepN','erpRewp','rERPRewp','blockType','outcome',...
    'trialResp_corr','trialResp_rt','rt_adjustment','rEEG','rEEG_time'};

writetable(data_driven_out,[dataFolder '/R/tbt_meta_data_driven.csv']);

%% get rt adjustment and EEG in each "on time" trial
data_fast_ontime = zeros(length(ps),350,200);
data_medium_ontime = zeros(length(ps),350,200);
data_slow_ontime = zeros(length(ps),350,200);
beh_fast_ontime= zeros(length(ps),350);
beh_medium_ontime = zeros(length(ps),350);
beh_slow_ontime = zeros(length(ps),350);
for p = 1:length(ps)
    this_data = data_driven_mat(p,:,11,:); %3 dim doesnt matter
    this_beh = data_driven_mat(p,:,10,1);%4 dim doesnt matter
    
    ind_fast_ontime = ((data_driven_mat(p,:,6,1) <= 4) & (data_driven_mat(p,:,7,1) == 1)); %the 4th dim doesnt matter here
    data_fast_ontime(p,1:sum(ind_fast_ontime),:) = squeeze(this_data(1,ind_fast_ontime,1,:));
    beh_fast_ontime(p,1:sum(ind_fast_ontime)) = this_beh(1,ind_fast_ontime);
    
    ind_medium_ontime = ((data_driven_mat(p,:,6,1) > 4) & (data_driven_mat(p,:,6,1) <= 8) & (data_driven_mat(p,:,7,1) == 1)); %the 4th dim doesnt matter here
    data_medium_ontime(p,1:sum(ind_medium_ontime),:) = squeeze(this_data(1,ind_medium_ontime,1,:));
    beh_medium_ontime(p,1:sum(ind_medium_ontime)) = this_beh(1,ind_medium_ontime);
    
    ind_slow_ontime = ((data_driven_mat(p,:,6,1) >8) & (data_driven_mat(p,:,7,1) == 1)); %the 4th dim doesnt matter here
    data_slow_ontime(p,1:sum(ind_slow_ontime),:) = squeeze(this_data(1,ind_slow_ontime,1,:));
    beh_slow_ontime(p,1:sum(ind_slow_ontime)) = this_beh(1,ind_slow_ontime);

end
%%
subplot(3,1,1);
for p = 1:length(ps)
    plot(times,squeeze(mean(data_fast_ontime(p,:,:),2,'omitnan'))); hold on;
end
plot(times,squeeze(mean(data_fast_ontime,[1,2],'omitnan')),'lineWidth',2,'Color','k'); hold on;
xline(times(times==0));

subplot(3,1,2);
for p = 1:length(ps)
    plot(times,squeeze(mean(data_medium_ontime(p,:,:),2,'omitnan'))); hold on;
end
plot(times,squeeze(mean(data_medium_ontime,[1,2],'omitnan')),'lineWidth',2,'Color','k'); hold on;
xline(times(times==0));

subplot(3,1,3);
for p = 1:length(ps)
    plot(times,squeeze(mean(data_slow_ontime(p,:,:),2,'omitnan'))); hold on;
end
plot(times,squeeze(mean(data_slow_ontime,[1,2],'omitnan')),'lineWidth',2,'Color','k'); hold on;
xline(times(times==0));
