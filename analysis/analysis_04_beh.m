%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

% Analyze behavioural data for the Drum Trainer project
%
% Other m-files required: 
% makefigure.m
% subtightplot.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

close all; clear all;

% Participants to include
ps = {'01','02','03','04','05','06','07','08','09','10','11','13','14','15','16','17','18','19','20','21'};

% Set data folder - change as needed
if ispc
    dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/data';
    outputFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2021_EEG_DrumTrainer_Hassall/analysis/output';
else
    dataFolder = '/Users/rh/Documents/ds004152/';
    outputFolder = '/Users/rh/Documents/ds004152/R/';
end

if ~exist(outputFolder)
    mkdir(outputFolder);
end

conditions = 1:12;
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

% Behavioural variables
rtByCond = [];
rtByOutcome = [];
rtByTempoAndOutcome = [];
rtAdjustByOutcome = [];
rtAdjustByTempoAndOutcome = [];
allMargins = [];
allSDs = [];
fbCountsByPattern = [];
fbCountsByTempo = [];
fbCounts = [];
allBonus = [];
incResponse = [];
corrByCond = [];
corrAll = [];
corrByPattern = [];

numPracticeTrials = 14; % Based on examination of the mean window for each condition, this is ~10% of trials in each condition

% Loop through participants
%figure();
for p = 1:length(ps)
    
    % Load data
    subName = ['sub-' ps{p}];
    rawFile = [subName '_task-drumtrainer_beh.tsv'];
    thisData = readtable(fullfile(dataFolder,subName,'beh',rawFile),'FileType','text');

    % Bonus payment
    numCor = sum((thisData.trialResp_corr==1) & (thisData.outcome == 1));
    numInc = sum((thisData.trialResp_corr==1) & (thisData.outcome ~= 1));
    bonus = 0.001 * 2 * numCor + 0.001 * 1 * numInc;
    
    % Trial types
    isTrial = ismember(thisData.blockType,1:12);
    isFastTrial = ismember(thisData.blockType,1:4);
    isMediumTrial = ismember(thisData.blockType,5:8);
    isSlowTrial = ismember(thisData.blockType,9:12);
    isPattern1 = ismember(thisData.blockType,[1,2,5,6,9,10]);
    isPattern2 = ismember(thisData.blockType,[3,4,7,8,11,12]);
    isEarly = thisData.outcome == 0;
    isOnTime = thisData.outcome == 1;
    isLate = thisData.outcome == 2;
    corrButton = thisData.trialResp_corr == 1;
    
    % Get correct flag for each tempo, pattern
    p1Data =thisData(isPattern1,:);
    p2Data = thisData(isPattern2,:);
    fData = thisData(isFastTrial,:);
    mData = thisData(isMediumTrial,:);
    sData = thisData(isSlowTrial,:);
    fCorrect  = fData.trialResp_corr;
    mCorrect  = mData.trialResp_corr;
    sCorrect  = sData.trialResp_corr;
    corrByCond(p,:,:) = [fCorrect(1:576) mCorrect(1:576) sCorrect(1:576)];
    corrAll(p,:) = thisData.trialResp_corr(1:1728);
    corrByPattern(p,:,:) = [p1Data.trialResp_corr(1:864) p2Data.trialResp_corr(1:864)];

    numX = sum(isTrial & ~corrButton);
    numCor = sum(isOnTime & corrButton);
    numInc = sum((isEarly | isLate) & corrButton);
    numPoints = 2*numCor - 1*numInc;
    if p == 1
        bonus = 0.5 * 0.001 * numPoints;
    else
        bonus = 0.001 * numPoints;
    end
    allBonus = [allBonus; bonus];
    
    % Basic RT thresholds (0.1 - 2 s)
    isOK = thisData.trialResp_rt < 2 & thisData.trialResp_rt > 0.1;
    
    % Conditions
    isPractice = zeros(size(thisData,1),1);
    isThisCond = [];
    for c = 1:length(conditions)
        isThisCond = ismember(thisData.blockType,conditions(c));
        thisCondI = find(isThisCond);
        thisCondPractice = thisCondI(1:numPracticeTrials);
        isPractice(thisCondPractice) = 1;
        theseMargins = thisData.margin(isThisCond==1);
        allMargins(p,c,:) = theseMargins(1:144);
        allSDs(p,c) = std(thisData.trialResp_rt(isOK & isThisCond==1));
    end
    theseMeanMargins = mean(allMargins(p,:,:),3);
    theseStdMargins = std(allMargins(p,:,:),[],3);
    marginThresh = theseMeanMargins + 0.5*theseStdMargins;
    
%     % Plot histograms of RT windows with threshold
%     for c = 1:12
%         subplot(21,12,(p-1)*12+c); histogram(allMargins(p,c,:)); xline(marginThresh(c));
%     end
    
    marginWithinThresh = nan(size(thisData,1),1);
    for c = 1:length(conditions)
        isThisCond = ismember(thisData.blockType,conditions(c));
        marginWithinThresh(isThisCond) = thisData.margin(isThisCond==1) < marginThresh(c);
    end
    
    % To include:
    % - hit the correct button AND
    % - is not a practice trial AND
    % - RT window was less than threshold defined above
    isOK = isOK & corrButton & ~isPractice & (marginWithinThresh==1);
        
    % Trial counts 
    fbCounts(p,1) = sum(isOK & isEarly);
    fbCounts(p,2) = sum(isOK & isOnTime);
    fbCounts(p,3) = sum(isOK & isLate);
    
    % Count the number of trial of eacc type (participant X pattern X outcome)
    fbCountsByPattern(p,1,1) = sum(isPattern1 & isEarly);
    fbCountsByPattern(p,1,2)  = sum(isPattern1 & isOnTime);
    fbCountsByPattern(p,1,3)  = sum(isPattern1 & isLate); 
    fbCountsByPattern(p,2,1) = sum(isPattern2 & isEarly);
    fbCountsByPattern(p,2,2)  = sum(isPattern2 & isOnTime);
    fbCountsByPattern(p,2,3)  = sum(isPattern2 & isLate);
    
    % Trial counts (participant X tempo X outcome)
    fbCountsByTempo(p,1,1) = sum(isFastTrial & isEarly);
    fbCountsByTempo(p,2,1) = sum(isMediumTrial & isEarly);
    fbCountsByTempo(p,3,1) = sum(isSlowTrial & isEarly);
    fbCountsByTempo(p,1,2) = sum(isFastTrial & isOnTime);
    fbCountsByTempo(p,2,2) = sum(isMediumTrial & isOnTime);
    fbCountsByTempo(p,3,2) = sum(isSlowTrial & isOnTime);
    fbCountsByTempo(p,1,3) = sum(isFastTrial & isLate);
    fbCountsByTempo(p,2,3) = sum(isMediumTrial & isLate);
    fbCountsByTempo(p,3,3) = sum(isSlowTrial & isLate);
    
    % RT by outcome (participants X outcome)
    rtByOutcome(p,1) = nanmean(thisData.trialResp_rt(isOK & isEarly));
    rtByOutcome(p,2) = nanmean(thisData.trialResp_rt(isOK & isOnTime));
    rtByOutcome(p,3) = nanmean(thisData.trialResp_rt(isOK & isLate));
    
    % RT by tempo and outcome (participant X tempo X outcome)
    % Not to be trusted since the trial count is low in some cases
    rtByTempoAndOutcome(p,1,1) = mean(thisData.trialResp_rt(isOK & isEarly & isFastTrial));
    rtByTempoAndOutcome(p,2,1) = mean(thisData.trialResp_rt(isOK & isEarly & isMediumTrial));
    rtByTempoAndOutcome(p,3,1) = mean(thisData.trialResp_rt(isOK & isEarly & isSlowTrial));
    rtByTempoAndOutcome(p,1,2) = mean(thisData.trialResp_rt(isOK & isOnTime & isFastTrial));
    rtByTempoAndOutcome(p,2,2) = mean(thisData.trialResp_rt(isOK & isOnTime & isMediumTrial));
    rtByTempoAndOutcome(p,3,2) = mean(thisData.trialResp_rt(isOK & isOnTime & isSlowTrial));
    rtByTempoAndOutcome(p,1,3) = mean(thisData.trialResp_rt(isOK & isLate & isFastTrial));
    rtByTempoAndOutcome(p,2,3) = mean(thisData.trialResp_rt(isOK & isLate & isMediumTrial));
    rtByTempoAndOutcome(p,3,3) = mean(thisData.trialResp_rt(isOK & isLate & isSlowTrial)); 
    
    for c = 1:12
        rtByCond(p,c) = nanmean(thisData.trialResp_rt(isOK & (thisData.blockType == conditions(c))));
    end
    
end

%% Compute bonus payment mean + ci
[~,~,CI,~] = ttest(allBonus);
disp(['mean bonus ' num2str(mean(allBonus))]);
disp(['ci ' num2str(CI')]);


%% Main behavioural figure

% Overwrite subplot with subtightplot for more control
subplot = @(m,n,iParticipant) subtightplot (m, n, iParticipant, [0.1 0.1], [0.2 0.15], [0.08 0.05]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

makefigure(18,5);
fontSize = 7;
lineWidth = 1;
condColours = brewermap(5,'Set1');
axs = {};

% Error count
axs{1} = subplot(1,3,1);
incAll = 1-corrAll;
incByPattern = 1-corrByPattern;

mmWin = 60;
mmIncorrect = movmean(incByPattern,mmWin,2);

meanErrorCount = squeeze(mean(mmIncorrect,1));
stdErrorCount = squeeze(std(meanErrorCount,[],1));
tval = abs(tinv(0.025,length(ps)-1));
ciErrorCount = tval * stdErrorCount / sqrt(length(ps));
% boundedline(1:1728, meanErrorCount, ciErrorCount);
% bar(1:length=,totalErrorAcrossPs');
[hls,hps] = boundedline(1:864,meanErrorCount,ciErrorCount,'cmap',condColours(4:5,:));
hls(1).LineWidth = lineWidth;
hls(2).LineWidth = lineWidth; 
% colororder(condColours(4:5,:));

% plot(squeeze(mean(mmIncorrect,1)),'LineWidth',lineWidth);
xlabel('Trial number');
ylabel('P(wrong key)');
% colororder(condColours(4:5,:));
l = legend(hls,{'AABA','AAABAA'});
l.Box = 'off';

% Response window
meanMargins = squeeze(mean(allMargins,1));

meanF = mean(meanMargins(1:4,:),1);
stdF = std(meanMargins(1:4,:),[],1);
ciF = tval * stdF / sqrt(length(ps));
meanM = mean(meanMargins(5:8,:),1);
stdM = std(meanMargins(5:8,:),[],1);
ciM = tval * stdM / sqrt(length(ps));
meanS = mean(meanMargins(9:12,:),1);
stdS = std(meanMargins(9:12,:),[],1);
ciS = tval * stdS / sqrt(length(ps));
axs{2} = subplot(1,3,2);
[hl1,~] = boundedline(1:length(meanF),meanF,ciF,'cmap',condColours(1,:),'alpha');
hl1.LineWidth = lineWidth;
hold on;
[hl2,~] = boundedline(1:length(meanM),meanM,ciM,'cmap',condColours(2,:),'alpha');
hl2.LineWidth = lineWidth;
[hl3,~] = boundedline(1:length(meanS),meanS,ciS,'cmap',condColours(3,:),'alpha');
hl3.LineWidth = lineWidth;
ylim([0.0,0.15]);
xlabel('Trial number');
ylabel('Response window (s)');
l = legend([hl1,hl2,hl3],{'Fast','Medium','Slow'},'Location','NorthWest');
l.Box = 'off';
arrowText = text(77,0.093,'\leftarrow New block','FontSize',fontSize,'FontAngle','italic');

% Mean RT
meanFRT = mean(rtByCond(:,1:4),2);
meanMRT = mean(rtByCond(:,5:8),2);
meanSRT = mean(rtByCond(:,9:12),2);
axs{3} = subplot(1,3,3);
nbp = notBoxPlot([meanFRT meanMRT meanSRT],'interval','tinterval');
formatNBP(nbp,condColours);
ylim([0,1]);
xlabel('Tempo');
xticklabels({'Fast','Medium','Slow'});
ylabel('Response time (s)');
xlim([0 4.5]);

% Target times
hold on;
targets = [0.4 0.6 1];
% targetNames = {'150 BPM','100 BPM','60 BPM'};
targetNames = {'Fast','Medium','Slow'};
for t = 1:length(targets)
    yline(targets(t),'LineStyle','-','LineWidth',lineWidth/2,'Color',condColours(t,:)); 
    text(3.5,targets(t)+0.03,targetNames{t},'FontSize',fontSize,'Color',condColours(t,:),'FontAngle','italic');
end


for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end
print(fullfile(outputFolder,'beh.tiff'),'-dtiff','-r600');

return;
%% FB count
makefigure();
nbp = notBoxPlot(fbCounts);
formatNBP(nbp);
ylabel('Feedback count');
set(gca,'XTickLabel',{'Early','OnTime','Late'});
%print(fullfile(outputFolder,'dt_powerpoint_fbcount.tiff'),'-dtiff','-r300');

%% FB count by pattern
temp = [];
for i = 1:2
    temp = [temp squeeze(fbCountsByPattern(:,i,:))];
end
makefigure();
nbp = notBoxPlot(temp);
formatNBP(nbp);
ylabel('Feedback count');
set(gca,'XTickLabel',{'P1Early','P1OnTime','P1Late','P2Early','P2OnTime','P2Late'});
%print(fullfile(outputFolder,'dt_powerpoint_fbcountbbytempo.tiff'),'-dtiff','-r300');

%% FB count by tempo
temp = [];
for i = 1:3
    temp = [temp squeeze(fbCountsByTempo(:,i,:))];
end
makefigure();
nbp = notBoxPlot(temp);
formatNBP(nbp);
ylabel('Feedback count');
set(gca,'XTickLabel',{'FastEarly','FastOnTime','FastLate','MediumEarly','MediumOnTime','MediumLate','SlowEarly','SlowOnTime','SlowLate'});
%print(fullfile(outputFolder,'dt_powerpoint_fbcountbbypattern.tiff'),'-dtiff','-r300');

%% RT by condtion
makefigure();
nbp = notBoxPlot(rtByCond);
formatNBP(nbp);
ylim([0 1.1]);
yline(0.4); yline(0.6); yline(1);
ylabel('Response time (s)');
set(gca,'XTickLabel',{'LP1F','RP1F','LP2F','RP2F','LP1M','RP1M','LP2M','RP2M','LP1S','RP1S','LP2S','RP2S'});
%print(fullfile(outputFolder,'dt_powerpoint_rt.tiff'),'-dtiff','-r300');

%% RT Window by condition
meanMargins = squeeze(mean(allMargins,1));
meanFP1 = mean(meanMargins(1:2,:),1);
meanFP2 = mean(meanMargins(3:4,:),1);
meanMP1 = mean(meanMargins(5:6,:),1);
meanMP2 = mean(meanMargins(7:8,:),1);
meanSP1 = mean(meanMargins(9:10,:),1);
meanSP2 = mean(meanMargins(11:12,:),1);
makefigure();
plot(meanFP1); hold on; 
plot(meanFP2); 
plot(meanMP1); 
plot(meanMP2); 
plot(meanSP1); 
plot(meanSP2);
legend('fp1','fp2','mp1','mp2','sp1','sp2');

%% Mean RT by prior outcome (early, on time, late);
makefigure(); 
subplot(1,2,1); 
nbp1 = notBoxPlot(rtByOutcome); 
formatNBP(nbp1); ax = gca;
ax.XTickLabel = {'Early','On Time','Late'};
ax.YLabel.String = 'Response time (s)';
earlyAdjust = abs(rtByOutcome(:,1) - rtByOutcome(:,2));
lateAdjust = abs(rtByOutcome(:,2) - rtByOutcome(:,3));
subplot(1,2,2); 
nbp2 = notBoxPlot([earlyAdjust lateAdjust]); 
formatNBP(nbp2); ax = gca;
ax.XTickLabel = {'Early','Late'};
ax.YLabel.String = 'Adjustment rel. to "on time" (s)';
[h,p] = ttest(earlyAdjust,lateAdjust)

%% Mean RT by temp (fast, med, slow) and prior outcome (too fast, on time, too slow)
tempos = {'Fast','Medium','Slow'};
feedbacks = {'early','on time','slow'};
makefigure();
for i = 1:3
    theseRTs = squeeze(rtByTempoAndOutcome(:,i,:));
    subplot(2,3,i); 
    nbp1 = notBoxPlot(theseRTs); 
    formatNBP(nbp1);
    ylim([0.3,1.1]);
    title(tempos{i});
    ax =gca;
    ax.XTickLabel = feedbacks;
    isEarlyDiff = abs(theseRTs(:,1)-theseRTs(:,2));
    isLateDiff = abs(theseRTs(:,2)-theseRTs(:,3));
    subplot(2,3,3+i); 
    nbp2 = notBoxPlot([isEarlyDiff isLateDiff]); 
    formatNBP(nbp2);
end