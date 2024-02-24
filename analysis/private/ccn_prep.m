function ccn_prep(dataFolder,subName,taskName,icaTriggers,icaWindow,filters,reference,badChannels,appendString)
%CCN_PREP Preprocess CCN lab data
%   dataFolder: high-level BIDS folder
%   subName: subject name, e.g. 'sub-01'
%   taskName: task name, e.g. 'cooltask'
%   icaTriggers: the triggers around which to train the ICA
%   icaWindow: window to train ICA
%   filters: low/high, e.g [0.1 30]
%   reference: e.g. {'TP9','TP10'};
%   badChannels: array of bad channels to remove, e.g.  {'T7','T8'};
%   appendString: = string to append to output file, e.g. 'test'

if nargin < 7
    filters = [0.1 30];
    reference = {'TP9','TP10'};
    badChannels = {};
    appendString = '';
end

%% Load data
rawFolder = fullfile(dataFolder,subName,'eeg');
headerFile = dir(fullfile(rawFolder, '*.vhdr'));
if length(headerFile) ~= 1
    error(['No header file (or more than one) in ' rawFolder]);
    return;
end
EEG = pop_loadbv(rawFolder,headerFile.name);

%% Add reference
EEG.data = [EEG.data; zeros(1,size(EEG.data,2))];
EEG.nbchan = EEG.nbchan + 1;
EEG.chanlocs(length(EEG.chanlocs)+1) =  EEG.chanlocs(end);
EEG.chanlocs(end).labels = 'Fz';

%% Channel locations
eegLocs = readlocs('ccnlabactichamp.locs');
EEG.chanlocs = eegLocs;

%% Downsample to 250 Hz.
EEG = pop_resample(EEG, 250);

%% Apply a bandpass filter
EEG = pop_eegfiltnew(EEG, filters(1), filters(2));

%% 50 Hz Notch filter
EEG = pop_eegfiltnew(EEG, 48, 52,[],1);

%% Re-reference to the average of the left and right mastoids and remove them from the analysis.
EEG = pop_reref(EEG, reference,'keepref', 'off');

%% Remove bad channels
fullLocs = EEG.chanlocs;
EEG = pop_select(EEG,'nochannel',badChannels);

%% Isolate some data on which to run the ICA
% Training data: a 3-second epoch around the initial fixation cross
icaEEG = pop_epoch(EEG,icaTriggers,icaWindow);

%% Remove bad trials from icaEEG (max - min sample > 500 uV)
isArtifactsCT = abs((max(icaEEG.data,[],2) -  min(icaEEG.data,[],2))) > 500; % size: channels X 1 X epochs
isArtifact = logical(squeeze(any(isArtifactsCT,1))); % size: epochs X 1
icaEEG = pop_select(icaEEG,'notrial',isArtifact);

%% Run ICA and get the results.
icaEEG = pop_runica(icaEEG,'runica'); % Possible algorithms: 'binica','fastica','runica'.
icaact = icaEEG.icaact;
icachansind = icaEEG.icachansind;
icasphere = icaEEG.icasphere;
icasplinefile = icaEEG.icasplinefile;
icaweights = icaEEG.icaweights;
icawinv = icaEEG.icawinv;

%% Transfer the ICA results from icaEEG to EEG
EEG.icaact = icaact;
EEG.icachansind = icachansind;
EEG.icasphere = icasphere;
EEG.icasplinefile = icasplinefile;
EEG.icaweights = icaweights;
EEG.icawinv = icawinv;

%% Perform IC rejection using the ICLabel EEGLAB extension.
EEG = iclabel(EEG, 'default');

%% Adding ocular correction here, CDH, 13 Jan 2022
eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
eyeI = EEG.etc.ic_classification.ICLabel.classifications(:,eyeLabel)  > 0.8;
whichOnes = find(eyeI);

% Remove ocular components
EEG = pop_subcomp(EEG,whichOnes,0);
EEG.numOcular = sum(eyeI);
disp(['removing ' num2str(EEG.numOcular) ' ocular components']);

% Interpolate
EEG = pop_interp(EEG,fullLocs);

%% Save preprocessed EEG
prepFile = [subName '_task-' taskName '_eegprep' appendString '.mat'];
prepFolder =  [dataFolder '/derivatives/eegprep/' subName];
if ~isfolder(prepFolder)
    mkdir(prepFolder);
end
save(fullfile(prepFolder,prepFile),'EEG');
end

