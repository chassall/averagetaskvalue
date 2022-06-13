% Load and preprocess raw data for the Average Task Value study
%
% Other required m-files: 
% EEGLAB with iclabel and bva-io extensions,
% find_artifacts.m
% Other required files: 
% erpclass2016.locs
% ccnlabactichamp.locs

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Initial call to EEGLAB to set path
eeglab(); clear all;

% Seed random number generator for consistency
rng(2016);

% Participants 11-2: pilot data, slightly different task
% Participants 2-26: Site 1 (UVic)
% Participants 27-38: Site 2 (Oxford)
participants = 1:38;
nParticipants = length(participants);

% Bad channels (as determined by the next script, analysis_02_checkpreprocessed)
badChannels = {{},{'FP1','FP2','F3','Fz','F4','F8','Pz'},{},{},{},{},{},{'FP2','T8'},{},{},{},{},{'FP1','FP2'},{},{},{},{},{'FP1','FP2','F7','F3','Fz','F4','F8','FC5','FC1','FCz','FC2','FC6'},{},{},{},{},{},{},{},{},{},{},{},{},{},{'TP10'},{},{},{},{},{},{}};

% S2, S20, S25 had some dead channels
badChannels{2} = [badChannels{2} {'F7','F3','FC2'}];
badChannels{20} = [badChannels{20} {'CP1'}];
badChannels{25} = [badChannels{25} {'T7'}];

% Participant loop
for iParticipant = 1:nParticipants
    
    % Specify location of raw data
    if participants(iParticipant) < 27
        recordingSite = 1;
        rawFolder = ['/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-1/sub-' num2str(participants(iParticipant),'%0.2i') '/eeg'];
    else
        recordingSite = 2;
        rawFolder = ['/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-2/sub-' num2str(participants(iParticipant),'%0.2i') '/eeg'];
    end

    %% Load the raw EEG.
    if recordingSite == 1
        pString = num2str(participants(iParticipant),'%0.04i');
        rawFile = ['erpclass2016_' pString '.vhdr']; % Not BIDS
    else
        rawFile = ['sub-' num2str(participants(iParticipant)) '_task-casinos_eeg.vhdr']; % BIDS
    end
    EEG = pop_loadbv(rawFolder, rawFile);
    
    %% Add reference
    if recordingSite == 2
        EEG.data = [EEG.data; zeros(1,size(EEG.data,2))];
        EEG.nbchan = EEG.nbchan + 1;
        EEG.chanlocs(length(EEG.chanlocs)+1) =  EEG.chanlocs(end);
        EEG.chanlocs(end).labels = 'Fz';
    end

    %% Remove some channels from P1, P2, P11 to align with the rest
    if ismember(participants(iParticipant),[1,2,11])
        EEG = pop_select(EEG, 'nochannel', [eeg_chaninds(EEG,{'FPz'}) eeg_chaninds(EEG,{'FT9'}) eeg_chaninds(EEG,{'TP7'}) eeg_chaninds(EEG,{'TP8'}) eeg_chaninds(EEG,{'FT10'}) eeg_chaninds(EEG,{'CPz'})] );
    end
    
    %% Load/set locations file
    if recordingSite == 1
        EEG.chanlocs = readlocs('site1channellocations.locs');
    else
        EEG.chanlocs = readlocs('site2channellocations.locs');
    end
    
    %% Downsample to 250 Hz.
    EEG = pop_resample(EEG, 250);
    
    %% Apply a bandpass filter (0.1-30 Hz).
    EEG = pop_eegfiltnew(EEG, 0.1, 30);
    
    %% Notch filter
    if recordingSite == 1
        EEG = pop_eegfiltnew(EEG, 58, 62,[],1); % 60 Hz
    else
        EEG = pop_eegfiltnew(EEG, 48, 52,[],1); % 50 Hz
    end
    
    %% Uncomment for manual inspection - click on "scroll data" then "OK".
    % EEG = pop_select(EEG);
    
    %%  Remove up to 4 bad channels (otherwise, exclude participant)
    theseBadChannels = badChannels{participants(iParticipant)};
    disp(['removing ' theseBadChannels{:}]);
    if length(theseBadChannels) <= 4
        EEG = pop_select(EEG, 'nochannel', theseBadChannels);
    end

    %% Re-reference to the average of the left and right mastoids.
    if recordingSite == 1
        EEG = pop_reref(EEG, {'M1','M2'});
    else
        if any(strcmp({EEG.chanlocs.labels},'TP9')) && any(strcmp({EEG.chanlocs.labels},'TP10'))
            EEG = pop_reref(EEG, {'TP9','TP10'});
        elseif any(strcmp({EEG.chanlocs.labels},'TP9'))
            EEG = pop_reref(EEG, {'TP9'});
        elseif any(strcmp({EEG.chanlocs.labels},'TP10'))
            EEG = pop_reref(EEG, {'TP10'});
        else
            error('no mastoid channels');
        end
    end
    
    %% Remove ocular channels if present (these are unused)
    if recordingSite == 1
        EEG = pop_select(EEG,'nochannel' ,{'LHEOG' 'RHEOG' 'VEOG'});
    end
      
    %% Isolate some data on which to run the ICA
    icaEEG = pop_epoch(EEG,{'S  1','S 11','S 21','S 31'},[0 3]);
    
    %% Remove bad trials from icaEEG (> 1000 uV change).
    badTrialIndex = find_artifacts(icaEEG, 500, 500, 40, 0.1);
    icaEEG = pop_select(icaEEG,'notrial',badTrialIndex);
    
    %% Run ICA and get the results.
    
    % Possible algorithms: 'binica','fastica','runica'.
    icaEEG = pop_runica(icaEEG,'runica');
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
    
    %% Save preprocessed EEG
    preprocessedFolder = ['/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data/derivatives/sub-' num2str(participants(iParticipant))]; 
    preprocessedFile = ['sub-' num2str(participants(iParticipant),'%0.2i') '_task-casinos_prep.mat'];

    if ~exist(preprocessedFolder,'dir')
        mkdir(preprocessedFolder);
    end
    save(fullfile(preprocessedFolder,preprocessedFile),'EEG');
end
