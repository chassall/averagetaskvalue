% Load and preprocess raw data for the Average Task Value study
%
% Other m-files required: EEGLAB with iclabel extension

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

participants = 1:26; % NB: 1,2,11 had a different montage
nParticipants = length(participants);
rawFolder = 'E:\tan_lab_raw\2016_EEG_PSYC574A_Hassall\Raw'; % Location of EEG data

% Bad channels (as determined by analysis_02_checkpreprocessed)
badChannels = {{},{'FP1','FP2','F3','Fz','F4','F8','Pz'},{},{},{},{},{},{'FP2','T8'},{},{},{},{},{'FP1','FP2'},{},{},{},{},{'FP1','FP2','F7','F3','Fz','F4','F8','FC5','FC1','FCz','FC2','FC6'},{},{},{},{},{},{},{},{}};

% Flat channels
badChannels{2} = [badChannels{2} {'F7','F3','FC2'}];
badChannels{20} = [badChannels{20} {'CP1'}];
badChannels{25} = [badChannels{25} {'T7'}];

% Participant loop
for iParticipant = 1:nParticipants
    
    %% Load the raw EEG.
    pString = num2str(participants(iParticipant),'%0.04i');
    rawFile = ['erpclass2016_' pString '.vhdr'];
    EEG = pop_loadbv(rawFolder, rawFile);
    
    %% Remove some channels from P1, P2, P11 to align with the rest
    if ismember(participants(iParticipant),[1,2,11])
        EEG = pop_select(EEG, 'nochannel', [eeg_chaninds(EEG,{'FPz'}) eeg_chaninds(EEG,{'FT9'}) eeg_chaninds(EEG,{'TP7'}) eeg_chaninds(EEG,{'TP8'}) eeg_chaninds(EEG,{'FT10'}) eeg_chaninds(EEG,{'CPz'})] );
    end
    
    %% Load/set locations file
    EEG.chanlocs = readlocs('erpclass2016.locs');
    
    %% Downsample to 250 Hz.
    EEG = pop_resample(EEG, 250);
    
    %% Apply a bandpass filter (0.1-20 Hz).
    EEG = pop_eegfiltnew(EEG, 0.1, 30);
    
    %% Notch filter
    EEG = pop_eegfiltnew(EEG, 58, 62,[],1);
    
    %% Uncomment for manual inspection - click on "scroll data" then "OK".
    % EEG = pop_select(EEG);
    
    %% Re-reference to the average of the left and right mastoids.
    EEG = pop_reref(EEG, {'M1','M2'});
    
    %% Remove ocular channels from main dataset
    EEG = pop_select(EEG,'nochannel' ,{'LHEOG' 'RHEOG' 'VEOG'});

    %% Make a copy of the full locations file
    fullLocs = EEG.chanlocs;
    
    %%  Remove up to 4 bad channels (otherwise, exclude participant)
    theseBadChannels = badChannels{iParticipant};
    if length(theseBadChannels) <= 4
        EEG = pop_select(EEG, 'nochannel', theseBadChannels);
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
    
    %% Interpolate missing channels (uncomment if needed).
    % EEG = pop_interp(EEG,fullLocs);
    
    %% Save preprocessed EEG
    preprocessedFolder = './preprocessed';
    preprocessedFile = ['erpclass2016_' pString '_prep'];
    save(fullfile(preprocessedFolder,preprocessedFile),'EEG','fullLocs');
end
