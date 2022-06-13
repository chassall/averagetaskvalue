% Check preprocessed data for the Average Task Value study
%
% Other m-files required: 
% EEGLAB with iclabel extension
% make_erp.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Initial call to EEGLAB to set path
eeglab(); clear all;

participants = 1:38;
nParticipants = length(participants);

% Set output folder
outputFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\analysis\output';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

allArtifacts = [];
removed = [];

% Participant loop
for iParticipant = 1:nParticipants

    % Set data folder
    if participants(iParticipant) <= 27
        dataFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\data_site-1';
    else
        dataFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\data_site-2';
    end

    pString = ['sub-' num2str(participants(iParticipant),'%0.02i')];
    preprocessedFolder = [dataFolder '/derivatives/' pString];
    preprocessedFile = [pString '_task-casinos_prep.mat'];
    eyeCorrFile = [pString '_task-casinos_prepoc.mat'];
    load(fullfile(preprocessedFolder,preprocessedFile),'EEG','fullLocs');
    
    % Remove ocular components using the results of ICLabel
    eyeLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Eye'));
    brainLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Brain'));
    muscleLabel = find(strcmp(EEG.etc.ic_classification.ICLabel.classes,'Muscle'));
    
    eyeI = EEG.etc.ic_classification.ICLabel.classifications(:,eyeLabel) > EEG.etc.ic_classification.ICLabel.classifications(:,brainLabel);
    whichOnes = find(eyeI);  
    removed(iParticipant) = sum(eyeI);
    disp(removed);
    
    % Remove ocular components
    EEG = pop_subcomp(EEG,whichOnes,0);
    
    % Interpolate as needed
    commonLocs = readlocs('common.locs');
    EEG = pop_interp(EEG,commonLocs);

    % Remove some channels if this is the second recording site
    if participants(iParticipant) >= 27
        EEG = pop_select(EEG,'nochannel',{'O1','O2','CPz',});
    end

    % Make ERPs
    [erp, allArtifacts(iParticipant,:),times] = make_erp(EEG,{'S  2','S 12','S 22','S 32','S  6','S 16','S 26','S 36','S  7','S 17','S 27','S 37'},[-0.2 .6],[-200 0]);
    save(fullfile(preprocessedFolder,eyeCorrFile),'EEG');

    % Test Plots
    [winERP,~,times] = make_erp(EEG,{'S  6','S 16','S 26','S 36'},[-0.2 .6],[-200 0]);
    [loseERP,~,times] = make_erp(EEG,{'S  7','S 17','S 27','S 37'},[-0.2 .6],[-200 0]);
    diff = winERP - loseERP;
    tWindow = 1000*[.256, .288];
    fbChannel = 'FCz'
    iWindow = dsearchn(times',tWindow');
    fbChannelI = eeg_chaninds(EEG,{fbChannel});
    figure();
    subplot(1,2,1);
    plot(times,winERP(fbChannelI,:)); hold on; plot(times,loseERP(fbChannelI,:)); legend('win','lose'); title(pString);
    diffTopo = mean(diff(:,iWindow(1):iWindow(2)),2);
    subplot(1,2,2); topoplot(diffTopo, EEG.chanlocs); hold off;
    pause();
end

% Save a copy of which ICs were removed
save(fullfile(outputFolder,'icsremoved.mat'),'removed');

%% Find bad channels

% Label a channel as bad if it would cause more than 20% of epochs to be
% rejected.
allLabels = {EEG.chanlocs.labels};
isBad = allArtifacts > 0.20;
badChannels = {};
for iParticipant = 1:nParticipants
    theseBadChannelIs = find(isBad(iParticipant,:));
    badChannels{iParticipant} = allLabels(theseBadChannelIs);
end

isBadString = ['{'];
for b = 1:length(badChannels)  
    
    isBadString = [isBadString '{'];
    
    for c = 1:length(badChannels{b})
        isBadString = [isBadString '''' badChannels{b}{c} ''','];
    end 
    if length(badChannels{b})
        isBadString(end) = [];
    end
    
    isBadString = [isBadString '},'];
end

isBadString(end) = '}';

% Save a copy of the bad channel list
save(fullfile(outputFolder,'badchannels.mat'),'isBadString');