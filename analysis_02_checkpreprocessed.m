% Check preprocessed data for the Average Task Value study
%
% Other m-files required: EEGLAB with iclabel extension

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

participants = 1:26;
nParticipants = length(participants);
preprocessedFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_PSYC574A\Analysis\preprocessed';
allArtifacts = [];
removed = [];

% Participant loop
for iParticipant = 1:nParticipants
    pString = num2str(participants(iParticipant),'%0.04i');
    preprocessedFile = ['erpclass2016_' pString '_prep.mat'];
    eyeCorrFile = ['erpclass2016_' pString '_eye.mat'];
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
    EEG = pop_interp(EEG,fullLocs);

    % Make ERPs
    [erp, allArtifacts(iParticipant,:),times] = make_erp(EEG,{'S  2','S 12','S 22','S 32','S  6','S 16','S 26','S 36','S  7','S 17','S 27','S 37'},[-0.2 .6],[-200 0]);
    save(fullfile(preprocessedFolder,eyeCorrFile),'EEG');
end
save('icsremoved.mat','removed');

%%
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