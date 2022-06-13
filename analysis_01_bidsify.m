% BIDS-ify data for participants 27-38 in the Average Task Value study
% (Ps 1-2 were pilot runs, Ps 3-26 did not consent to sharing their data)
%
% This script adds files and information to a BIDS folder that was created
% by Brain Visions's BV2BIDS command line tool. This only applies to the
% data collected at Site 2 (Oxford).
% 
% Other m-files required: 
% bids-matlab toolbox (https://github.com/bids-standard/bids-matlab)

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% specify BIDS folder
bidsFolder =  '/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-2';

% specify behavioural source folder (these are the raw behavioural files
% that will be bidsified below)
behSourceFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_source/beh';

%% read dataset_description.json, created by BV2BIDS, and modify
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
datasetFile = fullfile(bidsFolder, 'dataset_description.json');
dataset = bids.util.jsondecode(datasetFile);
dataset.DatasetType = 'raw';
dataset.ReferencesAndLinks = {'https://doi.org/10.1101/2021.09.16.460600'};
options.indent = '  '; % Adds white space, easier to read
bids.util.jsonwrite(datasetFile, dataset, options);

%% write participant tsv files, not created by BV2BIDS
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
participantFile = fullfile(bidsFolder, 'participants.tsv');
participants.participant_id = {'sub-27'; 'sub-28'; 'sub-29'; 'sub-30'; 'sub-31'; 'sub-32'; 'sub-33'; 'sub-34'; 'sub-35'; 'sub-36'; 'sub-37'; 'sub-38'};
participants.participant = [27; 28; 29; 30; 31; 32; 33; 34; 35; 36; 37; 38];
participants.date = {'18-Mar-2022 13:55:40'; '21-Mar-2022 14:57:27'; '31-Mar-2022 15:03:22'; '04-Apr-2022 15:47:01'; '06-Apr-2022 11:00:58'; '07-Apr-2022 13:33:53'; '08-Apr-2022 13:10:38'; '08-Apr-2022 15:52:25'; '13-Apr-2022 10:04:33'; '13-Apr-2022 12:21:42'; '25-Apr-2022 14:07:29'; '27-Apr-2022 14:06:19'};
participants.sex = {'F'; 'F'; 'M'; 'F'; 'F'; 'M'; 'F'; 'M'; 'M'; 'M'; 'F'; 'F'};
participants.handedness = {'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'L'; 'R'; 'R'; 'R'; 'R'; 'R'};
participants.age = [24; 30; 22; 24; 23; 52; 56; 32; 42; 77; 68; 46];
participants.order = {'50-50 mixed 80-20'; '80-20 mixed 50-50';	'80-20 50-50 mixed'; 'mixed 80-20 50-50'; 'mixed 50-50 80-20'; '50-50 80-20 mixed'; 'mixed 80-20 50-50'; '50-50 80-20 mixed'; '80-20 50-50 mixed'; 'mixed 80-20 50-50'; '50-50 80-20 mixed'; 'mixed 80-20 50-50'};
participants.stimuli = {'a lemon an orange'; 'a cherry an orange'; 'a cherry a lemon'; 'an orange a lemon'; 'a cherry a lemon'; 'a cherry an orange'; 'a lemon an orange'; 'an orange a lemon'; 'a cherry an orange'; 'a cherry a lemon'; 'a cherry an orange'; 'a cherry an orange'};
bids.util.tsvwrite(participantFile, participants);

%% write participant json file, not created by BV2BIDS
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
participantJSONFile = fullfile(bidsFolder, 'participants.json');
pInfoDesc.participant.Description = 'participant number as input by the tester';
pInfoDesc.date.Description = 'date and time at start of task';
pInfoDesc.sex.Description = 'self-reported sex of participant';
pInfoDesc.sex.Levels.M = 'male';
pInfoDesc.sex.Levels.F = 'female';
pInfoDesc.handedness.Description = 'self-reported handedness of participant';
pInfoDesc.handedness.Levels.L = 'left-handed';
pInfoDesc.handedness.Levels.R = 'right-handed';
pInfoDesc.handedness.Levels.LR = 'ambidextrous';
pInfoDesc.age.Description = 'self-reported age of participant';
pInfoDesc.age.Units = 'years';
pInfoDesc.order.Description = 'task order';
pInfoDesc.stimuli.Description = 'feedback stimuli (win, loss)';
options.indent = '  '; % Adds white space, easier to read
bids.util.jsonwrite(participantJSONFile,pInfoDesc,options);

%% read eeg json for each participant, created by BV2BIDS, and modify
% https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#electroencephalography
whichPs = 27:38;
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    eegJSONFile = fullfile(bidsFolder, pString,'eeg',[pString '_task-casinos_eeg.json']);
    eeg = bids.util.jsondecode(eegJSONFile);
    eeg.InstitutionName = 'University of Oxford';
    eeg.InstitutionAddress = 'Warneford Hospital, Oxford, OX3 7JX';
    eeg.InstitutionalDepartmentName = 'Department of Psychiatry';
    eeg.ManufacturersModelName = 'actiCHamp Plus';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(eegJSONFile, eeg, options);
end

%% write events json for each participant, not created by BV2BIDS
whichPs = 27:38;
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    eventsJSONFile = fullfile(bidsFolder, pString,'eeg',[pString '_task-casinos_events.json']);

    eInfoDesc.trial_type.Description = 'BrainVision event number and type (string)';
    eInfoDesc.value.Description = 'Event value (string)';
    eInfoDesc.value.Levels.S1 = 'Onset of fixation cross at start of trial (low-task, low-cue)';
    eInfoDesc.value.Levels.S11 = 'Onset of fixation cross at start of trial (mid-task, low-cue)';
    eInfoDesc.value.Levels.S21 = 'Onset of fixation cross at start of trial (mid-task, high-cue)';
    eInfoDesc.value.Levels.S31 = 'Onset of fixation cross at start of trial (high-task, high-cue)';
    eInfoDesc.value.Levels.S2 = 'Onset of cue (low-task, low-cue)';
    eInfoDesc.value.Levels.S13 = 'Onset of cue (mid-task, low-cue)';
    eInfoDesc.value.Levels.S23 = 'Onset of cue (mid-task, high-cue)';
    eInfoDesc.value.Levels.S33 = 'Onset of cue (high-task, high-cue)';
    eInfoDesc.value.Levels.S3 = 'Beep (low-task, low-cue)';
    eInfoDesc.value.Levels.S13 = 'Beep (mid-task, low-cue)';
    eInfoDesc.value.Levels.S23 = 'Beep (mid-task, high-cue)';
    eInfoDesc.value.Levels.S33 = 'Beep (high-task, high-cue)';
    eInfoDesc.value.Levels.S4 = 'Valid left response (low-task, low-cue)';
    eInfoDesc.value.Levels.S14 = 'Valid left response (mid-task, low-cue)';
    eInfoDesc.value.Levels.S24 = 'Valid left response (mid-task, high-cue)';
    eInfoDesc.value.Levels.S34 = 'Valid left response(high-task, high-cue)';
    eInfoDesc.value.Levels.S5 = 'Valid right response (low-task, low-cue)';
    eInfoDesc.value.Levels.S15 = 'Valid right response (mid-task, low-cue)';
    eInfoDesc.value.Levels.S25 = 'Valid right response (mid-task, high-cue)';
    eInfoDesc.value.Levels.S35 = 'Valid right response (high-task, high-cue)';
    eInfoDesc.value.Levels.S6 = 'Onset of win feedback (low-task, low-cue)';
    eInfoDesc.value.Levels.S16 = 'Onset of win feedback (mid-task, low-cue)';
    eInfoDesc.value.Levels.S26 = 'Onset of win feedback (mid-task, high-cue)';
    eInfoDesc.value.Levels.S36 = 'Onset of win feedback  (high-task, high-cue)';
    eInfoDesc.value.Levels.S7 = 'Onset of loss feedback (low-task, low-cue)';
    eInfoDesc.value.Levels.S17 = 'Onset of loss feedback (mid-task, low-cue)';
    eInfoDesc.value.Levels.S27 = 'Onset of loss feedback (mid-task, high-cue)';
    eInfoDesc.value.Levels.S37 = 'Onset of loss feedback (high-task, high-cue)';
    eInfoDesc.onset.Description = 'Event onset';
    eInfoDesc.onset.Units = 'milisecond';
    eInfoDesc.duration.Description = 'Event duration';
    eInfoDesc.duration.Units = 'milisecond';
    eInfoDesc.channel.Description = 'Channel number';
    eInfoDesc.channel.Levels.x0 = 'All channels';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(eventsJSONFile,eInfoDesc,options);
end

%% load behavioural data and save as tsv files
colNames = {'block','trial','task','cue','prob','r','g','b','shape','response','early','invalid','rt','outcome','optimal'};
whichPs = 27:38;
for p = 1:length(whichPs)
    pString = num2str(whichPs(p),'%0.02i');
    thisDir = dir(['/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_source/beh/casinos_*_' pString '.txt']);
    thisData = load(fullfile(behSourceFolder,thisDir.name));
    
    % Make a struct out of the behavioural data
    beh.block = thisData(:,1);
    beh.trial = thisData(:,2);
    beh.task = thisData(:,3);
    beh.cue = thisData(:,4);
    beh.prob = round(100*thisData(:,5));
    beh.red = thisData(:,6);
    beh.green = thisData(:,7);
    beh.blue = thisData(:,8);
    beh.shape = thisData(:,9);
    beh.response = thisData(:,10);
    beh.early = thisData(:,11);
    beh.invalid = thisData(:,12);
    beh.rt = thisData(:,13);
    beh.outcome = thisData(:,14);
    beh.optimal = thisData(:,15);
    
    pBIDSString = ['sub-' pString];
    behFolder = fullfile(bidsFolder,pBIDSString,'beh');
    behFile = [pBIDSString '_task-casinos_beh.tsv'];
    bids.util.tsvwrite(fullfile(behFolder,behFile), beh);
end

%% write beh json for each participant
whichPs = 27:38;
for p = 1:length(whichPs)
    pString = ['sub-' num2str(whichPs(p),'%0.02i')];
    behJSONFile = fullfile(bidsFolder, pString,'beh',[pString '_task-casinos_beh.json']);
    bInfoDesc.block.Description = 'Block number (integer)';
    bInfoDesc.trial.Description = 'Trial number (integer)';
    bInfoDesc.task.Description = 'Casino type (integer)';
    bInfoDesc.task.Levels.x1 = 'Low-value';
    bInfoDesc.task.Levels.x2 = 'Mid-value';
    bInfoDesc.task.Levels.x3 = 'High-value';
    bInfoDesc.cue.Description = 'Cue number (integer)';
    bInfoDesc.prob.Description = 'Feedback validity (integer)';
    bInfoDesc.prob.Levels.x50 = 'Outcomes are 50% likely';
    bInfoDesc.prob.Levels.x80 = 'Outcomes are 80% likely';
    bInfoDesc.red.Description = 'Red content of cue colour';
    bInfoDesc.green.Description = 'Green content of cue colour';
    bInfoDesc.blue.Description = 'Blue content of cue colour';
    bInfoDesc.shape.Description = 'Cue shape (integer)';
    bInfoDesc.shape.x1 = 'circle';
    bInfoDesc.shape.x2 = 'triangle';
    bInfoDesc.shape.x3 = 'diamond';
    bInfoDesc.shape.x4 = 'pentagon';
    bInfoDesc.shape.x5 = 'hexagon';
    bInfoDesc.shape.x6 = 'octagon';
    bInfoDesc.response.Description = 'Participant''s response (integer)';
    bInfoDesc.response.Levels.x1 = 'left';
    bInfoDesc.response.Levels.x2 = 'right';
    bInfoDesc.early.Description = 'Early response flag (boolean)';
    bInfoDesc.early.Levels.x0 = 'not too early';
    bInfoDesc.early.Levels.x1 = 'too early';
    bInfoDesc.invalid.Description = 'Invalid/late response flag (boolean)';
    bInfoDesc.invalid.Levels.x0 = 'not too late and valid response';
    bInfoDesc.invalid.Levels.x1 = 'too late or invalid response';
    bInfoDesc.rt.Description = 'Response time (milliseconds)';
    bInfoDesc.outcome.Description = 'Trial outcome (integer)';
    bInfoDesc.outcome.Levels.xneg1 = 'Invalid';
    bInfoDesc.outcome.Levels.x0 = 'Loss';
    bInfoDesc.outcome.Levels.x1 = 'Win';
    bInfoDesc.optimal.Description = 'Optimal response flag (boolean). Which response is optimal was randomly assigned if prob. was 50%';
    bInfoDesc.optimal.Levels.x0 = 'Non-optimal response';
    bInfoDesc.optimal.Levels.x1 = 'Optimal response';
    options.indent = '  '; % Adds white space, easier to read
    bids.util.jsonwrite(behJSONFile,bInfoDesc,options);
end