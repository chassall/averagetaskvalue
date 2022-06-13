% Make and analyze ERPs for the Average Task Value study
%
% Other m-files required: 
% EEGLAB with iclabel extension
% subtightplot.m https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
% swtest.m
% cbrewer.m
% makefigure.m
% rmANOVA.m 
% notBoxPlot.m (https://github.com/raacampbell/notBoxPlot)
% formatNBP.m
% makeMeans.m
% rmTTest.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Note that P1, P2 were pilot subjects and will not be included in any
% averages

eeglab; clear all;

% Set output folder
outputFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\analysis\output';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Participants
allPs = 1:38;
allIncludedPs = [3:17 19:38]; % 1-2: pilot, 18: noisy
learners = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26,27,28,31,34,35,36,37,38];
nonLearners = [3 5 9 12 20 24 25 29 30 32 33];

% numPs = length(allPs);
% numComponents = ones(numPs,1);
% getTwo = [3 4 14 24 25];
% getThree = [6 20];
% getFive = [2 13];
% numComponents(getTwo) = 2;
% numComponents(getThree) = 3;
% numComponents(getFive) = 5;
% numComponents = numComponents(allPs);

% EEG trigger index
% 1 low-value loss
% 2 low-value win
% 3 mid-value (50%) loss
% 4 mid-value (50%) win
% 5 mid-value (80%) loss
% 6 mid-value (80%) win
% 7 high-value loss
% 8 high-value win
triggers = {{'S  2'},{'S 12'},{'S 22'},{'S 32'},{'S  6'},{'S  7'},{'S 16'},{'S 17'},{'S 26'},{'S 27'},{'S 36'},{'S 37'},{'S  6','S 16','S 26','S 36'},{'S  7','S 17','S 27','S 37'}};

srate = 250;
nChannels = 30;
conditionNames = {'lowP-uniform','lowP-mixed','highP-mixed','highP-uniform','Low-Low Win','Low-Low Loss','Mid-Low Win','Mid-Low Loss','Mid-High Win','Mid-High Loss','High-High Win','High-High Loss','all win','all loss'};
numConditions = length(conditionNames);
eegInt = [-0.2 0.6];
times = eegInt(1):1/srate:eegInt(2)-1/srate;
erpLength = srate*(eegInt(2)-eegInt(1));
eegBl = [-200 0];
fbChannel = 'FCz';
allERPs = nan(length(allPs),numConditions,nChannels,erpLength);
artifactPs = [];

% Exclude first few trials of each task
numTrialsToExclude = 10;

% Participant loop
for pi = 1:length(allPs)

    % Set data folder and load data
    if allPs(pi) < 27
        dataFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\data_site-1';
    else
        dataFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\data_site-2';
    end
    pString = ['sub-' num2str( allPs(pi),'%0.02i')];
    preprocessedFolder = [dataFolder '/derivatives/' pString];
    eyeCorrFile = [pString '_task-casinos_prepoc.mat'];
    load(fullfile(preprocessedFolder, eyeCorrFile));

    fbChannelI = eeg_chaninds(EEG,{fbChannel});

    trialCount = [0 0 0];

    for i = 1:length(EEG.event)
        
        switch EEG.event(i).type


            case 'S  1'
                
                if trialCount(1) < numTrialsToExclude
                    
                    EEG.event(i).type = 'excluded';

                    % Rename all events associate with this trial
                    j = 1;
                    while ~any(strcmp(EEG.event(i+j).type,{'S  1','S 11','S 21','S 31'}))
                        EEG.event(i+j).type = 'excluded';
                        j = j + 1;
                    end

                end
                trialCount(1) = trialCount(1) + 1;

            case {'S 11','S 21'}

                if trialCount(2) < numTrialsToExclude
                    
                    EEG.event(i).type = 'excluded';

                    % Rename all events associate with this trial
                    j = 1;
                    while ~any(strcmp(EEG.event(i+j).type,{'S  1','S 11','S 21','S 31'}))
                        EEG.event(i+j).type = 'excluded';
                        j = j + 1;
                    end

                end
                trialCount(2) = trialCount(2) + 1;
                

            case 'S 31'

                if trialCount(3) < numTrialsToExclude
                    
                    EEG.event(i).type = 'excluded';

                    % Rename all events associate with this trial
                    j = 1;
                    while ~any(strcmp(EEG.event(i+j).type,{'S  1','S 11','S 21','S 31'}))
                        EEG.event(i+j).type = 'excluded';
                        j = j + 1;
                    end

                end
                trialCount(3) = trialCount(3) + 1;

        end

    end

    theseEEG = {};
    for ci = 1:numConditions
        theseEEG{ci} = pop_epoch(EEG,triggers{ci},eegInt);
        theseEEG{ci} = pop_rmbase(theseEEG{ci},eegBl);
    end
   
    % disp(toExclude);
    theseAr = {};
    for ci = 1:numConditions
        theseAr{ci} = find_artifacts(theseEEG{ci}, 150, 150, 40, 0.1);
        artifactPs(pi,ci) = mean(theseAr{ci});
    end
   

    % Remove bad trials
    for ci = 1:numConditions
        toRemove = theseAr{ci};
        if ~all(toRemove)
            theseEEG{ci} = pop_select(theseEEG{ci},'notrial',find(toRemove));
        else
            theseEEG{ci} = [];
        end
    end
    
    theseERPs = {};
    for ci = 1:numConditions
        if isempty(theseEEG{ci})
            theseERPs{ci} = nan(30,erpLength);
        else
            theseERPs{ci} = mean(theseEEG{ci}.data,3);
        end
    end
    
    for ci = 1:numConditions
        allERPs(pi,ci,:,:) = theseERPs{ci};
    end
    
end

% Save everything
save(fullfile(outputFolder,'erps.mat'));

%% Check artifacts, display artifact proportions

% Rough plot of artifact proportions for all participants and conditions
figure();
imagesc(artifactPs(learners,:)); colorbar();

disp('Cue-locked artifacts');
meanArtifacts = mean(mean(artifactPs(learners,1:4),2),1);
disp(meanArtifacts);
tval = abs(tinv(0.025,length(learners)-1));
ci = tval*std(mean(artifactPs(learners,1:4),2),1)/sqrt(length(learners));
disp(meanArtifacts-ci);
disp(meanArtifacts+ci);

disp('Feedback-locked artifacts');
meanArtifacts = mean(mean(artifactPs(learners,5:12),2),1);
disp(meanArtifacts);
tval = abs(tinv(0.025,length(learners)-1));
ci = tval*std(mean(artifactPs(learners,5:12),2),1)/sqrt(length(learners));
disp(meanArtifacts-ci);
disp(meanArtifacts+ci);

%% Is there a cue-locked RewP? Supplementary Figure 1.

% Electrode/time window
fbChannel = 'FCz';
tWindow = [0.240 0.340];
iWindow = dsearchn(times',tWindow');
fbChannelI = eeg_chaninds(EEG,{fbChannel});

% Make ERP plot
gap = [0.2,0.1];
marg_h = [0.2,0.18];
marg_w = [0.08,0.02];
plotColours = cbrewer('qual','Paired',8);
plotColours = plotColours([2 4 6 8],:);
nbpColours = cbrewer('qual','Dark2',4);
fontSize = 8;
lineWidth = 1.25;
axs = {};
labels = {};
legends = {};
makefigure(14,6);
axs{1} = subtightplot(1,3,[1:2],gap,marg_h,marg_w);
lowWave = squeeze(mean(mean(allERPs(learners,[1 2],:,:),1),2));
highWave = squeeze(mean(mean(allERPs(learners,[3 4],:,:),1),2));
diffWave = highWave - lowWave;
diffTopo = mean(diffWave(:,iWindow(1):iWindow(2)),2);
for ci = 1:4
    thisGA = squeeze(mean(allERPs(learners,ci,fbChannelI,:),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',plotColours(ci,:));
    hold on;
end
hold on;
axs{1}.FontSize = fontSize;
axs{1}.Box = 'off';
area(tWindow, [axs{1}.YLim(1) axs{1}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [axs{1}.YLim(2) axs{1}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
legends{1} = legend('Low Task, Low Cue','Mid Task, Low Cue','Mid Task, High Cue','High Task, High Cue','Location','NorthWest','Box','off');
xlabel('Time (s)');
ylabel('Voltage (\muV)');
text(axs{1}.XLim(1) + 0.015,axs{1}.YLim(2),fbChannel,'FontSize',fontSize);
labels{1} = title('(a)');
labels{1}.Units = 'normalized';
labels{1}.Position = [-0.06,1.1,0];

% Compute ERP scores and do stats
theseERPs = squeeze(mean(allERPs(learners,1:4,fbChannelI,iWindow(1):iWindow(2)),4));
rmANOVA(theseERPs);

% Make topo plot
thisTopo = squeeze(mean(mean(mean(allERPs(learners,1:4,:,iWindow(1):iWindow(2)),1),2),4));
[M,I] = max(thisTopo);
disp(EEG.chanlocs(I).labels);
axs{2} = subtightplot(1,3,3,gap,marg_h,marg_w);
topoColour = cbrewer('seq','OrRd',6);
numContours = 6;
tpLimits = [0 6];
tp = topoplot(thisTopo,EEG.chanlocs,'numcontour',numContours,'maplimits',tpLimits,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
colormap(axs{2},topoColour);
set(gca,'color','none');
c = colorbar('Location','South');
c.Position(2) = c.Position(2) - 0.25;
c.Label.String = 'Voltage(\muV)';
axs{2}.FontSize = fontSize;
axs{2}.Box = 'off';
labels{2} = title('(b)');
labels{2}.Units = 'normalized';
labels{2}.Position = [-0.06,1.1,0];

% Format all labels (a,b,c...)
for l = 1:length(labels)
   labels{l}.FontSize = 10;
   labels{l}.FontWeight = 'Normal';
end

% Format all legends
for i = 1:length(legends)
   legends{i}.FontSize = 6; 
end

print(fullfile(outputFolder,'cues_fcz.tiff'),'-dtiff','-r300');

%% Feedback plot (Difference Waves)

% Which participants?
toInclude  = learners; 

% Analysis electrode and window
fbChannels = {'FCz',}; % Sambrook and Goslin, 2015
fbChannelsI = eeg_chaninds(EEG,fbChannels);
tWindow = [.240, .340];  % Sambrook and Goslin, 2015
iWindow = dsearchn(times',tWindow');

% Plot conditional waveform for each condition
fontSize = 8;
lineWidth = 1.25;
gap = [0.2,0.1];
marg_h = [0.14,0.1];
marg_w = [0.08,0.02];
plotColours = cbrewer('qual','Paired',10);
reds = plotColours(5:6,:);
blues = plotColours(1:2,:);
purples = plotColours(9:10,:);
makefigure(19,11);
axs = {};
labels = {};
legends = {};
axs{1} = subtightplot(2,4,1:3,gap,marg_h,marg_w);
whichConditions = 5:12;
whichColours = [reds(1,:); blues(1,:); reds(1,:); blues(1,:); reds(2,:); blues(2,:); reds(2,:); blues(2,:)];
whichLines = {'-','-','--','--','--','--',':',':'};
for ci = 1:length(whichConditions)
    thisGA = squeeze(mean(mean(allERPs(toInclude,whichConditions(ci),fbChannelsI,:),3),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',whichColours(ci,:),'LineStyle',whichLines{ci});
    hold on;
end
ylim([-5 20]); xlabel('Time (s)'); ylabel('Voltage (\muV)');
legends{1} = legend({'Low-Low Win', 'Low-Low Loss', 'Mid-Low Win','Mid-Low Loss','Mid-High Win','Mid-High Loss','High-High Win','High-High Loss'},'Location','NorthWest','NumColumns',2,'Box','off');
text(axs{1}.XLim(1) + 0.015,axs{1}.YLim(2)+1,fbChannels,'FontSize',fontSize);
labels{1} = title('(a)');
labels{1}.Units = 'normalized';
labels{1}.Position = [-0.06,1.1,0];

% Compute difference waves
differenceWaves(:,1,:,:) = allERPs(:,5,:,:) - allERPs(:,6,:,:);
differenceWaves(:,2,:,:) = allERPs(:,7,:,:) - allERPs(:,8,:,:);
differenceWaves(:,3,:,:) = allERPs(:,9,:,:) - allERPs(:,10,:,:);
differenceWaves(:,4,:,:) = allERPs(:,11,:,:) - allERPs(:,12,:,:);

% Compute difference scores for each group (for correlation)
learnerDiffERPs = squeeze(max(mean(differenceWaves(learners,:,fbChannelsI,iWindow(1):iWindow(2)),3),[],4));
nonLearnerDiffERPs = squeeze(max(mean(differenceWaves(nonLearners,:,fbChannelsI,iWindow(1):iWindow(2)),3),[],4));

% Compute topo scores
topo1 = squeeze(mean(mean(mean(differenceWaves(toInclude,1,:,iWindow(1):iWindow(2)),1),2),4));
topo2 = squeeze(mean(mean(mean(differenceWaves(toInclude,2,:,iWindow(1):iWindow(2)),1),2),4));
topo3 = squeeze(mean(mean(mean(differenceWaves(toInclude,3,:,iWindow(1):iWindow(2)),1),2),4));
topo4 = squeeze(mean(mean(mean(differenceWaves(toInclude,4,:,iWindow(1):iWindow(2)),1),2),4));
diffTopo = squeeze(max(mean(mean(differenceWaves(toInclude,:,:,iWindow(1):iWindow(2)),1),2),[],4));

% Plot difference waves
axs{2} = subtightplot(2,4,[5 6 7],gap,marg_h,marg_w);
thisDiff1 = squeeze(mean(mean(allERPs(toInclude,5,fbChannelsI,:),3),1)) - squeeze(mean(mean(allERPs(toInclude,6,fbChannelsI,:),3),1));
plot(times,thisDiff1,'LineWidth',lineWidth,'Color',purples(1,:),'LineStyle','-'); hold on;
thisDiff2 = squeeze(mean(mean(allERPs(toInclude,7,fbChannelsI,:),3),1)) - squeeze(mean(mean(allERPs(toInclude,8,fbChannelsI,:),3),1));
plot(times,thisDiff2,'LineWidth',lineWidth,'Color',purples(1,:),'LineStyle','--'); 
thisDiff3 = squeeze(mean(mean(allERPs(toInclude,9,fbChannelsI,:),3),1)) - squeeze(mean(mean(allERPs(toInclude,10,fbChannelsI,:),3),1));
plot(times,thisDiff3,'LineWidth',lineWidth,'Color',purples(2,:),'LineStyle',':'); 
thisDiff4 = squeeze(mean(mean(allERPs(toInclude,11,fbChannelsI,:),3),1)) - squeeze(mean(mean(allERPs(toInclude,12,fbChannelsI,:),3),1));
plot(times,thisDiff4,'LineWidth',lineWidth,'Color',purples(2,:)); 
ylim([-5 10]); xlabel('Time (s)'); ylabel('Voltage (\muV)');
area(tWindow, [axs{1}.YLim(1) axs{1}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [axs{1}.YLim(2) axs{1}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
legends{2} = legend('Low Task, Low Cue','Mid Task, Low Cue','Mid Task, High Cue','High Task, High Cue','Location','NorthWest','NumColumns',1,'Box','off');
text(axs{2}.XLim(1) + 0.015,axs{2}.YLim(2)+1,fbChannels,'FontSize',fontSize);
labels{3} = title('(c)');
labels{3}.Units = 'normalized';
labels{3}.Position = [-0.06,1.1,0];

% Topo plot
topoColors = cbrewer('seq','OrRd',7);
numContour = 6;
tpLimits = [0 6];
axs{3} = subtightplot(2,4,4,gap,marg_h,marg_w);
tp = topoplot(diffTopo,EEG.chanlocs,'numcontour',numContour,'maplimits',tpLimits,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
colormap(axs{3},topoColors);
c = colorbar('Location','South');
c.Position(2) = c.Position(2) - 0.12;
c.Label.String = 'Voltage(\muV)';
labels{2} = title('(b)');
labels{2}.Units = 'normalized';
labels{2}.Position = [-0.4,1.1,0];

% Box plot of difference scores
axs{4} = subtightplot(2,4,8,gap,marg_h,marg_w);
nbpColours = plotColours([2 4 6 8],:);
nbp = notBoxPlot(learnerDiffERPs); formatNBP(nbp);
ylabel('Voltage (\muV)');
axs{4}.XTickLabel = {'Low-Low','Mid-Low','Mid-High','High-High'};
axs{4}.XTickLabelRotation = 45;
labels{4} = title('(d)');
labels{4}.Units = 'normalized';
labels{4}.Position = [-0.4,1.1,0];

% Format all axes
for a = 1:length(axs)
   axs{a}.FontSize = fontSize;
   axs{a}.Box = 'off';
end

% Format all labels (a,b,c...)
for l = 1:length(labels)
   labels{l}.FontSize = 10; 
   labels{l}.FontWeight = 'Normal';
end

% Format all legends
for i = 1:length(legends)
   legends{i}.FontSize = 6; 
end

% Means
makeMeans(learnerDiffERPs);

% T-Tests
rmTTest(learnerDiffERPs(:,2),learnerDiffERPs(:,1));
rmTTest(learnerDiffERPs(:,3),learnerDiffERPs(:,4));

print(fullfile(outputFolder,'rewp.tiff'),'-dtiff','-r600');

% Look at win-only or loss-only effects (for response to reviewers)
winWaves(:,1,:,:) = allERPs(:,5,:,:);
winWaves(:,2,:,:) = allERPs(:,7,:,:);
winWaves(:,3,:,:) = allERPs(:,9,:,:);
winWaves(:,4,:,:) = allERPs(:,11,:,:);
lossWaves(:,1,:,:) = allERPs(:,6,:,:);
lossWaves(:,2,:,:) = allERPs(:,8,:,:);
lossWaves(:,3,:,:) = allERPs(:,10,:,:);
lossWaves(:,4,:,:) = allERPs(:,12,:,:);
learnerWinERPs = squeeze(mean(mean(winWaves(learners,:,fbChannelsI,iWindow(1):iWindow(2)),3),4));
learnerLossERPs = squeeze(mean(mean(lossWaves(learners,:,fbChannelsI,iWindow(1):iWindow(2)),3),4));

% Win t-tests (unused)
disp('Win-only scores');
makeMeans(learnerWinERPs);
rmTTest(learnerWinERPs(:,2),learnerWinERPs(:,1));
rmTTest(learnerWinERPs(:,3),learnerWinERPs(:,4));

% Loss t-tests (unused)
disp('Loss-only scores');
makeMeans(learnerLossERPs);
rmTTest(learnerLossERPs(:,2),learnerLossERPs(:,1));
rmTTest(learnerLossERPs(:,3),learnerLossERPs(:,4));

%% Correlation between RewP and behaviour

% Compute relevant differences
nonlearnerERPDiff = nonLearnerDiffERPs(:,3) - nonLearnerDiffERPs(:,4);
learnerERPDiff = learnerDiffERPs(:,3) - learnerDiffERPs(:,4);

% Load average performance
load(fullfile(outputFolder,'avePerformance.mat'));
nonlearnerBeh = mean(avePerformance(nonLearners,:),2);
learnerBeh = mean(avePerformance(learners,:),2);

% Compute correlation
allBeh = [nonlearnerBeh; learnerBeh];
allERP = [nonlearnerERPDiff; learnerERPDiff];
[r,p] = corr(allBeh,allERP);
disp(r);
disp(p);

% Make the figure
fontSize = 8;
makefigure(9,7);
plot(allBeh,allERP,'LineStyle','none');
h = lsline(); hold on;
h.Color = 'k';
p1 = plot(nonlearnerBeh,nonlearnerERPDiff ,'k*'); hold on;
p2 = plot(learnerBeh,learnerERPDiff ,'ko');
legend([p1,p2],'Non-Learners','Learners','Location','NorthWest','Box','off');
xlabel('Performance (% correct)')
ylabel('\Delta RewP (\muV)');
ax = gca;
ax.Box = 'off';
ax.FontSize = fontSize;

print(fullfile(outputFolder,'rewpbehcorr.tiff'),'-dtiff','-r600');

%% Feedback plot (conditional waveforms) - Response to Reviewers

% Which participants?
toInclude = learners;

% Plot settings
nbpColours = cbrewer('qual','Dark2',4);
fontSize = 8;
lineWidth = 1.25;

% Electrode/window settings
fbChannel = 'POz';
tWindow = [.300, .380];
fbChannel = 'Pz';
tWindow = [.400, .600];
fbChannelI = eeg_chaninds(EEG,{fbChannel});
iWindow = dsearchn(times',tWindow');

% Make the figure
makefigure(19,16);
subplot(2,3,1:3);
plotColours = cbrewer('qual','Paired',8);
phs = {};
for ci = 5:12
    thisGA = squeeze(mean(allERPs(toInclude,ci,fbChannelI,:),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',plotColours(ci-4,:));
    hold on;
end
ylabel('Voltage (\muV)');
xlabel('Time (s)');
ax = gca;
ax.Box = 'off';
text(-0.18,15,fbChannel,'FontSize',fontSize);
hold on;
area(tWindow, [ax.YLim(1) ax.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [ax.YLim(2) ax.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
l = legend(conditionNames(5:12));
l.Box = 'off';
l.Location = 'EastOutside';

% ERP and topo scores
scores = squeeze(mean(allERPs(toInclude,5:12,fbChannelI,iWindow(1):iWindow(2)),4));
greatGrandAve = squeeze(mean(allERPs(toInclude,5:12,:,:),2));
greatGrandAveTopo = squeeze(mean(mean(greatGrandAve(:,:,iWindow(1):iWindow(2)),3),1));

% Boxplot
subplot(2,3,4:5);
nbpColours = cbrewer('qual','Paired',8);
nbp = notBoxPlot(scores); formatNBP(nbp,nbpColours);
xticklabels(conditionNames(5:12));
ylabel('P300 (\muV)');

% Topo
subplot(2,3,6)
topoColors = cbrewer('seq','OrRd',7);
numContour = 6;
tpLimits = [0 15];
tp = topoplot(greatGrandAveTopo,EEG.chanlocs,'numcontour',numContour,'maplimits',tpLimits,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
colormap(topoColors);

% Stats
between = array2table(scores,'VariableNames',{'V1','V2','V3','V4','V5','V6','V7','V8'});
conditions = [1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 2 1; 2 2 2; 2 1 1; 2 1 2];
pCond = {'1'; '1'; '1'; '1'; '2'; '2'; '2'; '2'};
envCond = {'1'; '1'; '2'; '2'; '2'; '2'; '3'; '3'};
fbCond = {'1'; '2'; '1'; '2'; '1'; '2'; '1'; '2'};
within = table(pCond,envCond,fbCond,'VariableNames',{'trialExpect','taskVal','valence'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V8 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','trialExpect*taskVal*valence');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

print(fullfile(outputFolder,'p300.tiff'),'-dtiff','-r600');