% Make ERPs for the Average Task Value study
%
% Other m-files required: EEGLAB with iclabel extension

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Note that P1, P2 were pilot subjects and will not be included in any
% averages

ps = 1:26;
numComponents = ones(26,1);
getTwo = [3 4 14 24 25];
getThree = [6 20];
getFive = [2 13];
numComponents(getTwo) = 2;
numComponents(getThree) = 3;
numComponents(getFive) = 5;
numComponents = numComponents(ps);

% EEG trigger index
% 1 low-value loss
% 2 low-value win
% 3 mid-value (50) loss
% 4 mid-value (50) win
% 5 mid-value (80) loss
% 6 mid-value (80) win
% 7 high-value loss
% 8 high-value win

markers = {{'S  2'},{'S 12'},{'S 22'},{'S 32'},{'S  6'},{'S  7'},{'S 16'},{'S 17'},{'S 26'},{'S 27'},{'S 36'},{'S 37'}};
srate = 250;
nChannels = 30;
numConditions = 12;
conditionNames = {'lowP-uniform','lowP-mixed','highP-mixed','highP-uniform','lv win','lv loss','mv win (50)','mv loss (50)','mv win (80)','mv loss (80)','hv win','hv loss'};
eegInt = [-0.2 0.6];
times = eegInt(1):1/srate:eegInt(2)-1/srate;
erpLength = srate*(eegInt(2)-eegInt(1));
eegBl = [-200 0];
fbChannel = 'FCz';
allERPs = nan(length(ps),numConditions,nChannels,erpLength);
artifactPs = [];

% Participant loop
for pi = 1:length(ps)
    
    whichP = ps(pi);
    load(['./preprocessed/erpclass2016_' num2str(whichP,'%0.04i') '.mat']);
    
    fbChannelI = eeg_chaninds(EEG,{fbChannel});
    
    theseEEG = {};
    for ci = 1:numConditions
        theseEEG{ci} = pop_epoch(EEG,markers{ci},eegInt);
        theseEEG{ci} = pop_rmbase(theseEEG{ci},eegBl);
    end
   
    theseAr = {};
    for ci = 1:numConditions
        theseAr{ci} = find_artifacts(theseEEG{ci}, 150, 150, 40, 0.1);
        artifactPs(pi,ci) = mean(theseAr{ci});
    end
   
    % Remove bad trials
    for ci = 1:numConditions
        if ~all(theseAr{ci})
            theseEEG{ci} = pop_select(theseEEG{ci},'notrial',find(theseAr{ci}));
        else
            theseEEG{ci} = [];
        end
    end
    
    theseERPs = {};
    for ci = 1:numConditions
        if isempty(theseEEG{ci})
            theseERPs{ci} = nan(30,erpLength); % TODO: might need to soft code this
        else
            theseERPs{ci} = mean(theseEEG{ci}.data,3);
        end
    end
    
    for ci = 1:numConditions
        allERPs(pi,ci,:,:) = theseERPs{ci};
    end
    
end

save('erps.mat');

return;

%% Check artifacts
toInclude = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only
figure();
imagesc(artifactPs(toInclude,:)); colorbar();

disp('Cue-locked artifacts');
meanArtifacts = mean(mean(artifactPs(toInclude,1:4),2),1);
disp(meanArtifacts);
tval = abs(tinv(0.025,length(toInclude)-1));
ci = tval*std(mean(artifactPs(toInclude,1:4),2),1)/sqrt(length(toInclude));
disp(meanArtifacts-ci);
disp(meanArtifacts+ci);

disp('Feedback-locked artifacts');
meanArtifacts = mean(mean(artifactPs(toInclude,5:12),2),1);
disp(meanArtifacts);
tval = abs(tinv(0.025,length(toInclude)-1));
ci = tval*std(mean(artifactPs(toInclude,5:12),2),1)/sqrt(length(toInclude));
disp(meanArtifacts-ci);
disp(meanArtifacts+ci);

%% Cue plot

gap = [0.2,0.1];
marg_h = [0.2,0.18];
marg_w = [0.08,0.02];

plotColours = cbrewer('qual','Paired',8);
plotColours = plotColours([2 4 6 8],:);

% toInclude = [3:17 19:26]; % All participants
toInclude = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only

nbpColours = cbrewer('qual','Dark2',4);
fontSize = 8;
lineWidth = 1.25;

tWindow = [.256, .288];
fbChannel = 'FCz'

iWindow = dsearchn(times',tWindow');
fbChannelI = eeg_chaninds(EEG,{fbChannel});

axs = {};
labels = {};
legends = {};
makefigure(14,6);
axs{1} = subtightplot(1,3,[1:2],gap,marg_h,marg_w);
lowWave = squeeze(mean(mean(allERPs(toInclude,[1 2],:,:),1),2));
highWave = squeeze(mean(mean(allERPs(toInclude,[3 4],:,:),1),2));
diffWave = highWave - lowWave;
diffTopo = mean(diffWave(:,iWindow(1):iWindow(2)),2);
for ci = 1:4
    thisGA = squeeze(mean(allERPs(toInclude,ci,fbChannelI,:),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',plotColours(ci,:));
    hold on;
end
hold on;

area(tWindow, [axs{1}.YLim(1) axs{1}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [axs{1}.YLim(2) axs{1}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
legends{1} = legend('Low Task, Low Cue','Mid Task, Low Cue','Mid Task, High Cue','High Task, High Cue','Location','NorthWest','Box','off');
xlabel('Time (s)');
ylabel('Voltage (\muV)');
text(axs{1}.XLim(1) + 0.015,axs{1}.YLim(2),fbChannel,'FontSize',fontSize);
labels{1} = title('(a)');
labels{1}.Units = 'normalized';
labels{1}.Position = [-0.06,1.1,0];

theseERPs = squeeze(mean(allERPs(toInclude,1:4,fbChannelI,iWindow(1):iWindow(2)),4));
between = array2table(theseERPs,'VariableNames',{'V1','V2','V3','V4'});
conditions = [1 1; 1 2; 2 2; 2 1];
pCond = {'1'; '1'; '2'; '2'};
envCond = {'1'; '2'; '2'; '1'};
within = table(pCond,envCond,'VariableNames',{'P','Environment'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P*Environment');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

theseERPs = squeeze(mean(allERPs(toInclude,1:4,fbChannelI,iWindow(1):iWindow(2)),4));
between = array2table(theseERPs,'VariableNames',{'V1','V2','V3','V4'});
conditions = [1; 2; 3; 4];
pCond = {'1'; '2'; '3'; '4'};
within = table(pCond,'VariableNames',{'P'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

rmANOVA(theseERPs);

thisTopo = squeeze(mean(mean(mean(allERPs(toInclude,1:4,:,iWindow(1):iWindow(2)),1),2),4));
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

labels{2} = title('(b)');
labels{2}.Units = 'normalized';
labels{2}.Position = [-0.06,1.1,0];

% format axes
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

print('cues_fcz.tiff','-dtiff','-r300');

%% Cue plot POz

%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 
gap = [0.2,0.1];
marg_h = [0.2,0.18];
marg_w = [0.08,0.02];

plotColours = cbrewer('qual','Paired',8);
plotColours = plotColours([2 4 6 8],:);

plotColours = cbrewer('qual','Paired',10);
reds = plotColours(5:6,:);
blues = plotColours(1:2,:);
purples = plotColours(9:10,:);
greens = plotColours(3:4,:);

% toInclude = [3:17 19:26]; % All participants
toInclude = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only

nbpColours = cbrewer('qual','Dark2',4);
fontSize = 8;
lineWidth = 1.25;

% Time window and electrode
tWindow = [.256, .288];
fbChannel = 'POz'

iWindow = dsearchn(times',tWindow');
fbChannelI = eeg_chaninds(EEG,{fbChannel});

axs = {};
labels = {};
legends = {};
makefigure(14,6);
axs{1} = subtightplot(1,3,[1:2],gap,marg_h,marg_w);
lowWave = squeeze(mean(mean(allERPs(toInclude,[1 2],:,:),1),2));
highWave = squeeze(mean(mean(allERPs(toInclude,[3 4],:,:),1),2));
diffWave = highWave - lowWave;

% Find min of difference wave
[M, minI] = min(diffWave(fbChannelI,:));
disp(times(minI));
peakP = 0.90;
peakV = peakP * M;
[B,I] = sort(abs(diffWave(fbChannelI,:)-peakV));
disp(times(I(1))); disp(times(I(2))); disp(times(I(3)));
whichColours = [greens(1,:); greens(1,:); greens(2,:); greens(2,:)];
whichLines = {'-','--','--',':'};

diffTopo = mean(diffWave(:,iWindow(1):iWindow(2)),2);
for ci = 1:4
    thisGA = squeeze(mean(allERPs(toInclude,ci,fbChannelI,:),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',whichColours(ci,:),'LineStyle',whichLines{ci});
    hold on;
end
plot(times,diffWave(fbChannelI,:),'Color',[.5,.5,.5],'LineWidth',1);
hold on;
area(tWindow, [axs{1}.YLim(1) axs{1}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [axs{1}.YLim(2) axs{1}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
legends{1} = legend('Low Task, Low Cue','Mid Task, Low Cue','Mid Task, High Cue','High Task, High Cue','High-Low Difference','Location','NorthWest','Box','off');
xlabel('Time (s)');
ylabel('Voltage (\muV)');
text(axs{1}.XLim(1) + 0.015,axs{1}.YLim(2),fbChannel,'FontSize',fontSize);
labels{1} = title('(a)');
labels{1}.Units = 'normalized';
labels{1}.Position = [-0.06,1.1,0];
ylim([-2 10]);

theseERPs = squeeze(mean(allERPs(toInclude,1:4,fbChannelI,iWindow(1):iWindow(2)),4));
between = array2table(theseERPs,'VariableNames',{'V1','V2','V3','V4'});
conditions = [1 1; 1 2; 2 2; 2 1];
pCond = {'1'; '1'; '2'; '2'};
envCond = {'1'; '2'; '2'; '1'};
within = table(pCond,envCond,'VariableNames',{'P','Environment'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P*Environment');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

theseERPs = squeeze(mean(allERPs(toInclude,1:4,fbChannelI,iWindow(1):iWindow(2)),4));
between = array2table(theseERPs,'VariableNames',{'V1','V2','V3','V4'});
conditions = [1; 2; 3; 4];
pCond = {'1'; '2'; '3'; '4'};
within = table(pCond,'VariableNames',{'P'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

rmANOVA(theseERPs);

[M,I] = min(diffTopo);
disp(EEG.chanlocs(I).labels);
axs{2} = subtightplot(1,3,3,gap,marg_h,marg_w);
topoColour = cbrewer('seq','OrRd',6);
numContours = 6;
tpLimits = [0 6];
tpLimits = [min(diffTopo) max(diffTopo)];
tpLimits = [-1.4, 1.4];
tp = topoplot(diffTopo,EEG.chanlocs,'numcontour',numContours,'maplimits',tpLimits,'electrodes','off','headrad','rim','shading','interp','whitebk','on','style','fill');
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
colormap(axs{2},topoColour);
set(gca,'color','none');
c = colorbar('Location','South');
c.Position(2) = c.Position(2) - 0.25;
c.Label.String = 'Voltage(\muV)';

labels{2} = title('(b)');
labels{2}.Units = 'normalized';
labels{2}.Position = [-0.06,1.1,0];

% format axes
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

print('cues_oz.tiff','-dtiff','-r300');

%% Cue Difference Waves (unused)
topo1 = squeeze(mean(mean(allERPs(toInclude,1,:,iWindow(1):iWindow(2)),4),1)) - squeeze(mean(mean(allERPs(toInclude,4,:,iWindow(1):iWindow(2)),4),1));
topo2 = squeeze(mean(mean(allERPs(toInclude,2,:,iWindow(1):iWindow(2)),4),1)) - squeeze(mean(mean(allERPs(toInclude,3,:,iWindow(1):iWindow(2)),4),1));

figure();
subplot(1,2,1);
topoplot(topo1,EEG.chanlocs);
subplot(1,2,2);
topoplot(topo2,EEG.chanlocs);

%% Feedback plot (Difference Waves)

fontSize = 8;
lineWidth = 1.25;

%toInclude = [3:17 19:26]; % All participants
%toInclude = [3 5 9 12 20 24 25]; % Non-learners
toInclude = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only

fbChannel = 'FCz';
tWindow = [.240, .276]; % 90% of peak p = .01

fbChannelI = eeg_chaninds(EEG,{fbChannel});
iWindow = dsearchn(times',tWindow');

% Subplot settings
gap = [0.2,0.1];
marg_h = [0.14,0.1];
marg_w = [0.08,0.02];

makefigure(19,11);
axs = {};
labels = {};
legends = {};
axs{1} = subtightplot(2,4,1:3,gap,marg_h,marg_w);
plotColours = cbrewer('qual','Paired',10);
reds = plotColours(5:6,:);
blues = plotColours(1:2,:);
purples = plotColours(9:10,:);

phs = {};
whichConditions = 5:numConditions;
whichColours = [reds(1,:); blues(1,:); reds(1,:); blues(1,:); reds(2,:); blues(2,:); reds(2,:); blues(2,:)];
whichLines = {'-','-','--','--','--','--',':',':'};
for ci = 1:length(whichConditions)
    thisGA = squeeze(mean(allERPs(toInclude,whichConditions(ci),fbChannelI,:),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',whichColours(ci,:),'LineStyle',whichLines{ci});
    hold on;
end
ylim([-5 20]); xlabel('Time (s)'); ylabel('Voltage (\muV)');
legends{1} = legend({'Low-Low Win', 'Low-Low Loss', 'Mid-Low Win','Mid-Low Loss','Mid-High Win','Mid-High Loss','High-High Win','High-High Loss'},'Location','NorthWest','NumColumns',2,'Box','off');
text(axs{1}.XLim(1) + 0.015,axs{1}.YLim(2),fbChannel,'FontSize',fontSize);

labels{1} = title('(a)');
labels{1}.Units = 'normalized';
labels{1}.Position = [-0.06,1.1,0];

differenceWaves(:,1,:,:) = allERPs(:,5,:,:) - allERPs(:,6,:,:);
differenceWaves(:,2,:,:) = allERPs(:,7,:,:) - allERPs(:,8,:,:);
differenceWaves(:,3,:,:) = allERPs(:,9,:,:) - allERPs(:,10,:,:);
differenceWaves(:,4,:,:) = allERPs(:,11,:,:) - allERPs(:,12,:,:);

highDiffDiff = differenceWaves(:,3,:,:) - differenceWaves(:,4,:,:);
highDiffDiffGA = squeeze(mean(mean(highDiffDiff(toInclude,:,fbChannelI,:),1),2));
highDiffDiffTopo = squeeze(mean(mean(mean(highDiffDiff(toInclude,:,:,iWindow(1):iWindow(2)),1),2),4));

topo1 = squeeze(mean(mean(mean(differenceWaves(toInclude,1,:,iWindow(1):iWindow(2)),1),2),4));
topo2 = squeeze(mean(mean(mean(differenceWaves(toInclude,2,:,iWindow(1):iWindow(2)),1),2),4));
topo3 = squeeze(mean(mean(mean(differenceWaves(toInclude,3,:,iWindow(1):iWindow(2)),1),2),4));
topo4 = squeeze(mean(mean(mean(differenceWaves(toInclude,4,:,iWindow(1):iWindow(2)),1),2),4));

diffGrandAve = squeeze(mean(mean(differenceWaves(toInclude,:,fbChannelI,:),1),2));
diffTopo = squeeze(mean(mean(mean(differenceWaves(toInclude,:,:,iWindow(1):iWindow(2)),1),2),4));
[M,I] = max(diffTopo)
EEG.chanlocs(I).labels



axs{2} = subtightplot(2,4,[5 6 7],gap,marg_h,marg_w);
thisDiff1 = squeeze(mean(allERPs(toInclude,5,fbChannelI,:),1)) - squeeze(mean(allERPs(toInclude,6,fbChannelI,:),1));
plot(times,thisDiff1,'LineWidth',lineWidth,'Color',purples(1,:),'LineStyle','-'); hold on;

thisDiff2 = squeeze(mean(allERPs(toInclude,7,fbChannelI,:),1)) - squeeze(mean(allERPs(toInclude,8,fbChannelI,:),1));
plot(times,thisDiff2,'LineWidth',lineWidth,'Color',purples(1,:),'LineStyle','--'); 

thisDiff3 = squeeze(mean(allERPs(toInclude,9,fbChannelI,:),1)) - squeeze(mean(allERPs(toInclude,10,fbChannelI,:),1));
plot(times,thisDiff3,'LineWidth',lineWidth,'Color',purples(2,:),'LineStyle',':'); 

thisDiff4 = squeeze(mean(allERPs(toInclude,11,fbChannelI,:),1)) - squeeze(mean(allERPs(toInclude,12,fbChannelI,:),1));

plot(times,thisDiff4,'LineWidth',lineWidth,'Color',purples(2,:)); 

[M,I] = max(diffGrandAve);
peakP = 0.90;
peakV = peakP * M;
[B,I] = sort(abs(diffGrandAve-peakV));
disp(['Max of great-grand difference wave ' num2str(M)]);
disp(['50% of max of great-grand difference wave ' num2str(.5*M)]);
times(I(1))
times(I(2))
times(I(3))

ylim([-5 10]); xlabel('Time (s)'); ylabel('Voltage (\muV)');
area(tWindow, [axs{1}.YLim(1) axs{1}.YLim(1)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');
area(tWindow, [axs{1}.YLim(2) axs{1}.YLim(2)],'FaceColor','k','FaceAlpha',0.06,'LineStyle','none');

legends{2} = legend('Low Task, Low Cue','Mid Task, Low Cue','Mid Task, High Cue','High Task, High Cue','Location','NorthWest','NumColumns',1,'Box','off');
text(axs{2}.XLim(1) + 0.015,axs{2}.YLim(2),fbChannel,'FontSize',fontSize);

labels{3} = title('(c)');
labels{3}.Units = 'normalized';
labels{3}.Position = [-0.06,1.1,0];

theseERPs = squeeze(max(differenceWaves(toInclude,:,fbChannelI,iWindow(1):iWindow(2)),[],4)); % Diff
topoColors = cbrewer('seq','OrRd',7);
numContour = 6;
tpLimits = [0 6];
axs{3} = subtightplot(2,4,4,gap,marg_h,marg_w)
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

axs{4} = subtightplot(2,4,8,gap,marg_h,marg_w);
nbpColours = plotColours([2 4 6 8],:);
nbp = notBoxPlot(theseERPs); formatNBP(nbp);
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

between = array2table(theseERPs,'VariableNames',{'V1','V2','V3','V4'});
conditions = [1 1; 1 2; 2 2; 2 1];
pCond = {'1'; '1'; '2'; '2'};
envCond = {'1'; '2'; '2'; '1'};
within = table(pCond,envCond,'VariableNames',{'P','Environment'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P*Environment');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

% Means
makeMeans(theseERPs);

% T-Tests
rmTTest(theseERPs(:,1),theseERPs(:,2));
rmTTest(theseERPs(:,3),theseERPs(:,4));

print('rewp.tiff','-dtiff','-r300');

% diff of diffs

%% Correlation
load('avePerformance.mat');

for i = 1:4
    [RHO,PVAL] = corr(mean(avePerformance(toInclude,:),2),theseERPs(:,i))
end

[RHO,PVAL] = corr(mean(avePerformance(toInclude,:),2),theseERPs(:,4)-theseERPs(:,1))

%% Feedback plot (conditional waveforms)
%toInclude = [3:17 19:26]; % All participants
%toInclude = [3 5 9 12 20 24 25]; % Non-learners
toInclude = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only

fbChannel = 'FCz';
tWindow = [.240, .336]; % 80% of "great grand ave diff wave"

fbChannelI = eeg_chaninds(EEG,{fbChannel});
iWindow = dsearchn(times',tWindow');

makefigure(19,10);
subplot(2,3,[1:3]);
plotColours = cbrewer('qual','Paired',8);
phs = {};
for ci = 5:numConditions
    thisGA = squeeze(mean(allERPs(toInclude,ci,fbChannelI,:),1));
    plot(times,thisGA,'LineWidth',lineWidth,'Color',plotColours(ci-4,:));
    hold on;
end
legend(conditionNames(5:12));

scores = squeeze(mean(allERPs(toInclude,5:numConditions,fbChannelI,iWindow(1):iWindow(2)),4));
greatGrandAve = squeeze(mean(allERPs(toInclude,5:numConditions,:,:),2));
greatGrandAveTopo = squeeze(mean(mean(greatGrandAve(:,:,iWindow(1):iWindow(2)),3),1));

subplot(2,3,[4:5]);
nbpColours = cbrewer('qual','Paired',8);
nbp = notBoxPlot(scores); formatNBP(nbp,nbpColours);

legend(diffConditionNames);

subplot(2,3,6)
topoplot(greatGrandAveTopo,EEG.chanlocs);

between = array2table(scores,'VariableNames',{'V1','V2','V3','V4','V5','V6','V7','V8'});
conditions = [1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 2 1; 2 2 2; 2 1 1; 2 1 2];
pCond = {'1'; '1'; '1'; '1'; '2'; '2'; '2'; '2'};
envCond = {'1'; '1'; '2'; '2'; '2'; '2'; '1'; '1'};
fbCond = {'1'; '2'; '1'; '2'; '1'; '2'; '1'; '2'};
within = table(pCond,envCond,fbCond,'VariableNames',{'P','Environment','Feedback'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V8 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P*Environment*Feedback');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));
