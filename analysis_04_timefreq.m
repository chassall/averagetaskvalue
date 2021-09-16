% Do time-frequency analysis for the Average Task Value study
%
% Other m-files required: EEGLAB with iclabel extension

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Note that P1, P2 were pilot subjects and will not be included in any
% averages

ps = 1:26;
markers = {{'S  2'},{'S 12'},{'S 22'},{'S 32'}};
numConditions = length(markers);
eegInt = [-0.6 1.2];
eegBl = [-200 0];
rawFolder = 'E:\OneDrive - Nexus365\Projects\RPS\data';
prepFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_PSYC574A\Analysis (New New)\preprocessed';
removedTrialCount = [];
for pi = 1:length(ps)
    whichP = ps(pi);
    prepFile = ['./preprocessed/erpclass2016_' num2str(whichP,'%0.04i') '_eye2.mat'];
    tfFile = ['./preprocessed/erpclass2016_' num2str(whichP,'%0.04i') '_tf_long.mat'];
    disp(prepFile);
    load(prepFile);
    
    % Epoching
    theseEEG = {};
    for ci = 1:numConditions
        theseEEG{ci} = pop_epoch(EEG,markers{ci},eegInt);
        theseEEG{ci} = pop_rmbase(theseEEG{ci},eegBl);
    end
    
    % Artifact detection
    theseAr = {};
    for ci = 1:numConditions
        [theseAr{ci},theseArCT{ci}] = find_artifacts(theseEEG{ci}, 150, 150, 40, 0.1);
        removedTrialCount(pi,ci) = mean(theseAr{ci});
    end
    
    % Remove bad trials
    for ci = 1:numConditions
        if ~all(theseAr{ci})
            theseEEG{ci} = pop_select(theseEEG{ci},'notrial',find(theseAr{ci}));
        %else
        %    theseEEG{ci} = [];
        end
    end
    
    % Do Wavelet
    theseWAV = {};
    for ci = 1:numConditions
        theseWAV{ci} = doWAV(theseEEG{ci},markers{ci},[],1,30,60,6);
    end
    
    % Save
    save(tfFile,'theseWAV');
end
save('removedTrialCountTFLong.mat','removedTrialCount');
return;

%% Artifact counts
load('removedTrialCountTFLong.mat');
whichPs = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only
whichArtifacts = mean(removedTrialCount(whichPs,:),2);
makeMeans(whichArtifacts);

%% Load preprocessed TF data
ps = 1:26;
markers = {{'S  2'},{'S 12'},{'S 22'},{'S 32'}};
numConditions = length(markers);
prepFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_PSYC574A\Analysis\preprocessed';
removedTrialCount = [];
allWav = [];
for pi = 1:length(ps)
    
    whichP = ps(pi);
    tfFile = ['./preprocessed/erpclass2016_' num2str(whichP,'%0.04i') '_logtf2.mat'];
    load(tfFile);
    
    for ci = 1:4
        allWav(pi,ci,:,:,:) = theseWAV{ci}.data;
        times = theseWAV{ci}.times;
        chanlocs = theseWAV{ci}.chanlocs;
        frequencies = theseWAV{ci}.frequencies;
    end
    
end

save('tf.mat');

%% Collapsed TF analysis

fontSize = 8;
lineWidth = 1.25;

% whichPs = [3:12 14:17 19:26]; % All Participants
% whichPs = [3 5 9 12 20 24 25]; % Non-learners
whichPs = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners only

tfBaselineS = [-400 -100];
tfBaselineP = dsearchn(times',tfBaselineS');
meanCue = squeeze(mean(allWav(:,:,:,:,:),2));
meanCueBaseline = mean(meanCue(:,:,:,tfBaselineP(1):tfBaselineP(2)),4);
allCuesBaselined = 10*log10( bsxfun(@rdivide, meanCue, meanCueBaseline) );

cueBaseline = mean(allWav(:,:,:,:,tfBaselineP(1):tfBaselineP(2)),5);
cueMeansBasedlined = 10*log10( bsxfun(@rdivide, allWav, cueBaseline) );

whichAllMeans = allCuesBaselined;
conditionMeansBaselined = cueMeansBasedlined;

% Theta
channelStr = 'FCz'; 
channelI = find(strcmp({chanlocs.labels},channelStr));
thisGreatGrandMean = squeeze(mean(whichAllMeans(whichPs,channelI,:,:),1));

allowableFHz = [3 8];
allowableTMs = [1 1000];
thresholdMap = thisGreatGrandMean >=  0;

% Set wavelet time limits (i.e. ignore start/end of interval)
wavXLimMS = [-200,1000];
wavXLimP = dsearchn(times',wavXLimMS');
thresholdMap(:,[1:wavXLimP(1)-1 wavXLimP(2)+1:end]) = 0;

% Set threshold frequency limits
allowableFP = dsearchn(frequencies',allowableFHz');
thresholdMap(1:allowableFP(1)-1,:) = 0;
thresholdMap(allowableFP(2)+1:end,:) = 0;

% Set threshold time limits
allowableTP = dsearchn(times',allowableTMs');
thresholdMap(:,1:allowableTP(1)-1) = 0;
thresholdMap(:,allowableTP(2)+1:end) = 0;

% Get freqs, times for map
[mapFreqsP,mapTimesP] = find(thresholdMap);
mapFreqs = frequencies(mapFreqsP);
mapTimes = times(mapTimesP);

% Mean in cluster for each participant
allPClusters = [];
for p = whichPs
    
    theseValues = [];
    for i = 1:length(mapFreqsP)
        theseValues = [theseValues; squeeze(whichAllMeans(p,:,mapFreqsP(i),mapTimesP(i)))];
    end
    allPClusters = [allPClusters; mean(theseValues,1)];
end
meanCluster = mean(allPClusters,1);

% Find max electrode
tfTopo = meanCluster;
[M,I] = max(tfTopo);
disp(['max ' chanlocs(I).labels]);
[M,I] = min(tfTopo);
disp(['min ' chanlocs(I).labels]);

% Find limits for wavelet plot
minWav = min(min(thisGreatGrandMean(:,wavXLimP(1):wavXLimP(2))));
maxWav = max(max(thisGreatGrandMean(:,wavXLimP(1):wavXLimP(2))));
wavLim = max(abs([minWav maxWav]));
wavLims = [-wavLim,wavLim];
disp(['Great-grand ave max ' num2str(maxWav)]);
disp(['Great-grand ave min ' num2str(minWav)]);

wavColours = cbrewer('div','RdGy',12); wavColours(wavColours > 1) = 1;
wavColours = flip(wavColours);

axs = {};
titles = {};

% Subplot settings
gap = [0.2,0.1];
marg_h = [0.2,0.18];
marg_w = [0.08,0.02];

makefigure(14,6);

axs{1} = subtightplot(1,3,[1:2],gap,marg_h,marg_w);
imagesc(times,frequencies,thisGreatGrandMean); hold on;
[I,J]= find(thresholdMap);
b = bwboundaries(thresholdMap);
for bi = 1:length(b)
    plot(times(b{bi}(:,2)),frequencies(b{bi}(:,1)),'k-','LineWidth',1.5);
end
set(gca,'YDir','normal');
set(gca, 'YScale', 'log');
set(gca,'YTick',[1 2 4 8 20 30]);
caxis(axs{1},wavLims);
axs{1}.XLim = wavXLimMS;
axs{1}.Colormap = wavColours;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
cb = colorbar();
cb.Label.String = 'Power (dB)';
titles{1} = title('(a)');
titles{1}.Units = 'normalized';
titles{1}.Position = [-0.06,1.1,0];

axs{2} = subtightplot(1,3,3,gap,marg_h,marg_w);
tp = topoplot(tfTopo,chanlocs,'style','both','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',wavColours);
tp.Parent.XLim = [-0.6 0.6];
tp.Parent.YLim = [-0.6 0.6];
set(gca,'color','none');
titles{2} = title('(b)');
titles{2}.Units = 'normalized';
titles{2}.Position = [-0.06,1.1,0];
c = colorbar('Location','South');
%c.Position(1) = c.Position(1) + 0.07;
c.Position(2) = c.Position(2) - 0.25;
c.Label.String = 'Power (dB)';

axs{1}.Colormap = wavColours;
axs{2}.Colormap = wavColours;

% Format all axes
for a = 1:length(axs)
   axs{a}.FontSize = fontSize;
   axs{a}.Box = 'off';
end
% 
% % Format all labels (a,b,c...)
for i = 1:length(titles)
   titles{i}.FontSize = 10;
   titles{i}.FontWeight = 'Normal';
end
print('fmtcollapsed.tiff','-dtiff','-r300');

%%
%
[mapFreqsP,mapTimesP] = find(thresholdMap);
mapFreqs = frequencies(mapFreqsP);
mapTimes = times(mapTimesP);

allMaxes = [];
allTopoMaxes = [];
axs = {};
tps = {};

condGAs = [];
condTopos = [];
condScores = [];

testMean = squeeze(mean(conditionMeansBaselined(whichPs,1,channelI,:,:),1));

allWavLims = [];
for whichCond = 1:4
    
    thisMean = squeeze(mean(conditionMeansBaselined(whichPs,whichCond,channelI,:,:),1));
    
    condGAs(whichCond,:,:)= thisMean;
        
    % Mean in cluster for each participant
    allPClusters = [];
    for p = whichPs

        theseValues = [];
        for i = 1:length(mapFreqsP)
            theseValues = [theseValues squeeze(conditionMeansBaselined(p,whichCond,:,mapFreqsP(i),mapTimesP(i)))];
        end
        allPClusters = [allPClusters mean(theseValues,2)];
    
    end
    condScores(:,whichCond) = allPClusters(channelI,:)';
    meanCluster = mean(allPClusters,2);
    
    tfTopo = meanCluster;
    condTopos(whichCond,:) = tfTopo;
    [M,I] = max(tfTopo)
    
    minWav = min(min(thisMean(:,wavXLimP(1):wavXLimP(2))));
    maxWav = max(max(thisMean(:,wavXLimP(1):wavXLimP(2))));
    wavLim = max(abs([minWav maxWav]));
    wavLims = [-wavLim,wavLim]
    allMaxes = [allMaxes wavLim];
    
    allTopoMaxes = [allTopoMaxes max(abs(tfTopo))];
    
end

%
between = array2table(condScores,'VariableNames',{'V1','V2','V3','V4'});
conditions = num2str([1 1; 1 2; 2 2; 2 1]);
pCond = {'1'; '1'; '2'; '2'};
envCond = {'1'; '2'; '2'; '1'};
within = table(pCond,envCond,'VariableNames',{'P','Environment'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no between-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P*Environment');
disp(ranovatbl);
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

between = array2table(condScores,'VariableNames',{'V1','V2','V3','V4'});
conditions = num2str([1; 2; 3; 4]);
pCond = {'1'; '2'; '3'; '4'};
within = table(pCond,'VariableNames',{'P'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no between-subject variable
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','P');
disp(ranovatbl);
% etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

% Means
makeMeans(condScores);

% T-Tests: within the mid-value
rmTTest(condScores(:,2),condScores(:,3));

% T-Tests: high versus low value
rmTTest(condScores(:,1),condScores(:,4));

% Non-learners
allPs = [3:12 14:17 19:26];
learners = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26]; % Learners
nonLearners = [3 5 9 12 20 24 25]; % Non-learners 
nlScores = [];
for whichCond = 1:4

    % Mean in cluster for each participant
    allPClusters = [];
    for p = nonLearners

        theseValues = [];
        for i = 1:length(mapFreqsP)
            theseValues = [theseValues squeeze(conditionMeansBaselined(p,whichCond,:,mapFreqsP(i),mapTimesP(i)))];
        end
        allPClusters = [allPClusters mean(theseValues,2)];
    
    end
    nlScores(:,whichCond) = allPClusters(channelI,:)';
end

% Correlations
disp('FMT-Performance Correlation');
load('avePerformance.mat');
orderedPerformance = [avePerformance(nonLearners,:); avePerformance(learners,:)];
orderedScores = [nlScores; condScores];
labels = [ones(1,length(nonLearners)) 2*ones(1,length(learners))];
[RHO,PVAL] = corr(orderedPerformance(:,1),orderedScores(:,2)) % p80m perf, p50m fmt
[RHO,PVAL] = corr(orderedPerformance(:,1),orderedScores(:,3)) % p80m perf, p80m fmt
[RHO,PVAL] = corr(orderedPerformance(:,2),orderedScores(:,4)) % p80 perf, p80 fmt

%% FMT plots

close all;

axs = {};
labels = {};

fontSize = 7;
lineWidth = 1.25;


% Subplot settings
gap = [0.10,0.08];
marg_h = [0.05,0.1];
marg_w = [0.06,0.14];

titleStrings = {{'Low Task','Low Cue'},{'Mid Task','Low Cue'},{'Mid Task','High Cue'},{'High Task','High Cue'}};
titles = {};
figA = makefigure(14,7);
for whichCond = 1:4
    axs{whichCond} = subtightplot(2,4,whichCond,gap,marg_h,marg_w);
    thisMean= squeeze(condGAs(whichCond,:,:));
    tfTopo = squeeze(condTopos(whichCond,:));
    wavColours = cbrewer('div','RdGy',11); wavColours(wavColours > 1) = 1;
    wavColours = flip(wavColours);
    imagesc(times,frequencies,thisMean); hold on;
    for bi = 1:length(b)
        plot(times(b{bi}(:,2)),frequencies(b{bi}(:,1)),'k-','LineWidth',1.5);
    end
    set(gca,'YDir','normal');
    set(gca, 'YScale', 'log');
    set(gca,'YTick',[1 2 4 8 20 30]);
    caxis(axs{whichCond},wavLims);
    axs{whichCond}.XLim = wavXLimMS;
    axs{whichCond}.Colormap = wavColours;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    titles{whichCond} = title(titleStrings{whichCond});
    if whichCond == 4
        c1 = colorbar('Position',[0.9216    0.6012    0.0168    0.3166]);
        c1.Label.String = 'Power (dB)';
        c1.FontSize = 6;
        c1.Label.FontSize = 6;
    end
 
    tps{whichCond} = subtightplot(2,4,4+whichCond,gap,marg_h,marg_w);
    tp = topoplot(tfTopo,chanlocs,'style','fill','electrodes','off','headrad','rim','shading','interp','whitebk','on','colormap',wavColours);
    tp.Parent.XLim = [-0.6 0.6];
    tp.Parent.YLim = [-0.6 0.6];
    set(gca,'color','none');
    axs{whichCond}.Colormap = wavColours;
    tps{whichCond}.Colormap = wavColours;
    
    if whichCond == 4
        c2 = colorbar('Position',[0.9216    0.0876    0.0168    0.3166]);
        c2.Label.String = 'Power (dB)';
        c2.FontSize = 6;
        c2.Label.FontSize = 6;
    end
    
end

for i = 1:length(titles)
   titles{i}.FontSize = 10;
   titles{i}.FontWeight = 'Normal';
end

maxMax = max(allMaxes);
wavLims = [-maxMax,maxMax];
for whichCond = 1:4
    caxis(axs{whichCond},wavLims);
end

maxTopoMax = max(allTopoMaxes);
topoLims = [-maxTopoMax,maxTopoMax];
for whichCond = 1:4
    caxis(tps{whichCond},topoLims);
end

% Box plots
figB = makefigure(7,6);
plotColours = cbrewer('qual','Paired',8);
nbpColours = plotColours([2 4 5 8],:);
nbp = notBoxPlot(condScores); formatNBP(nbp);
ylabel('FMT Power (dB)');
axs{end+1} = gca;
axs{end}.XTickLabel = {'Low Task\newlineLow Cue','Mid Task\newlineLow Cue','Mid Task\newlineHigh Cue','High Task\newlineHigh Cue'};

% Subplot settings
gap = [0.10,0.14];
marg_h = [0.14,0.2];
marg_w = [0.08,0.04];

figC = makefigure(19,6);
axs{end+1} = subtightplot(1,3,1,gap,marg_h,marg_w);
plot(orderedPerformance(:,1)*100,orderedScores(:,2),'LineStyle','none');
h = lsline(); hold on;
h.Color = 'k';
p1 = plot(avePerformance(learners,1)*100,condScores(:,2),'ko'); hold on;
p2 = plot(avePerformance(nonLearners,1)*100,nlScores(:,2),'k*'); hold on;
legend([p1,p2],'Learners','Non-Learners','Location','NorthWest','Box','off');
xlabel('Performance (% correct)');
ylabel('FMT (dB)');
ylim([-1 3]);

axs{end+1} = subtightplot(1,3,2,gap,marg_h,marg_w);
plot(orderedPerformance(:,1)*100,orderedScores(:,3),'LineStyle','none');
h = lsline(); hold on;
h.Color = 'k';
p1 = plot(avePerformance(learners,1)*100,condScores(:,3),'ko'); hold on;
p2 = plot(avePerformance(nonLearners,1)*100,nlScores(:,3),'k*'); hold on;
legend([p1,p2],'Learners','Non-Learners','Location','NorthWest','Box','off');
xlabel('Performance (% correct)');
ylabel('FMT power (dB)');
ylim([-1 3]);

axs{end+1} = subtightplot(1,3,3,gap,marg_h,marg_w);
plot(orderedPerformance(:,2)*100,orderedScores(:,4),'LineStyle','none');
h = lsline(); hold on;
h.Color = 'k';
p1 = plot(avePerformance(learners,2)*100,condScores(:,4),'ko'); hold on;
p2 = plot(avePerformance(nonLearners,2)*100,nlScores(:,4),'k*'); hold on;
legend([p1,p2],'Learners','Non-Learners','Location','NorthWest','Box','off');
xlabel('Performance (% correct)');
ylabel('FMT power (dB)');
ylim([-1 3]);

% Format all axes
for a = 1:length(axs)
   axs{a}.FontSize = 6;
   axs{a}.Box = 'off';
end

print(figA,'fmt.tiff','-dtiff','-r300');
print(figB,'fmtbar.tiff','-dtiff','-r300');
print(figC,'fmtcorr.tiff','-dtiff','-r300');