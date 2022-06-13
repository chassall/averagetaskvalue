% Do the behavioural analysis for the Average Task Value study
%
% Other m-files required: 
% EEGLAB with iclabel extension
% cbrewer.m
% notBoxPlot.m
% boundedline.m (https://github.com/kakearney/boundedline-pkg)

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

% Note that P1, P2 were pilot subjects and will not be included in any
% averages

close all; clear all;

% Set output folder
outputFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\analysis\output';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Participants
allPs = 1:38; % 1-2: pilot, 18: noisy
learners = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26,27,28,31,34,35,36,37,38];
nonLearners = [3 5 9 12 20 24 25 29 30 32 33]; % <60% in mid or high task

num_participants = length(allPs);
all_participant_acc = nan(num_participants,4,200); % allPs X conditions X trials
numTrialsMM = 200;
all_participant_acc_mm = nan(num_participants,2,numTrialsMM); % allPs X conditions X trials
all_participant_rsp = nan(num_participants,4,200);
acc = nan(num_participants,3,200);
win_size = 15;
newAcc = nan(num_participants,4,144);
aveRewards = [];
avePerformance = [];
aveChoice = [];

newNewAcc = nan(length(allPs),144,2);

numReward = [];

for p = 1:num_participants
    
    % Columns differed slightly between p1-2 (pilot) and p3-38
    if allPs(p) == 1 || allPs(p) == 2
        r_index = 9;
        e_index = 10;
        i_index = 11;
        o_index = 13;
        p_index = 14;
    else
        r_index = 10;
        e_index = 11;
        i_index = 12;
        o_index = 14;
        p_index = 15;
    end

    % Set data folder and load beh data
    % Colums: [b t this_block_type thi s_trial_stimulus this_trial_p this_trial_colour this_trial_shape chosen_side responded_early invalid_response response_time*1000 this_trial_a_winner this_trial_optimal];
    pString = ['sub-' num2str( allPs(p),'%0.02i')];
    if p <= 26
        dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-1';
        dString = [pString '_task-casinos_beh.txt'];
        this_data = load(fullfile(dataFolder,pString,'beh',dString));
    else
        dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-2';
        dString = [pString '_task-casinos_beh.tsv'];
        temp_data = tdfread(fullfile(dataFolder,pString,'beh',dString));
        this_data = [temp_data.block temp_data.trial temp_data.task temp_data.cue temp_data.prob temp_data.red temp_data.green temp_data.blue temp_data.shape temp_data.response temp_data.early temp_data.invalid temp_data.rt temp_data.outcome temp_data.optimal];
    end

    numReward(p) = sum(this_data(:,14));

    block_50_i = (this_data(:,3) == 1);
    block_50M_i = (this_data(:,3) == 2) & (this_data(:,4) <= 3 );
    block_80M_i = (this_data(:,3) == 2) & (this_data(:,4) >= 4);
    block_80_i = (this_data(:,3) == 3);
  
    responded_early_i = this_data(:,e_index);
    invalid_response_i = this_data(:,i_index);
    
    valid_i = ~responded_early_i & ~invalid_response_i; 
    
    tempData = this_data;
    tempData(block_50M_i | ~valid_i,p_index) = NaN;
    x =tempData(block_50M_i | block_80M_i,p_index);
    
    tempData = this_data;
    tempData(~valid_i,p_index) = NaN;
    y = tempData(block_80_i,p_index);
    
    if length(x) >= 144
        newNewAcc(p,:,1) = x(1:144);
    else
        newNewAcc(p,1:length(x),1) = x;
    end
    
    if length(y) >= 144
        newNewAcc(p,:,2) = y(1:144);
    else
        newNewAcc(p,1:length(y),2) = y;
    end
    
    data50 = this_data(block_50_i,p_index);
    dataM = this_data(block_50M_i | block_80M_i,p_index);
    data80 = this_data(block_80_i,p_index);
    
    data50(data50 == -1) = 0;
    dataM(dataM == -1) = 0;
    data80(data80 == -1) = 0;
    
    r_50 = this_data(block_50_i & valid_i,r_index);
    r_50M = this_data(block_50M_i & valid_i,r_index);
    r_80M = this_data(block_80M_i & valid_i,r_index);
    r_80 = this_data(block_80_i & valid_i,r_index);
    
    fb_50 = this_data(block_50_i & valid_i,o_index);
    fb_50M = this_data(block_50M_i & valid_i,o_index);
    fb_80M = this_data(block_80M_i & valid_i,o_index);
    fb_80 = this_data(block_80_i & valid_i,o_index);
    fb_M = this_data((block_50M_i | block_80M_i) & valid_i,o_index);
    
    acc_50 = this_data(block_50_i & valid_i,p_index);
    acc_50M = this_data(block_50M_i & valid_i,p_index);
    acc_80M = this_data(block_80M_i & valid_i,p_index);
    acc_80 = this_data(block_80_i & valid_i,p_index);
    
    rsp_50 = this_data(block_50_i & valid_i,r_index);
    rsp_50M = this_data(block_50M_i & valid_i,r_index);
    rsp_80M = this_data(block_80M_i & valid_i,r_index);
    rsp_80 = this_data(block_80_i & valid_i,r_index);
    
    v50 = block_50_i & valid_i;
    v50M = block_50M_i & valid_i;
    v80M = block_80M_i & valid_i;
    v80 = block_80_i & valid_i;
    
    acc(p,1,1:length(data50)) = data50;
    acc(p,2,1:length(dataM)) = dataM;
    acc(p,3,1:length(data80)) = data80;
    
    all_participant_acc(p,1,1:length(acc_50)) = acc_50;
    all_participant_acc(p,2,1:length(acc_50M)) = acc_50M;
    all_participant_acc(p,3,1:length(acc_80M)) = acc_80M;
    all_participant_acc(p,4,1:length(acc_80)) = acc_80;
    
    newAcc(p,1,1:length(acc_50)) = acc_50;
    newAcc(p,2,1:length(acc_50M)) = acc_50M;
    newAcc(p,3,1:length(acc_80)) = acc_80;
    newAcc(p,4,1:length(acc_80M)) = acc_80M;
    
    temp80m = movmean(acc_80M,[0 win_size-1]);
    temp80 = movmean(acc_80, [0 win_size-1]);
    temp80m = interp1(1:length(temp80m),temp80m,linspace(1,length(temp80m),numTrialsMM));
    temp80 = interp1(1:length(temp80),temp80,linspace(1,length(temp80),numTrialsMM));
    all_participant_acc_mm(p,1,:) = temp80m;
    all_participant_acc_mm(p,2,:) = temp80;
    
    all_participant_rsp(p,1,1:length(rsp_50)) = rsp_50;
    all_participant_rsp(p,2,1:length(rsp_50M)) = rsp_50M;
    all_participant_rsp(p,3,1:length(rsp_80M)) = rsp_80M;
    all_participant_rsp(p,4,1:length(rsp_80)) = rsp_80;
    
    aveRewards = [aveRewards; mean(fb_50) mean(fb_M) mean(fb_80)];
    avePerformance = [avePerformance; mean(acc_80M) mean(acc_80)];
    aveChoice = [aveChoice; mean(r_50) mean(r_50M) mean(r_80M) mean(r_80)];
    
end

% Compute the mean performance
ave_acc = nanmean(all_participant_acc,3);
gave_acc = nanmean(ave_acc,1);
gave_sem = std(ave_acc,1)./sqrt(num_participants);
t_value = tinv(0.975,num_participants-1);
gave_ci = gave_sem .* t_value;
payments = 0.018 * numReward';

%% Mean Performance plot

figure;
bar(avePerformance);
hline = refline(0,0.5);
hline.LineStyle = '--';
hline2 = refline(0,0.6);
hline2.LineStyle = '--';
hline2.Color = 'g';
hline3 = refline(0,0.80);
hline3.LineStyle = '--';
hline3.Color = 'b';
ax = gca;
ax.XTick = 1:length(allPs);
ax.XTickLabel = allPs;
save(fullfile(outputFolder,'avePerformance.mat'),'avePerformance');

%% Mean Rewards plot

figure;
bar(aveRewards);
hline = refline(0,0.5);
hline.LineStyle = '--';
hline2 = refline(0,0.6);
hline2.LineStyle = '--';
hline2.Color = 'r';
ax = gca;
ax.XTick = 1:length(allPs);
ax.XTickLabel = allPs;

%% Manuscript Figure 02: Behavioral data

axs = {};
labels = {};

fontSize = 8;
lineWidth = 1.25;

% Subplot settings
gap = [0.10,0.12];
marg_h = [0.24,0.24];
marg_w = [0.1,0.04];

% Set up figure for publication-ready image
makefigure(19,7);

axs{1} = subtightplot(1,3,1,gap,marg_h,marg_w);
[h2,stats2] = notBoxPlot(aveRewards(learners,:).*100,'interval','tInterval');
plotColours = cbrewer('qual','Dark2',3);
formatNBP(h2);
titles{1} = title('(a)');
titles{1}.Units = 'normalized';
titles{1}.Position = [-0.3,1.1,0];
xlabel('Task value');
xticklabels({'Low','Mid','High'});
ylabel('Reward (% wins)');

% Stats
rmANOVA(aveRewards(learners,:));

axs{2} = subtightplot(1,3,2,gap,marg_h,marg_w);
win_size = 30;
mmAcc = 100*movmean(newNewAcc,[0 win_size-1],2,'omitnan','Endpoints','discard');
meanAcc = squeeze(nanmean(mmAcc(learners,:,:),1));
stdAcc = squeeze(nanstd(mmAcc(learners,:,:),[],1));
tval = abs(tinv(0.025,length(learners)-1));
ciAcc = tval * stdAcc ./ sqrt(length(learners));
% blColours = cbrewer('qual','Set1',2);
blColours = plotColours(2:3,:);
[hl, hp] = boundedline(1:length(meanAcc),meanAcc,ciAcc,'alpha','cmap',blColours, 'transparency', 0.05);
for i = 1:length(hl)
    hl(i).LineWidth = 1.5;
end
l = legend(hl,'Mid-value task','High-value task','Box','off','Location','SouthEast');
xlabel('Trial');
ylabel('Performance (% correct)');
titles{2} = title('(b)');
titles{2}.Units = 'normalized';
titles{2}.Position = [-0.3,1.1,0];

axs{3} = subtightplot(1,3,3,gap,marg_h,marg_w);
[h1,stats1] = notBoxPlot(avePerformance(learners,:).*100,'interval','tInterval');
formatNBP(h1);
titles{3} = title('(c)');
titles{3}.Units = 'normalized';
titles{3}.Position = [-0.3,1.1,0];
xlabel('Task value');
xticklabels({'Mid','High'});
ylabel('Performance (% correct)');

% Stats
makeMeans(avePerformance);
rmTTest(avePerformance(learners,1),avePerformance(learners,2));

disp('Mean Rewards');
for i = 1:3
    disp([stats2(i).mu stats2(i).mu - stats2(i).interval stats2(i).mu + stats2(i).interval])
end

disp('Proportion Optimal');
for i = 1:2
    disp([stats1(i).mu stats1(i).mu - stats1(i).interval stats1(i).mu + stats1(i).interval])
end

[H,P,CI,STATS] = ttest(avePerformance(learners,1),avePerformance(learners,2));
disp(P);
disp(STATS);
disp(mean(avePerformance(learners,1)-avePerformance(learners,2))/std(avePerformance(learners,1)-avePerformance(learners,2)));

yLabels = {'Rewards (% wins)','Performance (% correct)'};
xLabels = {{'low','mid','high'},{'mid','high'}};

% Format all axes
for a = 1:length(axs)
   axs{a}.FontSize = fontSize;
   axs{a}.Box = 'off';
end
% 
% % Format all labels (a,b,c...)
for l = 1:length(titles)
   titles{l}.FontSize = 10;
   titles{l}.FontWeight = 'Normal';
end

print(fullfile(outputFolder,'behavioural'),'-dtiff','-r600');