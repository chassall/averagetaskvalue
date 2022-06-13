% Model choice behaviour in the Average Task Value study
%
% Other m-files required: 
% Optimization toolbox
% sigstar.m
% notboxplot.m
% makefigure.m
% rmTTest.m
% makeMeans.m
% runrl.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

close all; clear all;

% Seed random number generator for consistency
rng(2016);

% Set output folder
outputFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\analysis\output';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Participant lists (learners, non-learners, etc.)
allPs = 1:38;
allIncludedPs = [3:17 19:38]; % 1-2: pilot, 18: noisy
nonLearners = [3 5 9 12 20 24 25 29 30 32 33];
learners = [4,6,7,8,10,11,13,14,15,16,17,19,21,22,23,26,27,28,31,34,35,36,37,38];

% Set participants to be modelled
whichPs = allPs;

% FMINCON settings
options = optimoptions('fmincon','Display','off');
alphaLims = [0.1 1];
tauLims = [0.01 1];

% Modelled parameters and log-likelihoods
modelParams = []; % participants X model (rl, wsls) X dataset X parameter (tau, alpha)
modelLLs = []; % participants X model (rl, wsls) X dataset X parameter (p(win-stay), p(lost-shift))

% Tune the models for each participant
for p = 1:length(whichPs)

    % Load data
    pString = ['sub-' num2str( whichPs(p),'%0.02i')];
    if p <= 26
        dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-1';
        dString = [pString '_task-casinos_beh.txt'];
        data = load(fullfile(dataFolder,pString,'beh',dString));
    else
        dataFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2016_EEG_Casinos_Hassall/data_site-2';
        dString = [pString '_task-casinos_beh.tsv'];
        temp_data = tdfread(fullfile(dataFolder,pString,'beh',dString));
        data = [temp_data.block temp_data.trial temp_data.task temp_data.cue temp_data.prob temp_data.red temp_data.green temp_data.blue temp_data.shape temp_data.response temp_data.early temp_data.invalid temp_data.rt temp_data.outcome temp_data.optimal];
    end

    % Set indices  
    lowI = (data(:,3) == 1); % Low-value task
    midI = (data(:,3) == 2); % Mid-value task
    mid50I = (data(:,3) == 2) & (data(:,4) <= 3 ); % Mid-value task (low V stim only)
    mid80I = (data(:,3) == 2) & (data(:,4) >= 4); % Mid-value task (high V stim only)
    highI = (data(:,3) == 3); % High-value task
        
    % Data subsets to be modelled (i.e., run tuning for each set of indices)
    whichDatasets = [lowI midI mid50I mid80I highI];

    % Loop through data subsets
    for d = 1:size(whichDatasets,2)

        % This subset
        thisData = data(whichDatasets(:,d),:);

        % Tune RL model
        xMin = [alphaLims(1) tauLims(1)]; % alpha, tau
        xMax = [alphaLims(2) tauLims(2)];
        x0 = [0.8 0.1];
        f = @(x)runrl(x,thisData);
        [modelParams(p,d,:), modelLLs(p,d)] = fmincon(f,x0,[],[],[],[],xMin,xMax,[],options);

    end
end

% Save everything
save(fullfile(outputFolder,'modeloutput.mat'));

%% Plots
psToInclude = learners;

%% Neg. log-likelihood of 144 random choices
randomLL = -(144 - 6) * log(0.5); % Minus 6 because in the model fitting we don't count the first encounter with each stimulus

%% Summary of fit model parameters
fitAlpha = squeeze(modelParams(psToInclude,:,1));
fitTau = squeeze(modelParams(psToInclude,:,2));

[minAlpha, maxAlpha] = bounds(fitAlpha(:))
[minTau, maxTau] = bounds(fitTau(:))

%% Examine model fits and learning rates
subplot = @(m,n,i) subtightplot (m, n, i, [0.1 0.2], [0.2 0.14], [0.1 0.04]);
%     gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%              is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%              relatively large axis. 
%     marg_h  margins in height in normalized units (0...1)
%              or [lower uppper] for different lower and upper margins 
%     marg_w  margins in width in normalized units (0...1)
%              or [left right] for different left and right margins 

fontSize = 8;
labels = {};
makefigure(14,6);

subplot(1,2,1);
toCompare = [modelParams(psToInclude,1,1) modelParams(psToInclude,3,1) modelParams(psToInclude,4,1) modelParams(psToInclude,5,1)];
nbp2 = notBoxPlot(toCompare);
xticklabels({'low-low','mid-low','mid-high','high-high'});
formatNBP(nbp2);
ylabel('\alpha');
labels{1} = title('(a)');
labels{1}.Units = 'normalized';
labels{1}.Position = [-0.2,1.05,0];

disp('alpha, mid-low versus mid-high');
toCompare = [modelParams(psToInclude,3,1), modelParams(psToInclude,4,1)]; % mid-low, mid-high
makeMeans(toCompare);
p1 = rmTTest(toCompare);

disp('alpha, low-low versus high-high');
toCompare = [modelParams(psToInclude,1,1), modelParams(psToInclude,5,1)]; % low-low, high-high
makeMeans(toCompare);
p2 = rmTTest(toCompare);

H = sigstar([2 3], p1); set(H,'LineWidth',0.5);
% H = sigstar([1 4], p2); set(H,'LineWidth',0.5); % not significant
ax = gca;
ax.FontSize = fontSize;
ax.YLabel.FontSize = 12;
ax.YLabel.FontWeight = 'bold';
ax.XTickLabelRotation = 30;

subplot(1,2,2);
toCompare = [modelParams(psToInclude,1,2) modelParams(psToInclude,3,2) modelParams(psToInclude,4,2) modelParams(psToInclude,5,2)];
nbp3 =notBoxPlot(toCompare);
xticklabels({'low-low','mid-low','mid-high','high-high'});
formatNBP(nbp3);
ylabel('\tau');
labels{2} = title('(b)');
labels{2}.Units = 'normalized';
labels{2}.Position = [-0.2,1.05,0];

disp('tau, mid-low versus mid-high');
toCompare = [modelParams(psToInclude,3,2), modelParams(psToInclude,4,2)]; % mid-low, mid-high
makeMeans(toCompare);
p1 = rmTTest(toCompare);

disp('tau, low-low versus high-high');
toCompare = [modelParams(psToInclude,1,2), modelParams(psToInclude,5,2)]; % low-low, high-high
makeMeans(toCompare);
p2 = rmTTest(toCompare);

H = sigstar([2 3], p1); set(H,'LineWidth',0.5);
H = sigstar([1 4], p2); set(H,'LineWidth',0.5);

ax = gca;
ax.FontSize = fontSize;
ax.YLabel.FontSize = 12;
ax.YLabel.FontWeight = 'bold';
ax.XTickLabelRotation = 30;


% Format all labels (a,b,c...)
for l = 1:length(labels)
   labels{l}.FontSize = 10; 
   labels{l}.FontWeight = 'Normal';
end

print(fullfile(outputFolder,'fitlearningrateandtau'),'-dtiff','-r600');

%% Compare model fits between tasks

% Plot settings
fontSize = 8;
makefigure(9,6);

toCompare = [modelLLs(psToInclude,1) 2*modelLLs(psToInclude,3) 2*modelLLs(psToInclude,4) modelLLs(psToInclude,5)];
nbp1 = notBoxPlot(toCompare); hold on;
ylabel('-Log-Likelihood');
formatNBP(nbp1);
randLine = yline(randomLL);
randLine.LineStyle = ':';
t = text(4.5,1.05*randomLL,'Baseline','FontSize',fontSize,'FontAngle','italic'); 
xticklabels({'low-low','mid-low','mid-high','high-high'});

% Compare model fits
disp('model fits, mid-low versus mid-high');
toCompare = [modelLLs(psToInclude,3), modelLLs(psToInclude,4)]; % mid-low, mid-high
makeMeans(toCompare);
p1 = rmTTest(toCompare);

disp('model fits, low-low versus high-high');
toCompare = [modelLLs(psToInclude,1), modelLLs(psToInclude,5)]; % low-low, high-high
makeMeans(toCompare);
p2 = rmTTest(toCompare);

H = sigstar([2 3], p1); set(H,'LineWidth',0.5);
H = sigstar([1 4], p2); set(H,'LineWidth',0.5);
arrowText = text(5,15,'\leftarrow Better model fit','FontSize',8,'FontAngle','italic');
arrowText.Rotation = 90;

ax = gca;
ax.XTickLabelRotation = 30;
ax.FontSize = fontSize;

print(fullfile(outputFolder,'modelfits.tif'),'-dtiff','-r600');

%% Unused

% Scatter plots of -LL and tau
% Shows that a lower temperature is associated with a better model fit
figure();
whichSets = [1 2 5]; % low, mid, high tasks
for i = 1:length(whichSets)
    subplot(1,3,i);
    scatter(modelParams(psToInclude,whichSets(i),2), modelLLs(psToInclude,whichSets(i)));
    xlabel('\tau');
    ylabel('-LL');
    refline;
end