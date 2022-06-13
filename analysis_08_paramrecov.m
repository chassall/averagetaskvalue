% Do parameter and model recovery for the Average Task Value study
%
% Other m-files required:
% Optimization toolbox
% makefigure.m
% runrl.m
% sim_rl.m

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com

close all; clear all; clc; 

% Seed random number generator for consistency
rng(2016);

% Set output folder
outputFolder = 'E:\OneDrive - Nexus365\Projects\2016_EEG_Casinos_Hassall\analysis\output';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Model tuning (FMINCON) parameters
options = optimoptions('fmincon','Display','off');
alphaLims = [0.01 1];
tauLims = [0.01 1];
wsLims = [0 1];
lsLims = [0 1];

% Data generation parameters
alphas = [0.01 0.1:0.1:1];
taus = [0.01 0.1:0.1:1];
numPs = 100; % Note that running time is around 6 seconds per participant

% Test model performance
% [simData, simPerf] = sim_rl(100, 0.8,0.1);
% meanPerf = squeeze(mean(mean(simPerf,1),3));
% figure(); plot(meanPerf');

rlParamRecovParams = [];
rlParamRecovLLs = [];
rlModelRecovParams = [];
rlModelRecovLLs = [];

% Loop through different parameter combinations
for i = 1:length(alphas)
    for j = 1:length(taus)
        alpha = alphas(i);
        tau = taus(j);
        [simData, simPerf] = sim_rl(numPs, alpha, tau);
        for p = 1:numPs
            for b = 1:3

                % Generate data
                thisData = squeeze(simData(p,b,:,:));

                % Tune RL model
                xMin = [alphaLims(1) tauLims(1)];
                xMax = [alphaLims(2) tauLims(2)];
                x0 = [0.8 0.1]; % Initial alpha, tau
                f = @(x)runrl(x,thisData);
                [rlParamRecovParams(i,j,p,b,:), rlParamRecovLLs(i,j,p,b)] = fmincon(f,x0,[],[],[],[],xMin,xMax,[],options);

            end
        end
    end
end

save(fullfile(outputFolder,'modelandparamrecov.mat'));

%% Parameter recovery plot
clear subplot; % In case we've overloaded with subtightplot
makefigure(18,15);
fontSize = 8;
tasksStrings = {'Low-Value Task','Mid-Value Task','High-Value Task'};
axs = {};

% Alpha plot - collapse over taus
for b = 1:3 % Loop over tasks
    subplot(3,3,b);
    meanRecovAlphas = squeeze(mean(rlParamRecovParams(:,:,:,b,1),2));
    plot(10*alphas, mean(meanRecovAlphas,2),'LineStyle','none'); hold on;
    xlabel('Actual \alpha');
    ylabel('Recovered \alpha');
    lsl = lsline();
    lsl.LineWidth = 1.25;
    nbp = notBoxPlot(meanRecovAlphas', 10*alphas);
    axs{b} = gca;
    axs{b}.XTickLabels = alphas;
    formatNBP(nbp);
    title(tasksStrings{b});
end

% Tau plot - collapse over alphas
toCompare = [];
for b = 1:3 % Loop over tasks
    subplot(3,3,3+b);
    meanRecovTaus = squeeze(mean(rlParamRecovParams(:,:,:,b,2),1));
    plot(10*taus, mean(meanRecovTaus,2),'LineStyle','none'); hold on;
    xlabel('Actual \tau');
    ylabel('Recovered \tau');
    lsl = lsline();
    lsl.LineWidth = 1.25;
    toCompare(:,b) = meanRecovTaus(9,:)';
    nbp = notBoxPlot(meanRecovTaus',10*taus);
    axs{3+b} = gca;
    axs{3+b}.XTickLabels = taus;
    formatNBP(nbp);
end

% Common axis settings
for i = 1:length(axs)
    axs{i}.Box = 'off';
end

% Model recovery plot (RL versus random)

% Compute -LL for random model
randomLL = -(144 - 6) * log(0.5); % Minus 6 because in the model fitting we don't count the first encounter with each stimulus

rlModelProps = [];
for i = 1:length(alphas)
    for j = 1:length(taus)
        for b = 1:3
            theseRLLLs = squeeze(rlParamRecovLLs(i,j,:,b));
            rlModelProps(i,j,b) = mean(theseRLLLs < randomLL);
        end
    end
end

myColourMap = cbrewer('seq','Blues',127);
myColourMap = flip(myColourMap);
for b = 1:3

    subplot(3,3,6+b);
    imagesc(alphas,taus, rlModelProps(:,:,b),[0 1]);
    ylabel('\alpha');
    xlabel('\tau');
    colormap(myColourMap);
    c = colorbar();
    c.Label.String = 'Proportion LL_R_L > LL_R_a_n_d_o_m';
    set(gca,'YDir','normal');

    axs{6+b} = gca;
end

for i = 1:length(axs)
    axs{i}.Box = 'off';
    axs{i}.FontSize = fontSize;
end

print(fullfile(outputFolder,'modelparamrecov'),'-dtiff','-r600');

%% Unused

% Plot model fit for different alpha, tau combinations
% Shows that a lower temperature is associated with a better model fit
figure();
for b = 1:3
    meanLLs = squeeze(mean(rlParamRecovLLs(:,:,:,b),3));
    subplot(1,3,b);
    imagesc(alphas,taus, meanLLs,[0 100]);
    ylabel('\alpha');
    xlabel('\tau');
    c = colorbar();
    c.Label.String = '-LL';
    set(gca,'YDir','normal');
end