function [simData, simPerf] = sim_rl(numPs,alpha,tau)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

outcomePs = [0.5 0.5 0.5 0.5 0.5 0.5; 0.5 0.5 0.5 0.8 0.8 0.8; 0.8 0.8 0.8 0.8 0.8 0.8];
pWins(:,:,1) = outcomePs;
pWins(:,:,2) = 1-outcomePs;

% alpha = 9.956622622224070e-04;
% tau = 7.155307122748917e-04;

simPerf = [];
simData = [];
for p = 1:numPs
    for b = 1:3
        qVals = [0.5 0.5 0.5 0.5 0.5 0.5; 0.5 0.5 0.5 0.5 0.5 0.5];
        thisBData = [];
        for s = 1:6
            for t = 1:24
                qVal = qVals(:,s);
                softmaxVals = exp(qVal./tau) / sum(exp(qVal./tau));
                if rand < softmaxVals(1)
                    choice = 1;
                    simPerf(p,b,s,t) = 1;
                else
                    choice = 2;
                    simPerf(p,b,s,t) = 0;
                end

                % Rewards
                if rand < pWins(b,s,choice)
                    reward = 1;
                else
                    reward = 0;
                end

                % Update
                pe = reward - qVals(choice,s);
                qVals(choice,s) = qVals(choice,s) + alpha * pe;

                thisBData = [thisBData; nan(1,3) s nan(1,5) choice nan(1,3) reward];
            end

        end

        simData(p,b,:,:) = thisBData;
    end
end

end

