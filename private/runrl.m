function [ll, ps] = runrl(x,data)

alpha = x(1);
tau = x(2);
numTrials =size(data,1);
ll = 0;

verbose = 0;
ps = nan(numTrials,1); % block, trial, probability

% Keep track of last outcomes/actions as we do for WSLS
lastTrialOutcomes = nan(1,6);
lastTrialActions = nan(1,6);

qValues = 0.5 * ones(2,6);

for t = 1:numTrials
    
    thisStim = data(t,4);
    thisAction = data(t,10);
    thisLastOutcome = lastTrialOutcomes(thisStim);
    thisLastAction = lastTrialActions(thisStim);
    
    % Proceed if:
    % - thisAction is 1 or 2 (i.e. an action was taken)
    % - there was a valid "last action" and "last outcome"
    if thisAction && ~isnan(thisLastAction) && ~isnan(thisLastOutcome)

        thisOutcome = data(t,14);
        theseQValues = qValues(:,thisStim);

        % Assign action values using softmax
        softmaxValues = exp(theseQValues./tau) / sum(exp(theseQValues./tau));

        % Determine prob. of taken action
        ps(t) = softmaxValues(thisAction);
        ll = ll + log(ps(t));

        % Update q values
        pe = thisOutcome - theseQValues(thisAction);
        qValues(thisAction,thisStim) = qValues(thisAction,thisStim) + alpha * pe;
    end

    % Update last trial outcome and action as we do with WSLS
    lastTrialOutcomes(thisStim) = data(t,14);
    % Update last action if valid (thisAction == 1 or thisAction == 2)
    if thisAction
        lastTrialActions(thisStim) = thisAction;
    end
    
end


ll = -ll;

% if ll == Inf || isnan(ll)
%     disp('*');
% end

end