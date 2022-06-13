function P = rmTTest(x1,x2)

% % Means
% [~,~,CI] = ttest(x1);
% disp(['var 1 mean: ' num2str(mean(x1))]);
% disp(['var 1 CI: ' num2str(CI')]);
% 
% [~,~,CI] = ttest(x2);
% disp(['var 2 mean: ' num2str(mean(x2))]);
% disp(['var 2 CI: ' num2str(CI')]);

if nargin < 2
    x2 = x1(:,2);
    x1 = x1(:,1);
end

% Normality test
[H, ~, ~] = swtest(x1);
if H == 1
    disp('normality violated for var 1');
else
    disp('normality met for var 1');
end

[H, ~, ~] = swtest(x2);
if H == 1
    disp('normality violated for var 2');
else
    disp('normality met for var 2');
end

% T-test
[~,P,~,STATS] = ttest(x1-x2);
disp(['t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);

diff = x1 - x2;
mDiff = mean(diff);
sDiff = std(diff);
cohensD = mDiff/sDiff;
disp(['Cohen''s D = ' num2str(cohensD)]);

end

