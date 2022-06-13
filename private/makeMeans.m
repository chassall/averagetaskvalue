function makeMeans(x)
%MAKEMEANS Summary of this function goes here
%   Detailed explanation goes here


for i = 1:size(x,2)
    [~,~,CI] = ttest(x(:,i));
    disp(['mean: ' num2str(mean(x(:,i))) ', CI: ' num2str(CI')]);
end

end

