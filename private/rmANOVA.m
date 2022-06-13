function rmANOVA(x)

% Normality test
for i = 1:3
    [H, ~, ~] = swtest(x(1,:));
    if H == 1
        disp(['normality violated for var ' num2str(i)]);
    else
        disp(['normality met for var ' num2str(i)]);
    end
end

if size(x,2) == 3
between = array2table(x,'VariableNames',{'V1','V2','V3'});
within = table({'1';'2';'3'},'VariableNames',{'F1'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V3 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable
elseif size(x,2) == 4
between = array2table(x,'VariableNames',{'V1','V2','V3','V4'});
within = table({'1';'2';'3'; '4'},'VariableNames',{'F1'}); % Create a table reflecting the within subject factors 'Attention' and 'TMS' and their levels.
rm = fitrm(between,'V1-V4 ~ 1','WithinDesign',within); % Use ~1 since no betwee-subject variable   
end

[ranovatbl,~,~,~] = ranova(rm,'WithinModel','F1');
etap = ranovatbl.SumSq(3)/(ranovatbl.SumSq(3) + ranovatbl.SumSq(4));
etag = ranovatbl.SumSq(3)/(ranovatbl.SumSq(2) + ranovatbl.SumSq(3) + ranovatbl.SumSq(4));

disp(ranovatbl);
disp(['partial eta squared: ' num2str(etap)]);
disp(['generalized eta squared: ' num2str(etag)]);

end

