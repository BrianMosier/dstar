function [test_results,p_values] = roc3d_norm_test(x,D,alpha,plots)

if ~exist('plots','var')
   plots = 'no';
end

x1 = x(D==1);
x2 = x(D==2);
x3 = x(D==3);

[t1,p1] = swtest(x1,alpha);
[t2,p2] = swtest(x2,alpha);
[t3,p3] = swtest(x3,alpha);

p_values = [p1 p2 p3];
test_results = [t1 t2 t3];

for i = 1:3
    if test_results(i) == 0
         test_res(i) = append('Fail to reject normality for group ',string(i),': p-value = ',string(p_values(i)));
    elseif test_results(i) == 1
         test_res(i) = append('Reject normality for group ',string(i),': p-value = ',string(p_values(i)));
    end
end
disp('Results from the Shapiro-Wilk test for normality:')
disp(test_res')
switch plots
    case 'yes'
figure('Position',[0,0,1200,400])
subplot('Position',[.05 .15 .25 .7])
qqplot(x1)
title('Q-Q plot for Group 1')
subplot('Position',[.375 .15 .25 .7])
qqplot(x2)
title('Q-Q plot for Group 2')
subplot('Position',[.70 .15 .25 .7])
qqplot(x3)
title('Q-Q plot for Group 3')
set(gcf,'color','w');
end
end