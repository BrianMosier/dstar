function [test_results,p_values] = roc3d_norm_test_comb(x,D,alpha,plots)

if ~exist('plots','var')
   plots = 'no';
end

[d1,d2] = size(x);

if d1 < d2
   x = transpose(x);
end

x1 = x(D==1,:);
x2 = x(D==2,:);
x3 = x(D==3,:);


for i = 1:d2
    [t1(i),p1(i)] = swtest(x1(:,i),alpha);
    [t2(i),p2(i)] = swtest(x2(:,i),alpha);
    [t3(i),p3(i)] = swtest(x3(:,i),alpha);
end



p_values = [p1;p2;p3];
test_results = [t1;t2;t3];

for i = 1:3
    for j = 1:d2
         test_res(i+1,j+1) = append('p-value = ',string(p_values(i,j)));
    end
end

test_res(1,1) = 'Group\Marker';

for i = 2:4
    test_res(i,1) = append('Group ',string(i-1));
end

for i = 2:d2+1
    test_res(1,i) = append('Marker ',string(i-1));
end

disp('Results from the Shapiro-Wilk test for normality:')
disp(test_res)
switch plots
    case 'yes'
figure('Position',[0,0,900,800])
for i = 1:d2
    subplot(3,d2,i)
    qqplot(x1(:,i))
    title(append('Q-Q plot for Marker ',string(i),' Group 1'))
    subplot(3,d2,d2+i)
    qqplot(x2(:,i))
    title(append('Q-Q plot for Marker ',string(i),' Group 2'))
    subplot(3,d2,2*d2+i)
    qqplot(x3(:,i))
    title(append('Q-Q plot for Marker ',string(i),' Group 3'))
    set(gcf,'color','w');
end

end
end