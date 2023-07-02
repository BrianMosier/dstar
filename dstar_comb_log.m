function [beta,c1,c2,dstar,tcr1,tcr2,tcr3] = dstar_comb_log(x1,x2,x3)

[d1,d2] = size(x1);

if d1 < d2
   x1 = transpose(x1);
   x2 = transpose(x2);
   x3 = transpose(x3);
end

[~,nmarkers] = size(x1);

n1 = length(x1(:,1));
n2 = length(x2(:,1));
n3 = length(x3(:,1));

ind = [ones(n1,1);2*ones(n2,1);3*ones(n3,1)];
x = [x1;x2;x3];
    
[beta,~,~] = mnrfit(x,ind,'model','ordinal');

beta(1:2) = [];
beta = beta/max(abs(beta));

y1 = x1*beta;
y2 = x2*beta;
y3 = x3*beta;

mu1 = mean(y1);
mu2 = mean(y2);
mu3 = mean(y3);

if mu1 > mu2 && mu2 > mu3
    
beta = -beta;

y1 = x1*beta;
y2 = x2*beta;
y3 = x3*beta;

[dstar,c1,c2] =  Empirical_Dist_3D(y1,y2,y3);

tcr1 = sum(y1<=c1)/n1;
tcr2 = sum((y2<=c2).*(y2>c1))/n2;
tcr3 = sum(y3>c2)/n3;

end

