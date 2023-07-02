function [beta,c1,c2,dstar] = dstar_comb_ker(x1,x2,x3)

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

if nmarkers == 2
% Beta1 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[1 c(3)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[1 c(3)]') iqr(x1*[1 c(3)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[1 c(3)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3)]') iqr(x2*[1 c(3)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[1 c(3)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3)]') iqr(x2*[1 c(3)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[1 c(3)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[1 c(3)]') iqr(x3*[1 c(3)]')/1.34])*n3^(-0.2)))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_p1 = [1 c_p1(3)];      
       
       
% Beta1 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[-1 c(3)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[-1 c(3)]') iqr(x1*[-1 c(3)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[-1 c(3)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3)]') iqr(x2*[-1 c(3)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[-1 c(3)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3)]') iqr(x2*[-1 c(3)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[-1 c(3)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[-1 c(3)]') iqr(x3*[-1 c(3)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_n1 = [-1 c_n1(3)];         
       
% Beta2 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) 1]') iqr(x1*[c(3) 1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1]') iqr(x2*[c(3) 1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1]') iqr(x2*[c(3) 1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) 1]') iqr(x3*[c(3) 1]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_p2 = [c_p2(3) 1];

% Beta2 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) -1]') iqr(x1*[c(3) -1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1]') iqr(x2*[c(3) -1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1]') iqr(x2*[c(3) -1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) -1]') iqr(x3*[c(3) -1]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_n2 = [c_n2(3) -1];  
    
dstar_vec = [dstar_p1 dstar_n1 dstar_p2 dstar_n2];

dstar = min(dstar_vec);
if dstar_p1 == min(dstar_vec)
    c = c_p1;
    beta = beta_p1;
elseif dstar_n1 == min(dstar_vec)
    c = c_n1;
    beta = beta_n1;
elseif dstar_p2 == min(dstar_vec)
    c = c_p2;
    beta = beta_p2;
elseif dstar_n2 == min(dstar_vec)
    c = c_n2;
    beta = beta_n2;
end

c1 = c(1);
c2 = c(2);







elseif nmarkers == 3
% Beta1 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[1 c(3) c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[1 c(3) c(4)]') iqr(x1*[1 c(3) c(4)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[1 c(3) c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4)]') iqr(x2*[1 c(3) c(4)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[1 c(3) c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4)]') iqr(x2*[1 c(3) c(4)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[1 c(3) c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[1 c(3) c(4)]') iqr(x3*[1 c(3) c(4)]')/1.34])*n3^(-0.2)))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4)];      
       
       
% Beta1 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[-1 c(3) c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[-1 c(3) c(4)]') iqr(x1*[-1 c(3) c(4)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[-1 c(3) c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4)]') iqr(x2*[-1 c(3) c(4)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[-1 c(3) c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4)]') iqr(x2*[-1 c(3) c(4)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[-1 c(3) c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[-1 c(3) c(4)]') iqr(x3*[-1 c(3) c(4)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4)];              
       
% Beta2 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) 1 c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) 1 c(4)]') iqr(x1*[c(3) 1 c(4)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) 1 c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4)]') iqr(x2*[c(3) 1 c(4)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) 1 c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4)]') iqr(x2*[c(3) 1 c(4)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) 1 c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) 1 c(4)]') iqr(x3*[c(3) 1 c(4)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4)];

% Beta2 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) -1 c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) -1 c(4)]') iqr(x1*[c(3) -1 c(4)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) -1 c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4)]') iqr(x2*[c(3) -1 c(4)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) -1 c(4)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4)]') iqr(x2*[c(3) -1 c(4)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) -1 c(4)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) -1 c(4)]') iqr(x3*[c(3) -1 c(4)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4)];  

% Beta3 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) 1]') iqr(x1*[c(3) c(4) 1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1]') iqr(x2*[c(3) c(4) 1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1]') iqr(x2*[c(3) c(4) 1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) 1]') iqr(x3*[c(3) c(4) 1]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1];

% Beta3 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) -1]') iqr(x1*[c(3) c(4) -1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1]') iqr(x2*[c(3) c(4) -1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1]') iqr(x2*[c(3) c(4) -1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) -1]') iqr(x3*[c(3) c(4) -1]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_n3 = [c_n3(3) -1 c_n3(4)]; 
    
dstar_vec = [dstar_p1 dstar_n1 dstar_p2 dstar_n2 dstar_p3 dstar_n3];

dstar = min(dstar_vec);
if dstar_p1 == dstar
    c = c_p1;
    beta = beta_p1;
elseif dstar_n1 == dstar
    c = c_n1;
    beta = beta_n1;
elseif dstar_p2 == dstar
    c = c_p2;
    beta = beta_p2;
elseif dstar_n2 == dstar
    c = c_n2;
    beta = beta_n2;
elseif dstar_p3 == dstar
    c = c_p3;
    beta = beta_p3;
elseif dstar_n3 == dstar
    c = c_n3;
    beta = beta_n3;
end

c1 = c(1);
c2 = c(2);










elseif nmarkers == 4
% Beta1 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[1 c(3) c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[1 c(3) c(4) c(5)]') iqr(x1*[1 c(3) c(4) c(5)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[1 c(3) c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5)]') iqr(x2*[1 c(3) c(4) c(5)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[1 c(3) c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5)]') iqr(x2*[1 c(3) c(4) c(5)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[1 c(3) c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[1 c(3) c(4) c(5)]') iqr(x3*[1 c(3) c(4) c(5)]')/1.34])*n3^(-0.2)))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4) c_p1(5)];      
       
       
% Beta1 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[-1 c(3) c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[-1 c(3) c(4) c(5)]') iqr(x1*[-1 c(3) c(4) c(5)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[-1 c(3) c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5)]') iqr(x2*[-1 c(3) c(4) c(5)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[-1 c(3) c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5)]') iqr(x2*[-1 c(3) c(4) c(5)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[-1 c(3) c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[-1 c(3) c(4) c(5)]') iqr(x3*[-1 c(3) c(4) c(5)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4) c_n1(5)];              
       
% Beta2 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) 1 c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) 1 c(4) c(5)]') iqr(x1*[c(3) 1 c(4) c(5)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) 1 c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5)]') iqr(x2*[c(3) 1 c(4) c(5)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) 1 c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5)]') iqr(x2*[c(3) 1 c(4) c(5)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) 1 c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) 1 c(4) c(5)]') iqr(x3*[c(3) 1 c(4) c(5)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4) c_p2(5)];

% Beta2 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) -1 c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) -1 c(4) c(5)]') iqr(x1*[c(3) -1 c(4) c(5)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) -1 c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5)]') iqr(x2*[c(3) -1 c(4) c(5)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) -1 c(4) c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5)]') iqr(x2*[c(3) -1 c(4) c(5)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) -1 c(4) c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) -1 c(4) c(5)]') iqr(x3*[c(3) -1 c(4) c(5)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4) c_n2(5)];  

% Beta3 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) 1 c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) 1 c(5)]') iqr(x1*[c(3) c(4) 1 c(5)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) 1 c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5)]') iqr(x2*[c(3) c(4) 1 c(5)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) 1 c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5)]') iqr(x2*[c(3) c(4) 1 c(5)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) 1 c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) 1 c(5)]') iqr(x3*[c(3) c(4) 1 c(5)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1 c_p3(5)];

% Beta3 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) -1 c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) -1 c(5)]') iqr(x1*[c(3) c(4) -1 c(5)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) -1 c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5)]') iqr(x2*[c(3) c(4) -1 c(5)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) -1 c(5)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5)]') iqr(x2*[c(3) c(4) -1 c(5)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) -1 c(5)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) -1 c(5)]') iqr(x3*[c(3) c(4) -1 c(5)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1 c_n3(5)]; 

% Beta4 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) 1]') iqr(x1*[c(3) c(4) c(5) 1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1]') iqr(x2*[c(3) c(4) c(5) 1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1]') iqr(x2*[c(3) c(4) c(5) 1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) 1]') iqr(x3*[c(3) c(4) c(5) 1]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p4,dstar_p4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p4 = [c_p4(3) c_p4(4) c_p4(5) 1];

% Beta4 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) -1]') iqr(x1*[c(3) c(4) c(5) -1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1]') iqr(x2*[c(3) c(4) c(5) -1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1]') iqr(x2*[c(3) c(4) c(5) -1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) -1]') iqr(x3*[c(3) c(4) c(5) -1]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n4,dstar_n4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n4 = [c_n4(3) c_n4(4) c_n4(5) -1]; 
    
dstar_vec = [dstar_p1 dstar_n1 dstar_p2 dstar_n2 dstar_p3 dstar_n3 dstar_p4 dstar_n4];

dstar = min(dstar_vec);
if dstar_p1 == dstar
    c = c_p1;
    beta = beta_p1;
elseif dstar_n1 == dstar
    c = c_n1;
    beta = beta_n1;
elseif dstar_p2 == dstar
    c = c_p2;
    beta = beta_p2;
elseif dstar_n2 == dstar
    c = c_n2;
    beta = beta_n2;
elseif dstar_p3 == dstar
    c = c_p3;
    beta = beta_p3;
elseif dstar_n3 == dstar
    c = c_n3;
    beta = beta_n3;
elseif dstar_p4 == dstar
    c = c_p4;
    beta = beta_p4;
elseif dstar_n4 == dstar
    c = c_n4;
    beta = beta_n4;
end

c1 = c(1);
c2 = c(2);



















elseif nmarkers == 5
% Beta1 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[1 c(3) c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[1 c(3) c(4) c(5) c(6)]') iqr(x1*[1 c(3) c(4) c(5) c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[1 c(3) c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5) c(6)]') iqr(x2*[1 c(3) c(4) c(5) c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[1 c(3) c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5) c(6)]') iqr(x2*[1 c(3) c(4) c(5) c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[1 c(3) c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[1 c(3) c(4) c(5) c(6)]') iqr(x3*[1 c(3) c(4) c(5) c(6)]')/1.34])*n3^(-0.2)))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4) c_p1(5) c_p1(6)];      
       
       
% Beta1 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[-1 c(3) c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[-1 c(3) c(4) c(5) c(6)]') iqr(x1*[-1 c(3) c(4) c(5) c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[-1 c(3) c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5) c(6)]') iqr(x2*[-1 c(3) c(4) c(5) c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[-1 c(3) c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5) c(6)]') iqr(x2*[-1 c(3) c(4) c(5) c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[-1 c(3) c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[-1 c(3) c(4) c(5) c(6)]') iqr(x3*[-1 c(3) c(4) c(5) c(6)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4) c_n1(5) c_n1(6)];              
       
% Beta2 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) 1 c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) 1 c(4) c(5) c(6)]') iqr(x1*[c(3) 1 c(4) c(5) c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) 1 c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5) c(6)]') iqr(x2*[c(3) 1 c(4) c(5) c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) 1 c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5) c(6)]') iqr(x2*[c(3) 1 c(4) c(5) c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) 1 c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) 1 c(4) c(5) c(6)]') iqr(x3*[c(3) 1 c(4) c(5) c(6)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4) c_p2(5) c_p2(6)];

% Beta2 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) -1 c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) -1 c(4) c(5) c(6)]') iqr(x1*[c(3) -1 c(4) c(5) c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) -1 c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5) c(6)]') iqr(x2*[c(3) -1 c(4) c(5) c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) -1 c(4) c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5) c(6)]') iqr(x2*[c(3) -1 c(4) c(5) c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) -1 c(4) c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) -1 c(4) c(5) c(6)]') iqr(x3*[c(3) -1 c(4) c(5) c(6)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4) c_n2(5) c_n2(6)];  

% Beta3 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) 1 c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) 1 c(5) c(6)]') iqr(x1*[c(3) c(4) 1 c(5) c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) 1 c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5) c(6)]') iqr(x2*[c(3) c(4) 1 c(5) c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) 1 c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5) c(6)]') iqr(x2*[c(3) c(4) 1 c(5) c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) 1 c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) 1 c(5) c(6)]') iqr(x3*[c(3) c(4) 1 c(5) c(6)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1 c_p3(5) c_p3(6)];

% Beta3 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) -1 c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) -1 c(5) c(6)]') iqr(x1*[c(3) c(4) -1 c(5) c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) -1 c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5) c(6)]') iqr(x2*[c(3) c(4) -1 c(5) c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) -1 c(5) c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5) c(6)]') iqr(x2*[c(3) c(4) -1 c(5) c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) -1 c(5) c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) -1 c(5) c(6)]') iqr(x3*[c(3) c(4) -1 c(5) c(6)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1 c_n3(5) c_n3(6)]; 

% Beta4 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) 1 c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) 1 c(6)]') iqr(x1*[c(3) c(4) c(5) 1 c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) 1 c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1 c(6)]') iqr(x2*[c(3) c(4) c(5) 1 c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) 1 c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1 c(6)]') iqr(x2*[c(3) c(4) c(5) 1 c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) 1 c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) 1 c(6)]') iqr(x3*[c(3) c(4) c(5) 1 c(6)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p4,dstar_p4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p4 = [c_p4(3) c_p4(4) c_p4(5) 1 c_p4(6)];

% Beta4 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) -1 c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) -1 c(6)]') iqr(x1*[c(3) c(4) c(5) -1 c(6)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) -1 c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1 c(6)]') iqr(x2*[c(3) c(4) c(5) -1 c(6)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) -1 c(6)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1 c(6)]') iqr(x2*[c(3) c(4) c(5) -1 c(6)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) -1 c(6)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) -1 c(6)]') iqr(x3*[c(3) c(4) c(5) -1 c(6)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n4,dstar_n4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n4 = [c_n4(3) c_n4(4) c_n4(5) -1 c_n4(6)]; 

% Beta5 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) 1]') iqr(x1*[c(3) c(4) c(5) c(6) 1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) 1]') iqr(x2*[c(3) c(4) c(5) c(6) 1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) 1]') iqr(x2*[c(3) c(4) c(5) c(6) 1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) 1]') iqr(x3*[c(3) c(4) c(5) c(6) 1]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p5,dstar_p5]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p5 = [c_p5(3) c_p5(4) c_p5(5) c_p5(6) 1];

% Beta5 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) -1]') iqr(x1*[c(3) c(4) c(5) c(6) -1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) -1]') iqr(x2*[c(3) c(4) c(5) c(6) -1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) -1]') iqr(x2*[c(3) c(4) c(5) c(6) -1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) -1]') iqr(x3*[c(3) c(4) c(5) c(6) -1]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n5,dstar_n5]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n5 = [c_n5(3) c_n5(4) c_n5(5) c_n5(6) -1]; 
    
dstar_vec = [dstar_p1 dstar_n1 dstar_p2 dstar_n2 dstar_p3 dstar_n3 dstar_p4 dstar_n4 dstar_p5 dstar_n5];

dstar = min(dstar_vec);
if dstar_p1 == dstar
    c = c_p1;
    beta = beta_p1;
elseif dstar_n1 == dstar
    c = c_n1;
    beta = beta_n1;
elseif dstar_p2 == dstar
    c = c_p2;
    beta = beta_p2;
elseif dstar_n2 == dstar
    c = c_n2;
    beta = beta_n2;
elseif dstar_p3 == dstar
    c = c_p3;
    beta = beta_p3;
elseif dstar_n3 == dstar
    c = c_n3;
    beta = beta_n3;
elseif dstar_p4 == dstar
    c = c_p4;
    beta = beta_p4;
elseif dstar_n4 == dstar
    c = c_n4;
    beta = beta_n4;
elseif dstar_p5 == dstar
    c = c_p5;
    beta = beta_p5;
elseif dstar_n5 == dstar
    c = c_n5;
    beta = beta_n5;
end

c1 = c(1);
c2 = c(2);















elseif nmarkers == 6
% Beta1 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[1 c(3) c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[1 c(3) c(4) c(5) c(6) c(7)]') iqr(x1*[1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[1 c(3) c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5) c(6) c(7)]') iqr(x2*[1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[1 c(3) c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5) c(6) c(7)]') iqr(x2*[1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[1 c(3) c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[1 c(3) c(4) c(5) c(6) c(7)]') iqr(x3*[1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n3^(-0.2)))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4) c_p1(5) c_p1(6) c_p1(7)];      
       
       
% Beta1 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[-1 c(3) c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[-1 c(3) c(4) c(5) c(6) c(7)]') iqr(x1*[-1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[-1 c(3) c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5) c(6) c(7)]') iqr(x2*[-1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[-1 c(3) c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5) c(6) c(7)]') iqr(x2*[-1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[-1 c(3) c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[-1 c(3) c(4) c(5) c(6) c(7)]') iqr(x3*[-1 c(3) c(4) c(5) c(6) c(7)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4) c_n1(5) c_n1(6) c_n1(7)];              
       
% Beta2 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) 1 c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) 1 c(4) c(5) c(6) c(7)]') iqr(x1*[c(3) 1 c(4) c(5) c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) 1 c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5) c(6) c(7)]') iqr(x2*[c(3) 1 c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) 1 c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5) c(6) c(7)]') iqr(x2*[c(3) 1 c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) 1 c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) 1 c(4) c(5) c(6) c(7)]') iqr(x3*[c(3) 1 c(4) c(5) c(6) c(7)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4) c_p2(5) c_p2(6) c_p2(7)];

% Beta2 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) -1 c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) -1 c(4) c(5) c(6) c(7)]') iqr(x1*[c(3) -1 c(4) c(5) c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) -1 c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5) c(6) c(7)]') iqr(x2*[c(3) -1 c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) -1 c(4) c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5) c(6) c(7)]') iqr(x2*[c(3) -1 c(4) c(5) c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) -1 c(4) c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) -1 c(4) c(5) c(6) c(7)]') iqr(x3*[c(3) -1 c(4) c(5) c(6) c(7)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4) c_n2(5) c_n2(6) c_n2(7)];  

% Beta3 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) 1 c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) 1 c(5) c(6) c(7)]') iqr(x1*[c(3) c(4) 1 c(5) c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) 1 c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5) c(6) c(7)]') iqr(x2*[c(3) c(4) 1 c(5) c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) 1 c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5) c(6) c(7)]') iqr(x2*[c(3) c(4) 1 c(5) c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) 1 c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) 1 c(5) c(6) c(7)]') iqr(x3*[c(3) c(4) 1 c(5) c(6) c(7)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1 c_p3(5) c_p3(6) c_p3(7)];

% Beta3 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) -1 c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) -1 c(5) c(6) c(7)]') iqr(x1*[c(3) c(4) -1 c(5) c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) -1 c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5) c(6) c(7)]') iqr(x2*[c(3) c(4) -1 c(5) c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) -1 c(5) c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5) c(6) c(7)]') iqr(x2*[c(3) c(4) -1 c(5) c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) -1 c(5) c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) -1 c(5) c(6) c(7)]') iqr(x3*[c(3) c(4) -1 c(5) c(6) c(7)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1 c_n3(5) c_n3(6) c_n3(7)]; 

% Beta4 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) 1 c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) 1 c(6) c(7)]') iqr(x1*[c(3) c(4) c(5) 1 c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) 1 c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1 c(6) c(7)]') iqr(x2*[c(3) c(4) c(5) 1 c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) 1 c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1 c(6) c(7)]') iqr(x2*[c(3) c(4) c(5) 1 c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) 1 c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) 1 c(6) c(7)]') iqr(x3*[c(3) c(4) c(5) 1 c(6) c(7)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p4,dstar_p4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_p4 = [c_p4(3) c_p4(4) c_p4(5) 1 c_p4(6) c_p4(7)];

% Beta4 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) -1 c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) -1 c(6) c(7)]') iqr(x1*[c(3) c(4) c(5) -1 c(6) c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) -1 c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1 c(6) c(7)]') iqr(x2*[c(3) c(4) c(5) -1 c(6) c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) -1 c(6) c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1 c(6) c(7)]') iqr(x2*[c(3) c(4) c(5) -1 c(6) c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) -1 c(6) c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) -1 c(6) c(7)]') iqr(x3*[c(3) c(4) c(5) -1 c(6) c(7)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n4,dstar_n4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_n4 = [c_n4(3) c_n4(4) c_n4(5) -1 c_n4(6) c_n4(7)]; 

% Beta5 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) 1 c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) 1 c(7)]') iqr(x1*[c(3) c(4) c(5) c(6) 1 c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) 1 c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) 1 c(7)]') iqr(x2*[c(3) c(4) c(5) c(6) 1 c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) 1 c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) 1 c(7)]') iqr(x2*[c(3) c(4) c(5) c(6) 1 c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) 1 c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) 1 c(7)]') iqr(x3*[c(3) c(4) c(5) c(6) 1 c(7)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p5,dstar_p5]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_p5 = [c_p5(3) c_p5(4) c_p5(5) c_p5(6) 1 c_p5(7)];

% Beta5 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) -1 c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) -1 c(7)]') iqr(x1*[c(3) c(4) c(5) c(6) -1 c(7)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) -1 c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) -1 c(7)]') iqr(x2*[c(3) c(4) c(5) c(6) -1 c(7)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) -1 c(7)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) -1 c(7)]') iqr(x2*[c(3) c(4) c(5) c(6) -1 c(7)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) -1 c(7)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) -1 c(7)]') iqr(x3*[c(3) c(4) c(5) c(6) -1 c(7)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n5,dstar_n5]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_n5 = [c_n5(3) c_n5(4) c_n5(5) c_n5(6) -1 c_n5(7)]; 

% Beta6 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) c(7) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) c(7) 1]') iqr(x1*[c(3) c(4) c(5) c(6) c(7) 1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) 1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) 1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) 1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) 1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) c(7) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) c(7) 1]') iqr(x3*[c(3) c(4) c(5) c(6) c(7) 1]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p6,dstar_p6]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_p6 = [c_p6(3) c_p6(4) c_p6(5) c_p6(6) c_p6(7) 1];

% Beta6 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) c(7) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) c(7) -1]') iqr(x1*[c(3) c(4) c(5) c(6) c(7) -1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) -1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) -1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) -1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) -1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) c(7) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) c(7) -1]') iqr(x3*[c(3) c(4) c(5) c(6) c(7) -1]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n6,dstar_n6]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1],options);
beta_n6 = [c_n6(3) c_n6(4) c_n6(5) c_n6(6) c_n6(7) -1]; 
    
dstar_vec = [dstar_p1 dstar_n1 dstar_p2 dstar_n2 dstar_p3 dstar_n3 dstar_p4 dstar_n4 dstar_p5 dstar_n5 dstar_p6 dstar_n6];

dstar = min(dstar_vec);
if dstar_p1 == dstar
    c = c_p1;
    beta = beta_p1;
elseif dstar_n1 == dstar
    c = c_n1;
    beta = beta_n1;
elseif dstar_p2 == dstar
    c = c_p2;
    beta = beta_p2;
elseif dstar_n2 == dstar
    c = c_n2;
    beta = beta_n2;
elseif dstar_p3 == dstar
    c = c_p3;
    beta = beta_p3;
elseif dstar_n3 == dstar
    c = c_n3;
    beta = beta_n3;
elseif dstar_p4 == dstar
    c = c_p4;
    beta = beta_p4;
elseif dstar_n4 == dstar
    c = c_n4;
    beta = beta_n4;
elseif dstar_p5 == dstar
    c = c_p5;
    beta = beta_p5;
elseif dstar_n5 == dstar
    c = c_n5;
    beta = beta_n5;
elseif dstar_p6 == dstar
    c = c_p6;
    beta = beta_p6;
elseif dstar_n6 == dstar
    c = c_n6;
    beta = beta_n6;
end

c1 = c(1);
c2 = c(2);









elseif nmarkers == 7
% Beta1 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x1*[1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x3*[1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4) c_p1(5) c_p1(6) c_p1(7) c_p1(8)];      
       
       
% Beta1 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x1*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]') iqr(x3*[-1 c(3) c(4) c(5) c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4) c_n1(5) c_n1(6) c_n1(7) c_n1(8)];              
       
% Beta2 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]') iqr(x1*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]') iqr(x3*[c(3) 1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4) c_p2(5) c_p2(6) c_p2(7) c_p2(8)];

% Beta2 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]') iqr(x1*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]') iqr(x3*[c(3) -1 c(4) c(5) c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4) c_n2(5) c_n2(6) c_n2(7) c_n2(8)];  

% Beta3 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]') iqr(x1*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]') iqr(x3*[c(3) c(4) 1 c(5) c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1 c_p3(5) c_p3(6) c_p3(7) c_p3(8)];

% Beta3 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]') iqr(x1*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]') iqr(x3*[c(3) c(4) -1 c(5) c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1 c_n3(5) c_n3(6) c_n3(7) c_n3(8)]; 

% Beta4 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]') iqr(x1*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]') iqr(x3*[c(3) c(4) c(5) 1 c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p4,dstar_p4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p4 = [c_p4(3) c_p4(4) c_p4(5) 1 c_p4(6) c_p4(7) c_p4(8)];

% Beta4 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]') iqr(x1*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]') iqr(x3*[c(3) c(4) c(5) -1 c(6) c(7) c(8)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n4,dstar_n4]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n4 = [c_n4(3) c_n4(4) c_n4(5) -1 c_n4(6) c_n4(7) c_n4(8)]; 

% Beta5 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]') iqr(x1*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]') iqr(x3*[c(3) c(4) c(5) c(6) 1 c(7) c(8)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p5,dstar_p5]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p5 = [c_p5(3) c_p5(4) c_p5(5) c_p5(6) 1 c_p5(7) c_p5(8)];

% Beta5 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]') iqr(x1*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]') iqr(x3*[c(3) c(4) c(5) c(6) -1 c(7) c(8)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n5,dstar_n5]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n5 = [c_n5(3) c_n5(4) c_n5(5) c_n5(6) -1 c_n5(7) c_n5(8)]; 

% Beta6 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]') iqr(x1*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]') iqr(x3*[c(3) c(4) c(5) c(6) c(7) 1 c(8)]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p6,dstar_p6]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p6 = [c_p6(3) c_p6(4) c_p6(5) c_p6(6) c_p6(7) 1 c_p6(8)];

% Beta6 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]') iqr(x1*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]') iqr(x3*[c(3) c(4) c(5) c(6) c(7) -1 c(8)]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n6,dstar_n6]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n6 = [c_n6(3) c_n6(4) c_n6(5) c_n6(6) c_n6(7) -1 c_n6(8)]; 

% Beta7 = 1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) c(7) c(8) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) c(7) c(8) 1]') iqr(x1*[c(3) c(4) c(5) c(6) c(7) c(8) 1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) c(8) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) c(8) 1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) c(8) 1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) c(8) 1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) c(8) 1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) c(8) 1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) c(7) c(8) 1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) c(7) c(8) 1]') iqr(x3*[c(3) c(4) c(5) c(6) c(7) c(8) 1]')/1.34])*n3^(-0.2)))^2);

options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p7,dstar_p7]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_p7 = [c_p7(3) c_p7(4) c_p7(5) c_p7(6) c_p7(7) c_p7(8) 1];

% Beta7 = -1
dist = @(c) sqrt((1 -  ksdensity(x1*[c(3) c(4) c(5) c(6) c(7) c(8) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x1*[c(3) c(4) c(5) c(6) c(7) c(8) -1]') iqr(x1*[c(3) c(4) c(5) c(6) c(7) c(8) -1]')/1.34])*n1^(-0.2)))^2 + ...
                        (1 - (ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) c(8) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) c(8) -1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) c(8) -1]')/1.34])*n2^(-0.2)) - ...
                              ksdensity(x2*[c(3) c(4) c(5) c(6) c(7) c(8) -1]',c(1),'Function','cdf','Bandwidth',0.9*min([std(x2*[c(3) c(4) c(5) c(6) c(7) c(8) -1]') iqr(x2*[c(3) c(4) c(5) c(6) c(7) c(8) -1]')/1.34])*n2^(-0.2))))^2 + ...
                             (ksdensity(x3*[c(3) c(4) c(5) c(6) c(7) c(8) -1]',c(2),'Function','cdf','Bandwidth',0.9*min([std(x3*[c(3) c(4) c(5) c(6) c(7) c(8) -1]') iqr(x3*[c(3) c(4) c(5) c(6) c(7) c(8) -1]')/1.34])*n3^(-0.2)))^2);


options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n7,dstar_n7]= fminsearchbnd(dist,[mean(x1,'all'),mean(x2,'all'),0,0,0,0,0,0],[-inf,-inf,-1,-1,-1,-1,-1,-1],[inf,inf,1,1,1,1,1,1],options);
beta_n7 = [c_n7(3) c_n7(4) c_n7(5) c_n7(6) c_n7(7) c_n7(8) -1]; 
    
dstar_vec = [dstar_p1 dstar_n1 dstar_p2 dstar_n2 dstar_p3 dstar_n3 dstar_p4 dstar_n4 dstar_p5 dstar_n5 dstar_p6 dstar_n6 dstar_p7 dstar_n7];

dstar = min(dstar_vec);
if dstar_p1 == dstar
    c = c_p1;
    beta = beta_p1;
elseif dstar_n1 == dstar
    c = c_n1;
    beta = beta_n1;
elseif dstar_p2 == dstar
    c = c_p2;
    beta = beta_p2;
elseif dstar_n2 == dstar
    c = c_n2;
    beta = beta_n2;
elseif dstar_p3 == dstar
    c = c_p3;
    beta = beta_p3;
elseif dstar_n3 == dstar
    c = c_n3;
    beta = beta_n3;
elseif dstar_p4 == dstar
    c = c_p4;
    beta = beta_p4;
elseif dstar_n4 == dstar
    c = c_n4;
    beta = beta_n4;
elseif dstar_p5 == dstar
    c = c_p5;
    beta = beta_p5;
elseif dstar_n5 == dstar
    c = c_n5;
    beta = beta_n5;
elseif dstar_p6 == dstar
    c = c_p6;
    beta = beta_p6;
elseif dstar_n6 == dstar
    c = c_n6;
    beta = beta_n6;
elseif dstar_p7 == dstar
    c = c_p7;
    beta = beta_p7;
elseif dstar_n7 == dstar
    c = c_n7;
    beta = beta_n7;
end

c1 = c(1);
c2 = c(2);

end
end