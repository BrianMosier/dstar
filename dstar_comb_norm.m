function [beta,c1,c2,dstar] = dstar_comb_norm(x1,x2,x3)

[d1,d2] = size(x1);

if d1 < d2
   x1 = transpose(x1);
   x2 = transpose(x2);
   x3 = transpose(x3);
end

n1 = length(x1(:,1));
n2 = length(x2(:,1));
n3 = length(x3(:,1));

mu1 = mean(x1);
mu2 = mean(x2);
mu3 = mean(x3);

cov1 = cov(x1)*(n1-1)/n1;
cov2 = cov(x2)*(n2-1)/n2;
cov3 = cov(x3)*(n3-1)/n3;

[~,nmarkers] = size(x1);

if nmarkers == 2
    
% Beta1 = 1    
dist = @(c) sqrt((1 -  normcdf(c(1),[1 c(3)]*mu1',   sqrt([1 c(3)]*cov1*[1 c(3)]' )))^2 + ...
                 (1 - (normcdf(c(2),[1 c(3)]*mu2',   sqrt([1 c(3)]*cov2*[1 c(3)]')) -  ...
                       normcdf(c(1),[1 c(3)]*mu2',   sqrt([1 c(3)]*cov2*[1 c(3)]'))))^2 + ...
                      (normcdf(c(2),[1 c(3)]*mu3',   sqrt([1 c(3)]*cov3*[1 c(3)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_p1 = [1 c_p1(3)];

% Beta1 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[-1 c(3)]*mu1',   sqrt([-1 c(3)]*cov1*[-1 c(3)]' )))^2 + ...
                 (1 - (normcdf(c(2),[-1 c(3)]*mu2',   sqrt([-1 c(3)]*cov2*[-1 c(3)]')) -  ...
                       normcdf(c(1),[-1 c(3)]*mu2',   sqrt([-1 c(3)]*cov2*[-1 c(3)]'))))^2 + ...
                      (normcdf(c(2),[-1 c(3)]*mu3',   sqrt([-1 c(3)]*cov3*[-1 c(3)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_n1 = [-1 c_n1(3)];

% Beta2 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) 1]*mu1',   sqrt([c(3) 1]*cov1*[c(3) 1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) 1]*mu2',   sqrt([c(3) 1]*cov2*[c(3) 1]')) -  ...
                       normcdf(c(1),[c(3) 1]*mu2',   sqrt([c(3) 1]*cov2*[c(3) 1]'))))^2 + ...
                      (normcdf(c(2),[c(3) 1]*mu3',   sqrt([c(3) 1]*cov3*[c(3) 1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
beta_p2 = [c_p2(3) 1];

% Beta2 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) -1]*mu1',   sqrt([c(3) -1]*cov1*[c(3) -1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) -1]*mu2',   sqrt([c(3) -1]*cov2*[c(3) -1]')) -  ...
                       normcdf(c(1),[c(3) -1]*mu2',   sqrt([c(3) -1]*cov2*[c(3) -1]'))))^2 + ...
                      (normcdf(c(2),[c(3) -1]*mu3',   sqrt([c(3) -1]*cov3*[c(3) -1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0],[-inf,-inf,-1],[inf,inf,1],options);
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
dist = @(c) sqrt((1 -  normcdf(c(1),[1 c(3) c(4)]*mu1',   sqrt([1 c(3) c(4)]*cov1*[1 c(3) c(4)]' )))^2 + ...
                 (1 - (normcdf(c(2),[1 c(3) c(4)]*mu2',   sqrt([1 c(3) c(4)]*cov2*[1 c(3) c(4)]')) -  ...
                       normcdf(c(1),[1 c(3) c(4)]*mu2',   sqrt([1 c(3) c(4)]*cov2*[1 c(3) c(4)]'))))^2 + ...
                      (normcdf(c(2),[1 c(3) c(4)]*mu3',   sqrt([1 c(3) c(4)]*cov3*[1 c(3) c(4)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4)];

% Beta1 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[-1 c(3) c(4)]*mu1',   sqrt([-1 c(3) c(4)]*cov1*[-1 c(3) c(4)]' )))^2 + ...
                 (1 - (normcdf(c(2),[-1 c(3) c(4)]*mu2',   sqrt([-1 c(3) c(4)]*cov2*[-1 c(3) c(4)]')) -  ...
                       normcdf(c(1),[-1 c(3) c(4)]*mu2',   sqrt([-1 c(3) c(4)]*cov2*[-1 c(3) c(4)]'))))^2 + ...
                      (normcdf(c(2),[-1 c(3) c(4)]*mu3',   sqrt([-1 c(3) c(4)]*cov3*[-1 c(3) c(4)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4)];

% Beta2 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) 1 c(4)]*mu1',   sqrt([c(3) 1 c(4)]*cov1*[c(3) 1 c(4)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) 1 c(4)]*mu2',   sqrt([c(3) 1 c(4)]*cov2*[c(3) 1 c(4)]')) -  ...
                       normcdf(c(1),[c(3) 1 c(4)]*mu2',   sqrt([c(3) 1 c(4)]*cov2*[c(3) 1 c(4)]'))))^2 + ...
                      (normcdf(c(2),[c(3) 1 c(4)]*mu3',   sqrt([c(3) 1 c(4)]*cov3*[c(3) 1 c(4)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4)];

% Beta2 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) -1 c(4)]*mu1',   sqrt([c(3) -1 c(4)]*cov1*[c(3) -1 c(4)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) -1 c(4)]*mu2',   sqrt([c(3) -1 c(4)]*cov2*[c(3) -1 c(4)]')) -  ...
                       normcdf(c(1),[c(3) -1 c(4)]*mu2',   sqrt([c(3) -1 c(4)]*cov2*[c(3) -1 c(4)]'))))^2 + ...
                      (normcdf(c(2),[c(3) -1 c(4)]*mu3',   sqrt([c(3) -1 c(4)]*cov3*[c(3) -1 c(4)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4)];

% Beta3 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) 1]*mu1',   sqrt([c(3) c(4) 1]*cov1*[c(3) c(4) 1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) 1]*mu2',   sqrt([c(3) c(4) 1]*cov2*[c(3) c(4) 1]')) -  ...
                       normcdf(c(1),[c(3) c(4) 1]*mu2',   sqrt([c(3) c(4) 1]*cov2*[c(3) c(4) 1]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) 1]*mu3',   sqrt([c(3) c(4) 1]*cov3*[c(3) c(4) 1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1];

% Beta3 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) -1]*mu1',   sqrt([c(3) c(4) -1]*cov1*[c(3) c(4) -1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) -1]*mu2',   sqrt([c(3) c(4) -1]*cov2*[c(3) c(4) -1]')) -  ...
                       normcdf(c(1),[c(3) c(4) -1]*mu2',   sqrt([c(3) c(4) -1]*cov2*[c(3) c(4) -1]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) -1]*mu3',   sqrt([c(3) c(4) -1]*cov3*[c(3) c(4) -1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0],[-inf,-inf,-1,-1],[inf,inf,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1];

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
dist = @(c) sqrt((1 -  normcdf(c(1),[1 c(3) c(4) c(5)]*mu1',   sqrt([1 c(3) c(4) c(5)]*cov1*[1 c(3) c(4) c(5)]' )))^2 + ...
                 (1 - (normcdf(c(2),[1 c(3) c(4) c(5)]*mu2',   sqrt([1 c(3) c(4) c(5)]*cov2*[1 c(3) c(4) c(5)]')) -  ...
                       normcdf(c(1),[1 c(3) c(4) c(5)]*mu2',   sqrt([1 c(3) c(4) c(5)]*cov2*[1 c(3) c(4) c(5)]'))))^2 + ...
                      (normcdf(c(2),[1 c(3) c(4) c(5)]*mu3',   sqrt([1 c(3) c(4) c(5)]*cov3*[1 c(3) c(4) c(5)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4) c_p1(5)];

% Beta1 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[-1 c(3) c(4) c(5)]*mu1',   sqrt([-1 c(3) c(4) c(5)]*cov1*[-1 c(3) c(4) c(5)]' )))^2 + ...
                 (1 - (normcdf(c(2),[-1 c(3) c(4) c(5)]*mu2',   sqrt([-1 c(3) c(4) c(5)]*cov2*[-1 c(3) c(4) c(5)]')) -  ...
                       normcdf(c(1),[-1 c(3) c(4) c(5)]*mu2',   sqrt([-1 c(3) c(4) c(5)]*cov2*[-1 c(3) c(4) c(5)]'))))^2 + ...
                      (normcdf(c(2),[-1 c(3) c(4) c(5)]*mu3',   sqrt([-1 c(3) c(4) c(5)]*cov3*[-1 c(3) c(4) c(5)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4) c_n1(5)];

% Beta2 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) 1 c(4) c(5)]*mu1',   sqrt([c(3) 1 c(4) c(5)]*cov1*[c(3) 1 c(4) c(5)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) 1 c(4) c(5)]*mu2',   sqrt([c(3) 1 c(4) c(5)]*cov2*[c(3) 1 c(4) c(5)]')) -  ...
                       normcdf(c(1),[c(3) 1 c(4) c(5)]*mu2',   sqrt([c(3) 1 c(4) c(5)]*cov2*[c(3) 1 c(4) c(5)]'))))^2 + ...
                      (normcdf(c(2),[c(3) 1 c(4) c(5)]*mu3',   sqrt([c(3) 1 c(4) c(5)]*cov3*[c(3) 1 c(4) c(5)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4) c_p2(5)];

% Beta2 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) -1 c(4) c(5)]*mu1',   sqrt([c(3) -1 c(4) c(5)]*cov1*[c(3) -1 c(4) c(5)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) -1 c(4) c(5)]*mu2',   sqrt([c(3) -1 c(4) c(5)]*cov2*[c(3) -1 c(4) c(5)]')) -  ...
                       normcdf(c(1),[c(3) -1 c(4) c(5)]*mu2',   sqrt([c(3) -1 c(4) c(5)]*cov2*[c(3) -1 c(4) c(5)]'))))^2 + ...
                      (normcdf(c(2),[c(3) -1 c(4) c(5)]*mu3',   sqrt([c(3) -1 c(4) c(5)]*cov3*[c(3) -1 c(4) c(5)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4) c_n2(5)];

% Beta3 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) 1 c(5)]*mu1',   sqrt([c(3) c(4) 1 c(5)]*cov1*[c(3) c(4) 1 c(5)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) 1 c(5)]*mu2',   sqrt([c(3) c(4) 1 c(5)]*cov2*[c(3) c(4) 1 c(5)]')) -  ...
                       normcdf(c(1),[c(3) c(4) 1 c(5)]*mu2',   sqrt([c(3) c(4) 1 c(5)]*cov2*[c(3) c(4) 1 c(5)]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) 1 c(5)]*mu3',   sqrt([c(3) c(4) 1 c(5)]*cov3*[c(3) c(4) 1 c(5)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1 c_p3(5)];

% Beta3 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) -1 c(5)]*mu1',   sqrt([c(3) c(4) -1 c(5)]*cov1*[c(3) c(4) -1 c(5)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) -1 c(5)]*mu2',   sqrt([c(3) c(4) -1 c(5)]*cov2*[c(3) c(4) -1 c(5)]')) -  ...
                       normcdf(c(1),[c(3) c(4) -1 c(5)]*mu2',   sqrt([c(3) c(4) -1 c(5)]*cov2*[c(3) c(4) -1 c(5)]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) -1 c(5)]*mu3',   sqrt([c(3) c(4) -1 c(5)]*cov3*[c(3) c(4) -1 c(5)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1 c_n3(5)];

% Beta4 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) c(5) 1]*mu1',   sqrt([c(3) c(4) c(5) 1]*cov1*[c(3) c(4) c(5) 1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) c(5) 1]*mu2',   sqrt([c(3) c(4) c(5) 1]*cov2*[c(3) c(4) c(5) 1]')) -  ...
                       normcdf(c(1),[c(3) c(4) c(5) 1]*mu2',   sqrt([c(3) c(4) c(5) 1]*cov2*[c(3) c(4) c(5) 1]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) c(5) 1]*mu3',   sqrt([c(3) c(4) c(5) 1]*cov3*[c(3) c(4) c(5) 1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p4,dstar_p4]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
beta_p4 = [c_p4(3) c_p4(4) c_p4(5) 1];

% Beta4 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) c(5) -1]*mu1',   sqrt([c(3) c(4) c(5) -1]*cov1*[c(3) c(4) c(5) -1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) c(5) -1]*mu2',   sqrt([c(3) c(4) c(5) -1]*cov2*[c(3) c(4) c(5) -1]')) -  ...
                       normcdf(c(1),[c(3) c(4) c(5) -1]*mu2',   sqrt([c(3) c(4) c(5) -1]*cov2*[c(3) c(4) c(5) -1]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) c(5) -1]*mu3',   sqrt([c(3) c(4) c(5) -1]*cov3*[c(3) c(4) c(5) -1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n4,dstar_n4]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0],[-inf,-inf,-1,-1,-1],[inf,inf,1,1,1],options);
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
dist = @(c) sqrt((1 -  normcdf(c(1),[1 c(3) c(4) c(5) c(6)]*mu1',   sqrt([1 c(3) c(4) c(5) c(6)]*cov1*[1 c(3) c(4) c(5) c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[1 c(3) c(4) c(5) c(6)]*mu2',   sqrt([1 c(3) c(4) c(5) c(6)]*cov2*[1 c(3) c(4) c(5) c(6)]')) -  ...
                       normcdf(c(1),[1 c(3) c(4) c(5) c(6)]*mu2',   sqrt([1 c(3) c(4) c(5) c(6)]*cov2*[1 c(3) c(4) c(5) c(6)]'))))^2 + ...
                      (normcdf(c(2),[1 c(3) c(4) c(5) c(6)]*mu3',   sqrt([1 c(3) c(4) c(5) c(6)]*cov3*[1 c(3) c(4) c(5) c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p1,dstar_p1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p1 = [1 c_p1(3) c_p1(4) c_p1(5) c_p1(6)];

% Beta1 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[-1 c(3) c(4) c(5) c(6)]*mu1',   sqrt([-1 c(3) c(4) c(5) c(6)]*cov1*[-1 c(3) c(4) c(5) c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[-1 c(3) c(4) c(5) c(6)]*mu2',   sqrt([-1 c(3) c(4) c(5) c(6)]*cov2*[-1 c(3) c(4) c(5) c(6)]')) -  ...
                       normcdf(c(1),[-1 c(3) c(4) c(5) c(6)]*mu2',   sqrt([-1 c(3) c(4) c(5) c(6)]*cov2*[-1 c(3) c(4) c(5) c(6)]'))))^2 + ...
                      (normcdf(c(2),[-1 c(3) c(4) c(5) c(6)]*mu3',   sqrt([-1 c(3) c(4) c(5) c(6)]*cov3*[-1 c(3) c(4) c(5) c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n1,dstar_n1]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n1 = [-1 c_n1(3) c_n1(4) c_n1(5) c_n1(6)];

% Beta2 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) 1 c(4) c(5) c(6)]*mu1',   sqrt([c(3) 1 c(4) c(5) c(6)]*cov1*[c(3) 1 c(4) c(5) c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) 1 c(4) c(5) c(6)]*mu2',   sqrt([c(3) 1 c(4) c(5) c(6)]*cov2*[c(3) 1 c(4) c(5) c(6)]')) -  ...
                       normcdf(c(1),[c(3) 1 c(4) c(5) c(6)]*mu2',   sqrt([c(3) 1 c(4) c(5) c(6)]*cov2*[c(3) 1 c(4) c(5) c(6)]'))))^2 + ...
                      (normcdf(c(2),[c(3) 1 c(4) c(5) c(6)]*mu3',   sqrt([c(3) 1 c(4) c(5) c(6)]*cov3*[c(3) 1 c(4) c(5) c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p2,dstar_p2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p2 = [c_p2(3) 1 c_p2(4) c_p2(5) c_p2(6)];

% Beta2 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) -1 c(4) c(5) c(6)]*mu1',   sqrt([c(3) -1 c(4) c(5) c(6)]*cov1*[c(3) -1 c(4) c(5) c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) -1 c(4) c(5) c(6)]*mu2',   sqrt([c(3) -1 c(4) c(5) c(6)]*cov2*[c(3) -1 c(4) c(5) c(6)]')) -  ...
                       normcdf(c(1),[c(3) -1 c(4) c(5) c(6)]*mu2',   sqrt([c(3) -1 c(4) c(5) c(6)]*cov2*[c(3) -1 c(4) c(5) c(6)]'))))^2 + ...
                      (normcdf(c(2),[c(3) -1 c(4) c(5) c(6)]*mu3',   sqrt([c(3) -1 c(4) c(5) c(6)]*cov3*[c(3) -1 c(4) c(5) c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n2,dstar_n2]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n2 = [c_n2(3) -1 c_n2(4) c_n2(5) c_n2(6)];

% Beta3 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) 1 c(5) c(6)]*mu1',   sqrt([c(3) c(4) 1 c(5) c(6)]*cov1*[c(3) c(4) 1 c(5) c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) 1 c(5) c(6)]*mu2',   sqrt([c(3) c(4) 1 c(5) c(6)]*cov2*[c(3) c(4) 1 c(5) c(6)]')) -  ...
                       normcdf(c(1),[c(3) c(4) 1 c(5) c(6)]*mu2',   sqrt([c(3) c(4) 1 c(5) c(6)]*cov2*[c(3) c(4) 1 c(5) c(6)]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) 1 c(5) c(6)]*mu3',   sqrt([c(3) c(4) 1 c(5) c(6)]*cov3*[c(3) c(4) 1 c(5) c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p3,dstar_p3]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p3 = [c_p3(3) c_p3(4) 1 c_p3(5) c_p3(6)];

% Beta3 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) -1 c(5) c(6)]*mu1',   sqrt([c(3) c(4) -1 c(5) c(6)]*cov1*[c(3) c(4) -1 c(5) c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) -1 c(5) c(6)]*mu2',   sqrt([c(3) c(4) -1 c(5) c(6)]*cov2*[c(3) c(4) -1 c(5) c(6)]')) -  ...
                       normcdf(c(1),[c(3) c(4) -1 c(5) c(6)]*mu2',   sqrt([c(3) c(4) -1 c(5) c(6)]*cov2*[c(3) c(4) -1 c(5) c(6)]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) -1 c(5) c(6)]*mu3',   sqrt([c(3) c(4) -1 c(5) c(6)]*cov3*[c(3) c(4) -1 c(5) c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n3,dstar_n3]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n3 = [c_n3(3) c_n3(4) -1 c_n3(5) c_n3(6)];

% Beta4 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) c(5) 1 c(6)]*mu1',   sqrt([c(3) c(4) c(5) 1 c(6)]*cov1*[c(3) c(4) c(5) 1 c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) c(5) 1 c(6)]*mu2',   sqrt([c(3) c(4) c(5) 1 c(6)]*cov2*[c(3) c(4) c(5) 1 c(6)]')) -  ...
                       normcdf(c(1),[c(3) c(4) c(5) 1 c(6)]*mu2',   sqrt([c(3) c(4) c(5) 1 c(6)]*cov2*[c(3) c(4) c(5) 1 c(6)]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) c(5) 1 c(6)]*mu3',   sqrt([c(3) c(4) c(5) 1 c(6)]*cov3*[c(3) c(4) c(5) 1 c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p4,dstar_p4]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p4 = [c_p4(3) c_p4(4) c_p4(5) 1 c_p4(6)];

% Beta4 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) c(5) -1 c(6)]*mu1',   sqrt([c(3) c(4) c(5) -1 c(6)]*cov1*[c(3) c(4) c(5) -1 c(6)]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) c(5) -1 c(6)]*mu2',   sqrt([c(3) c(4) c(5) -1 c(6)]*cov2*[c(3) c(4) c(5) -1 c(6)]')) -  ...
                       normcdf(c(1),[c(3) c(4) c(5) -1 c(6)]*mu2',   sqrt([c(3) c(4) c(5) -1 c(6)]*cov2*[c(3) c(4) c(5) -1 c(6)]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) c(5) -1 c(6)]*mu3',   sqrt([c(3) c(4) c(5) -1 c(6)]*cov3*[c(3) c(4) c(5) -1 c(6)]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n4,dstar_n4]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_n4 = [c_n4(3) c_n4(4) c_n4(5) -1 c_n4(6)];

% Beta5 = 1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) c(5) c(6) 1]*mu1',   sqrt([c(3) c(4) c(5) c(6) 1]*cov1*[c(3) c(4) c(5) c(6) 1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) c(5) c(6) 1]*mu2',   sqrt([c(3) c(4) c(5) c(6) 1]*cov2*[c(3) c(4) c(5) c(6) 1]')) -  ...
                       normcdf(c(1),[c(3) c(4) c(5) c(6) 1]*mu2',   sqrt([c(3) c(4) c(5) c(6) 1]*cov2*[c(3) c(4) c(5) c(6) 1]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) c(5) c(6) 1]*mu3',   sqrt([c(3) c(4) c(5) c(6) 1]*cov3*[c(3) c(4) c(5) c(6) 1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_p5,dstar_p5]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
beta_p5 = [c_p5(3) c_p5(4) c_p5(5) c_p5(6) 1];

% Beta5 = -1
dist = @(c) sqrt((1 -  normcdf(c(1),[c(3) c(4) c(5) c(6) -1]*mu1',   sqrt([c(3) c(4) c(5) c(6) -1]*cov1*[c(3) c(4) c(5) c(6) -1]' )))^2 + ...
                 (1 - (normcdf(c(2),[c(3) c(4) c(5) c(6) -1]*mu2',   sqrt([c(3) c(4) c(5) c(6) -1]*cov2*[c(3) c(4) c(5) c(6) -1]')) -  ...
                       normcdf(c(1),[c(3) c(4) c(5) c(6) -1]*mu2',   sqrt([c(3) c(4) c(5) c(6) -1]*cov2*[c(3) c(4) c(5) c(6) -1]'))))^2 + ...
                      (normcdf(c(2),[c(3) c(4) c(5) c(6) -1]*mu3',   sqrt([c(3) c(4) c(5) c(6) -1]*cov3*[c(3) c(4) c(5) c(6) -1]')))^2);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-16);
[c_n5,dstar_n5]= fminsearchbnd(dist,[mean(x1,'all') mean(x3,'all'),0,0,0,0],[-inf,-inf,-1,-1,-1,-1],[inf,inf,1,1,1,1],options);
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
end
end
