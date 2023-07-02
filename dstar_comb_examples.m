close all
clear all
clc

rng(123)

n1 = 30;
n2 = 30;
n3 = 30;

n1_test = 100000;
n2_test = 100000;
n3_test = 100000;

mu1_vec = [0 0 0];
mu2_vec = [.9 1 1.1];
mu3_vec = [1.8 2 2.2];

sig1 = [1 .3 .3
       .3 1 .3
       .3 .3 1];
   
sig2 = [1 .3 .3
       .3 1 .3
       .3 .3 1];
   
sig3 = [1 .3 .3
       .3 1 .3
       .3 .3 1];
   
markers1 = mvnrnd(mu1_vec,sig1,n1);
markers2 = mvnrnd(mu2_vec,sig2,n2);
markers3 = mvnrnd(mu3_vec,sig3,n3);

ind = [ones(n1,1);2*ones(n2,1);3*ones(n3,1)];
markers = [markers1;markers2;markers3];

[beta, cutoffs, dstar, tcrs] = dstar_comb(markers,ind,'normal','yes','yes','yes');

%%

close all
clear all
clc

rng(123)

n1 = 50;
n2 = 50;
n3 = 50;

mus = [1.3 1.3 1.3
       1.5 1.8 1.5
       1.8 2   1.9];
   
sigmas = [.5 .5 .5
          .4 .4 .4
          .3 .3 .25];

rho = .32;
cor = [1    rho rho
       rho  1   rho
       rho  rho 1];
   
[~,nmarkers] = size(cor);
   
U1 = copularnd('gaussian',cor,n1);
U2 = copularnd('gaussian',cor,n2);
U3 = copularnd('gaussian',cor,n3);

for i = 1:nmarkers
    markers1(:,i) = logninv(U1(:,i),mus(1,i),sigmas(1,i));
    markers2(:,i) = logninv(U2(:,i),mus(2,i),sigmas(2,i));
    markers3(:,i) = logninv(U3(:,i),mus(3,i),sigmas(3,i));
end

ind = [ones(n1,1);2*ones(n2,1);3*ones(n3,1)];
markers = [markers1;markers2;markers3];

[beta, cutoffs, dstar, tcrs] = dstar_comb(markers,ind,'boxcox','yes','yes','yes');

%%


