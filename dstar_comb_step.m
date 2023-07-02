function [beta,c1,c2,dstar,tcr1,tcr2,tcr3] = dstar_comb_step(x1,x2,x3)

[d1,d2] = size(x1);

if d1 < d2
   x1 = transpose(x1);
   x2 = transpose(x2);
   x3 = transpose(x3);
end

[n1,nmarkers] = size(x1);
[n2,~] = size(x2);
[n3,~] = size(x3);

for i = 1:nmarkers
   dists(i) = Empirical_Dist_3D(x1(:,i),x2(:,i),x3(:,i));
end

marker = 1:nmarkers;
ind = [marker
        dists];
[~,order] = sort(ind(2,:),'ascend');
ord = ind(1,order);
x1_ord = x1(:,ord);
x2_ord = x2(:,ord);
x3_ord = x3(:,ord);

beta_vec = -1:.01:1;

pseudo_1 = x1_ord;
pseudo_2 = x2_ord;
pseudo_3 = x3_ord;

beta = ones(1,nmarkers);

for i = 1:(nmarkers-1)
    for j = 1:201
        pseudo_1_a = pseudo_1(:,i) + beta_vec(j)*pseudo_1(:,i+1);
        pseudo_2_a = pseudo_2(:,i) + beta_vec(j)*pseudo_2(:,i+1);
        pseudo_3_a = pseudo_3(:,i) + beta_vec(j)*pseudo_3(:,i+1);
        
        pseudo_1_b = beta_vec(j)*pseudo_1(:,i) + pseudo_1(:,i+1);
        pseudo_2_b = beta_vec(j)*pseudo_2(:,i) + pseudo_2(:,i+1);
        pseudo_3_b = beta_vec(j)*pseudo_3(:,i) + pseudo_3(:,i+1);
        
        [dist_a(j), c1_a(j), c2_a(j)] = Empirical_Dist_3D(pseudo_1_a,pseudo_2_a,pseudo_3_a);
        [dist_b(j), c1_b(j), c2_b(j)] = Empirical_Dist_3D(pseudo_1_b,pseudo_2_b,pseudo_3_b);

    end
    Dstar_a = min(dist_a);
    Dstar_b = min(dist_b);
    if Dstar_a <= Dstar_b
        Dstar = Dstar_a;
        tmp = beta_vec(Dstar_a==dist_a);
        beta_tmp = [1 tmp(1)];
        c_tmp1 = c1_a(Dstar_a==dist_a);
        c_tmp2 = c2_a(Dstar_a==dist_a);
        c1_tmp = c_tmp1(1);
        c2_tmp = c_tmp2(1);
        tmp = [];
    elseif Dstar_b < Dstar_a
        Dstar = Dstar_b;
        tmp = beta_vec(Dstar_b==dist_b);
        beta_tmp = [tmp(1) 1];
        c_tmp1 = c1_b(Dstar_b==dist_b);
        c_tmp2 = c2_b(Dstar_b==dist_b);
        c1_tmp = c_tmp1(1);
        c2_tmp = c_tmp2(1);
        tmp = [];
    end
    
    if i == 1
       beta(1) = beta_tmp(1);
       beta(2) = beta_tmp(2);
       beta;
    elseif i > 1
        for k = 1:i
            beta(k) = beta(k)*beta_tmp(1);
        end
        for l = i+1
            beta(l) = beta_tmp(2); 
        end
            beta;
    end
    pseudo_1(:,i+1) = x1_ord(:,1:i+1)*beta(1:i+1)';
    pseudo_2(:,i+1) = x2_ord(:,1:i+1)*beta(1:i+1)';
    pseudo_3(:,i+1) = x3_ord(:,1:i+1)*beta(1:i+1)';
    c1 = c1_tmp;
    c2 = c2_tmp;
    dstar = Dstar;
end

if max(abs(beta)<1)
   beta = beta/max(abs(beta));
   c1 = c1/max(abs(beta));
   c2 = c2/max(abs(beta));
end

beta_hld = [beta;ord]';
beta_hld = sortrows(beta_hld,2);
beta = beta_hld(:,1)';

y1 = x1*beta';
y2 = x2*beta';
y3 = x3*beta';

[~,c1,c2] = Empirical_Dist_3D(y1,y2,y3);

tcr1 = sum(y1<=c1)/n1;
tcr2 = sum((y2<=c2).*(y2>c1))/n2;
tcr3 = sum(y3>c2)/n3;

end


