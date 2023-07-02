function varargout = dstar_comb(y,D,method,normtest,qqplot,rocplots)

if nargin<3; error('The function ''dstar_comb'' requires at least 3 arguments');end
if ~exist('rocplots','var')
    rocplots = 'no';
end
if nargin>6; error('Too many input arguments in function ''dstar)comb'' ');end
if isnumeric(y)~=1;error('Biomarker scores must be numeric');end

listm={'normal','boxcox','kernel','stepwise','logistic'};
wm=sum(strncmp(method,listm,8));

if wm==0;error('Invalid value of ''function'' parameter. You must provide one of the following: ''normal'' or ''boxcox'' or ''kernel'' or ''stepwise'' or ''logistic'' ');end

list={'yes','no'};
w=sum(strncmp(rocplots,list,3));
if w==0;error('Invalid value of ''function'' parameter. You must provide one of the following: ''yes'' or ''no'' ');end

if ~exist('normtest','var')
    normtest = 'no';
end

[d1,d2] = size(y);

if d1 < d2
   y = transpose(y);
end

y1 = y(D==1,:);
y2 = y(D==2,:);
y3 = y(D==3,:);

[n1,nmarkers] = size(y1);
[n2,~] = size(y2);
[n3,~] = size(y3);

switch method
    case 'normal'
    switch normtest
        case 'yes'
            switch qqplot
                case 'yes'
                    roc3d_norm_test_comb(y,D,.05,'yes');
                case 'no'
                    roc3d_norm_test_comb(y,D,.05,'no');
            end
    end
    [beta,c1,c2,dstar] = dstar_comb_norm(y1,y2,y3);
    
    x1 = y1*beta';
    x2 = y2*beta';
    x3 = y3*beta';
    
    x = [x1;x2;x3];
    
    xbar1 = mean(x1);
    xbar2 = mean(x2);
    xbar3 = mean(x3);
    
    sd1 = std(x1)*(n1-1)/n1;
    sd2 = std(x2)*(n2-1)/n2;
    sd3 = std(x3)*(n3-1)/n3;
    
    TCR1ds = normcdf(c1,xbar1,sd1);
    TCR2ds = normcdf(c2,xbar2,sd2) - normcdf(c1,xbar2,sd2);
    TCR3ds = 1 - normcdf(c2,xbar3,sd3);
    
    for i = 1:nmarkers
        beta_vec1(i) = append('β',string(i));
        beta_vec2(i) = append(string(beta(i)));
    end
    
    
    
    disp('The combination coefficients are:')
    disp(append(string(beta_vec1),' = ',string(beta_vec2)))
    disp(' ')
    disp('The optimal cutoff values are:')
    disp(append('c1 = ',string(c1)))
    disp(append('c2 = ',string(c2)))
    disp(' ')
    disp('The Euclidean distance is:')
    disp(append('D* = ',string(dstar)))
    disp(' ')
    disp('The true classification rates are:')
    disp(append('[TCR1, TCR2, TCR3] = [',string(TCR1ds),', ',string(TCR2ds),', ',string(TCR3ds),']')) 
    
    varargout{1} = beta;
    varargout{2} = [c1 c2];
    varargout{3} = dstar;
    varargout{4} = [TCR1ds TCR2ds TCR3ds];
        
    switch rocplots
            case 'yes'
                [p1,p3] = meshgrid(0:0.01:1, 0:0.01:1);
                ROCsparam=(normcdf(norminv(1-p3,xbar3,sd3),xbar2,sd2)-normcdf(norminv(p1,xbar1,sd1),xbar2,sd2));
    
                figure('Position',[0 0 700 600])
                surf(p1,p3,ROCsparam) %plot the ROC surface
                title('ROC Surface for the Combined Scores')
                hold on
                view(220,10)
                xlim([0 1])
                ylim([0 1])
                zlim([0 1])
                line([1 TCR1ds], [1 TCR3ds], [1 TCR2ds],'Color','r','LineWidth',1.5)
                shading interp
                alpha 0.6
                xlabel('TCR1')
                ylabel('TCR3')
                zlabel('TCR2')
                set(gcf,'color','w');
                hold off
                
                figure('Position',[0 0 1000 600])
                h1 = histogram(x1,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 0 .7],'FaceColor',[0 0 .7],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                hold on
                h2 = histogram(x2,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 .6 0],'FaceColor',[0 .6 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                h3 = histogram(x3,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[.7 0 0],'FaceColor',[.7 0 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);

                gr = xbar1-3.5*sd1:.01:xbar3+3.5*sd3;
                plot(gr,normpdf(gr,xbar1,sd1),'Color',[0 0 .7])
                title('Distributions of the Combined Scores')
                hold on
                plot(gr,normpdf(gr,xbar2,sd2),'Color',[0 .7 0])
                plot(gr,normpdf(gr,xbar3,sd3),'Color',[.7 0 0])
                max1 = max(h1.Values);
                max2 = max(h2.Values);
                max3 = max(h3.Values);
                maxf = max([max1 max2 max3]);
                ylim([0 1.25*maxf])
                xlim([floor(xbar1-3.5*sd1) ceil(xbar3+3.5*sd3)])
                plot([c1 c1],[0 1.25*maxf],'--k')
                text(c1+.05,maxf*1.1,append('c_{1} = ',string(round(c1,4))))
                plot([c2 c2],[0 1.25*maxf],'--k')
                text(c2+.05,maxf*1.1,append('c_{2} = ',string(round(c2,4))))
                set(gcf,'color','w');
    end
    
    case 'boxcox'
    switch normtest
        case 'yes'
            switch qqplot
                case 'yes'
                    disp('For the original scores:')
                    roc3d_norm_test_comb(y,D,.05,'yes');
                    sgtitle('Q-Q Plots for the Original Scores')
                case 'no'
                    disp('For the original scores:')
                    roc3d_norm_test_comb(y,D,.05,'no');
            end
    end
    [beta,lambda,c1,c2,dstar] = dstar_comb_bc(y1,y2,y3);
    
    min_y = min(y,[],'all');

    if min_y < 0
        shift = abs(min_y)+1;
        y1 = y1 + shift;
        y2 = y2 + shift;
        y3 = y3 + shift;
    end
    
    for i = 1:nmarkers
    y1t(:,i) = (y1(:,i).^lambda(i)-1)/lambda(i);
    y2t(:,i) = (y2(:,i).^lambda(i)-1)/lambda(i);
    y3t(:,i) = (y3(:,i).^lambda(i)-1)/lambda(i);
    end
    yt = [y1t;y2t;y3t];
    switch normtest
        case 'yes'
            switch qqplot
                case 'yes'
                    disp('For the transformed scores:')
                    roc3d_norm_test_comb(yt,D,.05,'yes');
                    sgtitle('Q-Q Plots for the Transformed Scores')
                case 'no'
                    disp('For the transformed scores:')
                    roc3d_norm_test_comb(yt,D,.05,'no');
            end
    end
    
    x1 = y1t*beta';
    x2 = y2t*beta';
    x3 = y3t*beta';
    x = [x1;x2;x3];
    
    xbar1 = mean(x1);
    xbar2 = mean(x2);
    xbar3 = mean(x3);
    
    sd1 = std(x1)*(n1-1)/n1;
    sd2 = std(x2)*(n2-1)/n2;
    sd3 = std(x3)*(n3-1)/n3;
    
    TCR1ds = normcdf(c1,xbar1,sd1);
    TCR2ds = normcdf(c2,xbar2,sd2) - normcdf(c1,xbar2,sd2);
    TCR3ds = 1 - normcdf(c2,xbar3,sd3);
    
    for i = 1:nmarkers
        beta_vec1(i) = append('β',string(i));
        beta_vec2(i) = append(string(beta(i)));
    end
    
    disp('The combination coefficients are:')
    disp(append(string(beta_vec1),' = ',string(beta_vec2)))
    disp(' ')
    disp('The optimal cutoff values are:')
    disp(append('c1 = ',string(c1)))
    disp(append('c2 = ',string(c2)))
    disp(' ')
    disp('The Euclidean distance is:')
    disp(append('D* = ',string(dstar)))
    disp(' ')
    disp('The true classification rates are:')
    disp(append('[TCR1, TCR2, TCR3] = [',string(TCR1ds),', ',string(TCR2ds),', ',string(TCR3ds),']')) 
    disp(' ')
    disp('The transformation parameter, λ, is:')
    disp(append('λ = ',string(lambda)));
    
         
    varargout{1} = beta;
    varargout{2} = [c1 c2];
    varargout{3} = dstar;
    varargout{4} = [TCR1ds TCR2ds TCR3ds];
    varargout{5} = lambda;

    switch rocplots
            case 'yes'
                [p1,p3] = meshgrid(0:0.01:1, 0:0.01:1);
                ROCsparam=(normcdf(norminv(1-p3,xbar3,sd3),xbar2,sd2)-normcdf(norminv(p1,xbar1,sd1),xbar2,sd2));
    
                figure('Position',[0 0 700 600])
                surf(p1,p3,ROCsparam) %plot the ROC surface
                hold on
                view(220,10)
                xlim([0 1])
                ylim([0 1])
                zlim([0 1])
                line([1 TCR1ds], [1 TCR3ds], [1 TCR2ds],'Color','r','LineWidth',1.5)
                shading interp
                alpha 0.6
                xlabel('TCR1')
                ylabel('TCR3')
                zlabel('TCR2')
                set(gcf,'color','w');
                hold off
                
                figure('Position',[0 0 1000 600])
                h1 = histogram(x1,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 0 .7],'FaceColor',[0 0 .7],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                hold on
                h2 = histogram(x2,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 .6 0],'FaceColor',[0 .6 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                h3 = histogram(x3,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[.7 0 0],'FaceColor',[.7 0 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);

                gr = xbar1-3.5*sd1:.01:xbar3+3.5*sd3;
                plot(gr,normpdf(gr,xbar1,sd1),'Color',[0 0 .7])
                title('Distributions of the Combined Scores')
                hold on
                plot(gr,normpdf(gr,xbar2,sd2),'Color',[0 .7 0])
                plot(gr,normpdf(gr,xbar3,sd3),'Color',[.7 0 0])
                max1 = max(h1.Values);
                max2 = max(h2.Values);
                max3 = max(h3.Values);
                maxf = max([max1 max2 max3]);
                ylim([0 1.25*maxf])
                xlim([floor(xbar1-3.5*sd1) ceil(xbar3+3.5*sd3)])
                plot([c1 c1],[0 1.25*maxf],'--k')
                text(c1+.05,maxf*1.1,append('c_{1} = ',string(round(c1,4))))
                plot([c2 c2],[0 1.25*maxf],'--k')
                text(c2+.05,maxf*1.1,append('c_{2} = ',string(round(c2,4))))
                set(gcf,'color','w');
    end
    
    case 'kernel'
    
    if sum(strncmp(normtest,'yes',3)) == 1
    roc3d_norm_test_comb(y,D,.05,'no');
    end
        
    [beta,c1,c2,dstar] = dstar_comb_ker(y1,y2,y3);
    
    x1 = y1*beta';
    x2 = y2*beta';
    x3 = y3*beta';
    x = [x1;x2;x3];
    bw1 = 0.9*min([std(x1) iqr(x1)/1.34])*n1^(-0.2);
    bw2 = 0.9*min([std(x2) iqr(x2)/1.34])*n2^(-0.2);
    bw3 = 0.9*min([std(x3) iqr(x3)/1.34])*n3^(-0.2);
    
    TCR1ds = ksdensity(x1,c1,'Function','cdf','Bandwidth',bw1);
    TCR2ds = ksdensity(x2,c2,'Function','cdf','Bandwidth',bw2) - ...
             ksdensity(x2,c1,'Function','cdf','Bandwidth',bw2);
    TCR3ds = 1 - ksdensity(x3,c2,'Function','cdf','Bandwidth',bw3);    
    
    
    for i = 1:nmarkers
        beta_vec1(i) = append('β',string(i));
        beta_vec2(i) = append(string(beta(i)));
    end
    
    disp('The combination coefficients are:')
    disp(append(string(beta_vec1),' = ',string(beta_vec2)))
    disp(' ')
    disp('The optimal cutoff values are:')
    disp(append('c1 = ',string(c1)))
    disp(append('c2 = ',string(c2)))
    disp(' ')
    disp('The Euclidean distance is:')
    disp(append('D* = ',string(dstar)))
    disp(' ')
    disp('The true classification rates are:')
    disp(append('[TCR1, TCR2, TCR3] = [',string(TCR1ds),', ',string(TCR2ds),', ',string(TCR3ds),']')) 
    disp(' ')
    
    varargout{1} = beta;
    varargout{2} = c1;
    varargout{3} = c2;
    varargout{4} = dstar;
    
    bw1 = 0.9*min([std(x1) iqr(x1)/1.34])*length(x1)^(-0.2);
    bw2 = 0.9*min([std(x2) iqr(x2)/1.34])*length(x2)^(-0.2);
    bw3 = 0.9*min([std(x3) iqr(x3)/1.34])*length(x3)^(-0.2);
    
    TCR1ds = ksdensity(x1,c1,'Function','cdf','Bandwidth',bw1);
    TCR2ds = ksdensity(x2,c2,'Function','cdf','Bandwidth',bw2) - ksdensity(x2,c1,'Function','cdf','Bandwidth',bw2);
    TCR3ds = 1 - ksdensity(x3,c2,'Function','cdf','Bandwidth',bw3);
    
        switch rocplots
            case 'yes'      
               gr1 = 0:.02:1;
               gr2 = 0:.02:1;
               [pp1,pp3] = meshgrid(gr1,gr2);
               p1=pp1(1:end)';
               p3=pp3(1:end)';

                for i = 1:length(p1)
                        ROCsurf(i) = (ksdensity(x2,ksdensity(x3,1-p3(i),'Function','icdf','Bandwidth',bw3),'Function','cdf','Bandwidth',bw2) - ...
                         ksdensity(x2,ksdensity(x1,p1(i),'Function','icdf','Bandwidth',bw1),'Function','cdf','Bandwidth',bw2));
                end
    
                p1=reshape(p1,size(pp1,1),size(pp1,2));
                p3=reshape(p3,size(pp3,1),size(pp3,2));
    
                ROCsurf=reshape(ROCsurf,size(pp1,1),size(pp1,2));
    
                ROCsparam = @(p1,p3) (ksdensity(x3,1-p3,'Function','icdf','Bandwidth',bw3)>ksdensity(x1,p1,'Function','icdf','Bandwidth',bw1)).*(ksdensity(x2,ksdensity(x3,1-p3,'Function','icdf','Bandwidth',bw3),'Function','cdf','Bandwidth',bw2) - ...
                                     ksdensity(x2,ksdensity(x1,p1,'Function','icdf','Bandwidth',bw1),'Function','cdf','Bandwidth',bw2));
   
    
                figure('Position',[0 0 700 600])
                surf(p1,p3,ROCsurf)
                hold on
                view(220,10)
                xlim([0 1])
                ylim([0 1])
                zlim([0 1])
                line([1 TCR1ds], [1 TCR3ds], [1 TCR2ds],'Color','r','LineWidth',1.5)
                shading interp
                alpha 0.6
                xlabel('TCR1')
                ylabel('TCR3')
                zlabel('TCR2')
                set(gcf,'color','w');
                hold off
                
                figure('Position',[0 0 1000 600])
                
                [f1k,x1k] = ksdensity(x1);
                [f2k,x2k] = ksdensity(x2);
                [f3k,x3k] = ksdensity(x3);
                
                h1 = histogram(x1,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 0 .7],'FaceColor',[0 0 .7],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                hold on
                h2 = histogram(x2,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 .6 0],'FaceColor',[0 .6 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                h3 = histogram(x3,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[.7 0 0],'FaceColor',[.7 0 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);

                plot(x1k,f1k,'Color',[0 0 .7],'LineWidth',2)
                hold on
                plot(x2k,f2k,'Color',[0 .6 0],'LineWidth',2)
                plot(x3k,f3k,'Color',[.7 0 0],'LineWidth',2)
                max1 = max(h1.Values);
                max2 = max(h2.Values);
                max3 = max(h3.Values);
                maxf = max([max1 max2 max3]);
                ylim([0 1.25*maxf])
                plot([c1 c1],[0 1.25*maxf],'--k')
                text(c1+.05,maxf*1.1,append('c_{1} = ',string(round(c1,4))))
                plot([c2 c2],[0 1.25*maxf],'--k')
                text(c2+.05,maxf*1.1,append('c_{2} = ',string(round(c2,4))))
                set(gcf,'color','w');
        end
    
    
    case 'stepwise'
        
    if sum(strncmp(normtest,'yes',3)) == 1
    roc3d_norm_test_comb(y,D,.05,'no');
    end
    
    [beta,c1,c2,dstar,TCR1ds,TCR2ds,TCR3ds] = dstar_comb_step(y1,y2,y3);
    
    x1 = y1*beta';
    x2 = y2*beta';
    x3 = y3*beta';
    x = [x1;x2;x3];
    for i = 1:nmarkers
        beta_vec1(i) = append('β',string(i));
        beta_vec2(i) = append(string(beta(i)));
    end
    
    disp('The combination coefficients are:')
    disp(append(string(beta_vec1),' = ',string(beta_vec2)))
    disp(' ')
    disp('The optimal cutoff values are:')
    disp(append('c1 = ',string(c1)))
    disp(append('c2 = ',string(c2)))
    disp(' ')
    disp('The Euclidean distance is:')
    disp(append('D* = ',string(dstar)))
    disp(' ')
    disp('The true classification rates are:')
    disp(append('[TCR1, TCR2, TCR3] = [',string(TCR1ds),', ',string(TCR2ds),', ',string(TCR3ds),']')) 
    disp(' ')
    
    varargout{1} = beta;
    varargout{2} = c1;
    varargout{3} = c2;
    varargout{4} = dstar;
    
    switch rocplots
        case 'yes'
            figure('Position',[0 0 700 600])
            emp_dstar_plot(x1,x2,x3)
            view(220,10)
            
            figure('Position',[0 0 000 600])

                [f1k,x1k] = ksdensity(x1);
                [f2k,x2k] = ksdensity(x2);
                [f3k,x3k] = ksdensity(x3);
                
                h1 = histogram(x1,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 0 .7],'FaceColor',[0 0 .7],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                hold on
                h2 = histogram(x2,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 .6 0],'FaceColor',[0 .6 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                h3 = histogram(x3,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[.7 0 0],'FaceColor',[.7 0 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);

                plot(x1k,f1k,'Color',[0 0 .7],'LineWidth',2)
                hold on
                plot(x2k,f2k,'Color',[0 .6 0],'LineWidth',2)
                plot(x3k,f3k,'Color',[.7 0 0],'LineWidth',2)
                max1 = max(h1.Values);
                max2 = max(h2.Values);
                max3 = max(h3.Values);
                maxf = max([max1 max2 max3]);
                ylim([0 1.25*maxf])
                plot([c1 c1],[0 1.25*maxf],'--k')
                text(c1+.05,maxf*1.1,append('c_{1} = ',string(round(c1,4))))
                plot([c2 c2],[0 1.25*maxf],'--k')
                text(c2+.05,maxf*1.1,append('c_{2} = ',string(round(c2,4))))
                set(gcf,'color','w');
    end
    
    case 'logistic'
        
    if sum(strncmp(normtest,'yes',3)) == 1
    roc3d_norm_test_comb(y,D,.05,'no');
    end
    [beta,c1,c2,dstar,TCR1ds,TCR2ds,TCR3ds] = dstar_comb_log(y1,y2,y3);
    
    x1 = y1*beta;
    x2 = y2*beta;
    x3 = y3*beta;
    x = [x1;x2;x3];
    for i = 1:nmarkers
        beta_vec1(i) = append('β',string(i));
        beta_vec2(i) = append(string(beta(i)));
    end
    
    disp('The combination coefficients are:')
    disp(append(string(beta_vec1),' = ',string(beta_vec2)))
    disp(' ')
    disp('The optimal cutoff values are:')
    disp(append('c1 = ',string(c1)))
    disp(append('c2 = ',string(c2)))
    disp(' ')
    disp('The Euclidean distance is:')
    disp(append('D* = ',string(dstar)))
    disp(' ')
    disp('The true classification rates are:')
    disp(append('[TCR1, TCR2, TCR3] = [',string(TCR1ds),', ',string(TCR2ds),', ',string(TCR3ds),']')) 
    disp(' ')
    
    varargout{1} = beta;
    varargout{2} = [c1 c2];
    varargout{3} = dstar;
    varargout{4} = [TCR1ds TCR2ds TCR3ds];
    
    switch rocplots
        case 'yes'
            figure('Position',[0 0 700 600])
            emp_dstar_plot(x1,x2,x3)
            view(220,10)
            
            figure('Position',[0 0 1000 600])

                [f1k,x1k] = ksdensity(x1);
                [f2k,x2k] = ksdensity(x2);
                [f3k,x3k] = ksdensity(x3);
                
                h1 = histogram(x1,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 0 .7],'FaceColor',[0 0 .7],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                hold on
                h2 = histogram(x2,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[0 .6 0],'FaceColor',[0 .6 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);
                h3 = histogram(x3,'BinWidth',(max(x)-min(x))/15,'Normalization','pdf','EdgeColor',[.7 0 0],'FaceColor',[.7 0 0],'FaceAlpha',.3,'EdgeAlpha',.5,'LineWidth',1);

                plot(x1k,f1k,'Color',[0 0 .7],'LineWidth',2)
                hold on
                plot(x2k,f2k,'Color',[0 .6 0],'LineWidth',2)
                plot(x3k,f3k,'Color',[.7 0 0],'LineWidth',2)
                max1 = max(h1.Values);
                max2 = max(h2.Values);
                max3 = max(h3.Values);
                maxf = max([max1 max2 max3]);
                ylim([0 1.25*maxf])
                plot([c1 c1],[0 1.25*maxf],'--k')
                text(c1+.05,maxf*1.1,append('c_{1} = ',string(round(c1,4))))
                plot([c2 c2],[0 1.25*maxf],'--k')
                text(c2+.05,maxf*1.1,append('c_{2} = ',string(round(c2,4))))
                set(gcf,'color','w');
    end
end
end