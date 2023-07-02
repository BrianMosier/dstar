function [out]=likbox(x1,x2,x3,h)
%This calculates the value of the logL for a given h. The result is "out".
%Remember that x1,x2 and x3 are the original data.
n1=length(x1);
n2=length(x2);
n3=length(x3);
if h==0
    xh1=log(x1);
    xh2=log(x2);
    xh3=log(x3);
else 
    xh1=((x1.^h)-1)./h;
    xh2=((x2.^h)-1)./h;
    xh3=((x3.^h)-1)./h;
    
end

s1=sqrt(sum((xh1-sum(xh1)/n1).^2)./n1);%note that here I use the transformed data xh1,xh2,xh3 and NOT x1,x2,x3.
s2=sqrt(sum((xh2-sum(xh2)/n2).^2)./n2);
s3=sqrt(sum((xh3-sum(xh3)/n3).^2)./n3);

out1=(h-1).*(sum(log(x1))+sum(log(x2))+sum(log(x3))); %note that here I use the original data, x1,x2,x3.
out2=-n1/2.*(1+log(2*pi*s1.^2))-n2/2.*(1+log(2*pi*s2.^2))-n3/2.*(1+log(2*pi*s3.^2));
out=out1+out2;
end