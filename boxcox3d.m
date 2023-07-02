function [transx1, transx2, transx3, lam, exitflag]=boxcox3d(x1,x2,x3)
%Input: x1,x2,x3 i.e. the original data.
%Output: transx1, transx2, transx3 (i.e. the transformed data),
%"lam", which is the estimated BoxCox parameter 
%"exitflag", which informs us if convergence occured.

         logL=@(h) -likbox(x1,x2,x3,h); %this is minus the logL that depends on h (that is lambda)
         %the logL uses the likbox routine to calculate the value of logL
         %for a given h. i.e. it is a function of h, and will be used right
         %next to fminsearch for minimization:
         
         init=0;
         
         [lam , ~, exitflag]=fminsearch(logL,init); %minimize it with initial value 1.
         
         if lam==0 %which is never going to happen
             transx1=log(x1);
             transx2=log(x2);
             transx3=log(x3);
             
         else
             transx1=((x1.^lam)-1)./lam;
             transx2=((x2.^lam)-1)./lam;
             transx3=((x3.^lam)-1)./lam;
             
         end
             
end