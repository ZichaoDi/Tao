function [conv, flast1, ipivot1] = cnvtst (alpha, pnorm, xnorm, ...
		 dif, ftest, gnorm, gtp, f, flast, g, ipivot, accrcy)
global lambdaind 
global Cauchy 
%---------------------------------------------------------
% test for convergence
%---------------------------------------------------------
% set up
%---------------------------------------------------------
conv   = 0;
toleps = sqrt(accrcy) + sqrt(eps);
rtleps = accrcy + eps;
imax   = 0;
ltest  = (flast - f <= -0.5*gtp);
%---------------------------------------------------------
% if anti-zigzag test satisfied, test multipliers;
% if appropriate, modify active set, only reduce active set
%---------------------------------------------------------
if (~ltest)
%     disp('statisfy ipivot updated from cnvtst');
   ind   = find(ipivot ~= 0 & ipivot ~= 2);
   if (~isempty(ind));
      t = -ipivot(ind).*g(ind);
      [cmax, imax] = min(t);
      if (cmax >= 0); imax = 0; end;
   end;
end;
if (imax ~= 0)
    if(Cauchy==0)
 ipivot(ind(imax)) = 0;
    end
   flast = f;
else
   conv = (alpha*pnorm < toleps*(1 + xnorm) ...
           & abs(dif) < rtleps*ftest        ...
           & gnorm < accrcy^(1/3)*ftest)  ...
            | gnorm < .01*sqrt(accrcy)*ftest;
       if(conv)
           conv1=alpha*pnorm < toleps*(1 + xnorm)
           conv2=abs(dif) < rtleps*ftest 
           conv3=gnorm < accrcy^(1/3)*ftest 
           conv4=gnorm < .01*sqrt(accrcy)*ftest
       end
end;
flast1  = flast;
ipivot1 = ipivot;
lambdaind   = find((ipivot ~= 2) & (ipivot.*g>0));