function [conv, flast1, ipivot1] = cnvtst (alpha, pnorm, xnorm, ...
		 dif, ftest, gnorm, gtp, f, flast, g, ipivot, accrcy)
global lambdaind
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
% if appropriate, modify active set
%---------------------------------------------------------
if (~ltest)
   ind   = find(ipivot ~= 0 & ipivot ~= 2);
   if (length(ind)>0);
      t = -ipivot(ind).*g(ind);
      [cmax, imax] = min(t);
      if (cmax >= 0); imax = 0; end;
   end;
end;
if (imax ~= 0)
   ipivot(ind(imax)) = 0;
   flast = f;
else
   conv = (alpha*pnorm < toleps*(1 + xnorm) ...
           & abs(dif) < rtleps*ftest        ...
           & gnorm < accrcy^(1/3)*ftest)    ...
           | gnorm < .01*sqrt(accrcy)*ftest;
end;
flast1  = flast;
ipivot1 = ipivot;
lambdaind   = find((ipivot ~= 2) & (ipivot.*g>0));