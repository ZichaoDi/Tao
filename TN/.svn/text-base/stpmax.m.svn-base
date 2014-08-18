function spe = stpmax (stepmx, pe, x, p, ipivot, low, up)
%------------------------------------------------
% compute the maximum allowable step length
% (spe is the standard (unconstrained) max step)
%------------------------------------------------
spe  = stepmx / pe;
al   = spe;
au   = spe;
%------------------------------------------------
indl = find(ipivot==0 & p < 0);
if (length(indl) > 0);
   tl   = low(indl) - x(indl);
   [al,il]   = min(tl./p(indl));
end;
%------------------------------------------------
indu = find(ipivot==0 & p > 0);
if (length(indu) > 0);
   tu   = up (indu) - x(indu);
   au   = min(tu./p(indu));
end;
%------------------------------------------------
spe  = min([spe al au]);


