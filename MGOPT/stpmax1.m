function alpha = stpmax1 (v, e, v_low, v_up)
%------------------------------------------------
% compute the maximum allowable step length
%------------------------------------------------
% Usage: alpha = stpmax1 (v, e, v_low, v_up)
%------------------------------------------------
al = 1e20; % initialize in case "if" statements
au = 1e20; % are not executed
indl = find(e < 0);
if (~isempty(indl));
   tl   = v_low(indl) - v(indl);
   al   = min(tl./e(indl));

end;
%------------------------------------------------
indu = find(e > 0 );
if (~isempty(indu));
   tu   = v_up (indu) - v(indu);
   au   = min(tu./e(indu));
end;
%------------------------------------------------
alpha  = min([1 al au]);