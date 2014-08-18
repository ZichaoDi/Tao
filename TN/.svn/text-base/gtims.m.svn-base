function gv = gtims (v, x, g, accrcy, xnorm, sfun)
%%%%%%%%%%%%%%%%%##################
global bounds v_low v_up
%%%%%%%%%%%%%%%%%##################
%---------------------------------------------------------
% compute the product of the Hessian times the vector v;
% store result in the vector gv 
% (finite-difference version)
%---------------------------------------------------------
%%%%%%%%%%%%%%%%%##################
%%##  delta   = sqrt(accrcy)*(1 + xnorm);
delta   = sqrt(accrcy)*(1 + xnorm)/norm(v);
hg      = x + delta*v;
%%## if ( (max(hg) > max(v_up)+1) | (min(hg) < min(v_low)-1) )
%%##   format short e
%%##   more on
%%##   disp('GTIMS: hg')
%%##   accrcy
%%##   delta
%%##   nstep = norm(v)
%%##   vmat = [v_low hg v_up]
%%##   pause
%%## end;
%%## if ( (max(x) > max(v_up)+1) | (min(x) < min(v_low)-1) )
%%##   format short e
%%##   more on
%%##   disp('GTIMS: v')
%%##   vmat = [v_low x v_up]
%%##   pause
%%## end;
%%%%%%%%%%%%%%%%%##################
[f, gv] = feval (sfun, hg);
gv      = (gv - g)/delta;
