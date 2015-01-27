function [low,up] = get_bnd (v,step_bnd); 
%-----------------------------------------
global bounds v_low v_up
%-----------------------------------------------
% Generate upper and lower bounds for the
% step in the multigrid algorithm.
%-----------------------------------------------
e   = ones(size(v));
low = v - step_bnd * e;
up  = v + step_bnd * e;

if (bounds);
  low = max(low, v_low);
  up  = min( up, v_up);
end;