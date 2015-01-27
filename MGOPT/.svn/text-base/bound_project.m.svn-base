function v_out = bound_project(v_in,v_low,v_up)
%-------------------------------------------------------
% Project v_in onto the box defined by the bound 
% constraints v_low, v_up
%-------------------------------------------------------
% Usage: v_out = bound_project(v_in,v_low,v_up)
%-------------------------------------------------------
v_out = min(v_in,v_up);
v_out = max(v_out,v_low);