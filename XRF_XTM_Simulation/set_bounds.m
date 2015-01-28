function [v_low,v_up] = set_bounds(j,res_prob);
global NumElement

v_low=0*ones(N(j)^2*NumElement);
v_up=1e6*ones(N(j)^2*NumElement);
