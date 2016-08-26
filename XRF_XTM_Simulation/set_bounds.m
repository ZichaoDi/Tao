function [v_low,v_up] = set_bounds(j,res_prob);
global NumElement N 

ntot=N(j)^2*NumElement;
v_low=0*ones(ntot,1);
v_up=inf*ones(ntot,1);
