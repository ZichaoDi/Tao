function [v_low,v_up] = set_bounds(j,res_prob);
global NumElement N nTau numThetan

ntot=N(j)^2*NumElement*(nTau+1)*numThetan;
v_low=0*ones(ntot,1);
v_up=1e6*ones(ntot,1);
