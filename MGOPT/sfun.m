function [f, g, r] = sfun (W)
%--------------------------------------------------------------
global xtm_level I0 L_level current_n frame N

j = find(current_n==N);
[f, g, r] = sfun_radon(W,xtm_level{j},I0,L_level{j});



