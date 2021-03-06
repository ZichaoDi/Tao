function [F, G, F_XRF, F_XTM] = sfun_mg (v)
%--------------------------------------------------------------
% compute objective function and gradient [residual form]
%------------------------- -------------------------------------
global current_fnl 
%--------------------------------------------------------------
[F,G] = sfun(v);
F     = F - v'*current_fnl;
G     = G - current_fnl;
