function [F, G] = sfun_test(v)
%--------------------------------------------------------------
% compute objective function and gradient [residual form]
%--------------------------------------------------------------
load current_fnl
%--------------------------------------------------------------

[F,G] = sfun(v);
F     = F - v'*current_fnl;
G     = G - current_fnl;
