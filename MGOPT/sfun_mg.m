function [F, G, F_XRF, F_XTM] = sfun_mg (v)
%--------------------------------------------------------------
% compute objective function and gradient [residual form]
%------------------------- -------------------------------------
global current_fnl Joint
%--------------------------------------------------------------
if(Joint==1)
    [F,G, F_XRF, F_XTM] = sfun(v);
    F_XRF    = F_XRF - v'*current_fnl;
    F_XTM   = F_XTM - v'*current_fnl;
else
    [F,G] = sfun(v);
end
F     = F - v'*current_fnl;
G     = G - current_fnl;