function [f, g, shift_y, shift_yT] = sfun (W)
%--------------------------------------------------------------
global Joint

if(Joint==1)
    [f,g,shift_y, shift_yT]=sfun_J(W);
elseif(Joint==0)
    [f, g, shift_y] = sfun_R (W);
    shift_yT=[];
end