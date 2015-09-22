function [f, g, shift_y, shift_yT] = sfun (W)
%--------------------------------------------------------------
global Joint ReconAttenu

if(Joint==1)
    [f,g,shift_y, shift_yT]=sfun_J(W);
elseif(Joint==0)
    [f, g, shift_y] = sfun_R (W);
    shift_yT=[];
elseif(Joint==-1)
    if(ReconAttenu)
        [f, g, JJ, shift_yT] = sfun_T(W);
    else
        [f, g, shift_yT] = sfun_Tw(W);
    end
    shift_y=[];
end
