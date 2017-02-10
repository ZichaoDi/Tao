function [xnew, fnew, ConstNew, alpha] = ...
    lin3 (p, x, f, alphamax, sfun,ConstSub)
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
global icycle maxOut 
alpha  = alphamax;
xnew   = x;
fnew   = f;
ConstNew = ConstSub;
maxiTrial=2;
for itcnt = 1:maxiTrial;
    xt = x + alpha.*p;
    [ConstSub, ft] = feval (sfun, xt);
    if (ft < f )%| itcnt==1);
        ierror = 0;
        xnew   = xt;
        fnew   = ft;
        ConstNew = ConstSub;
        break;
    elseif(itcnt==maxiTrial)
        fprintf('no step taken\n')
        alpha=0;
    end;
    alpha = alpha./2;
end


