function [xnew, fnew, alpha] = ...
    lin3 (p, x, f, alphamax, sfun)
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
alpha  = alphamax;
xnew   = x;
fnew   = f;
for itcnt = 1:20;
    xt = x + alpha.*p;
    [~,~,~,~,ft] = feval (sfun, xt);
    if (ft < f);
        ierror = 0;
        xnew   = xt;
        fnew   = ft;
        break;
    elseif(itcnt==20)
        fprintf('no step taken\n')
        alpha=0;
    end;
    alpha = alpha./2;
end


