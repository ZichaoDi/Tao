function [xnew, fnew, InNew, OutNew, AttNew, alpha] = ...
    lin3 (p, x, f, alphamax, sfun,InTens, OutTens, AttenuM)
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
alpha  = alphamax;
xnew   = x;
fnew   = f;
InNew  = InTens;
OutNew = OutTens;
AttNew = AttenuM;
maxiTrial=5;
for itcnt = 1:maxiTrial;
    xt = x + alpha.*p;
    [InTens, OutTens, AttenuM, DW,ft] = feval (sfun, xt);
    if (ft < f);
        ierror = 0;
        xnew   = xt;
        fnew   = ft;
        InNew  = InTens;
        OutNew = OutTens;
        AttNew = AttenuM;
        break;
    elseif(itcnt==maxiTrial)
        fprintf('no step taken\n')
        alpha=0;
    end;
    alpha = alpha./2;
end


