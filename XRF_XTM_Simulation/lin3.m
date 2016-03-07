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
maxiTrial=5;
for itcnt = 1:maxiTrial;
    xt = x + alpha.*p;
    [ConstNew, I, O, ~, ~,ft] = feval (sfun, xt);
    if (ft < f | itcnt==1);
        ierror = 0;
        xnew   = xt;
        % if(icycle == maxOut) 
        Out= squeeze(sum(sum(O,1),2));
        xnew   = map1D(xt./map1D(Out(:),[1 80]),[min(xt),max(xt)]);
        % end
        fnew   = ft;
        ConstNew = ConstSub;
        break;
    elseif(itcnt==maxiTrial)
        fprintf('no step taken\n')
        alpha=0;
    end;
    alpha = alpha./2;
end


