function [xnew, fnew, gnew, nf1, ierror, alpha1, varargout] = ...
    lin_proj (p, x, f, g, alpha, sfun, low, up)
%---------------------------------------------------------
% line search: P(x+alpha*p)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
maxit  = 5;
%%%%%%%%%%%%%%%%%%%#######################################################
if (alpha == 0); ierror = 0; maxit = 1; end;
alpha1 = alpha;
trialLength=2;
if(alpha1<1)
    trialAlpha=linspace(1,alpha1,trialLength);
end
iproj=0;
q0=p'*g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
if(alpha1<1)%& alpha1>0
    for trial=1:length(trialAlpha)
        xt = x + trialAlpha(trial).*p;
        [~,~, xt] = crash (xt, low, up);
        [ft, gt] = feval (sfun, xt);
%         Armijo =ft<f+1e-4*trialAlpha(trial)*q0;
Armijo =ft<f+1e-4*g'*(xt-x);
%         Wolfe = abs(p'*gt)<0.25*abs(q0);
        if (Armijo);%
            fprintf('Armijo satisfied, trial= %d\n',trial);
            ierror = 0;
            iproj  = 1;
            xnew   = xt;
            fnew   = ft;
            gnew   = gt;
            alpha1=trialAlpha(trial);
            itcnt=trial;
            break;
        end;
    end
end
if (alpha1 == 0); ierror = 0; maxit = 1; end;
if(iproj==0)
    for itcnt = 1:maxit;
        xt = x + alpha1.*p;
        [ft, gt] = feval (sfun, xt);
        if (ft < f);
            ierror = 0;
            xnew   = xt;
            fnew   = ft;
            gnew   = gt;
            break;
        end;
        alpha1 = alpha1 ./ 2;
    end;
    itcnt=itcnt+trialLength-1;
end
if (ierror == 3); alpha1 = 0; end;
nf1 = itcnt;

if (nargout == 7)
    dfdp = dot(gt, p);
    varargout{1} = dfdp;
end