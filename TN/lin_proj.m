function [x, f, g, nf1, ierror, alpha1, ipivot,newcon,flast,f_xrf,f_xtm] = ...
    lin_proj (p, x, f, g, alpha, sfun, low, up,ipivot,flast,newcon)
global Joint
%---------------------------------------------------------
% line search: P(x+alpha*p)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
ierror = 3;
maxit  = 5;
%%%%%%%%%%%%%%%%%%%#######################################################
if (alpha == 0); ierror = 0; maxit = 1; end;
alpha1 = alpha;
trialLength=2;
if(alpha1<=1)
    trialAlpha=linspace(1,alpha1,trialLength);
end
iproj=0;
q0=p'*g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
f_xrf=0;
f_xtm=0;
if(alpha1<=1)%& alpha1>0
    for trial=1:length(trialAlpha)
        xt = x + trialAlpha(trial).*p;
        [ipivot1,~, xt] = crash (xt, low, up);
        
        if(Joint==1)
            [ft, gt,f_xrf,f_xtm] = feval (sfun, xt);
        else
            [ft, gt] = feval (sfun, xt);
        end
        Armijo =ft<f+1e-4*trialAlpha(trial)*q0;
%         Wolfe = abs(p'*gt)<0.25*abs(q0);%
        if (Armijo );%& Wolfe);
%             fprintf('Armijo and Wolfe satisfied, trial= %d\n',trial);
            ierror = 0;
            iproj  = 1;
            x   = xt;
            f   = ft;
            g   = gt;
            alpha1=trialAlpha(trial);
            itcnt=trial;
            ipivot=ipivot1;
            p = ztime(p,ipivot);
            newcon = 1;
            flast=f;
            break;
        end;
    end
    
%         if(ierror==0)
%             for trial=1:length(trialAlpha)
%                 xt = x + trialAlpha(trial).*p;
%                 [ipivot1,~, xt] = crash (xt, low, up);
%                 [ft, gt] = feval (sfun, xt);
%                 Armijo =ft<f+1e-4*trialAlpha(trial)*q0;
%                 Wolfe = abs(p'*gt)<0.25*abs(q0);%
%                 if (Armijo & Wolfe);
%                     fprintf('Projected Direction Armijo and Wolfe satisfied, trial= %d\n',trial);
%                     ierror = 0;
%                     iproj  = 1;
%                     x   = xt;
%                     f   = ft;
%                     g   = gt;
%                     alpha1=trialAlpha(trial);
%                     itcnt=trial;
%                     ipivot=ipivot1;
%                     p = ztime (p, ipivot);
%                     newcon = 1;
%                     flast=f;
%                     break;
%                 end;
%             end
%         end
end
if (alpha1 == 0); ierror = 0; maxit = 1; end;
if(iproj==0)
    for itcnt = 1:maxit;
        xt = x + alpha1.*p;
        [ft, gt] = feval (sfun, xt);
        if (ft < f);
            ierror = 0;
            x   = xt;
            f   = ft;
            g   = gt;
            [ipivot1,~, x] = crash (x, low, up);
%             if(norm(ipivot1-ipivot,1)~=0)
            newcon = 1;
            flast=f;
            ipivot=ipivot1;
%             end
            break;
        end;
        alpha1 = alpha1 ./ 2;
    end;
    itcnt=itcnt+trialLength;
end
if (ierror == 3); alpha1 = 0; end;
nf1 = itcnt;
