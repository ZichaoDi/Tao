function [xnew, fnew, gnew, nf1, ierror, alpha1, varargout] = ...
    lin1 (p, x, f, alpha, g, sfun)
global bounds
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
maxit  = 15;
%%%%%%%%%%%%%%%%%%%#####fix step length as 1 ##############################
% alpha=1;
%%%%%%%%%%%%%%%%%%%#######################################################
% disp('### In LINE SEARCH: lin1 ###')
% fprintf(' g''p = %e\n',g'*p);
% if(length(x)==1151)
%     format long e
%     alpha
% end
if (alpha == 0); ierror = 0; maxit = 1; end;
alpha1 = alpha;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
for itcnt = 1:maxit;
    xt = x + alpha1.*p;
%     if (length(x)==17)
%         load vll;
%         figure(124);plot(1:17,xt,'r.-',1:17,vll,'b.-',1:17,x,'g.-');pause
%     end
    %%%%%%%%%%############# keep binding in feasibel ########################
%     if(bounds)
% %     outup=find(xt>=5);
% %     xt(outup)=5;
%     outlow=find(xt<=-0.01);
%     xt(outlow)=-0.01;
%     end
    %%%%%%%%%%%#############################################################
    [ft, gt] = feval (sfun, xt);
    % fprintf('f(trial) = %e, f(current) = %e, alpha = %e\n',ft,f,alpha1);
    if (ft < f);
        ierror = 0;
        xnew   = xt;
        fnew   = ft;
        gnew   = gt;
        break;
    end;
    alpha1 = alpha1 ./ 2;
end;
if (ierror == 3); alpha1 = 0; end;
nf1 = itcnt;

if (nargout == 7)
    dfdp = dot(gt, p);
    varargout{1} = dfdp;
end