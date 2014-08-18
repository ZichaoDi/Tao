function [xnew, fnew, gnew, nf1, alpha, ierror, varargout] = ...
    lin3 (p, x, f, alphamax, g, sfun)
global lambdatn row it v_low v_up
%---------------------------------------------------------
% line search (modified unconstant step-length)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
alpha  = alphamax;
alpha
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
%=======================================================
%===========test merit function=========================

% alp=linspace(0,2,20);
% for i=1:length(alp)
%     vt=x+alp(i)*p;
%     [F,G]=feval(sfun,vt);
%     [fm,gm]=merit_tn(vt,F,G);
%     meritf(i)=fm;
%     objf(i)=F;
% end
% figure(22);%subplot(5,2,it);
% plot(alp,meritf,'r*-',alp,objf,'b.-'); title('red: merit function; blue: original function');%pause;
% %==========================
% format long e
% oldf=f
for itcnt = 1:20;
    xt = x + alpha.*p;
    %indl=find(xt<v_low)
    %     xt(indl)=v_low(indl);
    %          indu=find(xt>v_up);
    %     xt(indu)=v_up(indu);
    %     format long e
    %     if(~isempty(ind))
    %         alpha = alpha./2;
    %     else
    [ft, gt] = feval (sfun, xt);
    [ft,gt]=merit_tn(xt,ft,gt);
    if (ft < f);
        ierror = 0;
        xnew   = xt;
        fnew   = ft;
        gnew   = gt;
        fprintf('MG/Opt line search: alpha = %16.8e\n',alpha)
        break;
    elseif(itcnt==20)
        fprintf('no step taken')
    end;
    alpha = alpha./2;
    
    
    %end;
end


nf1 = itcnt;

if (nargout == 7)
    dfdp = dot(gt, p);
    varargout{1} = dfdp;
end
