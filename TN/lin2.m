function [xnew, fnew, gnew, nf1, alpha, ierror, varargout] = ...
	lin2 (p, x, f, alphamax, g, sfun)
%---------------------------------------------------------
% line search (naive)
%---------------------------------------------------------
% set up
%---------------------------------------------------------
ptn=p;
alpha  = alphamax;
ierror = 3;
xnew   = x;
fnew   = f;
gnew   = g;
%---------------------------------------------------------
% line search
%---------------------------------------------------------
for itcnt = 1:20;
   xt = x + alpha*p;
   [ft, gt] = feval (sfun, xt);
   if (ft < f);
      ierror = 0;
      xnew   = xt;
      fnew   = ft;
      gnew   = gt;
      fprintf('MG/Opt line search: alpha = %16.8e\n',alpha)
      break;
   end;
   alpha = alpha / 2;
end;

nf1 = itcnt;

if (nargout == 7)
  dfdp = dot(gt, p);
  varargout{1} = dfdp;
end
