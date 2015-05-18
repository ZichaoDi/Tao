function [x,fc,gc,eig_val,histout,costdata] = gaussn(x0,f,tol,maxit)
global WS
global NF N current_n
% C. T. Kelley, Dec 14, 1997
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,histout,costdata] = gaussn(x0,f)
%
% Damped Gauss-Newton with Armijo rule
% simple divide by 2 stepsize reduction
%
% Input: x0 = initial iterate
%        f = r^T r/2 = objective function,
%            the calling sequence for f should be
%            [fout,gout,jac]=f(x) where fout=f(x) is a scalar
%              gout = jac^T r = grad f(x) is a COLUMN vector
%              and jac = r' = Jacobian of r is an M x N matrix
%        tol = termination criterion norm(grad) < tol
%        maxit = maximum iterations (optional) default = 100
%
% Output: x = solution
%         histout = iteration history
%             Each row of histout is
%       [norm(grad), f, number of step length reductions, iteration count]
%         costdata = [num f, num grad, num hess] (for gaussn, num hess=0)
%
% At this stage all iteration parameters are hardwired in the code.

alp=1.d-4;
if nargin < 4
    maxit=100;
end
itc=1; xc=x0;
[fc,gc,jac]=feval(f,xc);
numf=1; numg=1; numh=0;
eig_val=[];
ithist=zeros(1,4);
ithist(1,1)=norm(gc); ithist(1,2) = fc; ithist(1,4)=itc-1; ithist(1,3)=0;
fprintf(1,'iter     LinIter          f               |g|                error\n');

while(norm(gc) > tol & itc <= maxit)
%     tic;
    dc= (jac'*jac)\gc;%pcg(jac'*jac,gc);%pinv(jac'*jac)*gc;%
%     eig_val{itc} = sort(eig(jac'*jac));
    xold=xc;
    lambda=1.0; xt=xc-lambda*dc;
    [ft,gt,jat]=feval(f,xt); numf=numf+1;
    iarm=0; itc=itc+1;
    fgoal= fc - alp*lambda*(gc'*dc);
    while(ft > fgoal)
        iarm=iarm+1;
        lambda=lambda/2;
        fgoal= fc - alp*lambda*(gc'*dc);
        xt=xc-lambda*dc;
        [ft,gt,jat]=feval(f,xt); numf=numf+1;
        if(iarm > 10)
            disp(' Armijo error in Gauss-Newton')
            x=xc; histout=ithist(1:itc-1,:);
            costdata=[numf, numg, numh];
            return; 
        end
    end
    xc=xt; fc=ft; gc=gt;jac=jat;numf=numf+1; numg=numg+1;
    ee(:,itc-1)=WS(:)-xc;%xc-xold;
    %%%-----------------------------------------------------------------
%     if(itc>2)
%     figure(15),plot(ee(:,itc-1)./ee(:,itc-2));%plot((WS(:)-xold)./(xt-xold),'ro-')
%     end
%     figure(16),subplot(1,3,1),plot(1:current_n^2,WS(:)-xold,'b.-',1:current_n^2,xt-xold,'ro-');
%     subplot(1,3,2),imagesc(reshape(WS(:)-xold,current_n,current_n));
%     subplot(1,3,3),imagesc(reshape(xt-xold,current_n,current_n));
%     drawnow;
%     %%%-----------------------------------------------------------------
%     save ee ee
%     gg(:,itc-1)=gc;
    ithist(itc,1)=norm(gc); ithist(itc,2) = fc;
    ithist(itc,3)=iarm; ithist(itc,4)=itc-1; 
    fprintf(1,'%3.0d       %3.1d       %d      %d      %d\n', ithist(itc,end:-1:1),norm(ee(:,itc-1)));
%     toc;
end

nind  = find(N==current_n);
        NF(1,nind) = NF(1,nind) + itc-1;
        NF(2,nind) = NF(2,nind) + iarm;
x=xc; histout=ithist(1:itc,:);
costdata=[numf, numg, numh];
