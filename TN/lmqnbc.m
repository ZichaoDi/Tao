function [xstar, f, g, ierror] = ...
    lmqnbc (x, sfun, low, up, maxit, maxfun, stepmx, accrcy)
%---------------------------------------------------------
% This is a bounds-constrained truncated-newton method.
% The truncated-newton method is preconditioned by a
% limited-memory quasi-newton method (computed by
% this routine) with a diagonal scaling (routine ndia3).
% For further details, see routine tnbc.
%---------------------------------------------------------
global sk yk sr yr yksk yrsr
global NF N current_n  fiter itertest
global ptest gv ipivot nit
global  i_cauchy
global W0 VarInd m
%---------------------------------------------------------
% check that initial x is feasible and that the bounds
% are consistent
%---------------------------------------------------------
[f, g] = feval (sfun, x);
oldf   = f;
fiter(1)=oldf;
[ipivot, ierror, x] = crash(x, low, up);
if (ierror ~= 0);
    disp('LMQNBC: termination 0')
    fprintf(1,'LMQNBC: terminating (no feasible point)');
    f = 0;
    g = zeros(size(x));
    xstar = x;
    return;
end;
%---------------------------------------------------------
% initialize variables, parameters, and constants
%---------------------------------------------------------
fprintf(1,'  it     nf     cg           f             |g|\n')
nind   = find(N==current_n);
upd1   = 1;
ncg    = 0;
conv   = 0;
xnorm  = norm(x,'inf');
ierror = 0;
if (stepmx < sqrt(accrcy) | maxfun < 1);
    disp('LMQNBC: termination 1')
    ierror = -1;
    xstar = x;
    return;
end;
%---------------------------------------------------------
% compute initial function value and related information
%---------------------------------------------------------
[f, g] = feval (sfun, x);
g0=g;
nf     = 1;
nit    = 0;
flast  = f;
itertest(1)=nf+ncg;
%---------------------------------------------------------
% Test if Lagrange multipliers are non-negative.
% Because the constraints are only bounds, the Lagrange
% multipliers are components of the gradient.
% Then form the projected gradient.
%---------------------------------------------------------
ind = find((ipivot ~= 2) & (ipivot.*g>0));
if (~isempty(ind));
    ipivot(ind) = zeros(length(ind),1);
end;
g = ztime (g, ipivot);
gnorm = norm(g,'inf');
fprintf(1,'%4i   %4i   %4i   % .8e   %.1e\n', ...
    nit, nf, ncg, f, gnorm)
%---------------------------------------------------------
% check if the initial point is a local minimum.
%---------------------------------------------------------
ftest = 1 + abs(f);
if (gnorm < .01*sqrt(eps)*ftest);
    disp('LMQNBC: termination 2')
    xstar = x;
    NF(1,nind) = NF(1,nind) + nit;
    NF(2,nind) = NF(2,nind) + nf;
    NF(3,nind) = NF(3,nind) + ncg;
    return;
end;
%---------------------------------------------------------
% set initial values to other parameters
%---------------------------------------------------------
n      = length(x);
icycle = n-1;
ireset = 0;
bounds = 1;
difnew = 0;
epsred = .05;
fkeep  = f;
i_cauchy=0;
d      = ones(n,1);
%---------------------------------------------------------
% ..........main iterative loop..........
%---------------------------------------------------------
% compute the search direction
%---------------------------------------------------------
argvec = [accrcy gnorm xnorm];
%%%%%%%%%%%%%%%%%%%%%##############################
% x_tran = x'
%%%%%%%%%%%%%%%%%%%%%##############################
[p, gtp, ncg1, d, eig_val] = ...
    modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun);
ncg = ncg + ncg1;
while (~conv);
    oldg = g;
    pnorm = norm(p, 'inf');
    oldf = f;
    fiter(nit+2)=oldf;
    itertest(nit+2)=nf+ncg;
    
    %---------------------------------------------------------
    % line search
    %---------------------------------------------------------
    if(i_cauchy~=0)
        p = ztime (p, ipivot);
    end
    pe = pnorm + eps;
    [spe] = stpmax (stepmx, pe, x, p, ipivot, low, up);
    alpha = step1 (f, gtp, spe);
    alpha0 = alpha;
    PieceLinear=1;
    if(PieceLinear)
        [x_new, f_new, g_new, nf1, ierror, alpha] = lin_proj (p, x, f, g, alpha0, sfun, low, up);
    else
        [x_new, f_new, g_new, nf1, ierror, alpha] = lin1 (p, x, f, alpha0, g, sfun);
    end
    Cauchy=0;
    %     if(gtp>0)
    %         Cauchy=1;
    %     end
    %---------------------------------------------------------
    save x_new x_new;
    if (alpha == 0 & alpha0 ~= 0 | ierror == 3);
        fprintf('Error in Line Search\n');
        fprintf('    ierror = %3i\n',    ierror);
        fprintf('    alpha  = %12.6f\n', alpha);
        fprintf('    alpha0 = %12.6f\n', alpha0);
        fprintf('    g''p    = %12.4e\n', gtp);
        %############################
        fprintf('    |g|     = %12.4e\n', norm(g));
        fprintf('    |p|     = %12.4e\n', norm(p));
        %         disp('Hit any key to continue')
        %         pause;
        if(gtp<0 & norm(pold-p)>0)
            [ipivot, ierror, x_new] = crash2(x_new, g_new, low, up);
        end  
    end;
    %#######################
    x   = x_new;
    f   = f_new;
    g   = g_new;
    g0=g;
    %#######################
    nf  = nf  + nf1;
    nit = nit +   1;
    %---------------------------------------------------------
    % update active set, if appropriate
    %---------------------------------------------------------
    newcon = 0;
    if (abs(alpha-spe) <= 10*eps & i_cauchy==0);%
        newcon = 1;
        ierror = 0;
        [ipivot, flast] = modz (x, p, ipivot, low, up, flast, f, alpha);
    end;
    
    if (ierror == 3);
        disp('LMQNBC: termination 3')
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    %---------------------------------------------------------
    % stop if more than maxfun evaluations have been made
    %---------------------------------------------------------
    if (nf > maxfun);
        disp('LMQNBC: termination 4')
        ierror = 2;
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    %---------------------------------------------------------
    % set up for convergence and resetting tests
    %---------------------------------------------------------
    difold = difnew;
    difnew = oldf - f;
    if (icycle == 1);
        if (difnew >  2*difold); epsred =  2*epsred; end;
        if (difnew < .5*difold); epsred = .5*epsred; end;
    end;
    gv    = ztime (g, ipivot);
    gnorm = norm(gv, 'inf');
    ftest = 1 + abs(f);
    xnorm = norm(x,'inf');
    fprintf(1,'%4i   %4i   %4i   % .8e   %.1e\n', ...
        nit, nf, ncg, f, gnorm)
    %--------------------------------- Error Plot
%     if(mod(nit,50)==0)
%         figure(100);
%         CurrentErr=abs(x-W0(VarInd));
%         subplot(1,2,1)
%         vs=sum(reshape(CurrentErr,m(1),m(2),length(VarInd)/prod(m)),3);
%         surf(vs);
%         %     [xx,yy] = meshgrid(1:.02:m(1),1:.02:m(2));
%         % [X,Y] = meshgrid(1:m(1),1:m(2));
%         % vv = interp2(X,Y,vs,xx,yy,'cubic');
%         % surfl(xx,yy,vv);
%         subplot(1,2,2)
%         plot(abs(x-W0(VarInd)),'ro-');
%     end
    %---------------------------------------------------------
    % test for convergence
    %---------------------------------------------------------
    [conv, flast, ipivot] = cnvtst (alpha, pnorm, xnorm, ...
        difnew, ftest, gnorm, gtp, f, flast, g, ...
        ipivot, accrcy);
    if (conv);
        disp('LMQNBC: termination 5')
        xstar = x;
        NF(1,nind) = NF(1,nind) + nit;
        NF(2,nind) = NF(2,nind) + nf;
        NF(3,nind) = NF(3,nind) + ncg;
        return;
    end;
    %%%######################## Update active set using Cauchy point
    if(Cauchy)
        [~,i_cauchy,alpha_cp,p_cp]=Cauchy_Point(x,g0, accrcy, xnorm, sfun,low,up);
        if(i_cauchy~=0)
            if(length(alpha_cp)==1)
                alpha_cp
            end
            [x, f, g, nf1, ierror] = lin_Cauchy (p_cp, x, f, g0, sfun, alpha_cp);
            if(ierror==3)
                disp('descent not found from Cauchy point')
            else
                nf = nf+nf1;
                [ipivot] = crash2 (x,g0, low, up);
            end
        else
            disp('Cauchy point not found')
        end
    end
    %------------------------------- Active set plot
%         figure(90);qpPlotAset(ipivot,nit,length(x));
    %%%===========================================================
    
    g = ztime (g, ipivot);
    %---------------------------------------------------------
    % modify data for LMQN preconditioner
    %---------------------------------------------------------
    if (~newcon);
        yk = g - oldg;
        sk = alpha*p;
        yksk = yk'*sk;
        ireset = (icycle == n-1 | difnew < epsred*(fkeep-f));
        if (~ireset);
            yrsr = yr'*sr;
            ireset = (yrsr <= 0);
        end;
        upd1 = (yksk <= 0);
    end;
    %---------------------------------------------------------
    % compute the search direction
    %---------------------------------------------------------
    argvec = [accrcy gnorm xnorm];
    %%%%%%%%%%%%%%%%%%%%%##############################
    % x_tran = x'
    %%%%%%%%%%%%%%%%%%%%%##############################
    pold=p;
    [p, gtp, ncg1, d, eig_val] = ...
        modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, argvec, sfun);
    ptest=p;
    ncg = ncg + ncg1;
    % %---------------------------------------------------------
    % update LMQN preconditioner
    %---------------------------------------------------------
    if (~newcon);
        if (ireset);
            sr     = sk;
            yr     = yk;
            fkeep  = f;
            icycle = 1;
        else
            sr     = sr + sk;
            yr     = yr + yk;
            icycle = icycle + 1;
        end
    end;
end;
