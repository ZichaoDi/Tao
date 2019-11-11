function [xstar, f1, g1, ierror, eig_val] = lBFGSm(x0,sfun, maxiter)
% lBFGS method with Armijo linesearch

% Inputs:
%     x       = initial guess
%     sfun    = objective function; function routine: [f,g] = sfun(x)

% Output:

% xstar  <-  computed solution.
% g      <-  final value of the gradient
% f      <-  final value of the objective function
% ierror <-  error code
%            ( 0 => normal return)
%            ( 2 => more than maxfun evaluations)
%            ( 3 => line search failed (may not be serious)
%            (-1 => error in input parameters)    
    global mbfgs
    global NF N current_n gnorm0
    global relGnormTol relGnormTolCoarsest relGnormTolFine
    global absGnormTol absGnormTolCoarsest absGnormTolFine
    global f0 onCoarsestGrid onFineGrid
    
    eig_val{1} = []; %hack
    if (~exist('m', 'var')); mbfgs = 5; end
    if (~exist('maxiter', 'var')); maxiter = 50; end

    format compact
    format short e
    fprintf(1,'iter\tlsit\tf\t\tf/f0\t\t|x-xOld|\t|g|\t\t|g|/|g0|\n')
    nind = find(N==current_n);
    
    %---------------------------------------------------------
    % create storage arrays for y_is and s_is
    %---------------------------------------------------------
    Ym = zeros(length(x0(:)), mbfgs);
    Sm = zeros(length(x0(:)), mbfgs);
    
    
    [finit, g0] = feval (sfun, x0);
    
    gnorm       = norm(g0);
    
    relGnorm     = gnorm/gnorm0;

    nf          = 0;
    nit         = 0;
    stepNorm    = 0; % Samy
    fprintf(1,'%d\t%d\t%1.3e\t%1.3e\t%1.2e\t%1.2e\t%1.2e\n', ...
    nit, nf, finit, finit/f0, stepNorm, gnorm, relGnorm)

    % choose tolerance according to current grid level
    if onCoarsestGrid 
        currentRelGnormTol = relGnormTolCoarsest;
        currentAbsGnormTol = absGnormTolCoarsest; 
    elseif onFineGrid
        currentRelGnormTol = relGnormTolFine;
        currentAbsGnormTol = absGnormTolFine;
    else
        currentRelGnormTol = relGnormTol;
        currentAbsGnormTol = absGnormTol;
    end
    
    % check if conv. crit. already satisfied
    if (relGnorm <= currentRelGnormTol ...
                || gnorm <= currentAbsGnormTol)
            fprintf('\n convergence criteria satisfied. No iterations required \n')
            ierror = 0;
            xstar = x0;
            f1 = finit; g1 = g0; 
            NF(1,nind) = NF(1,nind) + nit;
            NF(2,nind) = NF(2,nind) + nf;
            return;
    end

    conv   = 0;
    
    
    %---------------------------------------------------------
    % compute search direction
    %---------------------------------------------------------
    p = -g0;
    alpha0=1;
    [x1, f1, g1, nf1, ierror, alpha, stepNorm] = lin1 (p, x0, finit, alpha0, g0, sfun);
    
    
%%% ---------------- lBFGS Loop ------------------- %%%%% 
    while (~conv)
        
        nit = nit+1;


      s0 = x1-x0;
      y0 = g1-g0;
      hdiag = s0'*y0/(y0'*y0);
      p = zeros(length(g0),1);
      if (nit<=mbfgs)
        % update S,Y
        Sm(:,nit) = s0;
        Ym(:,nit) = y0;
        [p, flag] = getHg_lbfgs(g1,Sm(:,1:nit),Ym(:,1:nit),hdiag);
        % never forget the minus sign
        p = -p;
      elseif (nit>mbfgs)
        Sm(:,1:(mbfgs-1))=Sm(:,2:mbfgs);
        Ym(:,1:(mbfgs-1))=Ym(:,2:mbfgs);
        Sm(:,mbfgs) = s0;
        Ym(:,mbfgs) = y0;    
        [p, flag] = getHg_lbfgs(g1,Sm,Ym,hdiag);
        % never forget the minus sign
        p = -p;
      end 
    
        %---------------------------------------------------------
        % line search
        %---------------------------------------------------------
        
        if flag==0
            alpha0 = 1;
        else
            alpha0 = alpha*(g0'*g0)/(g1'*g1);
        end
        
        x0 = x1;
        g0 = g1;
        [x1, f1, g1, nf1, ierror, alpha, stepNorm] = lin1 (p, x1, f1, alpha0, g1, sfun);
        
        
        
        nf  = nf + nf1;

        gnorm    = norm(g1);
        relGnorm = gnorm/gnorm0;

        fprintf(1,'%d\t%d\t%1.3e\t%1.3e\t%1.2e\t%1.2e\t%1.2e\n', ...
        nit, nf, f1, f1/f0, stepNorm, gnorm, relGnorm)
        if (ierror == 3);
            xstar = x0;
            NF(1,nind) = NF(1,nind) + nit;
            NF(2,nind) = NF(2,nind) + nf;
            disp('ierror=3');
            return;
        end;
        
        %---------------------------------------------------------
        % convergence test
        %---------------------------------------------------------
        conv = (relGnorm <= currentRelGnormTol ...
                || gnorm <= currentAbsGnormTol ...
                || nit>=maxiter);
        

        if (conv);
            ierror = 0;
            xstar = x1;
            NF(1,nind) = NF(1,nind) + nit;
            NF(2,nind) = NF(2,nind) + nf;
            return;
        end
        
    end

end

%%%%%%%%%%%%%%%%
function [Hg, flag] = getHg_lbfgs(g,S,Y,hdiag)
% This function returns the approximate inverse Hessian multiplied by the gradient, H*g
% Input
%   S:    Memory matrix (n by k) , s{i}=x{i+1}-x{i}
%   Y:    Memory matrix (n by k) , df{i}=df{i+1}-df{i}
%   g:    gradient (n by 1)
%   hdiag value of initial Hessian diagonal elements (scalar)
% Output
%   Hg    the the approximate inverse Hessian multiplied by the gradient g
% Notice
% This funcion getHg_lbfgs is called by LBFGS_opt.m.
% Ref
%   Nocedal, J. (1980). "Updating Quasi-Newton Matrices with Limited Storage".
%   Wiki http://en.wikipedia.org/wiki/Limited-memory_BFGS
%   two loop recursion

% flag =-1 returns gradient descent step

    [n,k] = size(S);
    ro    = zeros(k,1);
    flag = 0;
    for i = 1:k
        ro(i) = 1/(Y(:,i)'*S(:,i));
        if ro(i)<0
            flag = -1;
        end
    end
    
    if flag==-1
        fprintf('Non-positive Hessian\n')
        Hg = g;
    else

        q = zeros(n,k+1);
    %     r = zeros(n,1);
        alpha =zeros(k,1);
        beta =zeros(k,1);

        % step 1
        q(:,k+1) = g;

        % first loop
        for i = k:-1:1
            alpha(i) = ro(i)*S(:,i)'*q(:,i+1);
            q(:,i) = q(:,i+1)-alpha(i)*Y(:,i);
        end

        % Multiply by Initial Hessian
        r = hdiag*q(:,1);

        % second loop
        for i = 1:k
            beta(i) = ro(i)*Y(:,i)'*r;
            r = r + S(:,i)*(alpha(i)-beta(i));
        end
        % 
        Hg=r;
    end
end % end of getHg_lbfgs