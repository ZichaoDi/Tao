%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [yc,His] = GaussNewton(fctn,yc,varargin)
%
% Gauss-Newton scheme with variable line search for minimizing J = fctn(yc)
%
% Input:
%   fctn        function handle
%   yc          starting guess
%   varargin    optional parameter, see below
%
% Output:
%   yc          numerical optimizer (current iterate)
%   his         iteration history
%
%==============================================================================

function [yc,His] = GaussNewton(fctn,yc,varargin)

global H yOld
set(gca,'NextPlot','replacechildren');

if nargin ==0, % help and minimal example
    help(mfilename);  E10_2Ddisc2C_hyperElastic;  yc = 'example finished'; return;
end;

% parameter initialization -----------------------------------------------
lineSearch   = @Armijo;         % default line search
maxIter      = 10000;              % maximum number of iterations
tolJ         = 1e-4;% 1e-3;%           % for stopping, objective function
tolY         = 1e-4;%1e-2;%            %   - " -     , current value
tolG         = 1e-2;%1;%            %   - " -     , norm of gradient
LSMaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
vecNorm      = @norm;           % norm to be used for dJ and dy
solver       = [];              % linear solver
yStop        = [];              % used for stopping in multi-level framework
J0        = [];              %
Plots        = @(iter,para) []; % for plots;
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if ~isa(Plots,'function_handle') && (Plots == 0 || strcmp(Plots,'off')),
    Plots        = @(iter,para) []; % for plots;
end;
%yStop=[];
if isempty(yStop), yStop  = trafo('w0'); % yc
end; % yStop used for stopping only
% -- end parameter setup   ----------------------------------------------

% some output
% FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
%   char(ones(1,10)*'-'),str,char(ones(1,60-length(str))*'-'));
FAIRmessage([mfilename '(JM 2011/01/03)']);
fprintf('[ maxIter=%s / tolJ=%s / tolY=%s / tolG=%s / length(yc)=%d ]\n',...
    num2str(maxIter),num2str(tolJ),num2str(tolY),num2str(tolG),length(yc));

% -- initialize  ---------------------------------------------------------
STOP = zeros(5,1);
if isempty(J0),
    % evaluate objective function for stopping values and plots
    [J0,para] = fctn(yStop); J0 = abs(J0) + (J0 == 0);
    Plots('stop',para);
end;
% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(yc);
% [Jc1,para1,dJ1,H1] = fctn(trafo('w0'));
Plots('start',para);
iter = 0; yOld = 0*yc; Jold = Jc; y0 = yc;

hisStr    = {'iter','J','Jold-J','|\nabla J|','|dy|','LS'};
his        = zeros(maxIter+2,6);
his(1,1:3) = [-1,J0,J0-Jc];
his(2,:)   = [0,Jc,J0-Jc,vecNorm(dJ),vecNorm(yc-yStop),0];

% some output
fprintf('%4s %-12s %-12s %-12s %-12s %4s\n%s\n',...
    hisStr{:},char(ones(1,64)*'-'));
dispHis = @(var) ...
    fprintf('%4d %-12.4e %-12.3e %-12.3e %-12.3e %4d\n',var);
dispHis(his(1,:));
% -- end initialization   ------------------------------------------------

%-- start the iteration --------------------------------------------------
while 1,
    % check stopping rules
    STOP(1) = (iter>0) && abs(Jold-Jc)  <= tolJ*(1+abs(J0));
    STOP(2) = (iter>0) && (norm(yc-yOld) <= tolY*(1+norm(y0)));
    STOP(3) = norm(dJ)      <= tolG*(1+abs(J0));
    STOP(4) = norm(dJ)      <= 1e6*eps;
    STOP(5) = (iter >= maxIter);
    if all(STOP(1:3)) || any(STOP(4:5)), break;  end;
    
    iter = iter + 1;
    % solve the Gauss-Newton System
    dy = solveGN(-dJ',H,solver);
    % check descent direction
    % note: descent is not granted if using an iterative solver
    descent =   dJ * dy;
    if descent > 0,
        warning('no descent direction, switch to -dy!') %#ok<WNTAG>
        dy      = -dy;
    end;
    
    % perform Armijo line-search
    [t,yt,LSiter] = lineSearch(fctn,yc,dy,Jc,dJ,...
        'LSMaxIter',LSMaxIter,'LSreduction',LSreduction,'m',para.m,'omega',para.omega);
    if (t == 0),
        break;
    end; % break if line-search fails
    % save old values and update
    yOld = yc; Jold = Jc; yc = yt;
    [Jc,para,dJ,H] = fctn(yc); % evalute objective function
    % some output
    his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(yc-yOld),LSiter];
    dispHis(his(iter+1,:));
    para.normdY = vecNorm(yc - yOld);
    Plots(iter,para);
    %
    %     Tomo1(iter) = getframe(gcf);
    %      save Tomo1 Tomo1;
    %    pause
end;%while; % end of iteration loop
%-------------------------------------------------------------------------
Plots(iter,para);

% clean up
His.str = hisStr;
His.his = his(1:iter+2,:);
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
    '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(J0)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
    '|yc-yOld|',norm(yc-yOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
    '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',tolG*(1+abs(J0)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
    'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
    'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%==============================================================================

function dy = solveGN(rhs,H,solver)
maxIterCG = 500; tolCG = 1e-1;

if isempty(solver)
    if isstruct(H),
        if isfield(H,'solver'),
            solver = H.solver;
        elseif isfield(H,'d2S') && isfield(H.d2S,'solver'),
            solver = H.d2S.solver;
        else
            error('solver has not been defined')
        end;
    else % no regularizer initialized, assuming PIR
        dy = H\rhs;
        return;
    end;
end;

if isa(solver, 'function_handle')
    dy = feval(solver, rhs, H, maxIterCG, tolCG);
    return
end
switch solver
    % matrix based
    % ------------
    case 'pcg'
        L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
        D   = diag(H); % L is lower, D is diagonal, U = L'
        SGS = @(x) L\(D.*(L'\x));
        [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,SGS);
    case 'ichol'
        L1   = cholinc(sparse(H),'0'); % Symmetric Gauss Seidel Preconditioning,
        [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,L1,L1');
    case 'jacobi-pcg'
        D   = diag(H); % D is diagonal
        PC = @(x) D.\x;
        [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,PC);
    case 'cg'
        [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
    case 'PIR-cg'
        if isstruct(H)
            [dy,flag,relres,iter] = pcg(H.operator,rhs,tolCG,maxIterCG);
        else
            [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
        end
        % matrix free
        % ------------
        % In all matrix free solvers the approximates Hessian is supplied by function handle.
        % Note, that the functionals have the common form
        %
        %   Jc = D + S (+ P) ==> d_2 Jc = d2D + d2S (+d2P)
        %
        % regularization and penalization (S,P) supply a struct stored in H.d2S or H.d2P
        % which already has a function handle describing the operator.
        % The  distance term (d2D) needs more care and the approximation to the Hessian
        % depends on the solver as well. Therefore the objective function does NOT supply
        % the function handle, but all variables needed to built such.
        % For multi-grid (available for elastic and curvature), only the diagonal of the distance
        % term is used
        % whereas in conjugate gradient (CG, available for elastic, curvature and hyperelastic) methods we can
        % model also the off-diagonals d2D(x) = (P'*dr'*d2psi*dr*P) * x
        %
        % The default solver for a chosen regularizer is parameterized by d2S.solver.
        % However, one can overload this setting using
        %    >> MLIR(..., 'solverNPIR','myFavoriteSolver', ... );
    case {'MG-elastic'}
        dy = MGsolver(rhs,H);
    case {'CG-elastic','CG-curvature'}
        M         = @(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn     = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m);
        [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG);
    case {'PCG-elastic','PCG-curvature'}
        M         = @(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn     = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m);
        % preconditioner
        Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
        D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.omega,H.m);
        PC = @(x) D.\x; % Jacobi preconditioner
        [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);
    case {'CG-hyperElastic'}
        M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
        [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG);
    case {'PCG-hyperElastic'}
        M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
        Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
        D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.d2S.yc);
        PC = @(x) D.\x; % Jacobi preconditioner
        [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);
    otherwise
        if isnumeric(H)
            % if H is a matrix, solve the linear system using MATLAB's backslash
            dy = H\rhs;
        else
            error(solver)
        end
end

if exist('flag','var')
    switch flag
        case 1
            fprintf('pcg iterated %d times but converged only to relres %e instead of %e\n',...
                iter,relres,tolCG);
        case 2
            fprintf('Preconditioner M was ill-conditioned.\n');
        case 3
            fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
        case 4
            fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n');
        otherwise
            fprintf('pcg success! %d iterations / relres= %1.2e / tolCG= %1.2e\n',iter,relres,tolCG);
    end
end

%{
  =======================================================================================
  FAIR: Flexible Algorithms for Image Registration, Version 2011
  Copyright (c): Jan Modersitzki
  Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
  Email: jan.modersitzki@mic.uni-luebeck.de
  URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
  =======================================================================================
  No part of this code may be reproduced, stored in a retrieval system,
  translated, transcribed, transmitted, or distributed in any form
  or by any means, means, manual, electric, electronic, electro-magnetic,
  mechanical, chemical, optical, photocopying, recording, or otherwise,
  without the prior explicit written permission of the authors or their
  designated proxies. In no event shall the above copyright notice be
  removed or altered in any way.

  This code is provided "as is", without any warranty of any kind, either
  expressed or implied, including but not limited to, any implied warranty
  of merchantibility or fitness for any purpose. In no event will any party
  who distributed the code be liable for damages or for any claim(s) by
  any other party, including but not limited to, any lost profits, lost
  monies, lost data or data rendered inaccurate, losses sustained by
  third parties, or any othe r special, incidental or consequential damages
  arrising out of the use or inability to use the program, even if the
  possibility of such damages has been advised against. The entire risk
  as to the quality, the performace, and the fitness of the program for any
  particular purpose lies with the party using the code.
  =======================================================================================
  Any use of this code constitutes acceptance of the terms of the above statements
  =======================================================================================
%}