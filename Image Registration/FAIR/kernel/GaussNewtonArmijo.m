%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% OBSOLETE! please use GaussNewton instead
%==============================================================================

function [yc,His] = GaussNewtonArmijo(fctn,yc,varargin)

yc = []; His= [];
warning('OBSOLETE! please use GaussNewton instead')
return;

if nargin ==0,
  help(mfilename);
end;

% parameter initialization -----------------------------------------------
maxIter      = 10;              % maximum number of iterations
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
LSMaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
vecNorm      = @norm;           % norm to be used for dJ and dy    
solver       = [];              % linear solver 
yStop        = [];              % used for stopping in multi-level framework
Jstop        = [];              % 
Plots        = @(iter,para) []; % for plots;
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(yStop), yStop  = yc; end; % yStop used for stopping only
% -- end parameter setup   ----------------------------------------------

% some output
FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
  char(ones(1,10)*'-'),str,char(ones(1,60-length(str))*'-'));
FAIRmessage([mfilename '(JM 2008/12/09)']);
fprintf('[ maxIter=%s / tolJ=%s / tolY=%s / tolG=%s / length(yc)=%d ]\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolY),num2str(tolG),length(yc));

% -- initialize  ---------------------------------------------------------                                        
STOP = zeros(5,1);

if isempty(Jstop),
  % evaluate objective function for stopping values and plots
  [Jstop,para] = fctn(yStop); Jstop = abs(Jstop) + (Jstop == 0);
  Plots('stop',para);
end;

% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(yc); 
Plots('start',para);
iter = 0; yOld = 0*yc; Jold = Jc; y0 = yc;

hisStr    = {'iter','J','Jold-J','|\nabla J|','|dy|','LS'};
his        = zeros(maxIter+2,6);
his(1,1:3) = [-1,Jstop,Jstop-Jc];
his(2,:)   = [0,Jc,Jstop-Jc,vecNorm(dJ),vecNorm(yc-yStop),0];

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
  STOP(1) = (iter>0) && abs(Jold-Jc)  <= tolJ*(1+abs(Jstop));
  STOP(2) = (iter>0) && (norm(yc-yOld) <= tolY*(1+norm(y0)));
  STOP(3) = norm(dJ)      <= tolG*(1+abs(Jstop));
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
    warning('no descent direction, switch to -dy!')
    dy      = -dy;
  end;

  % perform Armijo line-search
  [t,yt,LSiter] = Armijo(fctn,yc,dy,Jc,dJ,...
        'LSMaxIter',LSMaxIter,'LSreduction',LSreduction);
  if (t == 0), break; end; % break if line-search fails
  
  % save old values and update
  yOld = yc; Jold = Jc; yc = yt;    
  [Jc,para,dJ,H] = fctn(yc); % evalute objective function

  % some output
  his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(yc-yOld),LSiter];
  dispHis(his(iter+1,:));
  para.normdY = vecNorm(yc - yOld);
  Plots(iter,para);
% pause
end;%while; % end of iteration loop
%-------------------------------------------------------------------------
Plots(iter,para);

% clean up
His.str = hisStr;
His.his = his(1:iter+2,:);
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-yOld|',norm(yc-yOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%==============================================================================

function dy = solveGN(rhs,H,solver);
maxIterCG = 10; tolCG = 1e-2;

if isempty(solver) 
  if isstruct(H),
    solver = H.regularizer;
  else % no regularizer initialized, assuming PIR
    dy = H\rhs;
    return;
  end;
end;

switch solver,
  case 'mfElastic',     
    dy = MGsolver(rhs,H);
  case 'mfCurvature',    
    Afun     = @(dy) mfAy(dy,H); 
    [dy,FLAG] = pcg(Afun,rhs,tolCG,maxIterCG);  
  case 'pcg',  
    L   = tril(H); % Symmetric Gauss Seidel Preconditioning, 
    D   = diag(H); % L is lower, D is diagonal, U = L'
    SGS = @(x) L\(D.*(L'\x)); 
    [dy,FLAG] = pcg(H,rhs,tolCG,maxIterCG,SGS);

  otherwise,
   
    if isnumeric(H),
      % if H is a matrix, solve the linear system using MATLAB's backslash
      dy = H\rhs;
    else
      error(solver)
    end;
end;



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
  third parties, or any other special, incidental or consequential damages
  arrising out of the use or inability to use the program, even if the
  possibility of such damages has been advised against. The entire risk
  as to the quality, the performace, and the fitness of the program for any
  particular purpose lies with the party using the code.
  =======================================================================================
  Any use of this code constitutes acceptance of the terms of the above statements
  =======================================================================================
%}