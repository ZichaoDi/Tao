%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%==============================================================================
%
% function [Jc,para,dJ,H] = NPIRobjFctn(T,Rc,omega,m,yRef,yc)
%
% This is a template for an objective function to be supplied to an optimization scheme.
% 
% computes J(yc) = D(T(P*yc),Rc) + S(yc-yRef), where
% 
%  Tc       = T(yc) = inter(T,omega,P*yc), P cell-centered grid interpolation of yc
%  D(Tc,Rc) = distance(Tc,Rc,omega,m)
%  S(yc)    = regularizer(yc,omega,m) = 0.5*alpha*hd*norm(B*yc)^2 
%                    
% Modes
%   NPIRobjFctn     - displays the help
%   NPIRobjFctn([]) - reports the current setting of the distance and regularization
%   [Jc,para,dJ,H]  = NPIRobjFctn(T,Rc,omega,m,yRef,yc)
%                   - evaluates the objective function
%                   
% Input
%   T     - data for template image Tc = inter(T,omega,yc)
%   Rc    - reference image on grid Rc = inter(R,omega,xc)
%   omega - representation of computational domain
%   m     - discretization size
%   yRef  - reference for regularization, Sc = regularizer(yc-yRef,omega,m)
%   yc    - current grid
%   
% Output
%   Jc    - function value
%   para  - parameter for plots
%   dJ    - gradient
%   H     - GaussNewton style approximation to Hessian (either explicit or as operator)
% 
% see also NPIRFGSobjFctn (BFGS version), GaussNewton, E9_Hands_NPIRmb_GN
%==============================================================================

function [Jc,para,dJ,H] = NPIRobjFctn(T,Rc,omega,m,yRef,yc)

global regEnable
persistent P % cell-centered grid interpolation operator

if nargin == 0,
  help(mfilename)
  return
elseif ~exist('yc','var') || isempty(yc),
  if nargout == 1, Jc = 'NPIR';  return; end;
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  fprintf('Non-Parametric Image Registration\n'); 
  v = @(str) regularizer('get',str); % get regularization configuration like grid

  fprintf('  J(yc) = D(T(yc),R) + alpha*S(yc-yReg) != min\n');
  fprintf('  %20s : %s\n','INTERPOLATION',inter);
  fprintf('  %20s : %s\n','DISTANCE',distance);
  fprintf('  %20s : %s\n','REGULARIZER',regularizer);
  fprintf('  %20s : %s\n','alpha,mu,lambda',num2str([v('alpha'),v('mu'),v('lambda')]));
  fprintf('  %20s : %s\n','GRID',v('grid'));
  fprintf('  %20s : %s\n','MATRIX FREE',int2str(v('matrixFree')));
  fprintf('  %20s : %s\n','m',dimstr(m));
  fprintf('  %20s : %s\n','omega',dimstr(omega));
  return;
end;

% do the work ------------------------------------------------------------

% define interpolation for cell-centered grid
P = gridInterpolation(P,omega,m);

doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute interpolated image and derivative, formally: center(yc) = P*yc
[Tc,dT] = inter(T,omega,center(yc,m),'doDerivative',doDerivative);

% compute distance measure
[Dc,rc,dD,dres,d2psi] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

% compute regularizer
if regEnable,
    [Sc,dS,d2S] = regularizer(yc-yRef,omega,m,'doDerivative',doDerivative);
else
  Sc = 0;
  dS = 0;
  d2S  = 0;
end

% evaluate joint function and return if no derivatives need to be computed
Jc = Dc + Sc;

% store intermediates for outside visualization
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',center(yc,m),'Jc',Jc);

if ~doDerivative, return; end;

if isnumeric(P),
  % matrix based mode
  dr = dres*dT*P;
  dD = dD*dT*P;
  dJ = dD + dS;
  H  = dr'*d2psi*dr + d2S;
else

  % derivatives rather explicit
  dr = dres*dT;
  dD = dD*dT;  
  dJ = P(dD')' + dS;
  
  % approximation to d2D in matrix free mode
  % d2D   = P'*dr'*d2psi*dr*P 
  % P and P' are operators matrix free 
  H.omega     = omega;
  H.m         = m;
  H.d2D.how   = 'P''*dr''*d2psi*dr*P';
  H.d2D.P     = P;
  H.d2D.dr    = dr;
  H.d2D.d2psi = d2psi;

  H.d2S = d2S;  
end;

function P = gridInterpolation(P,omega,m)
switch regularizer,
  case 'mbElastic'
    if size(P,1) ~= length(omega)/2*prod(m), % update P
      P = stg2center(m); 
    end;
  case 'mfElastic', P = @(yc) stg2center(yc,m);
  case {'mbCurvature','mbTPS'},  P = 1;         % centered grid - matrix based
  case {'mfCurvature','mfTPS'},  P = @(yc) yc;  % centered grid - matrix free
  case 'mbHyperElastic',
    if size(P,1) ~= length(omega)/2 * prod(m)
      P = nodal2center(m);
    end
  case 'mfHyperElastic',
    P = @(y) nodal2center(y,m);
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