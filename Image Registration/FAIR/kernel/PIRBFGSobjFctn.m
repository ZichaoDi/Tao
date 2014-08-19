% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% function  [Jc,para,dJ,H] = PIRBFGSobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)
%
% Objective Function for Parametric Image Registration
%
% computes J(wC) = D(T(Y(wc)),R) + S(wc), where
%
% yc       = y(wc,xc) = trafo(wc,xc)
% Tc       = T(yc) = inter(T,omega,yc), 
% D(Tc,Rc) = distance(Tc,Rc,omega,m)
% S(wc)    = 0.5*(wc-wRef)'*M*(wc-wRef);
%                   
% Input:
%   T       coefficients for template image
%   Rc      sampled reference image
%   omega   spatial domain
%   m       number of discretization points
%   beta    adding beta*I to approximation of Hessian 
%   xc      discretization of Omega
%   wc      current parameters
% Output:
%  Jc       current function value J(wc)
%  para     struct {Tc=T(y(wc)), Rc, omega, m, yc=y(wc,xc), Jc}, for plots
%  dJ       gradient of J
%  H        approximation to Hessian of J
%
%==============================================================================

function [Jc,para,dJ,H] = PIRBFGSobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)

if nargin == 0,
  help(mfilename)
  return
elseif ~exist('wc','var') || isempty(wc),
  % if wc is not an input argument, reports status
  if nargout == 1, Jc = 'PIR';  return; end;
  % report current settings
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  wc      = trafo('w0');
  fprintf('Parametric Image Registration:');
  fprintf('    J(wc)=D(T(y(wc)),R) + (wc-wRef)''*M*(wc-wRef) != min\n');
  fprintf('  %20s : %s\n','m',dimstr(m));
  fprintf('  %20s : %s\n','omega',dimstr(omega));
  fprintf('  %20s : %s\n','INTERPOLATION',inter);
  fprintf('  %20s : %s\n','DISTANCE',distance);
  fprintf('  %20s : %s\n','TRAFO',trafo);
  fprintf('  %20s : %s\n','length(wc)',num2str(length(wc)));
  fprintf('  %20s : %s\n','M',sprintf('is %d-by-%d',size(M)));
  fprintf('  %20s : %s\n','beta',num2str(beta)); 
  Jc = wc; % return starting guess
  return;
end;

% do the work ------------------------------------------------------------
doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute transformation, distance, and regularization and combine these
[yc,dy]    = trafo(wc,center(xc,m),'doDerivative',doDerivative);
[Tc,dT]    = inter(T,omega,yc,'doDerivative',doDerivative);
[Dc,rc,dD] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

% add regularization
if isempty(M) || (all(size(M)==1) && (M==0)),
  Sc = 0;
  dS = 0;
  M  = 0;
else
  dS = (wc-wRef)'*M;
  Sc = 0.5*dS*(wc-wRef);
end;

Jc = Dc + Sc;                           

% collect variables for plots
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc);

if ~doDerivative, return; end;
dD = dD * dT;                    % multiply outer and inner derivatives
if size(dD,2) == size(dy,1),     % generic case, dy comes complete
  dJ = dD*dy + dS;     
else                             % tricky case,  dy comes sparse
  n  = size(dy{1},1);
  dJ = reshape(dy{1}'*reshape(dD*dT,n,[]),1,[]) + dS;
%   dr = [dr(:,1:n)*dy{1},dr(:,n+1:2*n)*dy{1},dr(:,2*n+1:3*n)*dy{1}];
end;
if nargout<4, return; end;

% approximation to Hessian
if M == 0,
  H = 1;
else
  H = M;
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