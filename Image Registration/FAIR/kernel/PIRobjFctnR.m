%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function  [Jc,para,dJ,H] = PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)
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
%
% Output:
%  Jc       current function value J(wc)
%  para     struct {Tc=T(y(wc)), Rc, omega, m, yc=y(wc,xc), Jc}, for plots
%  dJ       gradient of J
%  H        approximation to Hessian of J
%
% see also E6_HNSP_PIR_GN
%==============================================================================

function [Jc,para,dJ,H] = PIRobjFctnR(T,Rc,omega,omegat,m,mt,beta,M,wRef,xc,xt,wc)

global Ro
% do the work ------------------------------------------------------------
doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute transformation, distance, and regularization and combine these
xc=reshape(xc,m(1)*m(2),2);
% overlap=find((xc(:,1)-5).^2+(xc(:,2)-5).^2<=9);
ind=find(Ro);
overlap=ind;
  overlap=1:length(Ro);
[yc,dy] = trafo(wc,xt(:),'doDerivative',doDerivative);
% figure(111);plotGrid(yc,omegat,mt);pause;hold on; plotGrid(xt,omegat,mt,'color','r');hold off;
% pause;

[Tc,dT] = inter(T,omegat,yc,'doDerivative',doDerivative);
[Dc,rc,dD,dr,d2psi] = distance(Tc(overlap),Rc(overlap),omegat,mt,'doDerivative',doDerivative);

dTem=dT*dy;
if(doDerivative )
    xcLap=[overlap,size(dT,2)/2+overlap];
    dTem=dTem(overlap,:);
end
%==== add regularization
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
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'omegat',omegat,'m',m,'mt',mt,'yc',yc,'Jc',Jc);

if ~doDerivative, return; end;

% multiply outer and inner derivatives, note: dy might be sparse
if isnumeric(dy);
  dD = dD*dTem;
  dJ = dD + dS;
  if nargout<4, return; end;

  % approximation to Hessian
  dr = dr*dTem;
  H  = dr'*d2psi*dr + M + beta*speye(length(wc));
elseif isstruct(dy),
  dD = dy.Qadjoint((dD*dT)')';
  dJ = dD + dS;
  if nargout<4, return; end;
  
  % approximation to Hessian
  H.solver   = 'PIR-cg';
  H.operator = @(wc) dy.Qadjoint(dT'*(dr'*(d2psi* dr*(dT*dy.Q(wc))))) ...
    + M*wc + beta*wc;
elseif iscell(dy),
  dD = kronIdy(dD*dT,dy);
  dJ = dD + dS;
  if nargout<4, return; end;

  % approximation to Hessian
  dr = kronIdy(dr*dT,dy);
  H   = dr'*d2psi*dr + M + beta*speye(length(wc));
else
  error('nyi')
end;
% EIG=eig(H);
% figure(111);plot(1:length(EIG),EIG,'r.-');

function dd = kronIdy(Df,Dy)
% It is assumed that Q = I_d \otimes Dy{1},
%
% [df1,...dfd]*| Q       | = [df1*Q,...,dfd*Q]
%              |  \      |
%              |        Q|

m   = size(Df,2);
n   = size(Dy{1},1);
switch m/n,
  case 2, dd = [Df(:,1:n)*Dy{:},Df(:,n+1:end)*Dy{:}];
  case 3, dd = [...
      Df(:,1:n)*Dy{:},...
      Df(:,n+1:2*n)*Dy{:},...
      Df(:,2*n+1:end)*Dy{:}];
  otherwise,
    error('nyi')
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