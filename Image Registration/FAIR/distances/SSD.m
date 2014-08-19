% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [Dc,rc,dD,dr,d2psi] = SSD(Tc,Rc,omega,m,varargin);
% (c) Jan Modersitzki 2009/03/24, see FAIR.2 and FAIRcopyright.m.
% Sum of Squared Differences based distance measure 
% for usage in a general Gauss-Newton type framework
% computes D(Tc,Rc) = hd*psi(r(Tc)), r = Tc-Rc, psi = 0.5*r'*r 
% and derivatives, dr = dT, d2psi= hd*I, hd = prod(omega./m)
%
%   setupHandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Dc = SSD(Tc,Rc,omega,m);
%
% Input: 
%  Tc, Rc     template and reference
%  omega, m   represents domain and discretization
%  varargin   optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc         SSD(Tc,Rc)
%  rc         Tc-Rc
%  dD         dpsi*dr
%  dr         dT
%  d2psi      hd = prod(omega./m)
% 
% see also distances/contents

function [Dc,rc,dD,dr,d2psi] = SSD(Tc,Rc,omega,m,varargin);
global DrawError
if nargin == 0,
  help(mfilename);
  setupHandData;
  xc = getCellCenteredGrid(omega,m);
  Tc = linearInter(dataT,omega,xc);
  Rc = linearInter(dataR,omega,xc);
  D0 = feval(mfilename,Tc,Rc,omega,m);
  fprintf('%s: distance = %s\n',mfilename,num2str(D0))
  return;
end;
Dc  = []; dD = []; rc  = []; dr = []; d2psi = [];
doDerivative = (nargout > 2);

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration
rc = Tc - Rc;              		% the residual
if(DrawError)
    figure(221);
    subplot(1,3,1)
    xc = getCellCenteredGrid(omega,m);
    surf(reshape(xc(1:m(1)*m(2)),m(1),m(2)),reshape(xc(1+m(1)*m(2):end),m(1),m(2)),reshape(Tc,m(1),m(2)));
    colormap HSV
    subplot(1,3,2)
    surf(reshape(xc(1:m(1)*m(2)),m(1),m(2)),reshape(xc(1+m(1)*m(2):end),m(1),m(2)),reshape(Rc,m(1),m(2)));
    colormap HSV
     subplot(1,3,3)
    surf(reshape(xc(1:m(1)*m(2)),m(1),m(2)),reshape(xc(1+m(1)*m(2):end),m(1),m(2)),reshape(rc,m(1),m(2)));
    colormap HSV
end
Dc = 0.5*hd * (rc'*rc);       	% the SSD
% Dc1=Dc;
% Dc=Dc/(norm(Rc)*norm(Tc));
if ~doDerivative, return; end;
dr = 1; 						% or speye(length(rc),length(rc));
dD = hd * rc'*dr;
d2psi = hd;
% dA=(rc'*rc+2*rc'*Rc+Rc'*Rc)^(-3/2)*Tc';
% dD1=dD./(norm(Rc)*norm(Tc))-(Dc1/norm(Rc))*dA;
% dD=dD1;
% d2psi = d2psi/(norm(Rc)*norm(Tc))*eye(m(1)^2,m(1)^2)+hd/norm(Rc)*rc*dA- ...
% (-3*(rc'*rc+2*rc'*Rc+Rc'*Rc)^(-5/2)*Tc*Tc'+(rc'*rc+2*rc'*Rc+Rc'*Rc)^(-3/2)*eye(m(1)^2,m(1)^2));

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