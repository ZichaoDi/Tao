% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [Dc,rc,dD,dr,d2psi] = NCC(Tc,Rc,omega,m,varargin);
% 
% Normalized Cross Correlation based distance measure 
% using a general Gauss-Newton type framework
% computes D(Tc,Rc) =  1 - (Rc'*Tc(Y))^2/T2, Rc = Rc/norm(Rc); T2 = Tc(Y)'*Tc(Y); 
% or r = Tc, psi = 1-(r'*Rc)^2/(r'*r)
% and derivatives, dr = dT, d2psi = 2/(r'*r)
%
% Example: run('MIcc')
%   setupHandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Dc = NCC(Tc,Rc,omega,m);
%
% Input: 
%  Tc, Rc			template and reference
%  omega, m 	represents domain and discretization
%  varargin		optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc           NCC(Tc,Rc)
%  rc           Tc-Rc
%  dD           dpsi*dr
%  dr           dT
%  d2psi        2/(Tc'*Tc);
%
% see also distances/contents

function [Dc,rc,dD,dr,d2psi] = NCC(Tc,Rc,omega,m,varargin);

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
doDerivative = (nargout > 3);

for k=1:1:length(varargin)/2, % overwrite default parameter
  eval([varargin{2*k-1},'=varargin{',int2str(2*k),'};']);
end;

%hd  = prod(omega./m); not needed here
rc = Tc; Rc = Rc/norm(Rc); T2 = Tc'*Tc;
Dc = 1-(Rc'*Tc)^2/T2;
if ~doDerivative, return; end;
dr = 1; 
dD = (-2*(Rc'*Tc)/T2*Rc'+2*((Rc'*Tc)/T2)^2*Tc')*dr;
d2psi  = 2/T2;

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