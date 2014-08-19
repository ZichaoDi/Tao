% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% =======================================================================================
% function [Dc,rc,dD,dr,d2psi] = SSDweighted(T,R,omega,m,varargin);
% (c) Jan Modersitzki 2009/03/24, see FAIR.2 and FAIRcopyright.m.
% Sum of Squared Differences based distance measure 
% using a general Gauss-Newton type framework
% computes D(T,R) = hd*psi(r(T)), r = T-R, psi = 0.5*r'*r 
% and derivatives, dr = dT, d2psi= hd*I, hd = prod(omega./m)
%
%   setupHandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Wc = (Tc>100);
%   Dc = SSDweighted(Tc,Rc,omega,m,'W',Wc);
%
% Input: 
%  T, R			template and reference
%  omega, m 	represents domain and discretization
%  varargin		optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc           SSD(diag(W)*T,diag(W)*R)
%  rc           diag(W)*(T-R)
%  dD           dpsi*dr
%  dr           dT
%  d2psi        hd = prod(omega./m)
% 
% see also distances/contents

function [Dc,rc,dD,dr,d2psi] = SSDweighted(T,R,omega,m,varargin);

if nargin == 0,
  help(mfilename);
  setupHandData;
  xc = getCellCenteredGrid(omega,m);
  Tc = linearInter(dataT,omega,xc);
  Rc = linearInter(dataR,omega,xc);
  Wc = spdiags((Rc>100),0,length(Rc),length(Rc));
  FAIRfigure(1); clf; colormap(gray(256))
  subplot(2,2,1); viewImage2D(Tc,omega,m,'title','template');
  subplot(2,2,2); viewImage2D(Rc,omega,m,'title','reference');
  subplot(2,2,3); viewImage2Dsc(diag(Wc),omega,m,'title','weight');
  subplot(2,2,4); viewImage2D(128+Wc*(Tc-Rc)/2,omega,m,'title','distance');
  
  D0 = feval(mfilename,Tc,Rc,omega,m,'W',Wc);
  fprintf('%s: distance = %s\n',mfilename,num2str(D0))
  return;
end;

Dc  = []; dD = []; rc  = []; dr = []; d2psi = [];
dim = length(omega) / 2;
n   = prod(m);
W   = speye(length(T),length(T));

doDerivative = (nargout > 3);

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if size(W,1)~=n
    % In case of multi-level: Bring weighting mask to current resolution
    scalef = (size(W,1)/n)^(1/dim);
    W = spdiags(linearInter(reshape(full(diag(W)), round(m*scalef)), ...
        omega, getCellCenteredGrid(omega, m)), 0, n, n);
end

hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration
rc = T-R;                  		% the residual
Dc = 0.5*hd * rc'*W*rc;       	% the SSD
if ~doDerivative, return; end;
dr = 1; 						% or speye(length(rc),length(rc)); 
dD = hd * rc'*W*dr; 
d2psi = hd*W;

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