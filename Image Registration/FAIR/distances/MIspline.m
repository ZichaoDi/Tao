% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [Dc,rho,dD,drho,d2psi] = MIspline(Tc,Rc,omega,m,varargin);
%
% Mutual Information based distance measure  (using a matlab implementation for the joint 
% density, see MIspline for a C version)  using a general Gauss-Newton type framework.
%
% The MI distance measure is implemented using a Parzen-window estimator
% for the joint density rho. Given rho and drho, MI is computed as
% Dc    = psi(rho) = rho' * log(rho + tol) + ...
% dpsi  = log(rho + tol) + rho./(rho + tol) + ...
% d2psi = (rho + 2*tol)./(rho + tol)^2 + ...
% where the "..." indicates two more summands of similar type
%
% Example: run('MIcc')
%   setupHandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Dc = MIspline(Tc,Rc,omega,m);
%
% Input: 
%  Tc, Rc			template and reference
%  omega, m 	represents domain and discretization
%  varargin		optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc           MI(Tc,Rc)
%  rho          joint density estimator
%  dD           dpsi*drrho
%  drho         derivative of joint density estimator
%  d2psi        log2(length(Tc))*sqrt(nT*nR);
% 
% see also distances/contents

function [Dc,rho,dD,drho,d2psi] = MIspline(Tc,Rc,omega,m,varargin);

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

% initialize the output variables
Dc  = []; dD = []; rho  = []; drho = []; d2psi = [];

doDerivative = (nargout > 3);  % boolean for derivative computation

%setup default parameter for joint density estimator
tol  = 1e-6;  % tolerance for logarithm in entropy
minT = 0;     % (assumed) smallest value in Tc
maxT = 256;   % (assumed) largest value in Tc
nT   = 32;    % number of "bins" for Tc
minR = 0;     % (assumed) smallest value in Rc
maxR = 256;   % (assumed) largest value in Rc
nR   = 30;    % number of "bins" for Rc

% overwrite default parameter
for k=1:2:length(varargin), % overwrite defaults
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

widthT = (maxT-minT)/nT;
widthR = (maxR-minR)/nR;

% to ensure nothing is missed at the boundary, two artifical bins are added
minT  = minT-2*widthT;  maxT  = maxT+2*widthT;   nT = nT+4;
minR  = minR-2*widthR;  maxR  = maxR+2*widthR;   nR = nR+4;

% compute the Parzen-window estimator for the density and its derivative
% rho  is (nT*nR) x 1
% drho is (nT*nR) x length(Tc), large but sparse 
[rho,drho] = rhoSpline(Tc,Rc,minT,maxT,nT,minR,maxR,nR);

% reformat rho for efficient computation of marginal densities
% rhoT and rhoR and undo formating
% note, summation can be described via matrix muliplication:
% rhoT = ST*rho, rhoR = SR*rho, where
% ST   = kron( eye(nT,nT), ones(1,nR) ) and
% SR   = kron( ones(1,nT), eye(nR,nR) );
% this is used for the computation of the derivative 

rho  = reshape(rho,nT,nR);
rhoT = sum(rho,2);
rhoR = sum(rho,1)';
rho  = rho(:);

% compute MI
Dc = rhoT'*log(rhoT+tol) + rhoR'*log(rhoR+tol) - rho'*log(rho+tol);

if ~doDerivative, return; end;

% build the matrices ST, SR
ST    = sparse(kron(ones(1,nR),speye(nT,nT)));
SR    = sparse(kron(speye(nR,nR),ones(1,nT)));

dpsi  = (log(rhoT+tol)+rhoT./(rhoT+tol))'*ST ...
  +     (log(rhoR+tol)+rhoR./(rhoR+tol))'*SR ...
  -     (log(rho +tol)+rho ./(rho +tol))';

% apply chain rule for MI = psi(rho(Tc))
dD    = dpsi * drho;
d2psi = log2(length(Tc))*sqrt(nT*nR);

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