% =======================================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function varargout = distance(varargin)
%
% Main function for distance measures, possible calls:
%
% initialize: distance('reset','distance','SSD');
% use:        [Dc,rc,dD,drc,d2psi] = distance(Tc,Rc,omega,m);
% evaluates the distances between Tc and Rc, cell volume (pixel/voxel) 
% is based on omega and m;
% internals: residual rc and outer function psi, D(Tc) = psi(rc(Tc))
%
% administration, parameter setting:
%  for resetting, intitializing, setting, updating, clearing, displaying,
%  options - deals parameterization
%
% specific options:
%  [Dc,rc,dD,drc,d2psi] = distance(Tc,Rc,omega,m,specific{:});
%
% see also E9_Hands_MLIR_SSD_mbElas and BigTutorialDistance
% =======================================================================================

function varargout = distance(varargin);

persistent OPTN

if nargin == 0 && nargout == 0 && isempty(OPTN),
  help(mfilename);
  return;
end;

% handle options
[method,OPTN,task,stop] = options(OPTN,varargin{:});
if stop,
  varargout{1} = method;
  if nargout > 1, varargout{2} = OPTN;  end;
  return;
end

% do the work

Tc      = varargin{1};
Rc      = varargin{2};
omega   = varargin{3};
m       = varargin{4};
  
[method,optn]  = options(OPTN,'set','doDerivative',(nargout>3),varargin{5:end});
[Dc,rc,dD,dr,d2psi] =  feval(method,Tc,Rc,omega,m,optn{:});
  
varargout = {Dc,rc,dD,dr,d2psi};

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