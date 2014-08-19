% ==================================================================================
% (c) Jan Modersitzki, Fabian Gigengack,Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de
%
% function varargout = regularizer(varargin);
%
%  multipurpose regularizer for NPIR enabling matrix based and matrix free evaluation
% 
%  The regularization functional is either based on a linear differential operator B
%        S(Y) = 0.5*alpha*hd*|B*Y|^2
%
%  Or a general (non-linear) functional S(Y)
%
% Notes:
%  [Sc,dS,d2S] = regularizer(Y,omega,m)
%  alpha and yRef can be set using regularizer('set','alpha',1);
%  A = alpha*hd*B'*B is made persistent for eficiency
%  hd = prod(omega./m)
%  an option grid is introduce to indicate the appropriate discretization:
%  staggered for elastic and diffusive, 
%  cell-centered for curvature
%  nodal for TV, hyper elastic or diffusive EPI regularization
%
%
% call
%  [Sc,dS,d2S] = regularizer(Y,omega,m);
%  computes S(Y), derivative dS and Hessian
%
% for resetting, intitializing, setting, updating, clearing, displaying
% see options.m
%
% see also BigTutorialRegularizer for examples.
%=======================================================================================

function varargout = regularizer(varargin)

persistent OPTN A

if nargin == 0 && nargout == 0 && isempty(OPTN),
  help(mfilename);
  return;
end;

% --------------------------------------------------------------------------------------
% handle options
[method,OPTN,task,stop] = options(OPTN,varargin{:});

% check the grid according to the regularizer to be used
if strcmp(task,'set') || strcmp(task,'reset'),
  matrixFree = strcmp(method(1:2),'mf');
  fctn = method(3:end); fctn(1) = lower(fctn(1));
  A = [];
  switch fctn,
    case {'TV'},
      fprintf('tbi, TV is non-linear and a matrix based form cannot be provided\n');
      error('use option ''TV'' instead')
    case {'TV','hyperElastic','diffusionEPI'}
      grid = 'nodal';
    case {'elastic'},
      grid = 'staggered';
    case {'curvature','TPS','tPS'},
      grid = 'cell-centered';
    otherwise,
      fprintf('\n\n---> %s <---\n\n',method)
      error('nyi')
  end;
  [dummy,OPTN] = options(OPTN,'set','grid',grid,'matrixFree',matrixFree);
end;

% return, if no further work has to be handled
if stop,
  varargout{1} = method;
  if nargout > 1, varargout{2} = OPTN;  varargout{3} = A;  end;
  return;
end
% --------------------------------------------------------------------------------------
% do the work

% extract regularization parameters
alpha       = options(OPTN,'get','alpha');
matrixFree  = options(OPTN,'get','matrixFree');
fctn = method(3:end); fctn(1) = lower(fctn(1));

% extract variables
yc     = varargin{1};
omega  = varargin{2};
m      = varargin{3};
doDerivative = (nargout>1);
varargin = varargin(4:end);

for k=1:2:length(varargin), % overwrites defaults
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

[Sc,dS,d2S] = feval(fctn,yc,omega,m,'alpha',alpha,'matrixFree',matrixFree,...
    'doDerivative',doDerivative,OPTN{:},varargin{:});
varargout = {Sc,dS,d2S};
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
