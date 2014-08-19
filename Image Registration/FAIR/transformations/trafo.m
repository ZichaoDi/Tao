% =======================================================================================
% function varargout = trafo(varargin)
% (c) Jan Modersitzki 2009/03/24, see FAIR.2 and FAIRcopyright.m.
% 
% Wrapper for transformation model.
%
%  [y,dy] = trafo(w,X);
%  computes y = Q(X)*w, dy = Q(X), 
%  where Q(X) is a sample of the basis functions on the grid X,
%  Q is allocated persistent in the method to be used, X is used for proper size
%
%  initialize: trafo('reset','trafo','rigid2D');
%  usage:      w0 = trafo('w0'); % returns parameterization of identity
%              [yc,dy] = trafo(wc,xc;
%
% administration, parameter setting:
%  for resetting, intitializing, setting, updating, clearing, displaying,
%  see options.m
%
% starting value:
%  w0 = trafo('w0');
%  returns the parameterization for the identity transformation
%
% specific optn:
%  [y,dy] = trafo(w,X,specific{:}), uses specific for this call
%
% see also E4_US_trafo, transformations/contents.m, and BigTutorialTrafo.
% =======================================================================================

function varargout = trafo(varargin)

persistent OPTN 

if nargin == 0 && nargout == 0 && isempty(OPTN),
  help(mfilename);
  return;
end;

% handle options
[method,OPTN,task,stop] = options(OPTN,varargin{:});
if stop,
  varargout{1} = method;
  if any(strcmp(task,{'reset','set'})), trafo('w0'); end; % clear Q
  if nargout > 1, varargout{2} = OPTN;  end;
  return;
end

% return parameters for identity transformation
if strcmp(task, 'w0'),
  [~,w0] = feval(method,'w0',[],OPTN{:});
  varargout = {w0};
  return;
end;

% do the work
[method,optn] = options(OPTN,'set','doDerivative',(nargout>1),varargin{3:end});
doDerivative  = options(optn,'get','doDerivative');
w = varargin{1};
x = varargin{2};
if ~doDerivative,
  y  = feval(method,w,x,optn{:});
  varargout = {y,[]};
  return;
end;

[y,dy] = feval(method,w,x,optn{:});
if ~strcmp(options(optn,'get','debug'),'on'), 
  varargout = {y,dy};
  return; 
end;

% DEBUGGING mode (used for derivative checks),
% transform sparse presentation of dy to full
if iscell(dy)
  
  if length(dy) == 1,
    p = size(dy{1});
    dim = numel(w)/p(2);
    dy = kron(speye(dim),dy{1});
  else
    error('nyi')
  end;
end;
varargout = {y,dy};

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