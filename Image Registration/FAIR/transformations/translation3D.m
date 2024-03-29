% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [y,dy] = translation3D(w,x,varargin)
%
% computes y = Q*w = and the derivative wrt. w.
% x = reshape(x,[],2); 
% Q = [1 0 0; 0 1 0;0 0 1];  dy = Q:
% if no arguments are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m and BigTutorialTrafo

function [y,dy] = translation3D(w,x,varargin)

% the persistent variable Q stores the matrix 
% Q(x) = kron( I_2 , [x(:,1),x(:,2),1] );
persistent Q

if nargin == 0, 
  runMinimalExample;
  return;
else
  y = mfilename('fullfile'); 
  dy = [0;0;0];         % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

if isempty(w) || (size(Q,1) ~= numel(x)),
  n = length(x)/3;
  Q = sparse(kron(speye(3),ones(n,1)));
  if nargout == 0, return; end;
end;
y  = x + Q*w;
dy = Q;

function runMinimalExample
help(mfilename);
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,8,0,9]; m = [8,7,6]; 
c = (omega(2:2:end)-omega(1:2:end))'/2;
w = c/10
x = getCellCenteredGrid(omega,m);
z = feval(mfilename,w(:),x);
figure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;  

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