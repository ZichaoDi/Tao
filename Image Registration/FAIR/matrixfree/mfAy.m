% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [y,D] = mfAy(x,para)
%
% MatrixFree y =(M + alpha*h*B'*B)*y for elastic, diffusion, curvature
% where M is supposed to be diagonal and B is given via mfBu.m;
% the matrices M and B are represented by the struct para
% returns also D = diag(alpha*h*B'*B) 

function [y,D] = mfAy(x,para)

omega = para.omega; 
m     = para.m; 

if isfield(para,'d2D') && isfield(para.d2D,'M')
  if isnumeric(para.d2D.M),
    Mx    = para.d2D.M.*x;
  elseif isa(para.d2D.M, 'function_handle')
      Mx = para.d2D.M(x);
  else
    error('M needs to be specified')
  end;
end;
y = Mx + para.d2S.d2S(x,omega,m) ;

if nargout<2, return; end;
if isfield(para,'d2D') && isfield(para.d2D,'M'),
    M = para.d2D.M;
else
    error('M needs to be specified');
end
% only for preconditioning of Jacobi smoothing in multigrid
S = para.d2S.diag(omega,m);
n = numel(S);
D = spdiags(M+S,0,n,n);

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