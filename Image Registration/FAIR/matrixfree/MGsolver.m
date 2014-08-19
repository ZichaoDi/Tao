% (c) Jan Modersitzki 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function Y = MGsolver(rhs,H);
%
% prepares for multigrid solver for H*Y=rhs by initializing smoothing operator and such

function Y = MGsolver(rhs,H)

H.MGlevel      = log2(H.m(1))+1;%%(~strcmp(H.d2S.regularizer,'mfCurvature'));
H.MGcycle      = 1;
H.MGomega      = 0.5;   
H.MGsmoother   = 'mfJacobi';
H.MGpresmooth  = 3;
H.MGpostsmooth = 1;

d2D         = H.d2D;
H.d2D = [];


if isfield(d2D,'M'),
  M     = d2D.M;
elseif isfield(d2D,'dr'),
  M           = diag(d2D.dr'*d2D.d2psi*d2D.dr);
  M           = d2D.P(full(M));
else
  M = 0;
end;
% bring M into right format
if all(size(M)==1),
    M = M * ones(length(rhs),1);
elseif (any(size(M)) == 1)
    M = full(M(:));
else
    error('M needs to be diagonal and stored as a vector!');
end
H.d2D.M = M;
    
u0 = zeros(size(rhs));
Y  = mfvcycle(H,u0,rhs,1e-12,H.MGlevel,(max(H.m)>32));

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