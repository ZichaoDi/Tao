% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Landmark Based Thin-Plate-Spline Transformation
%
% function [Y,Z] = evalTPS(LM,c,X);
%
% evaluates the Thin-Plate-Spline parameterized by LM and c 
% at points X and the landmarks Z of the reference image.
%
% See also E5_Hands_TPS for an example

function [Y,Z] = evalTPS(LM,c,X)

if nargin == 0,
  help(mfilename); return;
end;

dim = size(LM,2)/2;

% prepare radial basis function
switch dim,
  case 2,  rho = @(r) r.^2.*log(r+(r==0));
  case 3,  rho = @(r) r;
end;

% L is number of landmarks, reshape X and lallocate memory for 
% Y=TPS(X) and Z = TPS(r)
L = size(LM,1); K = dim+1:2*dim; Z = zeros(size(LM,1),dim);
X = reshape(X,[],dim); x = @(i) X(:,i); Y = zeros(size(X)); 

% compute the linear part
for i=1:dim,
  Y(:,i) = c(L+1,i)+X*c(L+2:end,i);
  Z(:,i) = c(L+1,i)+LM(:,K)*c(L+2:end,i);
end;

% shortcut for ||x-r_j|| in vectorized version
r = @(X,j) sqrt(sum((X-ones(size(X))*diag(LM(j,K))).^2,2));

% add non-linear part
for j=1:L,
  rj = r(X,j); Rj = r(LM(:,K),j);
  for i=1:dim,
    Y(:,i) = Y(:,i) + c(j,i)*rho(rj);   Z(:,i) = Z(:,i) + c(j,i)*rho(Rj);
  end;    
end;
Y = Y(:);Z = Z(:);

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