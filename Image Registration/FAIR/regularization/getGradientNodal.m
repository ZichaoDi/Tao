% ==================================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% function B = getGradientNodal(omega,m);
% 
% Builds the gradient operator B for a domain defined by omega
% and a nodal grid discretization defined by m:
%
% 
% EXAMPLE FOR 2D:
%
% o--x--o--x--o--x--o    
% |     |     |     |   
% +     +     +     +  
% |     |     |     |    
% o--x--o--x--o--x--o    
% |     |     |     |       WHERE:
% +     +     +     +       'o'  - uc   	(nodal)
% |     |     |     |       'x'  - d_1 uc  	(staggered-2)
% o--x--o--x--o--x--o       '+'  - d_2 uc  	(staggered-1)
%
% Note: the matrices are simple, it is size that matters. 
% 
% see also hyperElastic and BigTutorialRegularizer
% ==================================================================================

function B = getGradientNodal(omega,m)
if nargin == 0,
  help(mfilename);
  runMinimalExample;
  return;
end;

h  	=  (omega(2:2:end)-omega(1:2:end))./m;
dim	=  size(omega,2)/2;
id 	=  @(i) speye(m(i)+1);										 % (m(i)+1) identity matrix
dx	=  @(i) spdiags(ones(m(i),1)*[-1,1],[0,1],m(i),m(i)+1)/h(i); % short differnce operator
switch dim 
    case 2
      B = [kron(id(2),dx(1));kron(dx(2),id(1))];
      z = sparse(size(B,1),size(B,2));
      B = sparse([B z; z B]);
    case 3
        B = [
        kron(id(3),kron(id(2),dx(1)))
        kron(id(3),kron(dx(2),id(1)))
        kron(dx(3),kron(id(2),id(1)))
        ];
        z = sparse(size(B,1),size(B,2));
        B = sparse([B z z; z B z; z z B]);
    otherwise
        error('Dimension must be either 2 or 3.')
end

function runMinimalExample
% 2D
  omega = [0,3,0,5]; m = [10,13];
  B2 = getGradientNodal(omega,m);
  figure(1); clf;
  subplot(1,2,1)
  spy(B2); 
  title(sprintf('%s matrix (2D)',mfilename));


% 3D
  omega = [0,1,0,2,0,3]; m = [4,5,6];
  B2 = getGradientNodal(omega,m);
  subplot(1,2,2);
  spy(B2);   
  title(sprintf('%s matrix (3D)',mfilename));
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
