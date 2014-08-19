% ==================================================================================
% (c) Jan Modersitzki 2010/12/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function B = getElasticMatrixNodal(omega,m,mu,lambda)
%
% Builds the elasticity matrix B for a domain defined by omega
% and a nodal grid discretization defined by m:
%
%     | a\nabla       0        0 | 
% B = |       0 a\nabla        0 |
%     |       0       0  a\nabla |
%     |           b \div         |,
%
% where \nabla and \div are weighted by a=\sqrt(mu) and b=\sqrt(\lambda+\mu),
% respectively; defaults: \mu = 1, \lambda = 0.
% Note: the matrices are simple, it is size that matters.
%
% see also regularizer and BigTutorialRegularizer.
% ===================================================================================
function B = getElasticMatrixNodal(omega,m,mu,lambda)
if nargin == 0,
  help(mfilename);
  runMinimalExample
  return;
end;

                                        % set defautls
if ~exist('mu','var'),     mu     = [];  end;
if ~exist('lambda','var'), lambda = [];  end;
if isempty(mu),     
  warning('mu has not yet been set! set mu=1');  
  mu = 1;  
end;
if isempty(lambda), 
  warning('lambda has not yet been set! lambda=0');  
  lambda = 0;  
end;

dim = length(omega)/2; 
h   = (omega(2:2:end)-omega(1:2:end))./m; % voxel size for integration
a   = sqrt(mu); 
b   = sqrt(mu+lambda);

% build the 2D elasticity operator
%
%      | a\nabla 0      |   |a D11  0     |
%  B = |                | = |a D12  0     |
%      | 0      a\nabla |   |0      a D21 |
%      |                |   |0      a D22 |
%      | b Div1 b Div2  |   |b D11  b D22 |

% Note: since the grid is nodal, dx(1) := D11 (= D21) and dx(2) := D12 (= D22)

dx = @(k) spdiags(ones(m(k),1)*[-1,1],0:1,m(k),m(k)+1)/h(k);
av = @(k) spdiags(ones(m(k),2)/2,0:1,m(k),m(k)+1);
  
switch dim,
  case 2,  
    D1   = kron(speye(m(2)+1),dx(1));
    D2   = kron(dx(2),speye(m(1)+1));
    A1   = kron( av(2) , speye(m(1)));
    A2   = kron( speye(m(2)),  av(1));
    [p1,p2] = size(D1);
    [p3,p4] = size(D2);
    
    B = [
      a*D1,sparse(p1,p2)
      a*D2,sparse(p3,p4)
      sparse(p1,p2),a*D1
      sparse(p3,p4),a*D2
      b*A1*D1,b*A2*D2
      ];
  case 3,
    D1   = kron(speye((m(3)+1)*(m(2)+1)),dx(1));
    D2   =kron(speye(m(3)+1),kron(dx(2),speye(m(1)+1)));
    D3   = kron(dx(3),speye((m(2)+1)*(m(1)+1)));
    A1   = kron( av(3), kron( av(2) , speye(m(1))));
    A2   = kron( av(3), kron( speye(m(2)),  av(1)));
    A3   = kron( speye(m(3)), kron( av(2),  av(1)));
    [p1,p2] = size(D1);
    [p3,p4] = size(D2);
    [p5,p6] = size(D3);
    
    B = [
      a*D1,sparse(p1,2*p2)
      a*D2,sparse(p3,2*p4)
      a*D3,sparse(p5,2*p6)
      sparse(p1,p2),a*D1,sparse(p1,p2)
      sparse(p3,p2),a*D2,sparse(p3,p2)
      sparse(p5,p2),a*D3,sparse(p5,p2)
      sparse(p1,2*p2),a*D1
      sparse(p3,2*p2),a*D2
      sparse(p5,2*p2),a*D3
      b*A1*D1,b*A2*D2,b*A3*D3
      ];
end;
%       sparse(p1,p2),a*D1;
%       sparse(p3,p4),a*D2;

function runMinimalExample
  omega = [0,1,0,2,0,3]; m = [4,5,6];
  B2 = getElasticMatrixNodal(omega(1:4),m(1:2),1,0);
  B3 = getElasticMatrixNodal(omega,m,1,0);
  figure(1); clf;
  subplot(1,2,1); spy(B2); title(sprintf('%s: B2',mfilename));
  subplot(1,2,2); spy(B3); title(sprintf('%s: B3',mfilename));
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

     