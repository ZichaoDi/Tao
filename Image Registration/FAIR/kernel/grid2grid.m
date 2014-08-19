%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%==============================================================================
%
% function [yc,y1,y2,y3] = grid2grid(yc,m,in,out)
%
% transfers input grid (c,s,n) to output grid (c,s,n) using inter- and extrapolation
% central operation are redce/expand, action each dimension individually
%
% Input:
%   yc        input grid
%   m         number of discretization points
%   in        type of input grid
%   out       type of output grid
% Output:
%   yc        output grid
%   y1,y2,y3  components of yc
%==============================================================================

function [yc,y1,y2,y3] = grid2grid(yc,m,in,out)

if nargin == 0, % help and minimal example
  help(mfilename);
  omega = [0,6,0,4]; m = [6 4];
  yC  = reshape(getCellCenteredGrid(omega,m),[],2);
  yN  = reshape(grid2grid(yC,m,'centered','nodal'),[],2);

  figure(1); clf; 
  plotGrid(yN,omega,m,'color','g'); hold on;  
  plot(yN(:,1),yN(:,2),'bs',yC(:,1),yC(:,2),'b*');
  title(mfilename); axis image
  return;
end;

dim = length(m);  y{1} = []; y{3} = []; y{3} = [];
if ~exist('out','var'), out = in;  end;

switch in,
  case {'centered','cell-centered'}, n = prod(m);
    % decompose Yin
    for j=1:dim, y{j} = reshape(yc((j-1)*n+(1:n)),m); end;
    switch out,
      case {'centered','cell-centered'},
      case 'staggered', for j=1:dim, y{j} = extend(y{j},j); end;
      case 'nodal',     
        for j=1:dim, 
          for k=1:dim, y{j} = extend(y{j},k); end; 
        end;
    end;

  case 'staggered', 
    % decompose Yin, get dimensions right
    E = eye(dim) + m'*ones(1,dim); n = cumsum([0,prod(E)]);
    for j=1:dim, y{j} = reshape(yc(n(j)+1:n(j+1)),E(:,j)'); end;
    switch out,
      case {'centered','cell-centered'}, for j=1:dim, y{j} = reduce(y{j},j); end;
      case 'staggered', 
      case 'nodal',    
        for j=1:dim, 
          for k=setdiff(1:dim,j), y{j} = extend(y{j},k); end; 
        end;
    end;

  case 'nodal', n = prod(m+1);
    % decompose Yin
    for j=1:dim, y{j} = reshape(yc((j-1)*n+(1:n)),m+1); end;
    switch out,
      case {'centered','cell-centered'},
        for j=1:dim, for k=1:dim, y{j} = reduce(y{j},k); end; end;
      case 'staggered',
        for j=1:dim, for k=setdiff(1:dim,j) y{j} = reduce(y{j},k); end; end;
      case 'nodal', 
    end;
    
end;
y1 = y{1}; y2 = y{2}; y3 = y{3}; yc = [y1(:);y2(:);y3(:)];

% -------------------------------------------------------------------------
% the following operators act on x=y{k} in the j-th direction,
% permuting x makes j the first direction, averaging (extend creates two 
% additional outside points based on linear BC) and permute back

function x = reduce(x,j)
m = size(x); J = [j,setdiff(1:length(size(x)),j)];
x = permute(x,J); x = (x(1:end-1,:,:)+x(2:end,:,:))/2; x = ipermute(x,J);

function x = extend(x,j)
m = size(x); J = [j,setdiff(1:length(size(x)),j)]; x = permute(x,J);
x = [1.5*x(1,:,:)-0.5*x(2,:,:);(x(1:end-1,:,:)+x(2:end,:,:))/2;...
     1.5*x(end,:,:)-0.5*x(end-1,:,:)];
x = ipermute(x,J);

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