% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% function varargout = viewSlices(T,omega,m,varargin)
%
% visualizes a 3D image as slides
%
% Input:
%   T           discretized image
%   omega		describing the domain
%	m 			number of discretization points
%   varargin    optional parameters like {'axis','off'}
%
% Output:
%   ih			image handle
%   B			the mosaic image
%   frames      number of frames for ij-directions

function varargout = viewSlices(T,omega,m,varargin)

if nargin==0
    runMinimalExample; return;
end

% set default parameter
s1 = round(m(1)/2);
s2 = round(m(2)/2);
s3 = round(m(3)/2);
% threshold = min(T(:))+0.1*(max(T(:))-min(T(:)));

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

T  = reshape(T,m);
h  = (omega(2:2:end)-omega(1:2:end))./m;

% show ortho-slices to first dimension
ih = []; h1 = [];
for i=1:length(s1), 
  x1 = (omega(1) + (s1(i)-.5) * h(1)) * ones(m(2:3));  % fix x1 coordinate
  X  = reshape(getCellCenteredGrid(omega(3:end),m(2:3)),[],2); % 2D grid
  x2 = reshape(X(:,1),m(2:3));
  x3 = reshape(X(:,2),m(2:3));

  Ti = squeeze(T(s1(i),:,:));
  
  sh(i) = surf(x1,x2,x3);
  axis(omega);
  hi = findobj('type','surface');
  h1(end+1) = hi(1);
  set(h1(end),'CData',Ti,'FaceColor','texturemap','edgecolor','none');
end;
ih = [ih;h1]; h1 = [];

% show ortho-slices to second dimension
for i=1:length(s2),
  x2 = (omega(3) + (s2(i)-.5) * h(2)) * ones(m([1 3]));
  X  = reshape(getCellCenteredGrid(omega([1 2 5 6]),m([1 3])),[],2);
  x1 = reshape(X(:,1),m([1 3]));
  x3 = reshape(X(:,2),m([1 3]));
  Ti = squeeze(T(:,s2(i),:));
  
  sh(i) = surf(x1,x2,x3);
  axis(omega);
  hi = findobj('type','surface');
  h1(end+1) = hi(1);
  set(h1(end),'CData',Ti,'FaceColor','texturemap','edgecolor','none');
end;
ih = [ih;h1]; h1 = [];

% show ortho-slices to third dimension
for i=1:length(s3),
  x3 = (omega(5) + (s3(i)-.5) * h(3)) * ones(m(1:2));
  X  = reshape(getCellCenteredGrid(omega(1:4),m(1:2)),[],2);
  x1 = reshape(X(:,1),m(1:2));
  x2 = reshape(X(:,2),m(1:2));
  Ti = squeeze(T(:,:,s3(i)));
  
  sh(i) = surf(x1,x2,x3);
  axis(omega);
  hi = findobj('type','surface');
  h1(end+1) = hi(1);
  set(h1(end),'CData',Ti,'FaceColor','texturemap','edgecolor','none');
end;
ih = [ih;h1];

if nargout == 1, varargout = {ih};  end;

function runMinimalExample
load mice3D; viewSlices(dataT,omega,m);
help(mfilename)

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