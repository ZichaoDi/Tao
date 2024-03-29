% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function pos = FAIRposition(varargin)
%
% computes the default position of a figure in FAIR

function pos = FAIRposition(varargin)

persistent theScreen

if isempty(theScreen),
  theScreen   = get(0,'screensize'); theScreen(4) = theScreen(4) - 95;
end;

fig       = getValue('fig',        1,varargin{:});
increment = getValue('increment', 20,varargin{:});
position  = getValue('position',  [],varargin{:});

if isempty(position),
  dx     = (fig-1)*increment;
  dy     = (fig-1)*increment;
  
  width  = getValue('width', min(floor(theScreen(3))/2,1200), varargin{:});
  width  = max(10,min(theScreen(3),width));
  height = getValue('height',min(floor(theScreen(4))/2, 800), varargin{:});
  height = max(10,min(theScreen(4),height));
  
  left   = getValue('left',  min(dx,theScreen(3)-width)+1,  varargin{:});
  left   = min(max(1,left),floor(0.95*theScreen(3)));
  bottom = getValue('bottom',theScreen(4)-height+1-min(dx,height),varargin{:});
  bottom = min(max(1,bottom),floor(0.95*theScreen(3)));
  
  pos = [left,bottom,width,height];
  return;
end;

if ~isnumeric(position), 
  if strcmp(position,'default'),
    pos = FAIRposition('fig',fig,...
      'width', 1200,...
      'height',800);
    return;
  end;
  
  error('17+4'); 
end;
  
if length(position) == 1,
  pos = FAIRposition('fig',fig,'height',position);
  return;
elseif length(position) == 4,
  pos = FAIRposition('fig',fig,...
    'left',  position(1),...
    'bottom',position(2),...
    'width', position(3),...
    'height',position(4));
  return;
end;

position
keyboard

function value = getValue(name,default,varargin)
value  = default;
j = strcmp(name,varargin);
if isempty(j), return; end;
if all(j==0),  return; end;
value = varargin{max(find(j))+1};

% f = 0;
% if isempty(p),
%   return;
% elseif strcmp(p,'default'),
% %   p=[800,1200]; f = fig;
%   p=[400,600]; f = fig;
% elseif strcmp(p,'print'),
%   p=[850,500,1600,1000]; f = fig;
% end;
%
%
% s   = get(0,'screensize');
% %pos = @(a,b) [s(3)+1-a-75-f*25,s(4)+1-b-100-f*25, a b];
% pos = @(a,b) [f*25,s(4)+1-b-100-f*25, a b];
% if p == 1;
%   p = pos(1200,800);
% elseif length(p) == 1,
%   p = pos(1200,p);
% elseif length(p) == 2;
%   p = pos(p(2),p(1));
% end;
% set(fig,'position',p);

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