% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function imgOverlay(T0,R0,omega,m)
%
% overlay of two images
%
% Input:
%   T0          discretized image
%   R0  		discretized image
%	omega 		describing the domain
%	m 			number of discretization points

function imgOverlay(T0,R0,omega,m)

if nargin==0
    runMinimalExample; return;
end

clf;

TT = rescale(T0,0,255,0.5);

RR = R0; RR(RR>220) = 220; RR(RR<50) = 5;
RR = rescale(RR,0,255,0.5);

RR = round(RR);
ZZ = zeros(length(RR),3);
grayMap = gray(256);

for j=0:255;
  J = find(RR==j);
  ZZ(J,:) = grayMap(j+1);
end;


TT = permute(reshape(TT*[1,0,0]/255,[m,3]),[2,1,3]);
RR = permute(reshape(ZZ,[m,3]),[2,1,3]);

alphadata  = 0.5*(reshape(T0,m)'>5);


ih=overlayImage2D(TT,RR,omega,m,'scale',0,'alphadata',alphadata); axis off

set(gca,'position',[0 0 1 1]);

function runMinimalExample
setupHandData;
imgOverlay(dataT(:),dataR(:),omega,m); %#ok<NODEF>
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