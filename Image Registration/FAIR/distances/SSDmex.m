% (c) Jan Modersitzki and Fabian Gigengack 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% http://www.uni-muenster.de/EIMI/
%
% function [Dc,rc,dD,dr,d2psi] = SSDmex(Tc,Rc,omega,m,varargin)
%
% Wrapper for C file for Sum of Squared Differences based distance measure,
% see also SSD for details 
% see also distances/contents

function [Dc,rc,dD,dr,d2psi] = SSDmex(Tc,Rc,omega,m,varargin)

if nargin == 0
    help(mfilename);
    setupHandData;
    xc = getCellCenteredGrid(omega,m);
    Tc = linearInter(dataT,omega,xc);
    Rc = linearInter(dataR,omega,xc);
    D0 = feval(mfilename,Tc,Rc,omega,m);
    fprintf('%s: distance = %s\n',mfilename,num2str(D0))
    return
end

dD = []; dr = []; d2psi = [];
doDerivative = (nargout > 3);

for k=1:2:length(varargin)      % overwrite default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration

if ~doDerivative
    [Dc,rc] = SSDmexC(Tc,Rc,hd);
else
    [Dc,rc,dD,dr,d2psi] = SSDmexC(Tc,Rc,hd);
end

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