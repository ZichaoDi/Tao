% (c) Fabian Gigengack 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI
%
% test MEX distance measures

% setupPETCTData;
% xc = getCellCenteredGrid(omega,m);
% Tc = linearInter(dataT,omega,xc);
% Rc = linearInter(dataR,omega,xc);

setup3DmiceData;
xc = getCellCenteredGrid(omega,m);
Tc = linearInter(dataT,omega,xc);
Rc = linearInter(dataR,omega,xc);

mfiles = {@NGFdot,@SSD,@NCC};
mexfiles = {@NGFmex,@SSDmex,@NCCmex};
tol = 1e-12;

for i=1:numel(mfiles)
    disp([func2str(mfiles{i}) ' vs. ' func2str(mexfiles{i})])
    % Without Derivatives
    tic;
    [Dc,rc] = mfiles{i}(Tc,Rc,omega,m);
    t = toc;
    tic;
    [Dcmex,rcmex] = mexfiles{i}(Tc,Rc,omega,m);
    tmex = toc;

    disp(['Speedup without derivatives: ' num2str(t/tmex)]);

    % With Derivatives
    tic;
    [Dc,rc,dD,drc,d2psi] = mfiles{i}(Tc,Rc,omega,m);
    t = toc;
    tic;
    [Dcmex,rcmex,dDmex,drcmex,d2psimex] = mexfiles{i}(Tc,Rc,omega,m);
    tmex = toc;

    disp(['   Speedup with derivatives: ' num2str(t/tmex)]);

    disp(['OK? ' num2str(abs((Dc-Dcmex)/min(Dc,Dcmex))<tol && max(abs(rc(:)-rcmex(:)))<tol ...
        && max(abs(dD(:)-dDmex(:)))<tol && full(max(abs(drc(:)-drcmex(:))))<tol ...
        && abs(d2psi-d2psimex)<tol)]);
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