% (c) Lars Ruthotto 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/
%
% test MEX distance measures and parallelization using OpenMP

mfiles = {'NGFdot','SSD','NCC'};
mexfiles = {'NGFmex','SSDmex','NCCmex'};
cfiles = {'NGFdotMex.cpp','SSDmexC.cpp','NCCmexC.cpp'};
cores = 2;
samples = 20;

% 1D
for int=1:size(mfiles,2),
    mfile = mfiles{int}; mexfile = mexfiles{int}; cfile=cfiles{int};
    for dim=1:3,
        fprintf('\nTESTING %s %d DIM interpolation \n',mfile,dim);
        switch dim
            case 1
                omega = [0,10];
                m     = 100;
                Tc    = rand(m,1);
                Rc    = rand(m,1);
            case 2
                omega = [0,10,0,8];
                m     = 512 * [1 1];
                Tc    = rand(m);
                Rc    = rand(m);
            case 3
                omega = [0,1,0,2,0,1];
                m =   4*[13,16,7];
                Tc    = rand(m);
                Rc    = rand(m);
        end
        Tc = Tc(:); Rc = Rc(:);
        tic;
        [Dc0,rc0,dD0,dr0,d2psi0] = feval(mfile,Tc,Rc,omega,m);
        MATtime = toc;
        dT0 = spdiags(dT0,prod(mf)*[0:dim-1]);
        FAIRmake(cfile);
        sMEXtime = zeros(1,samples);
        for i=1:samples,
            tic;
            [Dcs,rcs,dDs,drs,d2psis] = feval(mexfile,Tc,Rc,omega,m);
            sMEXtime(i) = toc;
        end
        dTs = spdiags(dTs,prod(mf)*[0:dim-1]);
        sMEXtime = mean(sMEXtime);
        FAIRmake(cfile,'cores',cores);
        pMEXtime = zeros(1,samples);
        for i=1:samples,
            tic;
            [Dcp,rcp,dDp,drp,d2psip] = feval(mexfile,Tc,Rc,omega,m);
            pMEXtime(i) = toc;
        end
        dTp = spdiags(dTp,prod(mf)*[0:dim-1]);
      
        pMEXtime = mean(pMEXtime);
        REs = norm(Dcs(:)-Dc0(:))/norm(Dc0(:));   REdDs = norm(dDs(:)-dD0(:))/norm(dD0(:));
        REp = norm(Dcp(:)-Dc0(:))/norm(Dc0(:));   REdDp = norm(dDp(:)-dD0(:))/norm(dD0(:));
        REp2 = norm(Dcs(:)-Dcp(:))/norm(Dcs(:));  REdDp2 = norm(dDs(:)-dDp(:))/norm(dDs(:)); 
        
        fprintf('Matlab vs. MEX(serial) \t: RE(Tc) %1.2e  RE(dT) %1.2e  speedup %2.2f \n',REs,REdDs,MATtime/sMEXtime);
        fprintf('Matlab vs. MEX(parallel): RE(Tc) %1.2e  RE(dT) %1.2e  speedup %2.2f \n',REp,REdDp,MATtime/pMEXtime);
        fprintf('MEX(s) vs. MEX(p) \t: RE(Tc) %1.2e  RE(dT) %1.2e  speedup %2.2f \n',REp2,REdDp2,sMEXtime/pMEXtime);
    end
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