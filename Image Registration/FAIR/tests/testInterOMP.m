% (c) Lars Ruthotto 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/
%
% test MEX interpolation schemes and parallelization using OpenMP

mfiles = {'linearInter','cubicInter','splineInter','linearInterSmooth'};
mexfiles = {'linearInterMex','cubicInterMex','splineInterMex','linearInterSmoothMex'};
cfiles = {'linearInterMexC.cpp','cubicInterMexC.cpp','splineInterMexC.cpp','linearInterSmoothMexC.cpp'};
cores = 2;
samples = 20;

% 1D
for int=1:size(mfiles,2),
    mfile = mfiles{int}; mexfile = mexfiles{int}; cfile=cfiles{int};
    for dim=1:3
        fprintf('\nTESTING %s %d DIM interpolation \n',mfile,dim);
        switch dim
            case 1
                omega = [0,10];
                m     = 100;
                TD    = rand(m,1);
                mf = 20001;
                xc    = linspace(-1,11,mf);
            case 2
                omega = [0,10,0,8];
                m     = 128 * [1 1];
                TD    = rand(m);
                mf = 5*m;
                xc    = getCellCenteredGrid(omega+[-1 1 -1 1],mf);
            case 3
                omega = [0,1,0,2,0,1];
                m =   [13,16,7];
                XD    = getCellCenteredGrid(omega,m);
                Y     = reshape(XD,[m,3]);
                TD    = (Y(:,:,:,1)-0.5).^2 + (Y(:,:,:,2)-0.75).^2 + (Y(:,:,:,3)-0.5).^2 <= 0.15;
                mf = 5*m;
                xc = getCellCenteredGrid(omega,mf);
        end
        
        tic;
        [T0,dT0] = feval(mfile,TD,omega,xc);
        MATtime = toc;
        dT0 = spdiags(dT0,prod(mf)*[0:dim-1]);
        FAIRmake(cfile,'verbose',-1);
        sMEXtime = zeros(1,samples);
        [Ts,dTs] = feval(mexfile,TD,omega,xc);
        for i=1:samples,
            tic;
            Ts = feval(mexfile,TD,omega,xc);
            sMEXtime(i) = toc;
        end
        dTs = spdiags(dTs,prod(mf)*[0:dim-1]);
        sMEXtime = mean(sMEXtime);
        FAIRmake(cfile,'cores',cores,'verbose',-1);
        pMEXtime = zeros(1,samples);
        [Tp,dTp] = feval(mexfile,TD,omega,xc);
        for i=1:samples,
            tic;
            Tp = feval(mexfile,TD,omega,xc);
            pMEXtime(i) = toc;
        end
        dTp = spdiags(dTp,prod(mf)*[0:dim-1]);
      
        pMEXtime = mean(pMEXtime);
        REs = norm(Ts(:)-T0(:))/norm(T0(:));   REdTs = norm(dTs(:)-dT0(:))/norm(dT0(:));
        REp = norm(Tp(:)-T0(:))/norm(T0(:));   REdTp = norm(dTp(:)-dT0(:))/norm(dT0(:));
        REp2 = norm(Ts(:)-Tp(:))/norm(Ts(:));  REdTp2 = norm(dTs(:)-dTp(:))/norm(dTs(:)); 
        
        fprintf('Matlab vs. MEX(serial) \t: RE(Tc) %1.2e  RE(dT) %1.2e  speedup %2.2f \n',REs,REdTs,MATtime/sMEXtime);
        fprintf('Matlab vs. MEX(parallel): RE(Tc) %1.2e  RE(dT) %1.2e  speedup %2.2f \n',REp,REdTp,MATtime/pMEXtime);
        fprintf('MEX(s) vs. MEX(p) \t: RE(Tc) %1.2e  RE(dT) %1.2e  speedup %2.2f \n',REp2,REdTp2,sMEXtime/pMEXtime);
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