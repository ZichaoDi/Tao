% (c) Lars Ruthotto 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/
%
% test MEX geometry schemes and parallelization using OpenMP

cfile = 'geometrymexC.cpp';

flags = {'V','A','Jac','dJacx','dJacadjx'};
cores = 2;
samples = 1;

for dim=2:3,
    FAIRmake(cfile,'verbose',-1);
    switch dim
        case 2
            omega = [0,10,0,8];
            m     = 512 * [1 1];
            xn    = getNodalGrid(omega,m);
            nTetra = 4;
          case 3
            omega = [0,1,0,2,0,1];
            m     = 4*[13,16,7];
            xn    = getNodalGrid(omega,m);
            nTriangle = 24;
            nTetra = 24;
    end
    xn = xn + 1e-1 * randn(size(xn));
    % get test vectors
    dJacxtest = randn(dim*prod(m+1),1);
    dJacadjxtest = randn(prod(m),1);
    d2Stest = randn(dim*prod(m+1),1);
    
    % get references
    h     = (omega(2:2:end) - omega(1:2:end)) ./ m;
    mOne  = ones(size(m));
    voxel = getNodalGrid(reshape([0;1]*h,1,[]),mOne);
   VRef = geometry(voxel,mOne,'V','matrixFree',1);
    FAIRmake(cfile,'clean',1);
    
    % TEST VOLUME
    flag = 'V';
    FAIRmake(cfile,'verbose',-1);
    sMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vs = geometrymexC(xn(:),m,'V');
        sMEXtime(i) = toc;
    end
    sMEXtime = mean(sMEXtime);
    FAIRmake(cfile,'cores',cores,'verbose',-1);
    pMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vp = geometrymexC(xn(:),m,'V');
        pMEXtime(i) = toc;
    end
    pMEXtime = mean(pMEXtime);
    REp = norm(Vs(:)-Vp(:))/norm(Vs(:));
    fprintf('\n%dD-Volume:\t RE(V) %1.2e   speedup %2.2f \n',dim,REp,sMEXtime/pMEXtime);
    FAIRmake(cfile,'clean',1);
    
    % TEST AREA
    if dim ==3
        ARef = geometry(voxel,mOne,'A','matrixFree',1);
     flag = 'A';
        FAIRmake(cfile,'verbose',-1);
        sMEXtime = zeros(1,samples);
        for i=1:samples,
            tic;
            Vs = geometrymexC(xn(:),m,'A');
            sMEXtime(i) = toc;
        end
        sMEXtime = mean(sMEXtime);
        FAIRmake(cfile,'cores',cores,'verbose',-1);
        pMEXtime = zeros(1,samples);
        for i=1:samples,
            tic;
            Vp = geometrymexC(xn(:),m,'A');
            pMEXtime(i) = toc;
        end
        pMEXtime = mean(pMEXtime);
        REp = norm(Vs(:)-Vp(:))/norm(Vs(:));
        fprintf('\n%dD-Area:\t RE(V) %1.2e   speedup %2.2f \n',dim,REp,sMEXtime/pMEXtime);
    end
    FAIRmake(cfile,'clean',1);
    
    % TEST Jac
    flag = 'Jac';
    FAIRmake(cfile,'verbose',-1);
    sMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vs = geometrymexC(xn(:),m,'Jac');
        sMEXtime(i) = toc;
    end
    sMEXtime = mean(sMEXtime);
    FAIRmake(cfile,'cores',cores,'verbose',-1);
    pMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vp = geometrymexC(xn(:),m,'Jac');
        pMEXtime(i) = toc;
    end
    pMEXtime = mean(pMEXtime);
    REp = norm(Vs(:)-Vp(:))/norm(Vs(:));
    fprintf('\n%dD-Jac:\t\t RE(V) %1.2e   speedup %2.2f \n',dim,REp,sMEXtime/pMEXtime);
    FAIRmake(cfile,'clean',1);
    
    % TEST dJacx
    flag = 'Jac';
    FAIRmake(cfile,'verbose',-1);
    sMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vs = geometrymexC(xn(:),m,'dJacx',dJacxtest);
        sMEXtime(i) = toc;
    end
    sMEXtime = mean(sMEXtime);
    FAIRmake(cfile,'cores',cores,'verbose',-1);
    pMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vp = geometrymexC(xn(:),m,'dJacx',dJacxtest);
        pMEXtime(i) = toc;
    end
    pMEXtime = mean(pMEXtime);
    REp = norm(Vs(:)-Vp(:))/norm(Vs(:));
    fprintf('\n%dD-dJacx:\t RE(V) %1.2e   speedup %2.2f \n',dim,REp,sMEXtime/pMEXtime);
    FAIRmake(cfile,'clean',1);
    % TEST dJacadjx
    flag = 'dJacadjx';
    FAIRmake(cfile,'verbose',-1);
    sMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vs = geometrymexC(xn(:),m,'dJacadjx',dJacadjxtest);
        sMEXtime(i) = toc;
    end
    sMEXtime = mean(sMEXtime);
    FAIRmake(cfile,'cores',cores,'verbose',-1);
    pMEXtime = zeros(1,samples);
    for i=1:samples,
        tic;
        Vp = geometrymexC(xn(:),m,'dJacadjx',dJacadjxtest);
        pMEXtime(i) = toc;
    end
    pMEXtime = mean(pMEXtime);
    REp = norm(Vs(:)-Vp(:))/norm(Vs(:));
    fprintf('\n%dD-dJacadjx:\t RE(V) %1.2e   speedup %2.2f \n',dim,REp,sMEXtime/pMEXtime);
    FAIRmake(cfile,'clean',1);
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