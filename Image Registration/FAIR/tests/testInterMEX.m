% (c) Fabian Gigengack 2011/04/27, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI
%
% Testfile for different interpolation schemes:
%   - linear interpolation
%   - linear smooth interpolation
%   - spline interpolation
%   - cubic interpolation
% 
% The C-Versions are compared to the MATLAB versions in terms of similarity
% and time.

%% -------------------------------- LINEAR --------------------------------
%% 1D linear
clc
clear

plots = false;
eps   = 1e-14;

omega = [0,10];
TD    = rand(1,20);
m     = length(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = linspace(-1,11,101);

[Tc,dT] = linearInter(TD,omega,xc);
[Tcmex,dTmex] = linearInterMex(TD,omega,xc);

if plots
    figure(1); clf;
    subplot(2,1,1); plot(xc,Tc,'b-',XD,TD,'ro');
    subplot(2,1,2); plot(xc,Tcmex,'b-',XD,TD,'ro');
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; linearInter(TD,omega,xc); t = toc;
tic; linearInterMex(TD,omega,xc); tmex = toc;

disp(' 1D:');
disp(['    Time for linearInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for linearInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                 Speedup: ' num2str(t/tmex)])


%% 2D linear
omega = [0,10,0,8];
TD    = rand(10,10);
m     = size(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
DD = reshape([XD;TD(:)],[],3);

[Tc,dT] = linearInter(TD,omega,xc);
Dc = reshape([xc;Tc],[5*m,3]);
[Tcmex,dTmex] = linearInterMex(TD,omega,xc);
Dcmex = reshape([xc;Tcmex],[5*m,3]);

if plots
    figure(1); clf;
    subplot(1,2,1); surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
    subplot(1,2,2); surf(Dcmex(:,:,1),Dcmex(:,:,2),Dcmex(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; linearInter(TD,omega,xc); t = toc;
tic; linearInterMex(TD,omega,xc); tmex = toc;

disp(' 2D:');
disp(['    Time for linearInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for linearInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                 Speedup: ' num2str(t/tmex)])


%% 3D linear
omega = [0,1,0,2,0,1];
TD    = rand(10,10,10);
m     = size(TD);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1 -1 1],5*m);

[Tc,dT] = linearInter(TD,omega,xc);
[Tcmex,dTmex] = linearInterMex(TD,omega,xc);

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; linearInter(TD,omega,xc); t = toc;
tic; linearInterMex(TD,omega,xc); tmex = toc;

disp(' 3D:');
disp(['    Time for linearInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for linearInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                 Speedup: ' num2str(t/tmex)])


%% ---------------------------- LINEAR SMOOTH -----------------------------
%% 1D linear smooth

omega = [0,10];
TD    = rand(1,20);
m     = length(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = linspace(-1,11,101);

[Tc,dT] = linearInterSmooth(TD,omega,xc);
[Tcmex,dTmex] = linearInterSmoothMex(TD,omega,xc);

if plots
    figure(1); clf;
    subplot(2,1,1); plot(xc,Tc,'b-',XD,TD,'ro');
    subplot(2,1,2); plot(xc,Tcmex,'b-',XD,TD,'ro');
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; linearInterSmooth(TD,omega,xc); t = toc;
tic; linearInterSmoothMex(TD,omega,xc); tmex = toc;

disp(' 1D:');
disp(['    Time for linearInterSmooth: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for linearInterSmoothMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                       Speedup: ' num2str(t/tmex)])

%% 2D linear smooth
omega = [0,10,0,8];
TD    = rand(10,10);
m     = size(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
DD = reshape([XD;TD(:)],[],3);

[Tc,dT] = linearInterSmooth(TD,omega,xc);
Dc = reshape([xc;Tc],[5*m,3]);
[Tcmex,dTmex] = linearInterSmoothMex(TD,omega,xc);
Dcmex = reshape([xc;Tcmex],[5*m,3]);

if plots
    figure(1); clf;
    subplot(1,2,1); surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
    subplot(1,2,2); surf(Dcmex(:,:,1),Dcmex(:,:,2),Dcmex(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; linearInterSmooth(TD,omega,xc); t = toc;
tic; linearInterSmoothMex(TD,omega,xc); tmex = toc;

disp(' 2D:');
disp(['    Time for linearInterSmooth: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for linearInterSmoothMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                       Speedup: ' num2str(t/tmex)])

%% 3D linear smooth
omega = [0,1,0,2,0,1];
TD    = rand(10,10,10);
m     = size(TD);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1 -1 1],5*m);

[Tc,dT] = linearInterSmooth(TD,omega,xc);
[Tcmex,dTmex] = linearInterSmoothMex(TD,omega,xc);

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; linearInterSmooth(TD,omega,xc); t = toc;
tic; linearInterSmoothMex(TD,omega,xc); tmex = toc;

disp(' 3D:');
disp(['    Time for linearInterSmooth: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for linearInterSmoothMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                       Speedup: ' num2str(t/tmex)])


%% -------------------------------- SPLINE --------------------------------
%% 1D spline

omega = [0,10];
TD    = rand(1,20);
m     = length(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = linspace(-1,11,101);

[Tc,dT] = splineInter(getSplineCoefficients(TD,'dim',1,'out',0),omega,xc);
[Tcmex,dTmex] = splineInterMex(getSplineCoefficients(TD,'dim',1,'out',0),omega,xc);

if plots
    figure(1); clf;
    subplot(2,1,1); plot(xc,Tc,'b-',XD,TD,'ro');
    subplot(2,1,2); plot(xc,Tcmex,'b-',XD,TD,'ro');
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
TDc = getSplineCoefficients(TD,'out',0);
tic; splineInter(TDc,omega,xc); t = toc;
tic; splineInterMex(TDc,omega,xc); tmex = toc;

disp(' 1D:');
disp(['    Time for splineInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for splineInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                 Speedup: ' num2str(t/tmex)])


%% 2D spline
omega = [0,10,0,8];
TD    = rand(10,10);
m     = size(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
DD = reshape([XD;TD(:)],[],3);

[Tc,dT] = splineInter(getSplineCoefficients(TD,'out',0),omega,xc);
Dc = reshape([xc;Tc],[5*m,3]);
[Tcmex,dTmex] = splineInterMex(getSplineCoefficients(TD,'out',0),omega,xc);
Dcmex = reshape([xc;Tcmex],[5*m,3]);

if plots
    figure(1); clf;
    subplot(1,2,1); surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
    subplot(1,2,2); surf(Dcmex(:,:,1),Dcmex(:,:,2),Dcmex(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
TDc = getSplineCoefficients(TD,'out',0);
tic; splineInter(TDc,omega,xc); t = toc;
tic; splineInterMex(TDc,omega,xc); tmex = toc;

disp(' 2D:');
disp(['    Time for splineInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for splineInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                 Speedup: ' num2str(t/tmex)])


%% 3D spline
omega = [0,1,0,2,0,1];
TD    = rand(10,10,10);
m     = size(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1 -1 1],5*m);

[Tc,dT] = splineInter(getSplineCoefficients(TD,'out',0),omega,xc);
[Tcmex,dTmex] = splineInterMex(getSplineCoefficients(TD,'out',0),omega,xc);

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
TDc = getSplineCoefficients(TD,'out',0);
tic; splineInter(TDc,omega,xc); t = toc;
tic; splineInterMex(TDc,omega,xc); tmex = toc;

disp(' 3D:');
disp(['    Time for splineInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for splineInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                 Speedup: ' num2str(t/tmex)])

%% -------------------------------- CUBIC ---------------------------------
%% 1D cubic

omega = [0,10];
TD    = rand(1,20);
m     = length(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = linspace(-1,11,101);

[Tc,dT] = cubicInter(TD,omega,xc);
[Tcmex,dTmex] = cubicInterMex(TD,omega,xc);

if plots
    figure(1); clf;
    subplot(2,1,1); plot(xc,Tc,'b-',XD,TD,'ro');
    subplot(2,1,2); plot(xc,Tcmex,'b-',XD,TD,'ro');
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; cubicInter(TD,omega,xc); t = toc;
tic; cubicInterMex(TD,omega,xc); tmex = toc;

disp(' 1D:');
disp(['    Time for cubicInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for cubicInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                Speedup: ' num2str(t/tmex)])


%% 2D cubic
omega = [0,10,0,8];
TD    = rand(10,10);
m     = size(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
DD = reshape([XD;TD(:)],[],3);

[Tc,dT] = cubicInter(TD,omega,xc);
Dc = reshape([xc;Tc],[5*m,3]);
[Tcmex,dTmex] = cubicInterMex(TD,omega,xc);
Dcmex = reshape([xc;Tcmex],[5*m,3]);

if plots
    figure(1); clf;
    subplot(1,2,1); surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
    subplot(1,2,2); surf(Dcmex(:,:,1),Dcmex(:,:,2),Dcmex(:,:,3));  hold on;
    plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40);  hold off;
end

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; cubicInter(TD,omega,xc); t = toc;
tic; cubicInterMex(TD,omega,xc); tmex = toc;

disp(' 2D:');
disp(['    Time for cubicInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for cubicInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                Speedup: ' num2str(t/tmex)])


%% 3D cubic
omega = [0,1,0,2,0,1];
TD    = rand(10,10,10);
m     = size(TD);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1 -1 1],5*m);

[Tc,dT] = cubicInter(TD,omega,xc);
[Tcmex,dTmex] = cubicInterMex(TD,omega,xc);

OK_T = max(abs(Tc(:)-Tcmex(:))) < eps;
OK_dT = full(max(abs(dT(:)-dTmex(:)))) < eps;

% Ohne Ableitung
tic; cubicInter(TD,omega,xc); t = toc;
tic; cubicInterMex(TD,omega,xc); tmex = toc;

disp(' 3D:');
disp(['    Time for cubicInter: ' num2str(t) ' - T OK? ' num2str(OK_T) '!']);
disp([' Time for cubicInterMex: ' num2str(tmex) ' - dT OK? ' num2str(OK_dT) '!']);
disp(['                Speedup: ' num2str(t/tmex)])

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