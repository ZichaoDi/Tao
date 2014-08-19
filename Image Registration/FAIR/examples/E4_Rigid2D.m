%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: rotation of an US image
%
% - load data                  (setupUSData)
% - transform                  (rigid2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData; fprintf('trafo=%s\n',rigid2D('[]'));
alpha = pi/6; R = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
center = (omega(2:2:end)+omega(1:2:end))'/2;
wc = [alpha;(eye(2)-R)*center]; 
xc = getCellCenteredGrid(omega,m);
yc = rigid2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
figure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
