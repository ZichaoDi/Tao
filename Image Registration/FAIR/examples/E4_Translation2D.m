%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: translation of an US image
%
% - load data                  (setupUSData)
% - transform                  (translation2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData;  
fprintf('trafo=%s\n',translation2D([]));
wc = [-50;0]; 
xc = getCellCenteredGrid(omega,m);  
yc = translation2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
figure(1); viewImage2D(dataT,[0 1000 0 1000],[513 386],'colormap','gray(256)'); 

figure(2); viewImage2D(Tc,[0 1000 0 1000],m,'colormap','gray(256)'); 

