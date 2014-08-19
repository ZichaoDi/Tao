%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: spline transformation of an US image
%
% - load data                  (setupUSData)
% - transform                  (splineTransformation2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData; p = [5,4]; 
xc = getCellCenteredGrid(omega,m);  
splineTransformation2D([],xc,'omega',omega,'m',m,'p',p);  
w1 = zeros(p); w2 = zeros(p);  w2(3,2) = 3; 
wc = [w1(:);w2(:)]; 
yc = splineTransformation2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
figure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
