%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: rotation of an US image, based on the book's version
%
% - load data                  (setupUSData)
% - transform                  (book_rigid2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData; c = (omega(2:2:end)+omega(1:2:end))'/2; alpha = pi/6; 
rot = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
wc  = [alpha;(eye(2)-rot)*c]; 
xc  = getCellCenteredGrid(omega,m);
yc  = book_rigid2D(wc,xc);
Tc  = linearInter(dataT,omega,yc);
figure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
 
