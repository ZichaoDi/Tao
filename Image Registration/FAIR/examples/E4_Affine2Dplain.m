%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 2D interpolation and transformations
%
% - load data                  (setupUSData)
% - generate affine linear map (book_affine2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setupUSData;
wc = [1 -0.2 50, 0, 0.75 50]'; 
xc = getCellCenteredGrid(omega,m);
yc = book_affine2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
figure(2); clf; viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
