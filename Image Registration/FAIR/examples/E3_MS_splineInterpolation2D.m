%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: multiscale spline interpolation in 2D
%
% - load data   (setupUSData)
% - interpolate (splineInter) enabling different theta's
% - visualize   (viewImage2D)
%==============================================================================

setupUSData; m = 128*[3,2]; xc = getCellCenteredGrid(omega,m);
T = getSplineCoefficients(dataT,'dim',2,'regularizer','gradient','theta',50);
figure(2); clf; viewImage2D(splineInter(T,omega,xc),omega,m); 
colormap(gray(256));
