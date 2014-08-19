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

setupUSData; m = 128*[3,2]; 
xc = getCellCenteredGrid(omega,m);
cc = @(reg,theta) getSplineCoefficients(dataT,'dim',2,'regularizer',reg,'theta',theta);

Tc   = @(theta) splineInter(cc('moments',theta),omega,xc);
name = @(theta) fullfile(FAIRpath,'temp',sprintf('US-MS-theta=%d.jpg',theta));

theta = 1e3
T = Tc(theta);
figure(2); clf; viewImage2D(T,omega,m); 
colormap(gray(256));
imwrite(uint8(flipud(reshape(T,m)')),name(theta))