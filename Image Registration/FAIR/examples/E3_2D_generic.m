%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 2D interpolation and visualization
%
% generic example for usage of interpolation and visualization
% loads data, initializes viewer and interpolator, shows some results
%
% - load data ('US.jpg')
% - display and visualize data  (viewImage2D)
% - view image in high res and low res
%==============================================================================

clear, close all, help(mfilename); 
echo on

fprintf('%s\n','generic interpolation and visualization')

fprintf('%s\n','load data, setup image viewer')
dataT = double(imread('US.jpg'));
m     = size(dataT);
omega = [0,size(dataT,1),0,size(dataT,2)];

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));

% shortcuts for labeling off plots
dimstr  = @(m) sprintf('%s=[%s]',inputname(1),sprintf(' %d',m));
lblstr  = @(str,m) sprintf('%s, %s',str,dimstr(m));

fprintf('%s\n','visualize original data')
FAIRfigure(1,'figname',mfilename); clf; 
viewImage(dataT,omega,m);
title(lblstr('highres',m),'fontsize',20); set(gca,'fontsize',20); pause;

fprintf('%s\n','setup interpolater')
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
T = getSplineCoefficients(dataT,'out',0);

fprintf('%s\n','generate some points and interpolate')
m  = [256,128];
xc = getCellCenteredGrid(omega,m);
Tc = inter(T,omega,xc);

FAIRfigure(2,'figname',mfilename); clf; 
viewImage(Tc,omega,m);
title(lblstr('lowres',m),'fontsize',20);  set(gca,'fontsize',20);

echo off