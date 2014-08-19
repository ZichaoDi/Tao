%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: 2D interpolation and transformations
%
% - load data ('USfair.jpg')
% - setup  viewer          (viewImage2D), 
%          interpolator    (linearInter), 
%          transformation  (rotation2D)
% - rotate image
%
% a non-trivial example for interpolation, transformtion, and visualization
% runs a loop over transformation parameters (rotation angles), 
%   computes the transformed grid, transformed image, and visualizes these
%==============================================================================

clear, close all, help(mfilename)

fprintf('%s\n','load data')
dataT=imread('LenaCropped.tiff');
m     = floor([512 512]./1);
omega = [0,512,0,512];  
xc    = getCellCenteredGrid(omega,m);

fprintf('%s\n','setup image viewer')
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));

fprintf('%s\n','setup interpolator')
inter('reset','inter','linearInter');


TC = inter('coefficients',dataT,[],omega,'regularizer','moments','theta',1e1);
viewImage(TC,[0 300 0 300],[300 300]);
center = (omega(2:2:end)-omega(1:2:end))'/2;
Tc = inter(TC,omega,xc);    % interpolate T on the grid
figure,vh = viewImage(Tc,omega,m);