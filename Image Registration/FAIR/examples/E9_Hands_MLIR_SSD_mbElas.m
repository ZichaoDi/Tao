% ===============================================================================
% Example for Non-Parametric Image Registration 
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - regularizer          mbElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

close all, help(mfilename)

% load data, initialize image viewer, interpolator, transformation, distance
setupHandData
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
trafo('reset','trafo','affine2D');
distance('reset','distance','SSD');
regularizer('reset','regularizer','mbElastic','alpha',1e3,'mu',1,'lambda',0);

[yc,wc,his] = MLIR(MLdata);
