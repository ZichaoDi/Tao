% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: MLPIR, MultiLevel Parametric Image Registration
% 
%   - data                 Hand, Omega=(0,20)x(0,25), level=3:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D, not regularized
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% load data, initialize image viewer, interpolator, transformation, distance
setupHandData
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
inter('reset','inter','splineInter','regularizer','moments','theta',1e-1);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');

% run Multilevel Parametric Image Registration 
wc = MLPIR(MLdata,'plotIter',0,'plotMLiter',0,'pause',5);
