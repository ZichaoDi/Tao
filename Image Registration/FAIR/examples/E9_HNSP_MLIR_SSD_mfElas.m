% ===============================================================================
% Tutorial for FAIR: Non-Parametric Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - transformation       affine2D
%   - regularizer          mfElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

close all, help(mfilename);

% load data, set viewer, interpolator, transformation, distance
setupHNSPData
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mfElastic','alpha',5e2,'mu',1,'lambda',0);

% run MLIR
[yc,wc,his] = MLIR(MLdata,'maxLevel',8);
