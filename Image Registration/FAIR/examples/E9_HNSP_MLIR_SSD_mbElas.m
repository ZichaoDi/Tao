% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2010/03/22, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), m=[ 512 256]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mbElastic
% ===============================================================================

close all, help(mfilename);

setupHNSPData
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbElastic','alpha',5e2,'mu',1,'lambda',0);
% regularizer('reset','regularizer','mfElastic','alpha',5e2,'mu',1,'lambda',0);
% regularizer('reset','regularizer','mbCurvature','alpha',1e1);
% regularizer('reset','regularizer','mfCurvature','alpha',1e1);

[yc,wc,his] = MLIR(MLdata,'maxLevel',7,'parametric',1);
