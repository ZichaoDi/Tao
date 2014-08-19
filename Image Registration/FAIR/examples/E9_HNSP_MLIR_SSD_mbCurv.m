% ===============================================================================
% Tutorial for FAIR: Multi-Level Non-Parametric Image Registration
% (C) 2008/08/13, Jan Modersitzki, see FAIR and FAIRcopyright.m.
%
% - load data (see setupHNSPData)
% - setup  viewer          (viewImage2D), 
%           interpolator    (splineInter), 
%           distance        (SSD), 
%           transformation  (affine2D, not regularized)
%           regularizer     (curvature,matrix based)
% - run optimization 
% ===============================================================================

clear, close all, help(mfilename);

% load data, set viewer, interpolator, transformation, distance
setupHNSPData
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbCurvature','alpha',1e1);

% run optimization
yc =  MLIR(MLdata,'maxLevel',7);

