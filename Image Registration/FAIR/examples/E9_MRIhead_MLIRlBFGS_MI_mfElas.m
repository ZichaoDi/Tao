% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% 
%   - data                 MRI (head), Omega=(0,128)x(0,128), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     rigid2D
%   - regularizer          mbElastic
%   - optimization         lBFGS
% ===============================================================================

close all, help(mfilename);

setupMRIData
inter('reset','inter','splineInter','regularizer','none','theta',1e-3);
distance('reset','distance','MI','nT',8,'nR',8);
trafo('reset','trafo','rigid2D');
regularizer('reset','regularizer','mfElastic','alpha',1e-4,'mu',1,'lambda',0);
[yc,wc,his] = MLIR(MLdata,...
  'PIRopt', @lBFGS,'PIRobj', @PIRBFGSobjFctn,...
  'NPIRopt',@lBFGS,'NPIRobj',@NPIRBFGSobjFctn,...
  'minLevel',4,'maxLevel',7,...
  'parametric',1,'plotMLiter',0);