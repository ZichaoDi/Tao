% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     rigid2D
%   - regularizer          mbCurvature
%   - optimization         lBFGS
% ===============================================================================

close all, help(mfilename);

setupPETCTdata
inter('reset','inter','splineInter','regularizer','none','theta',1e-3);
distance('reset','distance','MI','nT',32,'nR',32);
trafo('reset','trafo','rigid2D');
regularizer('reset','regularizer','mbCurvature','alpha',1e-1);
[yc,wc,his] = MLIR(MLdata,...
  'PIR',@lBFGS,'PIRobj',@PIRBFGSobjFctn,...
  'NPIR',@lBFGS,'NPIRobj',@NPIRBFGSobjFctn,...
  'minLevel',4,'maxLevel',7,'maxIterNPIR',25,'parametric',1,'plotMLiter',0);


