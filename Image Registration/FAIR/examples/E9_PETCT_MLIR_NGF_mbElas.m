% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             NGF
%   - pre-registration     rigid2D
%   - regularizer          mbElastic
%   - optimization         Gauss-Newton
% ===============================================================================

close all, help(mfilename);

setupPETCTdata
inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','NGF','edge',20);
trafo('reset','trafo','rigid2D');
regularizer('reset','regularizer','mbElastic','alpha',0.1,'mu',1,'lambda',0);
[yc,wc,his] = MLIR(MLdata,...
  'minLevel',4,'maxIterNPIR',25,'parametric',1,'plotMLiter',0);
