% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
% 
%   - data                 MRI (head), Omega=(0,128)x(0,128), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             NGF
%   - pre-registration     rigid2D
%   - regularizer          mbElastic
%   - optimization         Gauss-Newton
% ===============================================================================

close all, help(mfilename)

setupMRIData
theta = 1e-1; % something to play with
inter('reset','inter','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','NGF'); 
trafo('reset','trafo','rigid2D');
wStop = trafo('w0'); w0 = wStop;
regularizer('reset','regularizer','mbElastic','alpha',1e-2,'mu',1','lambda',0');
regularizer('disp');

fprintf(' - run MLIR using sufficient amount of details (minLevel=4)\n');
[yc,wc,his] =  MLIR(MLdata,'minLevel',4,'maxLevel',7,'plotIter',0,'plotMLiter',0);

  
