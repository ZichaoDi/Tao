% ===============================================================================
% Example for MLIR, MultiLevel Image Registration
% (c) Jan Modersitzki 2009/04/06, see FAIR.2 and FAIRcopyright.m.
% \url{http://www.cas.mcmaster.ca/~fair/index.shtml}
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             NGF
%   - pre-registration     affine2D
%   - regularizer          mfElastic
%   - optimization         Gauss-Newton
% ===============================================================================

clear; close all; help(mfilename)

setupPETCTdata
theta = 0; % something to play with
inter('reset','inter','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','NGF');
trafo('reset','trafo','affine2D'); wStop = trafo('w0'); w0 = wStop;
regularizer('reset','regularizer','mfElastic','alpha',0.05,'mu',1,'lambda',0);
regularizer('disp');

fprintf(' - run MLIR using sufficient amount of details (minLevel=4)\n');
yc =  MLIR(MLdata,'minLevel',4,'plotIter',0,'plotMLiter',0);

  
