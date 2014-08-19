% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Example for MLPIR, MultiLevel Parametric Image Registration
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     affine2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% load data, set viewer
setupPETCTdata;

% setup interpolator, distance
theta = 0; % something to play with
inter('reset','inter','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','MI'); 

% setup transformation, regularization 
trafo('reset','trafo','affine2D'); 
wc =  MLPIR(MLdata,'plotIter',0,'plotMLiter',0);
