% ===============================================================================
% (c) Jan Modersitzki 2011/01/02, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: MLPIR, MultiLevel Parametric Image Registration, regularized
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=3:7, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     rigid2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

setupHNSPData;                            % set up data
distance('reset','distance','SSD');       % specify distance measure
inter('reset','inter','splineInter');   % specify interpolator
trafo('reset','trafo','rigid2D');         % specify transformation
[wc,his] = MLPIR(MLdata,'minLevel',3,'plotMLiter',1);
