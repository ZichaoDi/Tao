%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR Distance Toolbox
%
% This folder provides several distance measures, wrappers and utilities.
%
% general purpose tools:
%   contents              - this file
%   distance              - the specific distance measure used in FAIR
%     initialize: distance('reset','distance','SSD');
%     usage:      D = distance(Tc,Rc,omega,m);
%   see also E9_Hands_MLIR_SSD_mbElas and BigTutorialDistance
%
% measures and wrappers:
%   MI          - wrapper for MIcc
%   MIcc        - Mutual Information based on rhoSplineC (C-implementation)
%   MIspline    - Mutual Information based on rhoSpline  (MATLAB-implementation)
%   rhoSpline   - MATLAB implementation of joint entropy estimator (spline kernel)  
%   rhoSplineC  - C based implementation of joint entropy estimator (spline kernel)  
%   NCC         - Normalized Cross Correlation
%   NCCmex      - Wrapper for NCCmexC
%   NCCmexC     - C version of NCC
%   NGF         - wrapper for NGFdot
%   NGFdot      - Normalized Gradient Fields, dot product based (THE implementation)
%   NGFcross    - Normalized Gradient Fields, cross product based (to be used with care)
%   NGFmex      - Wrapper for NGFdotMexC
%   NGFdotMexC  - C version of NGFdot
%   SSD         - Sum of Squared Differences 
%   SSDmex      - Wrapper for SSDmexC
%   SSDmexC     - C version of SSD
%   SSDweighted - enables a masked version
%  see also E9_Hands_MLIR_SSD_mbElas and BigTutorialDistance
%==============================================================================
help(mfilename)