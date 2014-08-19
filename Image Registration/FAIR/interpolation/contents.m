%==============================================================================
% (c) Jan Modersitzki 2011/02/11, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR INTERPOLATION TOOLBOX
%
% general purpose tools:
%   contents              - this file
%   inter                 - the specific interpolation scheme used in FAIR
%            initialize: inter('reset','inter','linearIOnter1D');
%            use:        [Tc,dT] = inter(T,omega,Y);
%   getSplineCoefficients - computes spline coefficients
%   motherSpline:         - 1D normalized basis spline and its derivative
%
%   schemes are available for dimension dim, where dim=1,2,3,
% 
%   nnInter               - next neighbor interpolation schemes (not recommended for optimization)
%   linearInter           - linear interpolation schemes
%   linearInterMatlab     - wrapper for MATLAB's interp(dim)
%   linearInterSmooth     - linear approximation schemes
%   splineInter           - spline interpolation schemes
%   cubicInter            - fast local cubic spline interpolation schemes
%
%   mex versions of selected schemes:
%   linearInterMex        - wrapper for C-code linearInterMexC
%   linearInterMexC       - CPP version of linear interpolation
%   linearInterSmoothMex  - wrapper for C-code linearInterSmoothMexC
%   linearInterSmoothMexC - CPP version of smooth linear approximation
%   splineInterMex        - wrapper for C-code splineInterMexC.c
%   splineInterMexC       - CPP version of spline interpolation
%   cubicInterMex         - wrapper for C-code cubicInterMexC.c
%   cubixInterMexC        - CPP version of local cubic spline interpolation
%
%  see also E9_Hands_MLIR_SSD_mbElas and BigTutorialInter
%==============================================================================
help(mfilename)