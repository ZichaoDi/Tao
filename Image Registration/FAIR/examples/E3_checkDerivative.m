%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: derivative check (use data from E3_splineInterpolation2D)
%
%==============================================================================

E3_splineInterpolation2D;
fctn = @(x) splineInter(T,omega,x);
[fig,ph,th] = checkDerivative(fctn,xf(:));
