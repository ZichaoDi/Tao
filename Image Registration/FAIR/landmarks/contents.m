%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR landmark based registration Toolbox
%
% contents           - this file
% LMreg		         - the landmark registration toolit, see E5_Hands_TPS for an example
% YTPS	             - computes the landmark based transformation Y
% evalTPS			 - evaluates the Thin-Plate-Spline solution 
%                      based on coefficients
% getLandmarks       - enables the determination of landmarks in 2D
% getTPScoefficients - computes the coefficients for the TPS solution
% plotLM             - convenient way for plotting landmarks
%   
%  see also E5_Hands_TPS and BigTutorialLandmarks
%==============================================================================
help(mfilename)