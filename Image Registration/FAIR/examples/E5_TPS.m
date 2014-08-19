%==============================================================================
% (c) Jan Modersitzki 2010/12/28, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR:  Thin plate spline LM registration for hand example
%
% - load data (see setupHandData)
% - setup  viewer (viewImage2D), interpolator (splineInter), 
% - setup landmarks (LM)
% - run affine
%==============================================================================

clear, close all, help(mfilename)

setupHandData; 
xc = getCellCenteredGrid(omega,m);
cc = getTPScoefficients(LM(:,1:4),'theta',0);
[yc,yLM] = evalTPS(LM,cc,xc); 
LM(:,[5,6]) = reshape(yLM,[],2); 
P5_LM; % for nice plots
