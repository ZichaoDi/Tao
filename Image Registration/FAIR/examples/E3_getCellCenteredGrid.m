%==============================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Tutorial for FAIR: Examples for grid generation
%
%==============================================================================

omega = [0,6,0,4,0,8], m = [3,2,2]
xc = getCellCenteredGrid(omega(1:2),m(1));   xc = reshape(xc,1,[])
xc = getCellCenteredGrid(omega(1:4),m(1:2)); xc = reshape(xc,[m(1:2),2])
xc = getCellCenteredGrid(omega(1:6),m(1:3)); xc = reshape(xc,[m(1:3),3])
